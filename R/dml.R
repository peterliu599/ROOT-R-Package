#' Cross fitted estimation of pseudo outcomes for two sample Double ML
#'
#' Trains nuisance models on each training fold and computes pseudo outcomes on
#' the corresponding test fold, then aggregates results. Returns the pseudo outcome table and aligned evaluation data.
#'
#' @param data A \code{data.frame} containing at least the outcome, treatment, and sample indicator columns.
#' @param outcome A \code{character(1)} name of the outcome column.
#' @param treatment A \code{character(1)} name of the treatment column with values \code{0} or \code{1}.
#' @param sample A \code{character(1)} name of the sample indicator column with values \code{0} or \code{1}.
#' @param crossfit An \code{integer(1)} number of folds for cross fitting where the value is at least \code{2}. Default \code{5}.
#'
#' @return A \code{list} with:
#' \item{df_v}{\code{data.frame} with one row per kept observation indexed by \code{primary_index}. Contains \code{te}, \code{a}, \code{b}, \code{te_sq}, \code{a_sq}. Only \code{S == 1} rows with finite values are kept.}
#' \item{data2}{\code{data.frame} subset of the original \code{data} corresponding to \code{df_v$primary_index}.}
#'
#' @note Rows with infinite or undefined weights are removed. Squared deviation columns are centered in the \code{S == 1} group.
estimate_dml <- function(data, outcome, treatment, sample, crossfit = 5) {
  # Input validation
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(data)) stop(sprintf("Column '%s' not found in data.", col), call. = FALSE)
  }

  # Allow NA in Y/Tr for S==0. Require:
  # - no NA in S
  # - no NA in covariates for all rows (used by pi_m)
  # - no NA in Y/Tr *within S==1* (used by a_i)
  if (any(is.na(data[[sample]]))) {
    stop("`", sample, "` contains NA.", call. = FALSE)
  }
  covariate_cols <- setdiff(names(data), c(outcome, treatment, sample))
  if (length(covariate_cols) > 0 && anyNA(data[, covariate_cols, drop = FALSE])) {
    stop("Covariates contain NA. Please impute/drop before training.", call. = FALSE)
  }
  s1 <- data[[sample]] == 1
  if (anyNA(data[s1, c(outcome, treatment), drop = FALSE])) {
    stop("`", outcome, "` or `", treatment, "` contains NA among S==1 rows.", call. = FALSE)
  }

  if (!is.numeric(crossfit) || length(crossfit) != 1 || crossfit < 2) {
    stop("`crossfit` must be an integer >= 2.", call. = FALSE)
  }
  crossfit <- as.integer(crossfit)
  if (all(data[[sample]] == 0) || all(data[[sample]] == 1)) {
    stop("Sample indicator `", sample, "` has no variation (all 0 or all 1). Cannot perform cross-fitting.", call. = FALSE)
  }

  n <- nrow(data)
  folds <- stratified_kfold(data[[sample]], K = crossfit)

  df_list <- vector("list", length(folds))

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)

    training_data <- data[train_idx, , drop = FALSE]
    testing_data <- data[test_idx, , drop = FALSE]

    # Train nuisances on training fold
    trn <- train(training_data, outcome = outcome, treatment = treatment, sample = sample)

    # Pseudo-outcomes on test fold
    est <- estimate(testing_data,
      outcome   = outcome,
      treatment = treatment,
      sample    = sample,
      pi        = trn$pi,
      pi_m      = trn$pi_m,
      e_m       = trn$e_m
    )

    df_list[[i]] <- data.frame(
      te            = as.numeric(est$v),
      primary_index = test_idx,
      a             = as.numeric(est$a),
      b             = as.numeric(est$b),
      row.names     = NULL
    )
  }

  # Combine folds
  df_v_all <- do.call(rbind, df_list)

  # Handle Inf/NA and drop bad rows (coerce to matrix so mocked is.infinite works)
  mat <- as.matrix(df_v_all)
  inf_mask <- is.infinite(mat)
  if (any(inf_mask)) mat[inf_mask] <- NA_real_
  df_v_all[] <- mat  # write back into the data.frame without changing its structure
  df_v_all <- stats::na.omit(df_v_all)

  # Means of te and a over S==1 group
  s_idx1 <- which(data[[sample]] == 1)
  te_mean_s1 <- mean(df_v_all$te[df_v_all$primary_index %in% s_idx1])
  a_mean_s1 <- mean(df_v_all$a[df_v_all$primary_index %in% s_idx1])

  # Add squared deviation columns
  df_v_all$te_sq <- (df_v_all$te - te_mean_s1)^2
  df_v_all$a_sq <- (df_v_all$a - a_mean_s1)^2

  # Aggregate (each index appears once with proper K-folds, but safe to aggregate)
  df_v_grp <- stats::aggregate(df_v_all[, c("te", "a", "b", "te_sq", "a_sq")],
    by = list(primary_index = df_v_all$primary_index),
    FUN = mean
  )

  # Keep only S==1 indices
  df_v_grp <- df_v_grp[df_v_grp$primary_index %in% s_idx1, , drop = FALSE]

  # Align original data rows
  data_s1 <- data[df_v_grp$primary_index, , drop = FALSE]

  list(df_v = df_v_grp, data2 = data_s1)
}

#' Train nuisance models for weighting
#'
#' Fits models to estimate sampling and treatment propensities on training data by logistic regression.
#'
#' @param training_data A \code{data.frame} that contains the training dataset.
#' @param outcome A \code{character(1)} name of the outcome column. Typically this is the observed outcome such as \code{"Yobs"}.
#' @param treatment A \code{character(1)} name of the treatment indicator column such as \code{"Tr"}.
#' @param sample A \code{character(1)} name of the sample inclusion indicator column such as \code{"S"}.
#'
#' @return A \code{list} with:
#' \item{pi}{\code{numeric(1)} prevalence of \code{S == 1} in the training data.}
#' \item{pi_m}{\code{glm} model with binomial family for \eqn{P(S=1 \mid X)}.}
#' \item{e_m}{\code{glm} model with binomial family for \eqn{P(Tr=1 \mid X, S=1)}.}
train <- function(training_data, outcome, treatment, sample) {
  # Helper function
  to01 <- function(x, allow_na = TRUE) {
    # Coerce common encodings to 0/1; NA handling controlled by allow_na
    if (is.numeric(x)) {
      out <- as.integer(x)
    } else if (is.logical(x)) {
      out <- as.integer(x)
    } else {
      if (is.factor(x)) x <- as.character(x)
      x_chr <- trimws(tolower(as.character(x)))
      map1 <- c("1", "yes", "treated", "t", "true")
      map0 <- c("0", "no", "control", "c", "false")
      out <- rep(NA_integer_, length(x_chr))
      out[x_chr %in% map1] <- 1L
      out[x_chr %in% map0] <- 0L
      # direct numeric strings
      suppressWarnings({
        need_num <- is.na(out) & !is.na(x_chr)
        out[need_num] <- as.integer(x_chr[need_num])
      })
    }
    if (!allow_na && any(is.na(out))) {
      bad <- unique(as.character(x)[is.na(out)])
      stop("Non 0/1 encodings found and NA not allowed: ", paste(bad, collapse = ", "), call. = FALSE)
    }
    out
  }

  # Input validation
  if (!is.data.frame(training_data)) {
    stop("`training_data` must be a data frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(training_data)) {
      stop(sprintf("Column '%s' not found in training_data.", col), call. = FALSE)
    }
  }

  # Partition covariates (everything except Y/Tr/S)
  covariate_cols <- setdiff(names(training_data), c(outcome, treatment, sample))

  if (any(is.na(training_data[[sample]]))) {
    stop("`", sample, "` contains NA.", call. = FALSE)
  }
  if (length(covariate_cols) > 0 && anyNA(training_data[, covariate_cols, drop = FALSE])) {
    stop("Covariates contain NA. Please impute/drop before training.", call. = FALSE)
  }

  # Coerce S and Tr to 0/1
  S_vec <- to01(training_data[[sample]], allow_na = FALSE)
  training_data[[sample]] <- S_vec

  Tr_vec <- training_data[[treatment]]
  # Allow NA in Tr where S==0; enforce 0/1 only on S==1 rows
  Tr_vec_s1 <- Tr_vec[S_vec == 1]
  if (!is.numeric(Tr_vec_s1)) Tr_vec_s1 <- to01(Tr_vec_s1, allow_na = FALSE)
  if (!all(Tr_vec_s1 %in% c(0L, 1L))) {
    stop("Treatment indicator column must be numeric 0/1 among S==1 observations.", call. = FALSE)
  }
  # Write back coerced values for S==1 rows
  if (!is.numeric(training_data[[treatment]])) {
    Tr_vec_full <- training_data[[treatment]]
    if (!is.numeric(Tr_vec_full)) {
      # Coerce entire Tr but allow NA for S==0
      Tr_all <- to01(Tr_vec_full, allow_na = TRUE)
      training_data[[treatment]] <- Tr_all
      Tr_vec <- Tr_all
    }
  }

  # Enforce no NA in Y/Tr among S==1
  if (anyNA(training_data[S_vec == 1, c(outcome, treatment), drop = FALSE])) {
    stop("`", outcome, "` or `", treatment, "` contains NA among S==1 rows.", call. = FALSE)
  }

  # Variation checks
  if (all(S_vec == 0) || all(S_vec == 1)) {
    stop("Cannot train model: sample indicator column `", sample, "` has no variation.", call. = FALSE)
  }
  Tr_s1_num <- as.integer(training_data[[treatment]][S_vec == 1])
  if (all(Tr_s1_num == 0) || all(Tr_s1_num == 1)) {
    stop("Cannot train model: treatment has no variation among S==1 observations.", call. = FALSE)
  }

  # Build design matrices
  X_df <- training_data[, covariate_cols, drop = FALSE]

  # Compute sample prevalence
  pi <- mean(S_vec)

  # Fit logistic model for P(S=1 | X)
  data_pi <- data.frame(S = S_vec, X_df, check.names = FALSE)
  pi_m <- stats::glm(formula = S ~ ., data = data_pi, family = stats::binomial())

  # Fit logistic model for P(Tr=1 | X, S=1) using only S==1 subset
  data_e <- data.frame(
    training_data[S_vec == 1, c(treatment, covariate_cols), drop = FALSE],
    check.names = FALSE
  )
  e_m <- stats::glm(
    formula = stats::as.formula(paste(treatment, "~ .")),
    data = data_e, family = stats::binomial()
  )

  list(pi = pi, pi_m = pi_m, e_m = e_m)
}


#' Compute pseudo outcome components \code{a} and \code{b} and their product \code{v}
#'
#' Uses nuisance model outputs to compute inverse probability weighted quantities for ATE in the trial sample.
#'
#' @param testing_data A \code{data.frame} that contains at least outcome, treatment, and sample indicator columns.
#' @param outcome A \code{character(1)} name of the outcome column in \code{testing_data}.
#' @param treatment A \code{character(1)} name of the treatment column in \code{testing_data} with values \code{0} or \code{1}.
#' @param sample A \code{character(1)} name of the sample indicator column in \code{testing_data} with values \code{0} or \code{1}.
#' @param pi A \code{numeric(1)} giving the estimated prevalence \eqn{P(S=1)} from the training data.
#' @param pi_m A fitted \code{glm} model for \eqn{P(S=1 \mid X)}.
#' @param e_m A fitted \code{glm} model for \eqn{P(Tr=1 \mid X, S=1)}.
#'
#' @return A \code{list} with \code{numeric} vectors of length \code{nrow(testing_data)}:
#' \item{v}{Pseudo outcome values.}
#' \item{a}{IPW adjusted outcome contrast.}
#' \item{b}{Overlap weight factors.}
#'
#' @details Predictions from \code{pi_m} and \code{e_m} are clamped to \code{[1e-8, 1 - 1e-8]} for stability.
estimate <- function(testing_data, outcome, treatment, sample, pi, pi_m, e_m) {
  # Input validation
  if (!is.data.frame(testing_data)) {
    stop("`testing_data` must be a data frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment, sample)) {
    if (!col %in% names(testing_data)) {
      stop(sprintf("Column '%s' not found in testing_data.", col), call. = FALSE)
    }
  }

  # Allow NA in Y/Tr for S==0. Require:
  # - no NA in S
  # - no NA in covariates for all rows (used by pi_m)
  # - no NA in Y/Tr *within S==1* (used by a_i)
  if (any(is.na(testing_data[[sample]]))) {
    stop("`", sample, "` contains NA in testing data.", call. = FALSE)
  }
  covariate_cols <- setdiff(names(testing_data), c(outcome, treatment, sample))
  if (length(covariate_cols) > 0 && anyNA(testing_data[, covariate_cols, drop = FALSE])) {
    stop("Covariates contain NA in testing data.", call. = FALSE)
  }
  s1 <- testing_data[[sample]] == 1
  if (anyNA(testing_data[s1, c(outcome, treatment), drop = FALSE])) {
    stop("`", outcome, "` or `", treatment, "` contains NA among S==1 in testing data.", call. = FALSE)
  }

  if (!is.numeric(pi) || length(pi) != 1) {
    stop("`pi` must be a numeric scalar.", call. = FALSE)
  }
  if (pi <= 0 || pi >= 1) {
    stop("Invalid `pi`: P(S=1) must be between 0 and 1 (exclusive).", call. = FALSE)
  }
  if (!inherits(pi_m, "glm") || !inherits(e_m, "glm")) {
    stop("`pi_m` and `e_m` must be model objects (e.g. from glm) for prediction.", call. = FALSE)
  }

  # Extract columns from testing_data
  S_vec <- testing_data[[sample]]
  Y_vec <- testing_data[[outcome]]
  Tr_vec <- testing_data[[treatment]]
  # Data frame of covariates (all other columns)
  covariate_cols <- setdiff(names(testing_data), c(outcome, treatment, sample))
  X_df <- testing_data[, covariate_cols, drop = FALSE]

  # Predict P(S=1 | X) and compute overlap factor lX
  p_s1x <- stats::predict(pi_m, newdata = X_df, type = "response")
  # Clamp predicted probabilities to avoid 0 or 1
  p_s1x <- pmin(pmax(as.numeric(p_s1x), 1e-8), 1 - 1e-8)
  lX <- (p_s1x / pi) / ((1 - p_s1x) / (1 - pi))

  # Predict P(Tr=1 | X, S=1) for all observations (treat those with S=0 similarly)
  p_t1x <- stats::predict(e_m, newdata = X_df, type = "response")
  p_t1x <- pmin(pmax(as.numeric(p_t1x), 1e-8), 1 - 1e-8)

  # Compute IPW components
  a_val <- ifelse(
    S_vec == 1,
    Tr_vec * Y_vec / p_t1x - (1 - Tr_vec) * Y_vec / (1 - p_t1x),
    0
  )
  b_val <- 1 / lX
  v_val <- a_val * b_val

  return(list(v = as.numeric(v_val), a = as.numeric(a_val), b = as.numeric(b_val)))
}

#' Stratified K fold index generator
#'
#' Splits indices into \code{K} folds while preserving the class distribution of a binary factor.
#' This mimics the behavior of stratified K fold allocation to keep the ratio of classes in each fold.
#'
#' @param S A \code{vector} or \code{factor} indicating class membership for stratification. Typical values are \code{0} or \code{1}.
#' @param K An \code{integer(1)} number of folds. If \code{K} exceeds the number of observations it is reduced to that number.
#'
#' @return A \code{list} of length \code{K} where each element is an \code{integer} vector of row indices assigned to that fold.
#' The union of all folds equals \code{seq_along(S)} and folds are close in size.
stratified_kfold <- function(S, K = 5) {
  # Input validation
  if (is.data.frame(S) || is.matrix(S)) {
    stop("`S` should be a vector or factor, not a data frame or matrix.", call. = FALSE)
  }
  if (!is.factor(S)) {
    # Coerce to factor; this works for numeric 0/1 as well
    S <- as.factor(S)
  }
  if (!is.numeric(K) || length(K) != 1 || K < 1) {
    stop("`K` must be a positive integer.", call. = FALSE)
  }
  K <- as.integer(K)
  if (K > length(S)) {
    warning("Requested K (", K, ") is greater than number of observations; using K = ", length(S))
    K <- length(S)
  }

  # Split indices by class
  idx_by_class <- split(seq_along(S), S)

  # For each class, distribute its indices across up to K parts in round-robin.
  parts_by_class <- lapply(idx_by_class, function(idx) {
    # If a class has fewer obs than K, it will have fewer parts.
    k_class <- min(K, max(1L, length(idx)))
    split(idx, rep_len(seq_len(k_class), length(idx)))
  })

  # Build K folds by taking the k-th part from each class when present.
  folds <- vector("list", K)
  for (k in seq_len(K)) {
    fold_k <- unlist(lapply(parts_by_class, function(p) {
      if (length(p) >= k) p[[k]] else integer(0)
    }), use.names = FALSE)
    folds[[k]] <- sort(fold_k)
  }

  folds
}

#' Train treatment propensity model for single sample mode
#'
#' Fits a logistic regression model for \eqn{P(T = 1 \mid X)} on the provided training data.
#' Used by the single sample Double ML path where no sample selection model is required.
#'
#' @param training_data A \code{data.frame} containing the outcome, treatment, and covariates.
#'   Only \code{treatment} and covariates are used for fitting.
#' @param outcome A \code{character(1)} name of the outcome column. Present for a consistent signature and not used here.
#' @param treatment A \code{character(1)} name of the binary treatment indicator column with values \code{0} or \code{1}.
#'
#' @return A \code{list} with:
#' \item{e_m}{A \code{glm} model with binomial family for the treatment propensity.}
train_single <- function(training_data, outcome, treatment) {
  covariate_cols <- setdiff(names(training_data), c(outcome, treatment))
  X_tr <- training_data[, covariate_cols, drop = FALSE]

  df <- data.frame(check.names = FALSE, row.names = rownames(training_data))
  df[[treatment]] <- training_data[[treatment]]   # <- give it the correct name
  df <- cbind(df, X_tr)

  e_m <- stats::glm(
    formula = stats::as.formula(paste(treatment, "~ .")),
    data = df,
    family = stats::binomial()
  )
  list(e_m = e_m)
}

#' Compute single sample pseudo outcomes
#'
#' Computes single sample pseudo outcome components for ATE style estimation:
#' \eqn{a_i = T_i Y_i / e_i - (1 - T_i) Y_i / (1 - e_i)} with \eqn{v_i = a_i} and \eqn{b_i \equiv 1}.
#'
#' @param testing_data A \code{data.frame} containing at least \code{outcome}, \code{treatment}, and covariates.
#' @param outcome A \code{character(1)} name of the outcome column.
#' @param treatment A \code{character(1)} name of the binary treatment indicator column with values \code{0} or \code{1}.
#' @param e_m A fitted \code{glm} model for \eqn{P(T = 1 \mid X)} such as the result of \code{\link{train_single}}.
#'
#' @return A \code{list} with \code{numeric} vectors of length \code{nrow(testing_data)}:
#' \item{v}{Pseudo outcome values which equal \code{a} in single sample mode.}
#' \item{a}{IPW adjusted outcome contrast.}
#' \item{b}{Vector of ones because there is no sample overlap weighting in single sample mode.}
estimate_single <- function(testing_data, outcome, treatment, e_m) {
  covariate_cols <- setdiff(names(testing_data), c(outcome, treatment))
  X_te <- testing_data[, covariate_cols, drop = FALSE]
  Y <- testing_data[[outcome]]
  Tr <- testing_data[[treatment]]

  p_t1x <- stats::predict(e_m, newdata = X_te, type = "response")
  p_t1x <- pmin(pmax(as.numeric(p_t1x), 1e-8), 1 - 1e-8)

  a <- Tr * Y / p_t1x - (1 - Tr) * Y / (1 - p_t1x)
  # Single-sample: b ≡ 1, v = a
  list(v = as.numeric(a), a = as.numeric(a), b = rep(1, nrow(testing_data)))
}

#' Cross fitted Double ML for single sample mode
#'
#' Runs K fold cross fitting to produce pseudo outcomes for ATE estimation when
#' no sample membership indicator is available or has no variation.
#'
#' @param data A \code{data.frame} containing \code{outcome}, \code{treatment}, and covariates.
#' @param outcome A \code{character(1)} name of the outcome column.
#' @param treatment A \code{character(1)} name of the binary treatment indicator column with values \code{0} or \code{1}.
#' @param crossfit An \code{integer(1)} number of folds for cross fitting where the value is at least \code{2}. Default \code{5}.
#'
#' @return A \code{list} with:
#' \item{df_v}{\code{data.frame} with one row per kept observation that contains \code{te}, \code{a}, \code{b}, \code{te_sq}, \code{a_sq}, and \code{primary_index}.}
#' \item{data2}{\code{data.frame} subset of \code{data} that corresponds to \code{df_v$primary_index}.}
estimate_dml_single <- function(data, outcome, treatment, crossfit = 5) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (col in c(outcome, treatment)) {
    if (!col %in% names(data)) stop(sprintf("Column '%s' not found in data.", col), call. = FALSE)
  }
  .check_no_na(data, colnames(data))
  if (!is.numeric(crossfit) || length(crossfit) != 1 || crossfit < 2) {
    stop("`crossfit` must be an integer >= 2.", call. = FALSE)
  }
  crossfit <- as.integer(crossfit)

  n <- nrow(data)
  # Reuse stratified_kfold for convenience; stratify on a constant is fine
  folds <- stratified_kfold(S = rep(1L, n), K = crossfit)

  df_list <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)
    trn <- train_single(data[train_idx, , drop = FALSE], outcome, treatment)
    est <- estimate_single(data[test_idx, , drop = FALSE], outcome, treatment, trn$e_m)

    df_list[[i]] <- data.frame(
      te            = est$v,
      primary_index = test_idx,
      a             = est$a,
      b             = est$b,
      row.names     = NULL
    )
  }

  df_v_all <- do.call(rbind, df_list)
  # Handle Inf/NA and drop bad rows (coerce to matrix so mocked is.infinite works)
  mat <- as.matrix(df_v_all)
  inf_mask <- is.infinite(mat)
  if (any(inf_mask)) mat[inf_mask] <- NA_real_
  df_v_all[] <- mat  # write back into the data.frame without changing its structure
  df_v_all <- stats::na.omit(df_v_all)

  te_mean <- mean(df_v_all$te)
  a_mean <- mean(df_v_all$a)
  df_v_all$te_sq <- (df_v_all$te - te_mean)^2
  df_v_all$a_sq <- (df_v_all$a - a_mean)^2

  df_v_grp <- stats::aggregate(df_v_all[, c("te", "a", "b", "te_sq", "a_sq")],
    by = list(primary_index = df_v_all$primary_index),
    FUN = mean
  )

  list(df_v = df_v_grp, data2 = data[df_v_grp$primary_index, , drop = FALSE])
}
