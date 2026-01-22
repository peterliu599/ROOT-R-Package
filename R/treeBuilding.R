#' Compute transport influence scores for generalization mode
#'
#' Internal helper used in the \code{generalization} path to construct
#' glm-based, IPW-style scores for transporting trial effects to a target
#' population without double machine learning.
#'
#' The function treats \code{data} as a stacked dataset with a sample
#' indicator \eqn{S} (\code{sample}) taking value 1 in the randomized trial
#' and 0 in the target sample. It proceeds in three steps:
#' \enumerate{
#'   \item Fit a sampling model \eqn{P(S = 1 | X)} using logistic regression
#'         on all rows of \code{data}.
#'   \item Within the trial subset \eqn{S = 1}, fit a treatment model
#'         \eqn{P(T = 1 | X, S = 1)} using logistic regression.
#'   \item For each row, form the density ratio
#'         \eqn{r(X) = P(S = 0 | X) / P(S = 1 | X)}
#'         and compute a Horvitz-Thompson-style transported score
#'         \deqn{v(X, T, Y) = r(X) \left[ \frac{T Y}{e(X)} - \frac{(1 - T) Y}{1 - e(X)} \right],}
#'         where \eqn{e(X) = P(T = 1 | X, S = 1)}.
#' }
#'
#' The resulting score vector \code{v} and its squared version \code{vsq}
#' can be used as pseudo-outcomes for tree-based search over transported
#' treatment effects.
#'
#' @param data A \code{data.frame} containing the outcome, treatment indicator,
#'   sample indicator, and covariates. All columns except \code{outcome},
#'   \code{treatment}, and \code{sample} are treated as covariates \eqn{X}
#'   and must be suitable for use in \code{stats::glm()} with
#'   \code{family = binomial()}.
#' @param outcome A length-1 character string giving the name of the outcome
#'   column in \code{data}. Must be numeric (e.g., a continuous outcome or
#'   a 0/1 indicator).
#' @param treatment A length-1 character string giving the name of the
#'   treatment indicator column in \code{data}. Must be coded 0/1.
#' @param sample A length-1 character string giving the name of the sample
#'   indicator column in \code{data}, with 1 for trial rows and 0 for target
#'   rows.
#'
#' @return A list with two numeric vectors of length \code{nrow(data)}:
#' \describe{
#'   \item{\code{v}}{The transported influence-style score.}
#'   \item{\code{vsq}}{The element-wise square of \code{v}, i.e., \code{v^2}.}
#' }
compute_transport_scores <- function(data, outcome, treatment, sample) {
  covars <- setdiff(colnames(data), c(outcome, treatment, sample))
  if (length(covars) == 0L) {
    stop("compute_transport_scores(): need at least one covariate.", call. = FALSE)
  }

  # P(S = 1 | X)
  formula_s <- stats::as.formula(paste(sample, "~", paste(covars, collapse = "+")))
  model_s   <- stats::glm(formula_s, data = data, family = stats::binomial())
  pi_s      <- stats::predict(model_s, type = "response")

  # P(Tr = 1 | X, S = 1)
  trial_data <- data[data[[sample]] == 1, , drop = FALSE]
  if (nrow(trial_data) == 0L) {
    stop("compute_transport_scores(): no S == 1 (trial) rows for treatment model.",
         call. = FALSE)
  }
  formula_t <- stats::as.formula(paste(treatment, "~", paste(covars, collapse = "+")))
  model_t   <- stats::glm(formula_t, data = trial_data, family = stats::binomial())
  e_x       <- stats::predict(model_t, newdata = data, type = "response")

  Y     <- data[[outcome]]
  T_ind <- data[[treatment]]

  # density ratio r(X) = P(S=0|X) / P(S=1|X)
  r_x <- (1 - pi_s) / pi_s

  # Horvitz–Thompson-style transported score
  v   <- r_x * ((T_ind * Y / e_x) - ((1 - T_ind) * Y / (1 - e_x)))
  vsq <- v^2

  list(v = v, vsq = vsq)
}

#' Recursive split builder for weighted tree
#'
#' Recursively builds a weighted decision tree to optimize a global objective,
#' using an exploration versus exploitation choice at leaves. Internal and used by \code{ROOT()}.
#'
#' @param split_feature A named \code{numeric} vector of feature selection probabilities. Must include the name \code{"leaf"}.
#' @param X A \code{data.frame} of current observations. Includes candidate split feature columns and may include a working copy of weights \code{w}.
#' @param D A \code{data.frame} representing the global state. Must include columns \code{w} and \code{vsq}. Row names must align to \code{X}.
#' @param parent_loss A \code{numeric(1)} giving the loss of the parent node. Used to decide if a split improves the objective.
#' @param depth An \code{integer(1)} giving the current tree depth.
#' @param explore_proba A \code{numeric(1)} between \code{0} and \code{1} for the probability of flipping the exploit choice at a leaf.
#' @param choose_feature_fn A \code{function} to choose the next feature. Default is \code{choose_feature}.
#' @param reduce_weight_fn A \code{function} to penalize the last tried feature on a rejected split. Default is \code{reduce_weight}.
#' @param global_objective_fn A \code{function} with signature \code{function(D) -> numeric} that scores the entire state.
#' @param max_depth An \code{integer(1)} giving the maximum depth. A node becomes a leaf at this depth.
#' @param min_leaf_n An \code{integer(1)} giving the minimum number of rows to attempt a split. Otherwise make a leaf.
#' @param log_fn A \code{function} for logging. Default is a function that performs no operation.
#' @param max_rejects_per_node An \code{integer(1)} giving the safety budget of rejected splits before forcing a leaf.
#'
#' @return A \code{list} representing the subtree. Includes updated \code{D} and a field named \code{"local objective"}.
#' @importFrom stats rbinom
split_node <- function(split_feature, X, D, parent_loss, depth,
                       explore_proba = 0.05,
                       choose_feature_fn = choose_feature,   # (renamed arg for clarity)
                       reduce_weight_fn = reduce_weight,
                       global_objective_fn = objective_default,
                       max_depth = 8,
                       min_leaf_n = 5,
                       log_fn = function(...) {},
                       max_rejects_per_node = 1000) {
  # Derive micro-evaluator: loss(val, indices, D) from global objective
  loss_fn <- loss_from_objective(global_objective_fn)

  # Helpers
  .log <- function(fmt, ...) log_fn(sprintf(fmt, ...))
  nearly_leq <- function(a,b,tol_abs=1e-12,tol_rel=1e-12) a <= b + tol_abs + tol_rel*max(1,abs(b))

  make_leaf <- function(X_sub, D_sub, depth_cur, reason, fj = NA_character_, cj = NA_real_) {
    losses <- c(loss_fn(0, rownames(X_sub), D_sub), loss_fn(1, rownames(X_sub), D_sub))
    w_exploit <- which.min(losses) - 1
    w_explore <- stats::rbinom(1, 1, 0.5)
    explore_flip <- stats::rbinom(1, 1, explore_proba)
    final_w <- if (explore_flip == 1) w_explore else w_exploit

    .log("[depth=%d] LEAF (%s): feature=%s, cut=%s, n=%d, losses={%.4f, %.4f}, w=%d",
         depth_cur, reason, as.character(fj),
         ifelse(is.na(cj), "NA", sprintf("%.4f", cj)),
         nrow(X_sub), losses[1], losses[2], final_w)

    if (nrow(X_sub) > 0) {
      D_sub[rownames(X_sub), "w"] <- final_w
      X_sub$w <- final_w
    }
    list(
      node = "leaf", w = final_w,
      `local objective` = min(losses),
      depth = depth_cur, D = D_sub,
      leaf_reason = reason, feature = fj, cut = cj
    )
  }

  # Stopping
  if (depth >= max_depth) return(make_leaf(X, D, depth, reason = "max-depth"))
  if (nrow(X) <= min_leaf_n) return(make_leaf(X, D, depth, reason = "min-leaf"))

  fj <- choose_feature_fn(split_feature, depth)
  if (identical(fj, "leaf")) return(make_leaf(X, D, depth, reason = "feature==leaf", fj = fj))

  cj <- midpoint(X[[fj]])
  X_left  <- X[X[[fj]] <= cj, , drop = FALSE]
  X_right <- X[X[[fj]] >  cj, , drop = FALSE]
  if (nrow(X_left) == 0 || nrow(X_right) == 0) {
    return(make_leaf(X, D, depth, reason = "empty-child", fj = fj, cj = cj))
  }

  # Evaluate best decision on children leaves
  loss_left  <- c(loss_fn(0, rownames(X_left),  D), loss_fn(1, rownames(X_left),  D))
  loss_right <- c(loss_fn(0, rownames(X_right), D), loss_fn(1, rownames(X_right), D))
  min_left   <- min(loss_left)
  min_right  <- min(loss_right)
  new_loss   <- (nrow(X_left) * min_left + nrow(X_right) * min_right) / nrow(X)

  if (nearly_leq(new_loss, parent_loss)) {
    .log("[depth=%d] SPLIT: feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.6f, parent_loss=%.6f",
         depth, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)
    .log("    left losses={%.6f, %.6f}, right losses={%.6f, %.6f}",
         loss_left[1], loss_left[2], loss_right[1], loss_right[2])

    w_left  <- which.min(loss_left)  - 1
    w_right <- which.min(loss_right) - 1
    D[rownames(X_left),  "w"] <- w_left;  X_left$w  <- w_left
    D[rownames(X_right), "w"] <- w_right; X_right$w <- w_right

    # Recurse (random order to break symmetry)
    if (stats::rbinom(1, 1, 0.5) == 1) {
      left_res  <- split_node(split_feature, X_left,  D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- left_res$D
      right_res <- split_node(split_feature, X_right, D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- right_res$D
    } else {
      right_res <- split_node(split_feature, X_right, D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- right_res$D
      left_res  <- split_node(split_feature, X_left,  D, new_loss, depth + 1,
                              explore_proba, choose_feature_fn, reduce_weight_fn,
                              global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
      D <- left_res$D
    }

    local_obj <- global_objective_fn(D); if (!is.finite(local_obj) || is.nan(local_obj)) local_obj <- Inf
    return(list(
      node = fj, split = cj,
      left_tree  = left_res, right_tree = right_res,
      `local objective` = local_obj, depth = depth, D = D
    ))
  }

  # Rejection loop with budget (stack-safe)
  .log("[depth=%d] REJECTED: feature=%s, cut=%.4f | new_loss=%.16f > parent_loss=%.16f",
       depth, fj, cj, new_loss, parent_loss)

  rej_sf <- split_feature
  for (attempt in seq_len(max_rejects_per_node)) {
    # Penalize the last tried feature
    rej_sf <- reduce_weight_fn(fj, rej_sf)

    # If "leaf" is now effectively certain, bail out as leaf
    if ("leaf" %in% names(rej_sf) && rej_sf["leaf"] >= 1 - 1e-8) {
      return(make_leaf(X, D, depth, reason = "forced-leaf-after-rejects", fj = fj, cj = cj))
    }

    # Choose a new feature and recompute children
    fj <- choose_feature_fn(rej_sf, depth)
    if (identical(fj, "leaf")) {
      return(make_leaf(X, D, depth, reason = "feature==leaf-after-rejects", fj = fj))
    }

    cj <- midpoint(X[[fj]])
    X_left  <- X[X[[fj]] <= cj, , drop = FALSE]
    X_right <- X[X[[fj]] >  cj, , drop = FALSE]
    if (nrow(X_left) == 0 || nrow(X_right) == 0) {
      return(make_leaf(X, D, depth, reason = "empty-child-after-rejects", fj = fj, cj = cj))
    }

    loss_left  <- c(loss_fn(0, rownames(X_left),  D), loss_fn(1, rownames(X_left),  D))
    loss_right <- c(loss_fn(0, rownames(X_right), D), loss_fn(1, rownames(X_right), D))
    min_left   <- min(loss_left); min_right <- min(loss_right)
    new_loss   <- (nrow(X_left) * min_left + nrow(X_right) * min_right) / nrow(X)

    if (nearly_leq(new_loss, parent_loss)) {
      .log("[depth=%d] SPLIT(after %d rejects): feature=%s, cut=%.4f | n_left=%d, n_right=%d | new_loss=%.6f, parent_loss=%.6f",
           depth, attempt, fj, cj, nrow(X_left), nrow(X_right), new_loss, parent_loss)

      w_left  <- which.min(loss_left)  - 1
      w_right <- which.min(loss_right) - 1
      D[rownames(X_left),  "w"] <- w_left;  X_left$w  <- w_left
      D[rownames(X_right), "w"] <- w_right; X_right$w <- w_right

      if (stats::rbinom(1, 1, 0.5) == 1) {
        left_res  <- split_node(rej_sf, X_left,  D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- left_res$D
        right_res <- split_node(rej_sf, X_right, D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- right_res$D
      } else {
        right_res <- split_node(rej_sf, X_right, D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- right_res$D
        left_res  <- split_node(rej_sf, X_left,  D, new_loss, depth + 1,
                                explore_proba, choose_feature_fn, reduce_weight_fn,
                                global_objective_fn, max_depth, min_leaf_n, log_fn, max_rejects_per_node)
        D <- left_res$D
      }

      local_obj <- global_objective_fn(D); if (!is.finite(local_obj) || is.nan(local_obj)) local_obj <- Inf
      return(list(
        node = fj, split = cj,
        left_tree  = left_res, right_tree = right_res,
        `local objective` = local_obj, depth = depth, D = D
      ))
    }

    .log("[depth=%d] REJECTED (%d/%d): feature=%s, cut=%.4f | new_loss=%.16f > parent_loss=%.16f",
         depth, attempt, max_rejects_per_node, fj, cj, new_loss, parent_loss)
  }

  # Budget exhausted -> become a leaf
  make_leaf(X, D, depth, reason = "max-rejects", fj = fj, cj = cj)
}

#' Randomly choose a split feature based on provided probabilities
#'
#' Selects one feature name at random according to a probability vector that may include a special \code{"leaf"} entry.
#'
#' @param split_feature A named \code{numeric} vector of feature selection probabilities. Names correspond to feature identifiers and may include \code{"leaf"}.
#' @param depth An \code{integer(1)} giving the current tree depth. Present for parity with the Python version and does not change probabilities here.
#'
#' @return A \code{character(1)} which is the chosen feature name or \code{"leaf"}.
#' @note The factor \eqn{2^{(0 \cdot \mathrm{depth} / 4)}} present in the code equals \code{1} and does not change the first element weight. All probabilities are normalized to sum to \code{1} before sampling.
choose_feature <- function(split_feature, depth) {
  # Input validation
  if (!is.numeric(split_feature) || length(split_feature) == 0) {
    stop("`split_feature` must be a non-empty numeric vector of probabilities.", call. = FALSE)
  }
  if (is.null(names(split_feature))) {
    stop("`split_feature` must have names for each feature (use 'leaf' for the special leaf option if applicable).", call. = FALSE)
  }
  if (!is.numeric(depth) || length(depth) != 1) {
    stop("`depth` must be a numeric scalar (current tree depth).", call. = FALSE)
  }

  # New strict probability checks
  if (anyNA(split_feature)) {
    stop("`split_feature` contains NA values.", call. = FALSE)
  }
  if (any(split_feature < 0)) {
    bad <- names(split_feature)[split_feature < 0]
    stop(sprintf(
      "All probabilities must be >= 0. Negative entries found in: %s",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  split_prob <- as.numeric(split_feature)
  # Adjust the first element probability by the (inactive) depth factor (kept for parity)
  if (length(split_prob) > 0) {
    split_prob[1] <- split_prob[1] * (2^(0 * depth / 4))
  }
  # Normalize to a probability distribution
  total <- sum(split_prob)
  if (!is.finite(total)) {
    stop("Invalid probabilities: sum is non-finite.", call. = FALSE)
  }
  split_prob <- split_prob / total

  # Perform random sampling
  choices <- names(split_feature)
  chosen <- sample(choices, size = 1, replace = TRUE, prob = split_prob)
  return(chosen)
}

#' Reduce a feature selection weight by one half and renormalize
#'
#' Lowers the probability weight of a given feature by one half and then renormalizes the full vector to sum to one.
#'
#' @param fj A \code{character(1)} feature name that must be present in \code{names(split_feature)}.
#' @param split_feature A named \code{numeric} vector of probabilities for features as used in splitting.
#'
#' @return A \code{numeric} vector of the same length as \code{split_feature} that sums to \code{1}.
#' @details This is used when a feature split was rejected. The feature probability is halved to reduce the chance of immediate reselection which encourages exploration of other features. If \code{fj} equals \code{"leaf"} its weight is also halved.
reduce_weight <- function(fj, split_feature) {
  # Input validation
  if (!is.character(fj) || length(fj) != 1) {
    stop("`fj` must be a single feature name (string).", call. = FALSE)
  }
  if (!is.numeric(split_feature) || is.null(names(split_feature))) {
    stop("`split_feature` must be a named numeric vector of probabilities.", call. = FALSE)
  }
  if (!(fj %in% names(split_feature))) {
    stop("Feature '", fj, "' not found in names(split_feature).", call. = FALSE)
  }
  # New strict probability checks
  if (anyNA(split_feature)) {
    stop("`split_feature` contains NA values.", call. = FALSE)
  }
  if (any(split_feature < 0)) {
    bad <- names(split_feature)[split_feature < 0]
    stop(sprintf(
      "All probabilities must be >= 0. Negative entries found in: %s",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  total_before <- sum(split_feature)
  if (!is.finite(total_before)) {
    stop("Invalid probabilities: sum is non-finite.", call. = FALSE)
  }

  # Reduce the specified feature's weight by half
  split_feature[fj] <- split_feature[fj] / 2
  # Renormalize to sum to 1
  split_feature <- split_feature / sum(split_feature)
  return(split_feature)
}

#' Compute the midpoint of a numeric vector
#'
#' Calculates the midpoint defined as \eqn{(\max(X) + \min(X)) / 2} while ignoring any \code{NA} values in \code{X}.
#'
#' @param X A \code{numeric} vector.
#'
#' @return A \code{numeric(1)} giving the midpoint of the finite values in \code{X}. Returns \code{NA_real_} when \code{X} is empty or has no finite values.
midpoint <- function(X) {
  if (!is.numeric(X)) {
    stop("`X` must be a numeric vector.", call. = FALSE)
  }
  if (length(X) == 0L) {
    warning("Empty vector provided to midpoint(). Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  # Work only with finite values to avoid base::range() warnings
  x <- X[is.finite(X)]
  if (length(x) == 0L) {
    warning("No finite values in `X`. Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  (min(x) + max(x)) / 2
}

#' Fit a shallow decision tree to characterize learned weights \code{w}
#'
#' Trains a classification tree on the covariates \code{X} to predict the binary membership \code{w}.
#' This provides an interpretable summary of how the weighted subgroup can be distinguished by \code{X}.
#'
#' @param X A \code{data.frame} of covariates.
#' @param w A \code{vector} of length \code{nrow(X)} that is binary. Accepts \code{0} and \code{1} or a \code{factor} with two levels.
#' @param max_depth An \code{integer(1)} giving the maximum tree depth. Default is \code{3}.
#'
#' @return An \code{rpart} object that represents the fitted classification tree.
#'
#' @details The tree uses the Gini index for classification and no pruning with complexity parameter \code{cp = 0}. Depth control is through \code{max_depth}. If \code{w} is not a factor it is converted internally. The resulting rules indicate which covariates and splits separate the two classes defined by \code{w}.
characterize_tree <- function(X, w, max_depth = 3) {
  # Input validation
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  .check_no_na(X, colnames(X))
  if (length(w) != nrow(X)) {
    stop("Length of `w` must equal the number of rows in `X`.", call. = FALSE)
  }
  # Coerce w to factor and check it has two levels
  w_factor <- as.factor(w)
  if (nlevels(w_factor) != 2) {
    stop("`w` must have exactly two classes (binary).", call. = FALSE)
  }

  # Prepare data for rpart
  df <- data.frame(w = w_factor, X)
  # Fit classification tree
  fit <- rpart::rpart(
    w ~ .,
    data = df,
    method = "class",
    parms = list(split = "gini"),
    control = rpart::rpart.control(maxdepth = max_depth, cp = 0.0, minsplit = 2, minbucket = 1),
    model = TRUE
  )
  return(fit)
}
