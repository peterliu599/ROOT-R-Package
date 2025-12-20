#' Ensemble of weighted trees for general optimization and Rashomon selection
#'
#' Builds multiple weighted trees, then identifies a "Rashomon set" of
#' top-performing trees and aggregates their weight assignments by majority vote.
#'
#' The function is framed as a general functional optimization routine:
#' given data D_n = {d_1,...,d_n} and a loss L(w, D_n), ROOT searches over
#' interpretable tree-based weight functions w(d) ∈ {0,1}.
#'
#' @param data A data.frame containing the dataset.
#'
#'   In *general optimization* mode (`generalizability_path = FALSE`), `data` can be
#'   any set of covariates and auxiliary columns. The user supplies a
#'   `global_objective_fn` that takes a data frame with a column `w`
#'   and returns a scalar loss.
#'
#'   In *generalizability_path* mode (`generalizability_path = TRUE`), `data` must contain
#'   columns `"Y"` (outcome), `"Tr"` (treatment indicator, 0/1),
#'   and `"S"` (sample indicator, 1 = trial, 0 = target). ROOT internally
#'   constructs transportability scores and, if no custom objective is given,
#'   uses a default variance-based loss.
#' @param global_objective_fn function with signature `function(D) -> numeric`
#'   scoring the entire state and minimized by ROOT. If NULL, a default
#'   variance-based objective is used (see `objective_default()`).
#' @param generalizability_path logical(1). If TRUE, use the built-in transportability
#'   objective based on (Y, Tr, S). If FALSE, treat `data` as arbitrary and rely
#'   on `global_objective_fn`. Default FALSE.
#' @param leaf_proba numeric(1) tuning parameter that increases the chance
#'   a node stops splitting by selecting a synthetic \code{"leaf"} feature.
#'   Internally, the probability of choosing \code{"leaf"} is
#'   \eqn{\text{leaf\_proba} / (1 + \text{leaf\_proba})} (assuming the
#'   covariate probabilities sum to 1). Default \code{0.25}.
#' @param seed optional numeric(1) seed for reproducibility.
#' @param num_trees integer(1) number of trees to grow. Default 10.
#' @param vote_threshold numeric(1) in (0.5, 1] giving the majority vote
#'   threshold for final `w = 1`. Default 2/3.
#' @param explore_proba numeric(1) giving the exploration probability at leaves.
#'   Default 0.05.
#' @param feature_est either a character(1) in c("Ridge", "GBM") or a
#'   function(X, y, ...) returning a named nonnegative numeric vector of
#'   importances with names matching columns of X. Used only to bias which
#'   covariates are chosen for splitting. If it fails, ROOT falls back to
#'   uniform feature sampling with a warning.
#' @param feature_est_args list of additional arguments passed to a user
#'   supplied `feature_est` function.
#' @param top_k_trees logical(1). If TRUE, select top `k` trees by objective;
#'   otherwise use `cutoff`. Default FALSE.
#' @param k integer(1) giving the number of top trees when
#'   `top_k_trees = TRUE`. Default 10.
#' @param cutoff numeric(1) or `"baseline"`. Used as the Rashomon cutoff when
#'   `top_k_trees = FALSE`. `"baseline"` uses the objective at w ≡ 1.
#' @param verbose logical(1). If TRUE, prints unweighted and (when available)
#'   weighted estimates and their standard errors in generalizability_path mode.
#'
#' @return An object of class "ROOT" (a list) with elements:
#'   - D_rash: data frame with Rashomon-set votes and w_opt.
#'   - D_forest: data frame with forest-level working columns.
#'   - w_forest: list of per-tree results from split_node().
#'   - rashomon_set: indices of selected trees.
#'   - global_objective_fn: the objective function used.
#'   - f: summary classifier (e.g., rpart tree) or NULL.
#'   - testing_data: data frame aligned to rows used to compute scores.
#'   - estimate: (only if generalizability_path = TRUE) list with unweighted and
#'       weighted estimands, SEs, and a note about the SE.
#'   - generalizability_path: logical flag.
#'
#' @export
ROOT <- function(data,
                 global_objective_fn = NULL,
                 generalizability_path = FALSE,
                 leaf_proba          = 0.25,
                 seed                = NULL,
                 num_trees           = 10,
                 vote_threshold      = 2 / 3,
                 explore_proba       = 0.05,
                 feature_est         = "Ridge",
                 feature_est_args    = list(),
                 top_k_trees         = FALSE,
                 k                   = 10,
                 cutoff              = "baseline",
                 verbose             = FALSE) {

  # ---- Small helpers ------------------------------------------------------
  coerce01 <- function(x, allow_na = TRUE) {
    if (is.numeric(x) || is.logical(x)) return(as.integer(x))
    if (is.factor(x)) x <- as.character(x)
    x_chr <- trimws(tolower(as.character(x)))
    map1 <- c("1", "yes", "treated", "t", "true")
    map0 <- c("0", "no", "control", "c", "false")
    out <- rep(NA_integer_, length(x_chr))
    out[x_chr %in% map1] <- 1L
    out[x_chr %in% map0] <- 0L
    if (!allow_na && any(is.na(out)))
      stop("Non 0/1 values found.", call. = FALSE)
    out
  }

  .norm_feat_prob <- function(imp, X_df) {
    if (!is.numeric(imp) || is.null(names(imp))) {
      stop("Custom feature importance must be a named numeric vector.", call. = FALSE)
    }
    imp <- imp[colnames(X_df)]
    if (any(is.na(imp))) stop("Importance missing for some X_df columns.", call. = FALSE)
    if (any(imp < 0))   stop("Importances must be non-negative.", call. = FALSE)
    s <- sum(imp)
    if (!is.finite(s) || s <= 0) {
      warning(
        "Feature importance was degenerate; now using a simple random probabilistic ",
        "setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- imp / s
    }
    feat_prob
  }

  compute_split_feature <- function(X_df, vsq_vals) {
    if (ncol(X_df) < 1L) {
      stop("compute_split_feature(): no features available for weighting.", call. = FALSE)
    }

    ## 1) User-supplied feature importance function
    if (is.function(feature_est)) {
      imp <- try(
        if (is.null(seed)) {
          do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
        } else {
          withr::with_seed(
            seed + 777L,
            do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
          )
        },
        silent = TRUE
      )
      if (inherits(imp, "try-error")) {
        warning(
          "feature_est() failed; now using a simple random probabilistic setting ",
          "which is uniform feature sampling.",
          call. = FALSE
        )
        feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
        if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
        return(feat_prob)
      }
      return(.norm_feat_prob(imp, X_df))
    }

    ## 2) Ridge-based importance
    if (tolower(feature_est) == "ridge") {
      use_ridge <- ncol(X_df) > 0 &&
        any(abs(vsq_vals) > .Machine$double.eps) &&
        is.finite(sum(abs(vsq_vals)))
      if (use_ridge) {
        ridge_fit <- try(
          MASS::lm.ridge(vsq_vals ~ ., data = data.frame(X_df), lambda = 1),
          silent = TRUE
        )
        if (!inherits(ridge_fit, "try-error")) {
          beta <- as.numeric(ridge_fit$coef)
          names(beta) <- colnames(X_df)
          if (!all(beta == 0) && is.finite(sum(abs(beta)))) {
            feat_prob <- abs(beta) / sum(abs(beta))
            return(feat_prob)
          }
        }
      }
      warning(
        "Ridge feature importance failed or was degenerate; now using a simple random ",
        "probabilistic setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
      return(feat_prob)
    }

    ## 3) GBM-based importance
    df_gbm <- data.frame(X_df, vsq = vsq_vals)
    fit_gbm <- function() {
      gbm::gbm(
        vsq ~ ., data = df_gbm, distribution = "gaussian",
        n.trees = 100, interaction.depth = 3, shrinkage = 0.1,
        bag.fraction = 0.5, train.fraction = 1.0,
        n.minobsinnode = 2, verbose = FALSE
      )
    }
    gbm_fit <- try(
      if (is.null(seed)) fit_gbm() else withr::with_seed(seed + 1000L, fit_gbm()),
      silent = TRUE
    )
    if (inherits(gbm_fit, "try-error")) {
      warning(
        "GBM feature importance failed; now using a simple random probabilistic ",
        "setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
      return(feat_prob)
    }
    rel_inf <- gbm::relative.influence(gbm_fit, n.trees = 100, sort. = FALSE)
    names(rel_inf) <- colnames(X_df)
    if (sum(rel_inf) <= 0 || !is.finite(sum(rel_inf))) {
      warning(
        "GBM feature importance was degenerate; now using a simple random probabilistic ",
        "setting which is uniform feature sampling.",
        call. = FALSE
      )
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- rel_inf / sum(rel_inf)
    }
    feat_prob
  }

  # ---- Argument checks ----------------------------------------------------
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.logical(generalizability_path) || length(generalizability_path) != 1L)
    stop("`generalizability_path` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(leaf_proba) || leaf_proba < 0 || leaf_proba > 1)
    stop("`leaf_proba` must be between 0 and 1.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L))
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  if (!is.numeric(num_trees) || num_trees < 1)
    stop("`num_trees` must be positive.", call. = FALSE)
  num_trees <- as.integer(num_trees)
  if (!is.numeric(vote_threshold) || vote_threshold <= 0 || vote_threshold > 1)
    stop("`vote_threshold` must be in (0, 1].", call. = FALSE)
  if (!is.numeric(explore_proba) || explore_proba < 0 || explore_proba > 1)
    stop("`explore_proba` must be between 0 and 1.", call. = FALSE)
  if (!(is.character(feature_est) || is.function(feature_est)))
    stop("`feature_est` must be \"Ridge\", \"GBM\", or a function.", call. = FALSE)
  if (!is.logical(top_k_trees) || length(top_k_trees) != 1L)
    stop("`top_k_trees` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(k) || k < 1) stop("`k` must be a positive integer.", call. = FALSE)
  k <- as.integer(k)
  if (!(is.numeric(cutoff) || (is.character(cutoff) && cutoff == "baseline")))
    stop("`cutoff` must be \"baseline\" or numeric.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L)
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  if (is.null(global_objective_fn)) global_objective_fn <- objective_default
  if (!is.function(global_objective_fn))
    stop("`global_objective_fn` must be a function(D) -> numeric.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # ---- Step 1: generalizability_path vs general optimization --------------------
  if (generalizability_path) {
    outcome_col   <- "Y"
    treatment_col <- "Tr"
    sample_col    <- "S"

    needed <- c(outcome_col, treatment_col, sample_col)
    missing_cols <- setdiff(needed, names(data))
    if (length(missing_cols) > 0L) {
      stop(
        "With generalizability_path = TRUE, `data` must contain columns ",
        "'Y' (outcome), 'Tr' (treatment indicator), and 'S' (sample indicator). ",
        "Missing: ", paste(missing_cols, collapse = ", "), ".",
        call. = FALSE
      )
    }

    if (!is.numeric(data[[outcome_col]])) {
      stop("In generalizability_path mode, `Y` must be numeric.", call. = FALSE)
    }

    data[[treatment_col]] <- coerce01(data[[treatment_col]], allow_na = FALSE)
    data[[sample_col]]    <- coerce01(data[[sample_col]], allow_na = FALSE)

    scores <- compute_transport_scores(
      data     = data,
      outcome  = outcome_col,
      treatment = treatment_col,
      sample   = sample_col
    )
    v_all   <- scores$v
    vsq_all <- scores$vsq
    S_full  <- data[[sample_col]]

    analysis_idx <- which(S_full == 1L)
    if (length(analysis_idx) == 0L)
      stop("No observations with S == 1 (trial units) found in `data`.", call. = FALSE)

    testing_data <- data[analysis_idx, , drop = FALSE]
    v_vals       <- v_all[analysis_idx]
    vsq_vals     <- vsq_all[analysis_idx]
    S_vec        <- rep(1L, length(analysis_idx))
  } else {
    testing_data <- data
    n <- nrow(testing_data)

    if ("v" %in% names(testing_data)) {
      v_vals <- testing_data[["v"]]
    } else {
      v_vals <- rep(0, n)
    }

    if ("vsq" %in% names(testing_data)) {
      vsq_vals <- testing_data[["vsq"]]
    } else if ("v" %in% names(testing_data) && is.numeric(v_vals)) {
      mu_v <- mean(v_vals, na.rm = TRUE)
      vsq_vals <- (v_vals - mu_v)^2
    } else {
      vsq_vals <- rep(1, n)
    }

    if ("S" %in% names(testing_data)) {
      S_vec <- coerce01(testing_data[["S"]], allow_na = TRUE)
    } else {
      S_vec <- rep(1L, n)
    }
  }

  n <- nrow(testing_data)

  # ---- Step 2: Prepare covariates for splitting --------------------------
  drop_cols <- intersect(
    names(testing_data),
    c("Y", "Tr", "S", "v", "vsq", "w", "lX")
  )
  X_df <- testing_data[, setdiff(names(testing_data), drop_cols), drop = FALSE]

  if (ncol(X_df) == 0L) {
    # no covariates to split on; construct leaf-only trees
    features      <- "leaf"
    split_feature <- stats::setNames(1, "leaf")
  } else {
    features      <- c("leaf", colnames(X_df))
    split_feature <- compute_split_feature(X_df, vsq_vals)
    split_feature <- c(leaf = leaf_proba, split_feature)
    split_feature <- split_feature / sum(split_feature)
    names(split_feature) <- features
  }

  # ---- Step 3: Build forest ----------------------------------------------
  D_forest <- data.frame(
    X_df,
    v   = v_vals,
    vsq = vsq_vals,
    S   = S_vec,
    stringsAsFactors = FALSE
  )

  w_forest <- vector("list", num_trees)
  for (t_idx in seq_len(num_trees)) {
    D_tree <- data.frame(
      X_df,
      v   = v_vals,
      vsq = vsq_vals,
      w   = rep(1, n),
      S   = S_vec,
      stringsAsFactors = FALSE
    )
    rownames(D_tree) <- as.character(seq_len(n))
    rownames(X_df)   <- as.character(seq_len(n))

    grow_one <- function() {
      split_node(
        split_feature       = split_feature,
        X                   = D_tree,
        D                   = D_tree,
        parent_loss         = Inf,
        depth               = 0L,
        explore_proba       = explore_proba,
        choose_feature_fn   = choose_feature,
        reduce_weight_fn    = reduce_weight,
        global_objective_fn = global_objective_fn
      )
    }

    tree_res <- if (is.null(seed)) grow_one() else withr::with_seed(seed + t_idx, grow_one())
    D_updated <- tree_res$D
    D_forest[[paste0("w_tree_", t_idx)]] <- D_updated$w
    w_forest[[t_idx]] <- tree_res
  }

  # ---- Step 4: Rashomon set selection ------------------------------------
  obj_values   <- vapply(w_forest, function(res) res[["local objective"]], numeric(1))
  tree_indices <- seq_along(w_forest)

  if (top_k_trees) {
    if (k > num_trees) {
      warning("k > num_trees; using k = num_trees.", call. = FALSE)
      k <- num_trees
    }
    ord <- order(obj_values, na.last = NA)
    rashomon_set <- utils::head(ord, k)
  } else {
    cutoff_val <- if (identical(cutoff, "baseline")) {
      D_base <- D_forest
      D_base$w <- 1
      global_objective_fn(D_base)
    } else as.numeric(cutoff)
    if (!is.finite(cutoff_val)) cutoff_val <- Inf
    rashomon_set <- tree_indices[obj_values < cutoff_val]
  }

  if (length(rashomon_set) == 0L) {
    if (identical(cutoff, "baseline")) {
      message("No trees selected into Rashomon set.")
    } else {
      warning("No trees selected into Rashomon set.", call. = FALSE)
    }
  }
  not_in_set <- setdiff(tree_indices, rashomon_set)

  D_rash <- D_forest
  if (length(not_in_set) > 0L) {
    drop_cols2 <- paste0("w_tree_", not_in_set)
    keep_cols  <- setdiff(names(D_rash), drop_cols2)
    D_rash     <- D_rash[, keep_cols, drop = FALSE]
  }

  weight_cols <- grep("^w_tree_", names(D_rash), value = TRUE)
  D_weights   <- if (length(weight_cols) > 0L) D_rash[, weight_cols, drop = FALSE] else data.frame()
  if (ncol(D_weights) > 0L) {
    row_means        <- rowMeans(D_weights)
    D_rash$w_opt     <- as.integer(row_means > vote_threshold)
    D_rash$vote_count <- rowSums(D_weights)
  } else {
    D_rash$w_opt      <- integer(nrow(D_rash))
    D_rash$vote_count <- integer(nrow(D_rash))
  }

  # ---- Step 5: Summary tree ----------------------------------------------
  final_classifier <- NULL
  if (ncol(X_df) > 0L) {
    w_classes <- unique(stats::na.omit(D_rash$w_opt))
    if (length(w_classes) == 2L && all(w_classes %in% c(0L, 1L))) {
      final_classifier <- characterize_tree(X_df, as.factor(D_rash$w_opt))
    } else {
      message("No summary tree available to plot (w_opt not binary with two classes or no covariates).")
    }
  }

  # ---- Step 6: Estimands (only for generalizability_path = TRUE) ----------------
  est_list <- NULL
  if (generalizability_path) {
    in_S  <- rep(TRUE, n)  # we already subset to S == 1
    v_sel <- v_vals
    w_sel <- if ("w_opt" %in% names(D_rash)) D_rash$w_opt[in_S] else rep(1L, length(v_sel))

    ## Unweighted SATE
    est_label_unw <- "SATE"
    ok_unw        <- !is.na(v_sel)
    n_eff_unw     <- sum(ok_unw)
    mu_unw        <- if (n_eff_unw > 0) mean(v_sel[ok_unw]) else NA_real_
    se_unw        <- if (n_eff_unw > 1) {
      diffs2 <- (v_sel[ok_unw] - mu_unw)^2
      sqrt(sum(diffs2) / (n_eff_unw * (n_eff_unw - 1)))
    } else NA_real_

    ## Weighted WTATE
    est_label_w  <- "WTATE"
    is_binary_w  <- all(is.na(w_sel) | (w_sel %in% c(0, 1)))
    den_w_any    <- sum(w_sel, na.rm = TRUE)
    mu_w_any     <- if (den_w_any > 0) sum(w_sel * v_sel, na.rm = TRUE) / den_w_any else NA_real_

    if (is_binary_w) {
      ok_w  <- ok_unw & !is.na(w_sel) & (w_sel == 1L)
      n_A   <- sum(ok_w)
      den_w <- sum(w_sel == 1L, na.rm = TRUE)
      if (n_A > 0) {
        mu_w <- mean(v_sel[ok_w])
        se_w <- if (n_A > 1) {
          diffs2_A <- (v_sel[ok_w] - mu_w)^2
          sqrt(sum(diffs2_A) / (n_A * (n_A - 1)))
        } else NA_real_

        se_w_note <- paste0(
          "Calculation of SE for WTATE uses sqrt( sum_{A} (v_i - vbar_A)^2 / ( n_A * (n_A - 1) ) ).\n",
          "  Here A = { i : w_i = 1 }, n_A = |A|, v_i are unit-level orthogonal scores, and vbar_A is their mean.\n"
        )
        if (!identical(global_objective_fn, objective_default)) {
          se_w_note <- paste0(
            se_w_note,
            "\nYou supplied a custom global_objective_fn; please verify this SE matches your ",
            "estimand, or perform additional variance analysis (e.g., bootstrap or ",
            "influence-function methods)."
          )
        }
      } else {
        mu_w     <- NA_real_
        se_w     <- NA_real_
        den_w    <- 0L
        n_A      <- 0L
        se_w_note <- "SE omitted: no kept observations (A = { i : w_i = 1 } is empty)."
      }
    } else {
      mu_w     <- mu_w_any
      den_w    <- den_w_any
      n_A      <- NA_integer_
      se_w     <- NA_real_
      se_w_note <- paste0(
        "SE omitted: non-binary weights detected; the subset-mean SE above is not appropriate.",
        " If you intentionally use a custom global_objective_fn with non-binary weights, ",
        "please conduct further variance analysis tailored to that objective."
      )
      warning("Non-binary w_opt detected; omitting weighted SE.", call. = FALSE)
    }

    if (verbose) {
      message(sprintf("%s (unweighted) = %.6f, SE = %.6f", est_label_unw, mu_unw, se_unw))
      if (is.finite(se_w) || (!is.na(se_w) && !is.nan(se_w))) {
        message(sprintf("%s (weighted)   = %.6f, SE = %.6f", est_label_w, mu_w, se_w))
        message(paste0("  Note: ", se_w_note))
      } else {
        message(sprintf("%s (weighted)   = %.6f", est_label_w, mu_w))
        message(paste0("  Note: ", se_w_note))
      }
    }

    est_list <- list(
      estimand_unweighted = est_label_unw,
      value_unweighted    = mu_unw,
      se_unweighted       = se_unw,
      estimand_weighted   = est_label_w,
      value_weighted      = mu_w,
      se_weighted         = se_w,
      se_weighted_note    = se_w_note,
      n_analysis          = sum(in_S),
      sum_w               = den_w,
      n_A                 = if (is_binary_w) n_A else NA_integer_
    )
  }

  res <- list(
    D_rash             = D_rash,
    D_forest           = D_forest,
    w_forest           = w_forest,
    rashomon_set       = rashomon_set,
    global_objective_fn = global_objective_fn,
    f                  = final_classifier,
    testing_data       = testing_data,
    estimate           = est_list,
    generalizability_path     = generalizability_path
  )
  class(res) <- c("ROOT", "list")
  res
}
