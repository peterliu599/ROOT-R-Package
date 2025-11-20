#' Ensemble of weighted trees (loss/objective-agnostic) and Rashomon selection
#'
#' Builds multiple weighted trees, then identifies a "Rashomon set" of
#' top-performing trees and aggregates their weight assignments by majority vote.
#'
#' @param data A data frame containing the dataset. Must include outcome, treatment, and sample indicator columns.
#' @param outcome A character string specifying the name of the outcome column in \code{data}.
#' @param treatment A character string specifying the name of the treatment indicator column (0/1) in \code{data}.
#' @param sample A character string specifying the name of the sample indicator column (0/1) in \code{data}. Use \code{NULL} for single-sample SATE mode.
#' @param leaf_proba A numeric value specifying the probability mass for the "leaf" option in each tree (default \code{0.25}).
#' @param seed An integer seed for reproducibility (default \code{NULL}).
#' @param num_trees An integer specifying the number of trees to grow in the forest (default \code{10}).
#' @param vote_threshold A numeric value in (0.5, 1] specifying the majority vote threshold for final \code{weight=1} (default \code{2/3}).
#' @param explore_proba A numeric value specifying the probability of exploration at leaves in each tree (default \code{0.05}).
#' @param feature_est A character string (\code{"Ridge"}, \code{"GBM"}) or a function \code{function(X, y, ...)} returning a named, non-negative vector of importances.
#' @param feature_est_args A named list of extra arguments for a user-supplied \code{feature_est} function.
#' @param top_k_trees A logical value. If \code{TRUE}, select top-k trees by objective; else use cutoff (default \code{FALSE}).
#' @param k An integer specifying the number of top trees if \code{top_k_trees=TRUE} (default \code{10}).
#' @param cutoff A numeric value or character string "baseline". If \code{top_k_trees=FALSE}, this defines the Rashomon set cutoff.
#' @param verbose A logical value. If \code{TRUE}, prints 2 lines with (unweighted and weighted) estimate + SE. Default \code{FALSE}.
#' @param global_objective_fn A function scoring the entire state (minimize).
#'
#' @return An S3 object of class \code{"ROOT"}. This object is a list with the following elements:
#'   \item{D_rash}{The data frame with weights from the Rashomon set.}
#'   \item{f}{The summary \code{rpart} tree object.}
#'   \item{estimate}{A list containing the unweighted and weighted estimates.}
#'   \item{...}{Other elements containing the internal forest structures.}
#'
#' @seealso
#'   \code{\link{summary.ROOT}} for summarizing results,
#'   \code{\link{plot.ROOT}} for visualizing the decision tree.
#'
#' @examples
#' \dontrun{
#' data(diabetes_data)
#'
#' # Run ROOT
#' res <- ROOT(data = diabetes_data, outcome = "Y", treatment = "Tr", sample = "S")
#'
#' # Summary of results
#' summary(res)
#'
#' # Plot the characterization tree
#' plot(res)
#' }
#' @export
ROOT <- function(data,
                 outcome,
                 treatment,
                 sample,
                 leaf_proba = 0.25,
                 seed = NULL,
                 num_trees = 10,
                 vote_threshold = 2 / 3,
                 explore_proba = 0.05,
                 feature_est = "Ridge",
                 feature_est_args = list(),
                 top_k_trees = FALSE,
                 k = 10,
                 cutoff = "baseline",
                 verbose = FALSE,
                 global_objective_fn = objective_default) {
  # Helpers
  coerce01 <- function(x, allow_na = TRUE) {
    if (is.numeric(x) || is.logical(x)) return(as.integer(x))
    if (is.factor(x)) x <- as.character(x)
    x_chr <- trimws(tolower(as.character(x)))
    map1 <- c("1","yes","treated","t","true")
    map0 <- c("0","no","control","c","false")
    out <- rep(NA_integer_, length(x_chr))
    out[x_chr %in% map1] <- 1L; out[x_chr %in% map0] <- 0L
    if (!allow_na && any(is.na(out))) stop("Non 0/1 values found.")
    out
  }
  .norm_feat_prob <- function(imp, X_df) {
    if (!is.numeric(imp) || is.null(names(imp))) {
      stop("Custom feature importance must be a *named* numeric vector.", call. = FALSE)
    }
    imp <- imp[colnames(X_df)]
    if (any(is.na(imp))) stop("Importance missing for some X_df columns.", call. = FALSE)
    if (any(imp < 0)) stop("Importances must be non-negative.", call. = FALSE)
    s <- sum(imp)
    if (!is.finite(s) || s <= 0) {
      feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
      if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
    } else {
      feat_prob <- imp / s
    }
    feat_prob
  }

  if (!is.function(global_objective_fn)) stop("`global_objective_fn` must be a function(D)->numeric.")

  covariate_cols <- setdiff(names(data), c(outcome, treatment, sample))
  if (length(covariate_cols) == 0L) {
    stop("ROOT(): no covariate columns found (need at least one feature besides outcome/treatment/sample).",
         call. = FALSE)
  }

  # Basic tests
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(outcome)   || length(outcome)   != 1L || !(outcome   %in% names(data)))
    stop("`outcome` must be a single column name present in `data`.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L || !(treatment %in% names(data)))
    stop("`treatment` must be a single column name present in `data`.", call. = FALSE)
  if (!is.null(sample)) {
    if (!is.character(sample) || length(sample) != 1L)
      stop("`sample` must be NULL or a single column name string.", call. = FALSE)
    if (!(sample %in% names(data)))
      stop("`sample` column not found; pass `sample = NULL` to run single-sample mode.", call. = FALSE)
  }
  if (!is.numeric(leaf_proba) || leaf_proba < 0 || leaf_proba > 1)
    stop("`leaf_proba` must be between 0 and 1.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L))
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  if (!is.numeric(num_trees) || num_trees < 1) stop("`num_trees` must be positive.", call. = FALSE)
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

  # Coerce treatment/sample to 0/1
  data[[treatment]] <- coerce01(data[[treatment]], allow_na = TRUE)
  if (!is.null(sample)) data[[sample]] <- coerce01(data[[sample]], allow_na = FALSE)
  Tr_check <- data[[treatment]]
  if (!all(Tr_check[!is.na(Tr_check)] %in% c(0L,1L)))
    stop("`treatment` column must be 0/1 (after coercion).", call. = FALSE)
  if (!is.null(sample)) {
    S_check <- data[[sample]]
    if (!all(S_check %in% c(0L,1L)))
      stop("`sample` column must be 0/1.", call. = FALSE)
  }

  # Single vs two-sample mode
  single_sample_mode <- is.null(sample) ||
    (all(data[[sample]] %in% c(0,1)) && (all(data[[sample]] == 1) || all(data[[sample]] == 0)))
  if (!single_sample_mode && (all(data[[sample]] == 0) || all(data[[sample]] == 1)))
    stop("Sample indicator `", sample, "` has no variation (all 0 or all 1).", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # Estimate orthogonal scores (v)
  if (single_sample_mode) {
    dml <- estimate_dml_single(data, outcome = outcome, treatment = treatment, crossfit = 5)
    df_v <- dml$df_v
    testing_data <- dml$data2
    S_vec <- rep(1L, nrow(testing_data))      # all target in single-sample path
    lX_col <- NA_real_                        # no overlap factor
  } else {
    dml <- estimate_dml(data, outcome = outcome, treatment = treatment, sample = sample, crossfit = 5)
    df_v <- dml$df_v
    testing_data <- dml$data2
    S_vec <- testing_data[[sample]]
    lX_col <- 1 / df_v$b                      # overlap factor = 1/b
  }

  Y_vec <- testing_data[[outcome]]
  Tr_vec <- testing_data[[treatment]]
  drop_cols <- c(outcome, treatment, if (!is.null(sample) && (sample %in% names(testing_data))) sample)
  X_df <- testing_data[, setdiff(names(testing_data), drop_cols), drop = FALSE]

  n <- nrow(testing_data)
  v_vals <- df_v$te
  vsq_vals <- df_v$te_sq

  # Feature probabilities for splitting
  no_cov <- (ncol(X_df) == 0)
  features <- c("leaf", colnames(X_df))

  compute_split_feature <- function() {
    if (ncol(X_df) < 1L) {
      stop("compute_split_feature(): no features available for weighting.", call. = FALSE)
    }
    if (is.function(feature_est)) {
      imp <- if (is.null(seed)) {
        do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
      } else {
        withr::with_seed(seed + 777L,
                         do.call(feature_est, c(list(X = X_df, y = vsq_vals), feature_est_args))
        )
      }
      feat_prob <- .norm_feat_prob(imp, X_df)
    } else if (tolower(feature_est) == "ridge") {
      if (ncol(X_df) == 0 || all(abs(vsq_vals) < .Machine$double.eps) || !is.finite(sum(abs(vsq_vals)))) {
        feat_prob <- rep(1 / max(1, ncol(X_df)), ncol(X_df))
        if (ncol(X_df) > 0) names(feat_prob) <- colnames(X_df)
      } else {
        ridge_fit <- MASS::lm.ridge(vsq_vals ~ ., data = data.frame(X_df), lambda = 1)
        beta <- as.numeric(ridge_fit$coef); names(beta) <- colnames(X_df)
        if (all(beta == 0) || !is.finite(sum(abs(beta)))) {
          feat_prob <- rep(1 / ncol(X_df), ncol(X_df)); names(feat_prob) <- colnames(X_df)
        } else {
          feat_prob <- abs(beta) / sum(abs(beta))
        }
      }
    } else { # GBM
      df_gbm <- data.frame(X_df, vsq = vsq_vals)
      fit_gbm <- function() {
        gbm::gbm(vsq ~ ., data = df_gbm, distribution = "gaussian",
                 n.trees = 100, interaction.depth = 3, shrinkage = 0.1,
                 bag.fraction = 0.5, train.fraction = 1.0, n.minobsinnode = 2, verbose = FALSE)
      }
      gbm_fit <- if (is.null(seed)) fit_gbm() else withr::with_seed(seed + 1000L, fit_gbm())
      rel_inf <- gbm::relative.influence(gbm_fit, n.trees = 100, sort. = FALSE)
      names(rel_inf) <- colnames(X_df)
      if (sum(rel_inf) <= 0 || !is.finite(sum(rel_inf))) {
        feat_prob <- rep(1 / ncol(X_df), ncol(X_df)); names(feat_prob) <- colnames(X_df)
      } else {
        feat_prob <- rel_inf / sum(rel_inf)
      }
    }
    proba <- c(leaf_proba, feat_prob); proba <- proba / sum(proba)
    stats::setNames(proba, features)
  }

  split_feature <- if (no_cov) stats::setNames(1, "leaf") else compute_split_feature()

  # Build forest
  D_forest <- data.frame(X_df, v = v_vals, vsq = vsq_vals, S = S_vec, stringsAsFactors = FALSE)
  D_forest$lX <- lX_col

  w_forest <- vector("list", num_trees)
  for (t_idx in seq_len(num_trees)) {
    D_tree <- data.frame(X_df, v = v_vals, vsq = vsq_vals, w = rep(1, n), S = S_vec, stringsAsFactors = FALSE)
    rownames(D_tree) <- as.character(seq_len(n)); rownames(X_df) <- as.character(seq_len(n))
    grow_one <- function() {
      split_node(
        split_feature = split_feature,
        X = D_tree, D = D_tree, parent_loss = Inf, depth = 0L,
        explore_proba = explore_proba,
        choose_feature_fn = choose_feature,
        reduce_weight_fn = reduce_weight,
        global_objective_fn = global_objective_fn
      )
    }
    tree_res <- if (is.null(seed)) grow_one() else withr::with_seed(seed + t_idx, grow_one())
    D_updated <- tree_res$D
    D_forest[[paste0("w_tree_", t_idx)]] <- D_updated$w
    w_forest[[t_idx]] <- tree_res
  }

  # Rashomon set selection
  obj_values <- vapply(w_forest, function(res) res[["local objective"]], numeric(1))
  tree_indices <- seq_along(w_forest)
  if (top_k_trees) {
    if (k > num_trees) { warning("k > num_trees; using k = num_trees.", call. = FALSE); k <- num_trees }
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
  if (length(rashomon_set) == 0) warning("No trees selected into Rashomon set.", call. = FALSE)
  not_in_set <- setdiff(tree_indices, rashomon_set)

  # Keep selected trees and votes
  D_rash <- D_forest
  if (length(not_in_set) > 0) {
    drop_cols <- paste0("w_tree_", not_in_set)
    keep_cols <- setdiff(names(D_rash), drop_cols)
    D_rash <- D_rash[, keep_cols, drop = FALSE]
  }
  weight_cols <- grep("^w_tree_", names(D_rash), value = TRUE)
  D_weights <- if (length(weight_cols) > 0) D_rash[, weight_cols, drop = FALSE] else data.frame()
  if (ncol(D_weights) > 0) {
    row_means <- rowMeans(D_weights)
    D_rash$w_opt <- as.integer(row_means > vote_threshold)
    D_rash$vote_count <- rowSums(D_weights)
  } else {
    D_rash$w_opt <- integer(nrow(D_rash))
    D_rash$vote_count <- integer(nrow(D_rash))
  }

  # Final summarized tree
  final_classifier <- if (no_cov) NULL else characterize_tree(X_df, as.factor(D_rash$w_opt))

  # --- CHANGE: Removed plotting logic here entirely ---

  # Estimands: unweighted vs weighted (with SEs only)
  in_S  <- if (single_sample_mode) rep(TRUE, n) else (S_vec == 1L)
  v_sel <- v_vals[in_S]
  w_sel <- if ("w_opt" %in% names(D_rash)) D_rash$w_opt[in_S] else rep(1L, length(v_sel))

  ## Unweighted (ATE in RCT or TATE)
  est_label_unw <- if (single_sample_mode) "SATE" else "TATE"
  ok_unw   <- !is.na(v_sel)
  n_eff_unw <- sum(ok_unw)
  mu_unw   <- if (n_eff_unw > 0) mean(v_sel[ok_unw]) else NA_real_
  # SE(tbar) = sqrt( 1 / ( n * (n - 1) ) * sum_i (t_i - tbar)^2 )
  se_unw <- if (n_eff_unw > 1) {
    diffs2 <- (v_sel[ok_unw] - mu_unw)^2
    sqrt( sum(diffs2) / (n_eff_unw * (n_eff_unw - 1)) )
  } else NA_real_

  ## Weighted (WATE in RCT or WTATE)
  est_label_w <- if (single_sample_mode) "WATE" else "WTATE"

  # Is the final weight vector binary? (subset-mean WTATE is appropriate only then)
  is_binary_w <- all(is.na(w_sel) | (w_sel %in% c(0, 1)))
  # For completeness, still compute weighted mean even if non-binary
  den_w_any <- sum(w_sel, na.rm = TRUE)
  mu_w_any  <- if (den_w_any > 0) sum(w_sel * v_sel, na.rm = TRUE) / den_w_any else NA_real_

  if (is_binary_w) {
    # Subset mean over A = {w_opt == 1}
    ok_w <- ok_unw & !is.na(w_sel) & (w_sel == 1L)
    n_A  <- sum(ok_w)
    den_w <- sum(w_sel == 1L, na.rm = TRUE)  # equals n_A for binary w
    if (n_A > 0) {
      mu_w <- mean(v_sel[ok_w])
      # SE(vbar_A) = sqrt( 1 / { n_A * (n_A - 1) } * sum_{i in A} (v_i - vbar_A)^2 )
      se_w <- if (n_A > 1) {
        diffs2_A <- (v_sel[ok_w] - mu_w)^2
        sqrt( sum(diffs2_A) / (n_A * (n_A - 1)) )
      } else NA_real_

      # Compose the explanatory note that will always print under WTATE
      se_w_note <- paste0(
        "Calculation of SE for WTATE uses sqrt( sum_{A} (v_i - vbar_A)^2 / ( n_A * (n_A - 1) ) ).\n",
        "  Here A = { i : w_i = 1 }, n_A = |A|, v_i are unit-level orthogonal scores, and vbar_A is their mean. \n"
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
      mu_w <- NA_real_; se_w <- NA_real_
      den_w <- 0L
      se_w_note <- "SE omitted: no kept observations (A = { i : w_i = 1 } is empty)."
    }
  } else {
    # Non-binary weights: subset-mean SE is not appropriate.
    mu_w <- mu_w_any
    den_w <- den_w_any
    n_A   <- NA_integer_
    se_w  <- NA_real_
    se_w_note <- paste0(
      "SE omitted: non-binary weights detected; the subset-mean SE above is not appropriate.",
      " If you intentionally use a custom global_objective_fn with non-binary weights, ",
      " please conduct further variance analysis tailored to that objective."
    )
    warning("Non-binary w_opt detected; omitting weighted SE.", call. = FALSE)
  }

  # Verbose printing
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


  # Assemble result (Remove tree_plot)
  results <- list(
    D_rash = D_rash,
    D_forest = D_forest,
    w_forest = w_forest,
    rashomon_set = rashomon_set,
    global_objective_fn = global_objective_fn,
    f = final_classifier,
    testing_data = testing_data,
    # tree_plot = tree_plot,  <-- REMOVED

    estimate = list(
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
  )
  results$single_sample_mode <- single_sample_mode
  class(results) <- c("ROOT", "list")
  return(results)
}
