#' Characterize under-represented subgroups (wraps ROOT)
#'
#' Combines an RCT (S=1) and a target dataset (S=0), then calls \code{ROOT()} to
#' learn a weighted tree that identifies subgroups with different representation
#' in the target population compared to the trial.
#'
#' @param DataRCT A data frame containing the randomized clinical trial data.
#'   Must include treatment, outcome, and covariate columns.
#' @param covariateColName_RCT A character vector of covariate column names in \code{DataRCT}.
#' @param trtColName_RCT A character string specifying the treatment column name in \code{DataRCT} (0/1).
#' @param outcomeColName_RCT A character string specifying the outcome column name in \code{DataRCT}.
#' @param DataTarget A data frame containing the target population data (covariates only).
#' @param covariateColName_TargetData A character vector of covariate column names in \code{DataTarget}.
#' @param leaf_proba A numeric value for the "leaf" probability in the ROOT tree growth (default 0.25).
#' @param seed An integer seed for reproducibility (default 123).
#' @param num_trees An integer specifying the number of trees to grow (default 10).
#' @param vote_threshold A numeric value (0.5, 1] for the majority vote threshold (default 2/3).
#' @param explore_proba A numeric value for exploration probability (default 0.05).
#' @param feature_est A string ("Ridge", "GBM") or function for feature importance estimation.
#' @param feature_est_args A named list of arguments passed to the feature estimator.
#' @param top_k_trees Logical; if TRUE, selects the top k trees instead of using a cutoff.
#' @param k Integer; number of trees to select if \code{top_k_trees} is TRUE.
#' @param cutoff A numeric value or "baseline" to determine the Rashomon set cutoff.
#' @param verbose Logical; if TRUE, prints progress and estimand summaries.
#' @param global_objective_fn A function(D) -> numeric to minimize (default \code{objective_default}).
#' @param keep_threshold Unused; kept for backward compatibility.
#' @param lX_threshold Unused; kept for backward compatibility.
#'
#' @return A \code{characterizing_underrep} object (S3 list) containing:
#'   \item{root}{The resulting \code{ROOT} object.}
#'   \item{combined}{The combined data frame (RCT + Target) used for analysis.}
#'   \item{leaf_summary}{A data frame summarizing the terminal nodes (rules, counts, and labels).}
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(diabetes_data)
#'
#' # Split into Trial (S=1) and Target (S=0)
#' trial  <- subset(diabetes_data, S == 1)
#' target <- subset(diabetes_data, S == 0)
#'
#' # Run characterization
#' res <- characterizing_underrep(
#'   DataRCT = trial,
#'   covariateColName_RCT = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
#'   trtColName_RCT = "Tr",
#'   outcomeColName_RCT = "Y",
#'   DataTarget = target,
#'   covariateColName_TargetData = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
#'   seed = 123
#' )
#'
#' # View Summary
#' summary(res)
#'
#' # Plot the annotated tree
#' plot(res)
#' }
#' @export
characterizing_underrep <- function(
    DataRCT,
    covariateColName_RCT,
    trtColName_RCT,
    outcomeColName_RCT,
    DataTarget,
    covariateColName_TargetData,
    leaf_proba = 0.25,
    seed = 123,
    num_trees = 10,
    vote_threshold = 2 / 3,
    explore_proba = 0.05,
    feature_est = "Ridge",
    feature_est_args = list(),
    top_k_trees = FALSE,
    k = 10,
    cutoff = "baseline",
    verbose = FALSE,
    global_objective_fn = objective_default,
    keep_threshold = 0.50,
    lX_threshold = NULL
) {
  # 1. Input Validation
  if (!is.data.frame(DataRCT) || !is.data.frame(DataTarget)) {
    stop("DataRCT and DataTarget must be data.frames.", call. = FALSE)
  }
  if (!all(c(trtColName_RCT, outcomeColName_RCT) %in% names(DataRCT))) {
    stop("trtColName_RCT and/or outcomeColName_RCT not found in DataRCT.", call. = FALSE)
  }
  if (!all(covariateColName_RCT %in% names(DataRCT))) {
    miss <- setdiff(covariateColName_RCT, names(DataRCT))
    stop("Missing RCT covariates in DataRCT: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!all(covariateColName_TargetData %in% names(DataTarget))) {
    miss <- setdiff(covariateColName_TargetData, names(DataTarget))
    stop("Missing target covariates in DataTarget: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # 2. Harmonize Covariates
  # We take the intersection of covariates available in both datasets
  common_covs <- intersect(covariateColName_RCT, covariateColName_TargetData)
  if (length(common_covs) == 0L) {
    stop("No overlapping covariate names between RCT and Target.", call. = FALSE)
  }

  # Warn if some requested covariates are being dropped
  if (!identical(sort(common_covs), sort(covariateColName_RCT)) ||
      !identical(sort(common_covs), sort(covariateColName_TargetData))) {
    warning("Using intersection of covariates: ", paste(common_covs, collapse = ", "))
  }

  X_rct <- DataRCT[, common_covs, drop = FALSE]
  X_tgt <- DataTarget[, common_covs, drop = FALSE]

  # 3. Construct Combined Dataset
  # S = 1 for RCT, S = 0 for Target
  combined_rct <- data.frame(
    X_rct,
    Y  = DataRCT[[outcomeColName_RCT]],
    Tr = DataRCT[[trtColName_RCT]],
    S  = 1L
  )
  combined_tgt <- data.frame(
    X_tgt,
    Y  = NA_real_,
    Tr = NA_integer_,
    S  = 0L
  )
  combined <- rbind(combined_rct, combined_tgt)
  rownames(combined) <- NULL

  # 4. Call ROOT (Logic Only)
  # We strip away all plotting arguments here.
  root_out <- ROOT(
    data = combined,
    outcome = "Y",
    treatment = "Tr",
    sample = "S",
    leaf_proba = leaf_proba,
    seed = seed,
    num_trees = num_trees,
    vote_threshold = vote_threshold,
    explore_proba = explore_proba,
    feature_est = feature_est,
    feature_est_args = feature_est_args,
    top_k_trees = top_k_trees,
    k = k,
    cutoff = cutoff,
    verbose = verbose,
    global_objective_fn = global_objective_fn
  )

  # 5. Summarize Terminal Nodes
  # Extract rules from the summary tree (if it exists)
  leaf_summary <- NULL
  if (!is.null(root_out$f)) {
    f <- root_out$f
    frm <- f$frame
    is_leaf <- frm$var == "<leaf>"

    # Get paths (rules)
    # path.rpart returns a list of paths for the specified nodes
    paths <- rpart::path.rpart(f, nodes = as.numeric(rownames(frm)[is_leaf]), print.it = FALSE)
    path_text <- vapply(paths, function(p) paste(p, collapse = " & "), character(1))

    # Predicted class label per leaf
    ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)
    pred_class <- ylv[frm$yval[is_leaf]]

    n_in_leaf <- frm$n[is_leaf]
    total_n <- sum(n_in_leaf)
    pct <- if (total_n > 0) round(100 * n_in_leaf / total_n, 1) else NA_real_

    # Human-readable label
    label <- ifelse(pred_class %in% c("0", "0.0", "0L", "FALSE"),
                    "Under-represented (drop, w=0)",
                    "Represented (keep, w=1)")

    leaf_summary <- data.frame(
      leaf_id     = rownames(frm)[is_leaf],
      rule        = path_text,
      predicted_w = pred_class,
      n           = n_in_leaf,
      pct         = pct,
      label       = label,
      row.names   = NULL,
      check.names = FALSE
    )
  }

  # 6. Construct and Return S3 Object
  res <- list(
    root         = root_out,
    combined     = combined,
    leaf_summary = leaf_summary
  )
  class(res) <- c("characterizing_underrep", "list")
  return(res)
}
