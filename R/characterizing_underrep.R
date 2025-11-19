#' Characterize under-represented subgroups (wraps ROOT)
#'
#' Combines an RCT (S=1) and a target dataset (S=0), calls \code{ROOT()} to learn
#' a binary selector \code{w(x)}.
#'
#' @param DataRCT A data frame with trial data.
#' @param covariateColName_RCT A character vector of covariate column names in \code{DataRCT}.
#' @param trtColName_RCT A character string: treatment column name in \code{DataRCT} (0/1).
#' @param outcomeColName_RCT A character string: outcome column name in \code{DataRCT}.
#' @param DataTarget A data frame with target-population covariates.
#' @param covariateColName_TargetData A character vector of covariate column names in \code{DataTarget}.
#' @param leaf_proba,seed,num_trees,vote_threshold,explore_proba,feature_est,feature_est_args,top_k_trees,k,cutoff,verbose Passed to \code{ROOT()}.
#' @param global_objective_fn Global objective/loss function for ROOT.
#' @param keep_threshold,lX_threshold Kept for API compatibility (unused).
#'
#' @return A \code{characterizing_underrep} object.
#' @examples
#' \dontrun{
#' data(diabetes_data)
#' trial  <- subset(diabetes_data, S == 1)
#' target <- subset(diabetes_data, S == 0)
#'
#' res <- characterizing_underrep(
#'   DataRCT = trial,
#'   covariateColName_RCT = c("Race_Black", "Sex_Male", "DietYes", "Age45"),
#'   trtColName_RCT = "Tr",
#'   outcomeColName_RCT = "Y",
#'   DataTarget = target,
#'   covariateColName_TargetData = c("Race_Black", "Sex_Male", "DietYes", "Age45")
#' )
#' summary(res)
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
  # Basic tests
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

  # Harmonize covariates
  common_covs <- intersect(covariateColName_RCT, covariateColName_TargetData)
  if (length(common_covs) == 0L) stop("No overlapping covariate names between RCT and Target.", call. = FALSE)
  if (!identical(sort(common_covs), sort(covariateColName_RCT)) ||
      !identical(sort(common_covs), sort(covariateColName_TargetData))) {
    warning("Using intersection of covariates: ", paste(common_covs, collapse = ", "))
  }

  X_rct <- DataRCT[, common_covs, drop = FALSE]
  X_tgt <- DataTarget[, common_covs, drop = FALSE]

  # 1 = RCT (has Y, Tr), 0 = Target (covariates only)
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

  # Call ROOT on the combined frame (Remove plotting args)
  root_args <- list(
    data = combined,
    outcome = "Y",
    treatment = "Tr",
    sample = "S"
  )
  dots <- list(
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

  root_out <- do.call(ROOT, c(root_args, dots))

  # Summarize terminal nodes (from the fitted subset)
  leaf_summary <- NULL
  if (!is.null(root_out$f)) {
    f <- root_out$f
    frm <- f$frame
    is_leaf <- frm$var == "<leaf>"

    # Get paths (rules)
    paths <- rpart::path.rpart(f, nodes = as.numeric(rownames(frm)[is_leaf]), print.it = FALSE)
    path_text <- vapply(paths, function(p) paste(p, collapse = " & "), character(1))

    # Predicted class label per leaf (prefer f$ylevels if present)
    ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)
    pred_class <- ylv[frm$yval[is_leaf]]

    n_in_leaf <- frm$n[is_leaf]
    total_n <- sum(n_in_leaf)
    pct <- if (total_n > 0) round(100 * n_in_leaf / total_n, 1) else NA_real_
    label <- ifelse(pred_class %in% c("0", "0.0", "0L"),
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

  # Return (Remove plot objects)
  res <- list(
    root               = root_out,
    combined           = combined,
    leaf_summary       = leaf_summary
  )
  class(res) <- c("characterizing_underrep", "list")
  return(res)
}
