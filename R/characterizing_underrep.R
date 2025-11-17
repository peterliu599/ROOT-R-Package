#' Characterize under-represented subgroups (wraps ROOT)
#'
#' Combines an RCT (S=1) and a target dataset (S=0), calls \code{ROOT()} to learn
#' a binary selector \code{w(x)}, and (optionally) renders an annotated tree that
#' highlights represented (w=1) vs underrepresented (w=0) leaves.
#'
#' @param DataRCT data.frame with trial data; must include \code{trtColName_RCT}, \code{outcomeColName_RCT}, and the covariates named in \code{covariateColName_RCT}.
#' @param covariateColName_RCT character vector of covariate column names in \code{DataRCT}.
#' @param trtColName_RCT single string: treatment column name in \code{DataRCT} (0/1).
#' @param outcomeColName_RCT single string: outcome column name in \code{DataRCT}.
#' @param DataTarget data.frame with target-population covariates (no treatment/outcome required).
#' @param covariateColName_TargetData character vector of covariate column names in \code{DataTarget}.
#' @param leaf_proba,seed,num_trees,vote_threshold,explore_proba,feature_est,feature_est_args,top_k_trees,k,cutoff,verbose Passed to \code{ROOT()}.
#' @param global_objective_fn Global objective/loss function for ROOT (default \code{objective_default}).
#' @param root_plot_tree Logical; pass-through to \code{ROOT(plot_tree=...)}.
#' @param root_plot_args Optional list passed to \code{ROOT(plot_tree_args=...)}.
#' @param plot_underrep Logical; if \code{TRUE}, draws an annotated represented/underrepresented tree.
#' @param keep_threshold,lX_threshold Kept for API compatibility (unused).
#' @param plot_main Title for the annotated plot.
#'
#' @return A \code{characterizing_underrep} object with
#'   \item{root}{the fitted \code{ROOT} object}
#'   \item{combined}{stacked RCT+Target data used for fitting}
#'   \item{leaf_summary}{data.frame with terminal rules and sizes (if a summary tree exists)}
#'   \item{tree_plot_root}{recorded plot of the ROOT summary tree (if produced)}
#'   \item{tree_plot_underrep}{recorded plot of the annotated underrep tree (if produced)}
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
    root_plot_tree = TRUE,
    root_plot_args = list(
      type = 2, extra = 109, under = TRUE, faclen = 0, tweak = 1.1,
      fallen.leaves = TRUE, box.palette = c("pink", "palegreen3"),
      shadow.col = c("gray"), branch.lty = 3,
      main = "Final Characterized Tree from Rashomon Set"
    ),
    plot_underrep = TRUE,
    keep_threshold = 0.50,
    lX_threshold = NULL,
    plot_main = "Underrepresented Population Characterization Tree"
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

  # Call ROOT on the combined frame
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
    plot_tree = root_plot_tree,
    global_objective_fn = global_objective_fn,   # <-- pass user-supplied objective
    plot_tree_args = root_plot_args
  )
  if (!is.null(root_plot_args) && length(root_plot_args)) {
    dots$plot_tree_args <- root_plot_args
  }
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

  # Annotated underrepresented plotter
  .plot_underrep_tree <- function(fo_res, main) {
    if (is.null(fo_res$f)) return(NULL)
    if (!requireNamespace("rpart.plot", quietly = TRUE)) {
      warning("Package 'rpart.plot' is not installed; skipping annotated tree plot.", call. = FALSE)
      return(NULL)
    }
    f <- fo_res$f
    frm <- f$frame
    isLeaf <- frm$var == "<leaf>"
    node_ids <- as.integer(row.names(frm))
    ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)

    node_pred_label <- rep(NA_character_, nrow(frm))
    node_pred_label[isLeaf] <- if (length(ylv) >= 2) ylv[frm$yval[isLeaf]] else as.character(frm$yval[isLeaf])

    is_keep_label <- function(lbl) lbl %in% c("1", "1.0", "1L")
    col_keep <- "#4E79A7"; col_drop <- "#F28E2B"

    box_col <- rep(NA, nrow(frm))
    if (any(isLeaf)) {
      lab <- node_pred_label[isLeaf]
      box_col[isLeaf] <- ifelse(is_keep_label(lab), col_keep, col_drop)
    }

    leaf_ids <- node_ids[isLeaf]
    leaf_count <- frm$n[isLeaf]
    total_cnt <- sum(leaf_count)
    leaf_stats <- data.frame(
      leaf = leaf_ids,
      n    = leaf_count,
      pct  = if (total_cnt > 0) leaf_count / total_cnt else 0
    )
    stat_map <- split(leaf_stats, as.character(leaf_stats$leaf))

    node_fun <- function(x, labs, digits, varlen) {
      fr <- x$frame
      ids <- as.integer(row.names(fr))
      out <- character(length(ids))
      for (i in seq_along(ids)) {
        id <- ids[i]
        if (fr$var[i] == "<leaf>") {
          st <- stat_map[[as.character(id)]]
          size_pct <- if (!is.null(st) && !is.na(st$pct)) 100 * st$pct else 0
          lab <- node_pred_label[which(node_ids == id)]
          out[i] <- if (is_keep_label(lab)) sprintf("%.0f%%", size_pct) else sprintf("UNDERREP\n%.0f%%", size_pct)
        } else {
          out[i] <- labs[i]
        }
      }
      out
    }

    rp <- NULL
    try({
      rpart.plot::prp(
        f,
        main             = main,
        type             = 2,
        extra            = 0,
        under            = TRUE,
        faclen           = 0,
        tweak            = 1.1,
        fallen.leaves    = TRUE,
        branch.lty       = 3,
        shadow.col       = "gray",
        box.col          = box_col,
        node.fun         = node_fun,
        split.box.col    = NA,
        split.border.col = "gray40",
        branch.col       = "gray40",
        roundint         = FALSE
      )
      old_par <- graphics::par(xpd = NA); on.exit(graphics::par(old_par), add = TRUE)
      graphics::legend("topleft",
                       legend = c("w(x) = 1 Represented", "w(x) = 0 Underrepresented"),
                       fill = c(col_keep, col_drop),
                       border = NA, bty = "n"
      )
      rp <- grDevices::recordPlot()
    }, silent = TRUE)

    rp
  }

  # Optional annotated plot
  tree_plot_underrep <- NULL
  if (plot_underrep && !is.null(root_out$f)) {
    tree_plot_underrep <- .plot_underrep_tree(root_out, main = plot_main)
  }

  # Return
  res <- list(
    root               = root_out,
    combined           = combined,
    leaf_summary       = leaf_summary,
    tree_plot_root     = root_out$tree_plot,
    tree_plot_underrep = tree_plot_underrep
  )
  class(res) <- c("characterizing_underrep", "list")
  return(res)
}
