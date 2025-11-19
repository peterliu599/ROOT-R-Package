#' Summarize a characterizing_underrep fit
#'
#' Prints the \code{ROOT} summary (un/weighted estimates with standard errors; the
#' \emph{weighted} SE is omitted when a custom \code{global_objective_fn} was used in \code{ROOT()})
#' and a brief overview of terminal rules from the annotated summary tree, if available.
#'
#' @param object A \code{characterizing_underrep} object.
#' @param ... Unused; included for S3 compatibility.
#'
#' @return \code{object}, invisibly.
#'
#' @details Delegates core statistics to \code{summary(object$root)}; previews up to
#'   ten terminal rules when a summary tree exists, and reports plot availability.
#'
#' @method summary characterizing_underrep
#' @export
summary.characterizing_underrep <- function(object, ...) {
  if (!inherits(object, "characterizing_underrep")) stop("Not a characterizing_underrep object.")

  cat("characterizing_underrep object\n")
  cat("  --- ROOT summary ---\n")
  summary(object$root)  # uses summary.ROOT() defined above

  # Leaf summary
  if (is.null(object$leaf_summary)) {
    cat("  Leaf summary:    none (no summarized tree)\n")
  } else {
    cat("  Leaf summary:    ", nrow(object$leaf_summary), " terminal nodes\n", sep = "")
    prev <- utils::head(object$leaf_summary, 10)
    # Trim very long rule strings for console
    if ("rule" %in% names(prev)) prev$rule <- substr(prev$rule, 1, 60)
    print(prev, row.names = FALSE)
    if (nrow(object$leaf_summary) > 10) cat("  ...\n")
  }

  # Plots availability
  cat("  Plots available:\n")
  cat("    * ROOT tree:          ", if (!is.null(object$tree_plot_root))     "yes" else "no", "\n", sep = "")
  cat("    * Underrep tree:      ", if (!is.null(object$tree_plot_underrep)) "yes" else "no", "\n", sep = "")

  invisible(object)
}



#' Plot Under-represented Population Characterization
#'
#' Visualizes the decision tree derived from the ROOT analysis, highlighting
#' which subgroups are represented (w=1) versus underrepresented (w=0).
#'
#' @param x A \code{characterizing_underrep} object.
#' @param ... Additional arguments passed to \code{rpart.plot::prp()} (e.g., \code{main}, \code{tweak}).
#'
#' @return No return value; the plot is drawn to the active graphics device.
#' @importFrom rpart.plot prp
#' @importFrom graphics par legend
#' @export
plot.characterizing_underrep <- function(x, ...) {
  # 1. Safety Check
  if (is.null(x$root$f)) {
    message("No summary tree available to plot (possibly no covariates or tree failed to grow).")
    return(invisible(NULL))
  }

  f <- x$root$f
  frm <- f$frame
  is_leaf <- frm$var == "<leaf>"

  # 2. Determine Class Labels (Keep vs Drop)
  # Try to get explicit levels, otherwise assume factor levels of the response
  ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)

  # Helper to check if a label means "Keep" (w=1)
  is_keep_label <- function(lbl) {
    lbl %in% c("1", "1.0", "1L", "TRUE", "True", "yes", "Yes")
  }

  # Map every node to its predicted label
  node_pred_label <- rep(NA_character_, nrow(frm))
  if (length(ylv) >= 2) {
    node_pred_label[is_leaf] <- ylv[frm$yval[is_leaf]]
  } else {
    node_pred_label[is_leaf] <- as.character(frm$yval[is_leaf])
  }

  # 3. Define Colors
  col_keep <- "#4E79A7" # Blue-ish
  col_drop <- "#F28E2B" # Orange-ish

  box_col <- rep(NA, nrow(frm))
  if (any(is_leaf)) {
    labs <- node_pred_label[is_leaf]
    box_col[is_leaf] <- ifelse(is_keep_label(labs), col_keep, col_drop)
  }

  # 4. Define Node Text Function
  # This ensures underrepresented nodes get explicitly labeled "UNDERREP"
  total_n <- frm$n[1] # Root node count

  node_fun <- function(tree_obj, labs, digits, varlen) {
    fr <- tree_obj$frame
    # Calculate percentage of total data in this node
    pct <- if (total_n > 0) fr$n / total_n else 0

    out <- character(nrow(fr))
    for (i in seq_len(nrow(fr))) {
      if (fr$var[i] == "<leaf>") {
        # It's a leaf: formatting logic
        lbl <- node_pred_label[i]
        pct_str <- sprintf("%.0f%%", 100 * pct[i])

        if (is_keep_label(lbl)) {
          # Represented: Just show %
          out[i] <- pct_str
        } else {
          # Underrepresented: Show Warning + %
          out[i] <- paste0("UNDERREP\n", pct_str)
        }
      } else {
        # Internal node: use standard split label
        out[i] <- labs[i]
      }
    }
    out
  }

  # 5. Prepare Arguments for rpart.plot
  # We set defaults, but allow '...' to override them
  args <- list(...)

  # Defaults
  if (is.null(args$type))          args$type <- 2
  if (is.null(args$extra))         args$extra <- 0  # We generate text manually
  if (is.null(args$under))         args$under <- TRUE
  if (is.null(args$faclen))        args$faclen <- 0
  if (is.null(args$tweak))         args$tweak <- 1.1
  if (is.null(args$fallen.leaves)) args$fallen.leaves <- TRUE
  if (is.null(args$shadow.col))    args$shadow.col <- "gray"
  if (is.null(args$branch.lty))    args$branch.lty <- 3
  if (is.null(args$split.border.col)) args$split.border.col <- "gray40"
  if (is.null(args$branch.col))    args$branch.col <- "gray40"
  if (is.null(args$roundint))      args$roundint <- FALSE
  if (is.null(args$main))          args$main <- "Underrepresented Population Characterization Tree"

  # Hard overrides (Logic depends on these)
  args$x <- f
  args$box.col <- box_col
  args$node.fun <- node_fun
  args$split.box.col <- NA # Transparent split boxes

  # 6. Draw Plot
  do.call(rpart.plot::prp, args)

  # 7. Draw Legend
  # We use par(xpd=NA) to allow drawing outside plot margins if needed
  old_par <- graphics::par(xpd = NA)
  on.exit(graphics::par(old_par))

  graphics::legend(
    "topleft",
    legend = c("w(x) = 1 Represented", "w(x) = 0 Underrepresented"),
    fill = c(col_keep, col_drop),
    border = NA,
    bty = "n",
    cex = 0.9
  )
}
