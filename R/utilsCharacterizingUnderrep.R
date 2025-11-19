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
#' @param ... Additional arguments passed to \code{rpart.plot::prp()}.
#' @return No return value; draws a plot.
#' @importFrom rpart.plot prp
#' @importFrom graphics par legend
#' @export
plot.characterizing_underrep <- function(x, ...) {
  if (is.null(x$root$f)) {
    message("No summary tree available to plot (possibly no covariates or tree failed to grow).")
    return(invisible(NULL))
  }

  f   <- x$root$f
  frm <- f$frame
  is_leaf <- frm$var == "<leaf>"

  # --- colors
  col_keep <- "#4E79A7" # represented (w=1)
  col_drop <- "#F28E2B" # underrepresented (w=0)

  # helper to interpret labels as "keep"
  is_keep_label <- function(lbl) {
    lbl %in% c("1", "1.0", "1L", "TRUE", "True", "yes", "Yes", "w=1", "keep")
  }

  # Determine predicted class/value at leaves and map to keep/drop
  keep_leaf <- rep(FALSE, nrow(frm))
  if (!is.null(f$ylevels)) {
    # classification tree: yval indexes into ylevels
    ylv <- f$ylevels
    leaf_lbl <- character(nrow(frm))
    leaf_lbl[is_leaf] <- ylv[frm$yval[is_leaf]]
    keep_leaf[is_leaf] <- is_keep_label(leaf_lbl[is_leaf])
  } else {
    # regression/probability tree: treat >= 0.5 as keep
    keep_leaf[is_leaf] <- frm$yval[is_leaf] >= 0.5
  }

  box_col <- rep(NA_character_, nrow(frm))
  box_col[is_leaf] <- ifelse(keep_leaf[is_leaf], col_keep, col_drop)

  total_n <- frm$n[1]

  # --- IMPORTANT: the formal args must be (x, labs, digits, varlen)
  node_fun <- function(x, labs, digits, varlen) {
    fr <- x$frame
    is_leaf <- fr$var == "<leaf>"
    # percent of total in each node
    pct <- if (total_n > 0) fr$n / total_n else 0

    out <- character(nrow(fr))
    for (i in seq_len(nrow(fr))) {
      if (is_leaf[i]) {
        pct_str <- sprintf("%.0f%%", 100 * pct[i])
        if (keep_leaf[i]) {
          out[i] <- pct_str
        } else {
          out[i] <- paste0("UNDERREP\n", pct_str)
        }
      } else {
        out[i] <- labs[i]  # use prp's default split text for internal nodes
      }
    }
    out
  }

  args <- list(...)
  if (is.null(args$type))             args$type <- 2
  if (is.null(args$extra))            args$extra <- 0
  if (is.null(args$under))            args$under <- TRUE
  if (is.null(args$faclen))           args$faclen <- 0
  if (is.null(args$tweak))            args$tweak <- 1.1
  if (is.null(args$fallen.leaves))    args$fallen.leaves <- TRUE
  if (is.null(args$shadow.col))       args$shadow.col <- "gray"
  if (is.null(args$branch.lty))       args$branch.lty <- 3
  if (is.null(args$split.border.col)) args$split.border.col <- "gray40"
  if (is.null(args$branch.col))       args$branch.col <- "gray40"
  if (is.null(args$roundint))         args$roundint <- FALSE
  if (is.null(args$main))             args$main <- "Underrepresented Population Characterization Tree"

  args$x <- f
  args$box.col <- box_col
  args$node.fun <- node_fun
  args$split.box.col <- NA

  do.call(rpart.plot::prp, args)

  old_par <- graphics::par(xpd = NA)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::legend(
    "topleft",
    legend = c("w(x) = 1 Represented", "w(x) = 0 Underrepresented"),
    fill = c(col_keep, col_drop),
    border = NA,
    bty = "n",
    cex = 0.9
  )
}
