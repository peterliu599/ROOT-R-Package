#' Summarize a ROOT fit
#'
#' Provides a human-readable summary of a \code{ROOT} object, including:
#' \enumerate{
#'   \item the summary characterization tree \code{f},
#'   \item the first few rows of \code{testing_data},
#'   \item the \code{global_objective_fn} used during optimization, and
#'   \item in generalization mode (\code{generalization = TRUE}), the
#'         unweighted and weighted estimands with their standard errors
#'         and an explanatory note for the weighted SE.
#' }
#'
#' When \code{generalization = TRUE}, the unweighted estimand corresponds
#' to a SATE-type quantity and the weighted estimand to a WTATE-type
#' quantity for the transported target population. When \code{generalization = FALSE},
#' ROOT is used for general functional optimization and no causal labels
#' are imposed; the summary focuses on the tree and diagnostics.
#'
#' @section Diagnostics:
#' The summary also reports:
#' \itemize{
#'   \item the number of trees grown,
#'   \item the size of the Rashomon set,
#'   \item the percentage of observations with ensemble vote \code{w_opt == 1}.
#' }
#'
#' @param object A \code{"ROOT"} S3 object returned by \code{ROOT()}.
#' @param ... Currently unused and included for S3 compatibility.
#'
#' @return \code{object} returned invisibly. Printed output is for inspection.
#'
#' @method summary ROOT
#' @examples
#' \dontrun{
#' ROOT.output = ROOT(diabetes_data,generalizability_path = TRUE, seed = 123)
#' summary(ROOT.output)
#' }
#' @export
summary.ROOT <- function(object, ...) {
  x <- object

  cat("ROOT object\n")
  cat("  Generalization mode:", isTRUE(x$generalization), "\n\n")

  ## Summary tree
  cat("Summary classifier (f):\n")
  if (!is.null(x$f)) {
    print(x$f)
  } else {
    cat("  <no summary tree available>\n")
  }
  cat("\n")

  ## Testing data
  #cat("Testing data (head):\n")
  #if (!is.null(x$testing_data)) {
  #  print(utils::head(x$testing_data))
  #} else {
  #  cat("  <no testing_data stored>\n")
  #}
  #cat("\n")

  ## Objective function
  cat("Global objective function:\n")
  if (!is.null(x$global_objective_fn)) {
    print(x$global_objective_fn)
  } else {
    cat("  <no global_objective_fn stored>\n")
  }
  cat("\n")

  ## Estimands (only in generalization mode)
  if (isTRUE(x$generalization_path) && !is.null(x$estimate)) {
    est <- x$estimate
    cat("Estimand summary (generalization mode):\n")
    cat("  Unweighted ", est$estimand_unweighted,
        " = ", est$value_unweighted,
        ", SE = ", est$se_unweighted, "\n", sep = "")
    cat("  Weighted   ", est$estimand_weighted,
        " = ", est$value_weighted,
        ", SE = ", est$se_weighted, "\n", sep = "")
    cat("  Note: ", est$se_weighted_note, "\n", sep = "")
    cat("\n")
  }

  ## Core diagnostics
  cat("Diagnostics:\n")
  # number of trees
  if (!is.null(x$w_forest)) {
    cat("  Number of trees grown: ", length(x$w_forest), "\n", sep = "")
  } else {
    cat("  Number of trees grown: <unknown>\n")
  }

  # Rashomon size
  if (!is.null(x$rashomon_set)) {
    cat("  Rashomon set size: ", length(x$rashomon_set), "\n", sep = "")
  } else {
    cat("  Rashomon set size: <unknown>\n")
  }

  # % with w_opt == 1
  if (!is.null(x$D_rash) && "w_opt" %in% names(x$D_rash)) {
    pct_keep <- mean(x$D_rash$w_opt == 1, na.rm = TRUE) * 100
    cat(sprintf("  %% observations with w_opt == 1: %.1f%%\n", pct_keep))
  } else {
    cat("  % observations with w_opt == 1: <not available>\n")
  }

  invisible(x)
}


#' Plot the ROOT summary tree
#'
#' Visualizes the decision tree that characterizes the weighted subgroup
#' (the weight function \eqn{w(d)} in \code{\{0,1\}}) identified by \code{ROOT()},
#' using \code{rpart.plot::prp()}.
#'
#' @param x A \code{"ROOT"} S3 object returned by \code{ROOT()} with
#'   \code{x$f} an \code{rpart} object representing the summary / characterization tree.
#' @param ... Additional arguments passed to \code{rpart.plot::prp()}.
#'
#' @return No return value; the plot is drawn to the active graphics device.
#' @importFrom rpart.plot prp
#'
#' @examples
#' \dontrun{
#' ROOT.output = ROOT(diabetes_data,generalizability_path = TRUE, seed = 123)
#' plot(ROOT.output)
#' }
#' @export
plot.ROOT <- function(x, ...) {
  if (is.null(x$f) || !inherits(x$f, "rpart")) {
    message("No summary tree available to plot (x$f is NULL or not an 'rpart' object).")
    return(invisible(NULL))
  }

  # Default arguments (can be overridden by ...)
  args <- list(...)
  if (!"type" %in% names(args))          args$type <- 2
  if (!"extra" %in% names(args))         args$extra <- 109
  if (!"under" %in% names(args))         args$under <- TRUE
  if (!"faclen" %in% names(args))        args$faclen <- 0
  if (!"tweak" %in% names(args))         args$tweak <- 1.1
  if (!"fallen.leaves" %in% names(args)) args$fallen.leaves <- TRUE
  if (!"box.palette" %in% names(args))   args$box.palette <- c("pink", "palegreen3")
  if (!"shadow.col" %in% names(args))    args$shadow.col <- "gray"
  if (!"branch.lty" %in% names(args))    args$branch.lty <- 3
  if (!"main" %in% names(args))          args$main <- "Final Characterized Tree from Rashomon Set"

  args$x <- x$f
  do.call(rpart.plot::prp, args)
}
