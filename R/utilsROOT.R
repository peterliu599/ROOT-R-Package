#' Detect single-sample mode in a ROOT fit
#'
#' Returns \code{TRUE} if a \code{ROOT} object represents single-sample (trial-only)
#' mode (ATE in RCT/Weighted ATE in RCT), and \code{FALSE} otherwise (TATE/WTATE).
#'
#' @param object A \code{ROOT} object (list with component \code{D_forest}).
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.root_is_single_sample <- function(object) {
  # Prefer the explicit flag if present
  if (!is.null(object$single_sample_mode)) return(isTRUE(object$single_sample_mode))
  # Fallback: infer from lX only (do NOT check S=1, since testing_data is S=1 by design)
  lx <- if ("lX" %in% names(object$D_forest)) object$D_forest$lX else NA_real_
  all(is.na(lx))
}

#' Extract covariate names used by ROOT
#'
#' Returns the names of covariate columns in \code{object$D_forest}, excluding internal
#' columns such as \code{v}, \code{vsq}, \code{S}, \code{lX}, and any \code{w_tree_*}.
#'
#' @param object A \code{ROOT} object.
#' @return Character vector of column names.
#' @keywords internal
#' @noRd
.root_covariate_names <- function(object) {
  wt_cols <- grep("^w_tree_", names(object$D_forest), value = TRUE)
  setdiff(names(object$D_forest), c("v","vsq","S","lX", wt_cols))
}

#' Compute the baseline loss for a ROOT fit
#'
#' Computes the baseline loss (no selection) as \eqn{\sqrt{\sum vsq}/n^2} using
#' \code{object$D_forest$vsq} and the number of rows \code{n}.
#'
#' @param object A \code{ROOT} object.
#' @return Numeric scalar (baseline loss).
#' @keywords internal
#' @noRd
.root_baseline_loss <- function(object) {
  n <- nrow(object$D_forest)
  sqrt(sum(object$D_forest$vsq, na.rm = TRUE) / (n^2))
}

#' Get objectives of the selected Rashomon trees
#'
#' Returns the vector of objective values for trees included in the Rashomon set.
#' If none are selected, returns \code{numeric(0)}.
#'
#' @param object A \code{ROOT} object.
#' @return Numeric vector of objectives (possibly empty), named by tree index.
#' @keywords internal
#' @noRd
.root_selected_objectives <- function(object) {
  if (is.null(object$w_forest) || !length(object$w_forest)) return(numeric(0))
  vals <- vapply(
    object$w_forest,
    function(t) {
      v <- t[["local objective"]]
      if (is.null(v)) NA_real_ else unname(v)
    },
    numeric(1)
  )
  if (length(object$rashomon_set) == 0) return(numeric(0))
  out <- vals[object$rashomon_set]
  names(out) <- paste0("w_tree_", object$rashomon_set)
  out
}

#' Summarize a ROOT fit
#'
#' Summarizes a \code{ROOT} object by reporting the primary estimands and key
#' model diagnostics. The first lines report:
#' \enumerate{
#'   \item the \strong{unweighted} estimate (ATE in RCT for single-sample or TATE for two-sample)
#'         and its \strong{standard error (SE)};
#'   \item the \strong{weighted} estimate (WATE in RCT or WTATE using \code{w_opt}) and its
#'         \strong{SE} whenever \code{w_opt} is effectively \emph{binary} (subset-mean SE);
#'         if \code{w_opt} is non-binary, the SE is omitted with a note.
#' }
#' Subsequent lines describe the estimand type, number of trees, size of the Rashomon set,
#' presence of a summary tree, covariate count, observation count, baseline loss,
#' selected-tree losses, and the proportion kept by \code{w_opt}.
#'
#' @param object A \code{ROOT} object returned by \code{ROOT()}.
#' @param ... Unused; included for S3 compatibility.
#'
#' @return The input \code{object}, invisibly. Printed output is a human-readable summary.
#'
#' @details
#' This method prefers pre-computed estimates. If unavailable, it recomputes:
#' \itemize{
#'   \item Unweighted effect as \eqn{\bar v}{mean(v)} over the analysis set
#'         (all rows in single-sample; \code{S = 1} in two-sample);
#'   \item Unweighted SE as \eqn{\sqrt{ \frac{1}{n(n-1)} \sum (v_i - \bar v)^2 }}
#'         {\code{sqrt( sum((v - vbar)^2) / ( n * (n - 1) ) )}};
#'   \item Weighted effect when \code{w_opt} is binary, with
#'         \eqn{A = \{ i : w_i = 1 \}}{\code{A <- which(w == 1)}}, i.e.
#'         \eqn{\bar v_A = \frac{1}{n_A}\sum_{i \in A} v_i}{\code{mean(v[w == 1])}};
#'   \item Weighted SE (WTATE/WATE), for binary \code{w_opt}, as
#'         \eqn{\sqrt{ \frac{1}{n_A(n_A-1)} \sum_{i \in A} (v_i - \bar v_A)^2 }}
#'         {\code{sqrt( sum( (v[w == 1] - vbarA)^2 ) / ( nA * (nA - 1) ) )}}.
#' }
#'
#' @method summary ROOT
#' @export
summary.ROOT <- function(object, ...) {
  if (!inherits(object, "ROOT")) stop("Not a ROOT object.")

  # Analysis set (single-sample: all rows; two-sample: S=1)
  single <- .root_is_single_sample(object)
  in_S <- if (single) rep(TRUE, nrow(object$D_forest)) else (object$D_forest$S == 1L)

  # Prefer precomputed estimates from ROOT()
  if (!is.null(object$estimate)) {
    eu <- object$estimate

    # Unweighted line
    if (!is.null(eu$estimand_unweighted) && !is.null(eu$value_unweighted) && !is.null(eu$se_unweighted)) {
      cat(sprintf("%s (unweighted) = %.6f, SE = %.6f\n",
                  eu$estimand_unweighted, eu$value_unweighted, eu$se_unweighted))
    }

    # Weighted line + always print the explanatory note (if present)
    if (!is.null(eu$estimand_weighted) && !is.null(eu$value_weighted)) {
      if (!is.null(eu$se_weighted) && is.finite(eu$se_weighted)) {
        cat(sprintf("%s (weighted)   = %.6f, SE = %.6f\n",
                    eu$estimand_weighted, eu$value_weighted, eu$se_weighted))
      } else {
        cat(sprintf("%s (weighted)   = %.6f\n", eu$estimand_weighted, eu$value_weighted))
      }
      if (!is.null(eu$se_weighted_note) && nzchar(eu$se_weighted_note)) {
        cat(sprintf("  Note: %s\n", eu$se_weighted_note))
      }
    }

  } else {
    # Fallback recompute (SE-only) with updated rules
    label_unw <- if (single) "ATE in RCT" else "TATE"
    label_w   <- if (single) "Weighted ATE in RCT" else "WTATE"

    v <- object$D_forest$v[in_S]
    w <- if ("w_opt" %in% names(object$D_rash)) object$D_rash$w_opt[in_S] else rep(1L, length(v))

    # Unweighted
    ok_unw <- !is.na(v)
    n_unw  <- sum(ok_unw)
    mu_unw <- if (n_unw > 0) mean(v[ok_unw]) else NA_real_
    se_unw <- if (n_unw > 1) sqrt(sum((v[ok_unw] - mu_unw)^2) / (n_unw * (n_unw - 1))) else NA_real_
    cat(sprintf("%s (unweighted) = %.6f, SE = %.6f\n", label_unw, mu_unw, se_unw))

    # Weighted
    is_binary_w <- all(is.na(w) | (w %in% c(0L, 1L)))
    if (is_binary_w) {
      ok_w <- ok_unw & (w == 1L) & !is.na(w)
      n_A  <- sum(ok_w)
      if (n_A > 0) {
        mu_w <- mean(v[ok_w])
        se_w <- if (n_A > 1) sqrt(sum((v[ok_w] - mu_w)^2) / (n_A * (n_A - 1))) else NA_real_
        cat(sprintf("%s (weighted)   = %.6f, SE = %.6f\n", label_w, mu_w, se_w))

        # Note block (always shown for binary w)
        se_note <- paste0(
          "Calculation of SE for WTATE uses sqrt( sum_{A} (v_i - vbar_A)^2 / ( n_A * (n_A - 1) ) ).\n",
          "  Here A = { i : w_i = 1 }, n_A = |A|, v_i are unit-level orthogonal scores, and vbar_A is their mean. \n"
        )
        if (exists("global_objective_fn") && !identical(object$global_objective_fn, objective_default)) {
          se_note <- paste0(
            se_note,
            "\nYou supplied a custom global_objective_fn; please verify this SE matches your ",
            "estimand, or perform additional variance analysis (e.g., bootstrap or ",
            "influence-function methods)."
          )
        }
        cat(se_note)
      } else {
        cat(sprintf("%s (weighted)   = NA (no kept observations)\n", label_w))
        cat("Note on SE:\n  SE omitted because no kept observations (A = { i : w_i = 1 } is empty)..\n")
      }
    } else {
      mu_w <- sum(w * v, na.rm = TRUE) / sum(w, na.rm = TRUE)
      cat(sprintf("%s (weighted)   = %.6f\n", label_w, mu_w))

      # Non-binary weights: explain why SE is omitted
      se_note <- paste0(
        "Note on SE:\n",
        "  SE omitted: non-binary w_opt detected; the subset-mean SE above is not appropriate.\n",
        "  Consider bootstrap or influence-function methods tailored to your objective.\n"
      )
      if (exists("global_objective_fn") && !identical(object$global_objective_fn, objective_default)) {
        se_note <- paste0(
          se_note,
          "\nYou supplied a custom global_objective_fn; please ensure your variance method ",
          "matches that estimand.\n"
        )
      }
      cat(se_note)
    }
  }


  # Diagnostics
  cat("ROOT object\n")
  estimand <- if (single) "ATE in RCT (single-sample)" else "TATE/PATE (transported ATE)"
  cat("  Estimand:       ", estimand, "\n", sep = "")

  wt_cols <- grep("^w_tree_", names(object$D_forest), value = TRUE)
  cat("  Trees grown:    ", length(wt_cols), "\n", sep = "")
  cat("  Rashomon size:  ", length(object$rashomon_set), "\n", sep = "")

  if (!is.null(object$f)) cat("  Summary tree:    present (rpart)\n") else cat("  Summary tree:    none\n")

  covs <- .root_covariate_names(object)
  cat("  Covariates:     ", length(covs), " (", paste(utils::head(covs, 4), collapse = ", "),
      if (length(covs) > 4) ", ..." else "", ")\n", sep = "")
  cat("  Observations:   ", nrow(object$testing_data), " (analysis set for trees)\n", sep = "")

  base <- .root_baseline_loss(object)
  sel  <- .root_selected_objectives(object)
  cat("  Baseline loss:  ", sprintf("%.5f", base), "\n", sep = "")
  if (length(sel)) {
    cat("  Selected loss:  min/median = ",
        sprintf("%.5f/%.5f", min(sel, na.rm = TRUE), stats::median(sel, na.rm = TRUE)), "\n", sep = "")
  } else {
    cat("  Selected loss:  (no trees selected)\n")
  }

  if ("w_opt" %in% names(object$D_rash)) {
    keep_rate <- mean(object$D_rash$w_opt[in_S] %in% c(1L, 1), na.rm = TRUE)
    cat("  Kept (w_opt=1): ", sprintf("%.1f%%", 100 * keep_rate), "\n", sep = "")
  }

  invisible(object)
}

#' Plot the ROOT Summary Tree
#'
#' Visualizes the decision tree that characterizes the weighted subgroup identified by ROOT.
#'
#' @param x A \code{ROOT} object returned by \code{ROOT()}.
#' @param ... Additional arguments passed to \code{rpart.plot::prp()}.
#'
#' @return No return value; the plot is drawn to the active graphics device.
#' @importFrom rpart.plot prp
#' @export
plot.ROOT <- function(x, ...) {
  if (is.null(x$f)) {
    message("No summary tree available to plot (possibly no covariates).")
    return(invisible(NULL))
  }

  # Default arguments (can be overridden by ...)
  args <- list(...)
  if (!"type" %in% names(args)) args$type <- 2
  if (!"extra" %in% names(args)) args$extra <- 109
  if (!"under" %in% names(args)) args$under <- TRUE
  if (!"faclen" %in% names(args)) args$faclen <- 0
  if (!"tweak" %in% names(args)) args$tweak <- 1.1
  if (!"fallen.leaves" %in% names(args)) args$fallen.leaves <- TRUE
  if (!"box.palette" %in% names(args)) args$box.palette <- c("pink", "palegreen3")
  if (!"shadow.col" %in% names(args)) args$shadow.col <- "gray"
  if (!"branch.lty" %in% names(args)) args$branch.lty <- 3
  if (!"main" %in% names(args)) args$main <- "Final Characterized Tree from Rashomon Set"

  # Add the model object
  args$x <- x$f

  # Call rpart.plot::prp
  do.call(rpart.plot::prp, args)
}

