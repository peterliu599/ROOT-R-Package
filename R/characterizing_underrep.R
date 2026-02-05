#' Characterize under-represented subgroups (wrapper around ROOT)
#'
#' A convenience wrapper around \code{ROOT()} for under-representation analyses.
#' It takes a single \code{data} set and calls \code{ROOT()} either in
#' generalizability_path mode (when \code{generalizability_path = TRUE}) or in
#' general optimization mode (\code{generalizability_path = FALSE}).
#'
#' When \code{generalizability_path = TRUE}, \code{data} must contain
#' standardized columns:
#' \itemize{
#'   \item \code{Y}: outcome,
#'   \item \code{Tr}: treatment indicator (0/1),
#'   \item \code{S}: sample indicator (1 = trial, 0 = target).
#' }
#'
#' @param data A \code{data.frame} containing covariates and, in generalizability_path
#'   mode, also columns \code{Y}, \code{Tr}, and \code{S}.
#' @param global_objective_fn function with signature `function(D) -> numeric`
#'   scoring the entire state and minimized by ROOT. If NULL, a default
#'   variance-based objective is used (see `objective_default()`).
#' @param generalizability_path Logical. If \code{TRUE}, calls
#'   \code{ROOT()} with \code{generalizability_path = TRUE} and expects columns
#'   \code{Y}, \code{Tr}, and \code{S} in \code{data}. If \code{FALSE}, calls
#'   \code{ROOT()} in general optimization mode. Default \code{FALSE}.
#' @param leaf_proba A \code{numeric(1)} tuning parameter that increases the chance
#'   a node stops splitting by selecting a synthetic \code{"leaf"} feature.
#'   Internally, the probability of choosing \code{"leaf"} is
#'   \code{leaf_proba / (1 + leaf_proba)} (assuming the
#'   covariate probabilities sum to 1). Default \code{0.25}.
#' @param seed Random seed for reproducibility.
#' @param num_trees Number of trees to grow.
#' @param vote_threshold Majority vote threshold used for \code{w_opt}.
#' @param explore_proba Exploration probability in tree growth.
#' @param feature_est Either \code{"Ridge"}, \code{"GBM"}, or a custom
#'   feature importance function.
#' @param feature_est_args List of extra arguments passed to \code{feature_est}
#'   when it is a function.
#' @param top_k_trees Logical; if \code{TRUE}, uses top \code{k} trees by
#'   objective, otherwise a cutoff rule.
#' @param k Number of trees when \code{top_k_trees = TRUE}.
#' @param cutoff Numeric or \code{"baseline"} Rashomon cutoff.
#' @param verbose Logical; if \code{TRUE}, prints progress/estimands from
#'   \code{ROOT()}.
#'
#' @references
#' Parikh H, Ross RK, Stuart E, Rudolph KE (2025).
#' "Who Are We Missing?: A Principled Approach to Characterizing the Underrepresented Population."
#' *Journal of the American Statistical Association*.
#' doi:10.1080/01621459.2025.2495319
#'
#' @return A \code{characterizing_underrep} S3 object (a \code{list}) with:
#'   \item{root}{The \code{ROOT} object returned by \code{ROOT()}.}
#'   \item{combined}{The input \code{data} (for continuity with prior API).}
#'   \item{leaf_summary}{Data frame of terminal node rules and labels, or \code{NULL}.}
#'
#' @examples
#' \dontrun{
#' char.output = characterizing_underrep(diabetes_data,generalizability_path = TRUE, seed = 123)
#' }
#' @export
characterizing_underrep <- function(data,
                                    global_objective_fn   = NULL,
                                    generalizability_path = FALSE,
                                    leaf_proba            = 0.25,
                                    seed                  = 123,
                                    num_trees             = 10,
                                    vote_threshold        = 2 / 3,
                                    explore_proba         = 0.05,
                                    feature_est           = "Ridge",
                                    feature_est_args      = list(),
                                    top_k_trees           = FALSE,
                                    k                     = 10,
                                    cutoff                = "baseline",
                                    verbose               = FALSE) {
  # Data frame check
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.logical(generalizability_path) || length(generalizability_path) != 1L) {
    stop("`generalizability_path` must be TRUE or FALSE.", call. = FALSE)
  }

  # If we are in generalizability_path mode, enforce Y / Tr / S presence
  if (isTRUE(generalizability_path)) {
    needed <- c("Y", "Tr", "S")
    missing_cols <- setdiff(needed, names(data))
    if (length(missing_cols) > 0L) {
      stop(
        "For generalizability_path = TRUE, `data` must contain columns: ",
        paste(needed, collapse = ", "),
        ". Missing: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # 1. Call ROOT with the requested path
  root_out <- ROOT(
    data                = data,
    generalizability_path = generalizability_path,
    leaf_proba          = leaf_proba,
    seed                = seed,
    num_trees           = num_trees,
    vote_threshold      = vote_threshold,
    explore_proba       = explore_proba,
    feature_est         = feature_est,
    feature_est_args    = feature_est_args,
    top_k_trees         = top_k_trees,
    k                   = k,
    cutoff              = cutoff,
    verbose             = verbose,
    global_objective_fn = global_objective_fn
  )

  # 2. Summarize terminal nodes (if characterization tree exists)
  leaf_summary <- NULL
  if (!is.null(root_out$f)) {
    f   <- root_out$f
    frm <- f$frame
    is_leaf <- frm$var == "<leaf>"

    paths <- rpart::path.rpart(
      f,
      nodes    = as.numeric(rownames(frm)[is_leaf]),
      print.it = FALSE
    )
    path_text <- vapply(paths, function(p) paste(p, collapse = " & "), character(1))

    ylv <- if (!is.null(f$ylevels)) f$ylevels else levels(stats::model.frame(f)$w)
    pred_class <- ylv[frm$yval[is_leaf]]

    n_in_leaf <- frm$n[is_leaf]
    total_n   <- sum(n_in_leaf)
    pct       <- if (total_n > 0) round(100 * n_in_leaf / total_n, 1) else NA_real_

    label <- ifelse(
      pred_class %in% c("0", "0.0", "0L", "FALSE"),
      "Under-represented (drop, w = 0)",
      "Represented (keep, w = 1)"
    )

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

  # 3. Construct and return S3 object
  res <- list(
    root         = root_out,
    combined     = data,
    leaf_summary = leaf_summary
  )
  class(res) <- c("characterizing_underrep", "list")
  res
}
