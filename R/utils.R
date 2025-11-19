#' Check for missing values in training data
#'
#' Ensures there are no NA values in any of the relevant columns of the dataset.
#'
#' @param data A data frame.
#' @param cols Character vector of column names to check.
#' @return Invisibly returns TRUE if no NA found; otherwise throws an error.
#' @keywords internal
check_no_na <- function(data, cols) {
  data_name <- deparse(substitute(data))
  for (col in cols) {
    if (anyNA(data[[col]])) {
      stop(sprintf("Data `%s` column '%s' contains missing values. Please handle NA before training.", data_name, col),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Generic objective interface (source of truth)
#'
#' Default objective: SE proxy of (W)TATE/PATE
#'
#' Computes \eqn{\sqrt{\sum_i vsq_i * w_i / (\sum_i w_i)^2}}.
#' Requires columns `vsq` and `w` in `D`. Minimize this.
#' Supply your own function(D) -> scalar to use a different objective.
#' @param D data.frame with at least numeric columns `vsq` and `w`.
#' @return numeric scalar objective value; `Inf` if undefined.
objective_default <- function(D) {
  stopifnot(is.data.frame(D), all(c("vsq","w") %in% names(D)))
  num <- sum(D$vsq * D$w, na.rm = TRUE)
  den <- (sum(D$w, na.rm = TRUE))^2
  out <- sqrt(num / den)
  if (!is.finite(out) || is.nan(out)) Inf else out
}

#' Helper: evaluate objective after a hypothetical local change
#'
#' @param val 0/1 assignment to apply
#' @param indices integer or rownames to receive `val`
#' @param D data.frame used by `global_objective_fn`
#' @param global_objective_fn function(D)->scalar
#' @return numeric scalar objective after the hypothetical change
objective_if <- function(val, indices, D, global_objective_fn) {
  stopifnot(is.function(global_objective_fn), length(val) == 1, val %in% c(0,1))
  rows <- integer(0)
  if (length(indices)) {
    rows <- if (is.numeric(indices)) as.integer(indices)
    else which(rownames(D) %in% as.character(indices))
  }
  Dtmp <- D
  if (length(rows)) Dtmp[rows, "w"] <- val
  global_objective_fn(Dtmp)
}

#' Backward/fast-path micro-evaluator adaptor
#'
#' Wrap a global objective \code{global_objective_fn(D)} into a splitter-compatible
#' loss function \code{loss_fn(val, indices, D)} by evaluating
#' \code{\link{objective_if}} on a temporary copy of \code{D}.
#'
#' @param global_objective_fn Function of one argument \code{D} returning a numeric
#'   scalar to be minimized (e.g., \code{\link{objective_default}}).
#'
#' @return A function \code{loss_fn(val, indices, D)} suitable for use in
#'   \code{\link{ROOT}} and \code{\link{split_node}}. It sets \code{w = val}
#'   on \code{indices} (non-mutating), then returns \code{global_objective_fn(D)}.
loss_from_objective <- function(global_objective_fn) {
  force(global_objective_fn)
  function(val, indices, D) objective_if(val, indices, D, global_objective_fn)
}
