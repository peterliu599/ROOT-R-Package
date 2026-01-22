#' Check for missing values in training data
#'
#' Ensures there are no \code{NA} values in any of the relevant columns of the dataset.
#'
#' @param data A \code{data.frame}.
#' @param cols A \code{character} vector of column names to check in \code{data}.
#'
#' @return Invisibly returns \code{TRUE} if no \code{NA} is found. Otherwise an error is thrown.
#' @keywords internal
.check_no_na <- function(data, cols) {
  data_name <- deparse(substitute(data))
  for (col in cols) {
    if (anyNA(data[[col]])) {
      stop(sprintf("Data `%s` column '%s' contains missing values. Please handle NA before training.", data_name, col),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Generic objective interface
#'
#' Default objective that serves as a proxy for Standard Error in Weighted Transported
#' Average Treatment Effect and Population Average Treatment Effect.
#'
#' Computes \code{sqrt(sum(vsq_i * w_i)) / (sum(w_i))^2}.
#' Requires columns \code{vsq} and \code{w} in \code{D}. The goal is to minimize the value.
#' Supply your own \code{function(D) -> numeric} to use a different objective.
#'
#' @section Abbreviations:
#' ATE means Average Treatment Effect. SE means Standard Error. TATE means Transported ATE.
#' WTATE means Weighted TATE. WATE means Weighted ATE. PATE means Population ATE.
#'
#' @param D A \code{data.frame} with at least numeric columns \code{vsq} and \code{w}.
#'
#' @return A \code{numeric(1)} objective value. Returns \code{Inf} when undefined.
#' @export
objective_default <- function(D) {
  if (!("w" %in% names(D))) {
    stop("objective_default() expects a column `w` in D.", call. = FALSE)
  }
  w <- D$w
  if (all(is.na(w)) || sum(w, na.rm = TRUE) <= 0) return(Inf)

  if ("vsq" %in% names(D) && is.numeric(D$vsq) && any(is.finite(D$vsq))) {
    vsq <- D$vsq
  } else if ("v" %in% names(D) && is.numeric(D$v) && any(is.finite(D$v))) {
    v  <- D$v
    mu <- stats::weighted.mean(v, w = w, na.rm = TRUE)
    vsq <- (v - mu)^2
  } else {
    # Fallback when no variance proxy is available:
    # SE of mean of a Bernoulli with p ≈ mean(w>0) and n = #kept.
    n_keep <- sum(w > 0, na.rm = TRUE)
    if (n_keep <= 1) return(Inf)
    p <- mean(w > 0, na.rm = TRUE)
    return(sqrt(p * (1 - p) / n_keep))
  }

  num <- sum(w * vsq, na.rm = TRUE)
  den <- sum(w, na.rm = TRUE)^2
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(Inf)
  sqrt(num / den)
}

#' Helper to evaluate the objective after a hypothetical local change
#'
#' Evaluates \code{global_objective_fn} on a temporary copy of \code{D} after setting
#' \code{w = val} for the rows selected by \code{indices}.
#'
#' @param val A \code{numeric(1)} that must be either \code{0} or \code{1}.
#' @param indices An \code{integer} vector of row indices or a \code{character} vector of row names that receive \code{val}.
#' @param D A \code{data.frame} used by \code{global_objective_fn}. Must contain columns \code{w} and \code{vsq}.
#' @param global_objective_fn A \code{function} with signature \code{function(D) -> numeric}.
#'
#' @return A \code{numeric(1)} objective value after the hypothetical change.
#' @export
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

#' Backward and fast path micro evaluator adaptor
#'
#' Wraps a global objective \code{global_objective_fn(D)} into a splitter compatible
#' loss function \code{loss_fn(val, indices, D)} by calling \code{\link{objective_if}}
#' on a temporary copy of \code{D}.
#'
#' @param global_objective_fn A \code{function} with signature \code{function(D) -> numeric}
#'   that returns a scalar to be minimized. For example \code{\link{objective_default}}.
#'
#' @return A \code{function} \code{loss_fn(val, indices, D)} suitable for use in
#'   \code{\link{ROOT}} and \code{\link{split_node}}. It sets \code{w = val} on \code{indices}
#'   without mutation of the original \code{D} and then returns \code{global_objective_fn(D)}.
#' @export
loss_from_objective <- function(global_objective_fn) {
  force(global_objective_fn)
  function(val, indices, D) objective_if(val, indices, D, global_objective_fn)
}
