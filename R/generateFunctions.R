#' Generate covariates X and potential outcomes (Y0, Y1)
#'
#' Simulates a regression problem based on Friedman number one and defines a
#' heterogeneous treatment effect. Uses \code{mlbench.friedman1} to generate
#' feature matrix \code{X} and a baseline outcome \code{Y0}. The treatment
#' potential outcome \code{Y1} is defined as \code{Y1 = Y0 + log(Y0 + 1)}.
#'
#' @param n A \code{numeric(1)} or \code{integer(1)} giving the number of
#'   observations to simulate. Must be positive.
#' @param seed Optional \code{numeric(1)} for the random number generator seed.
#'   If \code{NULL} no seed is set.
#'
#' @return A \code{list} with two components:
#'   \item{X}{\code{data.frame} of simulated covariates with columns
#'     \code{X0, X1, ...}.}
#'   \item{Y}{\code{data.frame} with columns \code{Y0} and \code{Y1}.}
#'
#' @details The \code{mlbench.friedman1} generator creates ten continuous
#'   features and a baseline outcome \code{Y0} with additive noise. The
#'   potential outcome \code{Y1} adds a nonlinear term \code{log(Y0 + 1)} to
#'   \code{Y0}.
#'
#' @keywords internal
#' @noRd
gen_XY <- function(n = 1000, seed = NULL) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a positive numeric value.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }
  if (!requireNamespace("mlbench", quietly = TRUE)) {
    stop(
      "Package 'mlbench' is required for gen_XY(); please install it ",
      "or use a custom generator.",
      call. = FALSE
    )
  }

  # Generate Friedman #1 data
  friedman_data <- if (is.null(seed)) {
    mlbench::mlbench.friedman1(n = n, sd = 1)
  } else {
    withr::with_seed(seed, mlbench::mlbench.friedman1(n = n, sd = 1))
  }

  X  <- as.data.frame(friedman_data$x)
  p  <- ncol(X)
  colnames(X) <- paste0("X", seq_len(p) - 1L)

  Y0 <- as.numeric(friedman_data$y)
  Y1 <- Y0 + log(Y0 + 1)
  Y  <- data.frame(Y0 = Y0, Y1 = Y1)

  list(X = X, Y = Y)
}

#' Generate sample indicator S drawn from a Bernoulli distribution
#'
#' Generates a binary sample inclusion indicator \code{S} for each observation,
#' using a logistic model influenced by a rectangular region in \code{X0} and
#' \code{X1}. In the ROOT generalizability path, \code{S = 1} is interpreted as
#' belonging to the trial sample and \code{S = 0} as belonging to the target
#' sample.
#'
#' @param X A \code{data.frame} of covariates that contains at least
#'   \code{X0} and \code{X1}.
#' @param seed Optional \code{numeric(1)} seed. If \code{NULL} no seed is set.
#'
#' @return A \code{data.frame} with one column \code{S} of values in \code{0}
#'   or \code{1}.
#'
#' @details The inclusion probability is \eqn{p = plogis(a)} where
#'   \eqn{a = 0.25 - 2 * I(X0 in (0.5,1) and X1 in (0.5,1))}.
#'
#' @keywords internal
#' @noRd
gen_S <- function(X, seed = NULL) {
  # Input validation
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  if (!all(c("X0", "X1") %in% names(X))) {
    stop(
      "`X` must contain columns 'X0' and 'X1' for computing inclusion probabilities.",
      call. = FALSE
    )
  }
  .check_no_na(X, c("X0", "X1"))
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }

  # Linear predictor 'a' defining inclusion odds
  a <- 0.25 - 2 * (X$X0 > 0.5 & X$X0 < 1 & X$X1 > 0.5 & X$X1 < 1)
  # Inclusion probability p = logit^{-1}(a)
  p <- stats::plogis(a)

  # Draw S ~ Bernoulli(p) for each observation
  S_vec <- if (is.null(seed)) {
    stats::rbinom(n = nrow(X), size = 1, prob = p)
  } else {
    withr::with_seed(seed + 1L, stats::rbinom(n = nrow(X), size = 1, prob = p))
  }

  data.frame(S = as.integer(S_vec))
}

#' Generate treatment indicator Tr drawn from a Bernoulli distribution
#'
#' Assigns treatment indicators combining a randomized design for \code{S == 1}
#' (trial sample) and an observational assignment driven by \code{X0} for
#' \code{S == 0} (target sample).
#'
#' @param X A \code{data.frame} of covariates.
#' @param S A \code{data.frame} with column \code{S} in \code{0} or \code{1}.
#'   This should be the output of \code{gen_S()}.
#' @param seed Optional \code{numeric(1)} seed. If \code{NULL} no seed is set.
#'
#' @return A \code{list} with:
#'   \item{Tr}{\code{data.frame} with one column \code{Tr} in \code{0} or \code{1}.}
#'   \item{pi}{\code{numeric} vector of assignment probabilities per observation.}
#'
#' @details For \code{S == 1} the treatment probability is \code{0.5}. For
#'   \code{S == 0} the treatment probability is \code{plogis(X0)}. The combined
#'   probability is \eqn{pi_i = S_i * 0.5 + (1 - S_i) * plogis(X0_i)}.
#'
#' @keywords internal
#' @noRd
gen_T <- function(X, S, seed = NULL) {
  # Input validation
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  if (!is.data.frame(S) || !"S" %in% names(S)) {
    stop("`S` must be a data frame with a column named 'S'.", call. = FALSE)
  }
  if (nrow(S) != nrow(X)) {
    stop("`X` and `S` must have the same number of rows.", call. = FALSE)
  }
  .check_no_na(S, colnames(S))
  .check_no_na(X, colnames(X))
  if (!all(S$S %in% c(0, 1))) {
    stop("`S$S` must contain only 0 or 1 values.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }

  # Define probabilities for treatment assignment
  pi_exp <- 0.5                     # experimental probability for S == 1
  pi_obs <- stats::plogis(X$X0)     # observational probability for S == 0

  # Combined assignment probability for each observation:
  pi_vec <- S$S * pi_exp + (1 - S$S) * pi_obs

  # Draw treatment assignments
  Tr_vec <- if (is.null(seed)) {
    stats::rbinom(n = nrow(X), size = 1, prob = pi_vec)
  } else {
    withr::with_seed(seed + 2L, stats::rbinom(n = nrow(X), size = 1, prob = pi_vec))
  }

  list(
    Tr = data.frame(Tr = as.integer(Tr_vec)),
    pi = pi_vec
  )
}


#' Convenience wrapper to generate a full simulated dataset
#'
#' Generates covariates, sample inclusion, treatment assignments, and observed
#' outcomes for a given sample size by calling \code{gen_XY()}, \code{gen_S()},
#' and \code{gen_T()}. The returned \code{data} element is ready for use with
#' \code{ROOT(data = ..., generalizability_path = TRUE)}, since it contains
#' columns \code{Y}, \code{Tr}, and \code{S} in the expected format.
#'
#' @param n A \code{numeric(1)} or \code{integer(1)} sample size.
#' @param seed Optional \code{numeric(1)} base seed. If provided, internal
#'   generators use simple offsets for reproducibility.
#'
#' @return A \code{list} with:
#'   \item{data}{\code{data.frame} of length \code{n} with covariates
#'     \code{X0, X1, ...}, observed outcome \code{Y}, sample indicator
#'     \code{S}, and treatment indicator \code{Tr}. This is directly usable by
#'     \code{ROOT} in \code{generalizability_path = TRUE} mode.}
#'   \item{Y}{\code{data.frame} of potential outcomes with columns \code{Y0} and
#'     \code{Y1}.}
#'
#' @examples
#' sim <- get_data(n = 100, seed = 599)
#' dim(sim$data)
#' head(sim$data$Y)
#' head(sim$Y)
#'
#' @keywords internal
#' @noRd
get_data <- function(n = 1000, seed = NULL) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a positive numeric value.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }

  # Generate covariates and potential outcomes
  xy <- gen_XY(n = n, seed = seed)
  X  <- xy$X
  Y  <- xy$Y

  # Generate sample inclusion and treatment assignment
  S_df <- gen_S(X, seed = seed)
  T_lst <- gen_T(X, S_df, seed = seed)

  Tr_vec <- as.numeric(T_lst$Tr$Tr)

  # Construct observed outcome Y (this is what ROOT uses as "Y")
  Y_obs <- Tr_vec * Y$Y1 + (1 - Tr_vec) * Y$Y0

  # Combine into one data frame (no Y0/Y1 here so ROOT won't treat them as covariates)
  full_data <- cbind(
    X,
    Y  = Y_obs,
    S  = as.integer(S_df$S),
    Tr = as.integer(Tr_vec)
  )

  list(data = full_data, Y = Y)
}
