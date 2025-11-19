#' Generate covariates `X` and potential outcomes (`Y0`, `Y1`)
#'
#' Simulates a regression problem (Friedman #1) and defines a treatment effect.
#' Uses `mlbench.friedman1` to generate `X` features and a baseline outcome `Y0`.
#' The treatment potential outcome `Y1` is defined as `Y1 = Y0 + log(Y0 + 1)`,
#' introducing a heterogeneous treatment effect.
#'
#' @param n Integer or numeric. Number of observations to simulate (must be positive).
#' @param seed Optional. Single numeric value for RNG seed. If provided, a global
#'   seed is set for reproducibility. If `NULL` (default), no seed is set (results
#'   will vary on each run).
#' @return A list with two components:
#' \item{X}{A data frame of simulated covariates with columns `X0, X1, ...` up to `X(p-1)`.}
#' \item{Y}{A data frame of potential outcomes with columns `Y0` (baseline outcome)
#'   and `Y1` (outcome under treatment).}
#' @details The `mlbench.friedman1` function from the \pkg{mlbench} package is used
#'   to generate 10 independent continuous features and a baseline outcome `Y0` with additive noise.
#'   The treatment outcome `Y1` is defined by adding a non-linear term `log(Y0 + 1)`
#'   to the baseline. If a seed is specified, the random number generator state is
#'   reset at the start of the function (which affects other random operations).
gen_XY <- function(n = 1000, seed = NULL) {
  # Input validation
  if (!requireNamespace("mlbench", quietly = TRUE)) {
    stop("Package 'mlbench' is required for gen_XY(); please install it or use a custom generator.",
         call. = FALSE)
  }
  friedman_data <- if (is.null(seed)) {
    mlbench::mlbench.friedman1(n = n, sd = 1)
  } else {
    withr::with_seed(seed, mlbench::mlbench.friedman1(n = n, sd = 1))
  }
  X <- friedman_data$x  # matrix of features
  Y0 <- friedman_data$y  # baseline outcome
  Y1 <- Y0 + log(Y0 + 1)  # define treatment outcome

  # Convert X to data frame and name the columns X0, X1, ..., X(p-1)
  X_df <- as.data.frame(X)
  p <- ncol(X_df)
  colnames(X_df) <- paste0("X", seq_len(p) - 1)

  # Assemble Y data frame
  Y_df <- data.frame(Y0 = Y0, Y1 = Y1)

  return(list(X = X_df, Y = Y_df))
}

#' Generate sample indicator `S` ~ Bernoulli(plogis(a))
#'
#' Generates a binary sample inclusion indicator `S` for each observation,
#' using a logistic model influenced by a rectangular region in the first two covariates (`X0` and `X1`).
#'
#' @param X A data frame of covariates (must contain at least columns `X0` and `X1`).
#' @param seed Optional numeric seed for RNG. If provided, `set.seed(seed + 1)`
#'   is invoked for reproducibility. If `NULL` (default), no specific seed is set.
#' @return A data frame with a single column `S` of 0/1 values indicating inclusion (1) or exclusion (0).
#' @details The inclusion probability is defined as \eqn{p = \mathrm{plogis}(a)},
#'   where \eqn{a = 0.25 - 2 * I\{X0, X1 \text{ in region } (0.5,1)\}}.
#'   In other words, observations for which both `X0` and `X1` lie in (0.5, 1)
#'   have a lower odds of being included (due to a negative contribution in the linear predictor).
#'   This mirrors a scenario where a specific region in feature space is under-sampled.
#'   If a seed is set, it uses `seed + 1` to differentiate from other generators.
gen_S <- function(X, seed = NULL) {
  # Input validation
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame of covariates.", call. = FALSE)
  }
  if (!all(c("X0", "X1") %in% names(X))) {
    stop("`X` must contain columns 'X0' and 'X1' for computing the inclusion probabilities.", call. = FALSE)
  }
  check_no_na(X, c("X0", "X1"))
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
    withr::with_seed(seed + 1, stats::rbinom(n = nrow(X), size = 1, prob = p))
  }

  return(data.frame(S = S_vec))
}

#' Generate treatment indicator `Tr` ~ Bernoulli(pi)
#'
#' Assigns a treatment indicator for each observation, combining an experimental design for included samples (S==1)
#' and an observational assignment for excluded samples (S==0).
#'
#' @param X A data frame of covariates.
#' @param S A data frame with column `S` (0/1 indicating sample inclusion for each observation).
#' @param seed Optional numeric seed for RNG. If provided, `set.seed(seed - 1)` is used. Default `NULL` means no explicit seeding.
#' @return A list with two elements:
#' \item{Tr}{A data frame with a single column `Tr` (treatment assignments 0/1 for each observation).}
#' \item{pi}{A numeric vector of length equal to number of observations, giving the treatment probability used for each observation.}
#' @details For observations with `S==1` (in sample), treatment is assigned with probability `0.5` (mimicking a randomized experiment).
#'   For those with `S==0` (out of sample), treatment probability is \eqn{\mathrm{plogis}(X0)}, i.e., it increases with the value of covariate `X0`.
#'   The overall assignment probability for each observation is \eqn{\pi_i = S_i * 0.5 + (1 - S_i) * \mathrm{plogis}(X0_i)}.
#'   If a seed is provided, an offset `seed - 1` is used to differentiate from other generation steps.
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
  check_no_na(S, colnames(S))
  check_no_na(X, colnames(X))
  if (!all(S$S %in% c(0, 1))) {
    stop("`S$S` must contain only 0 or 1 values.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }


  # Define probabilities for treatment assignment
  pi_exp <- 0.5                       # experimental probability for S==1
  pi_obs <- stats::plogis(X$X0)       # observational probability for S==0 (depends on X0)
  # Combined assignment probability for each observation:
  pi_vec <- S$S * pi_exp + (1 - S$S) * pi_obs

  # Draw treatment assignments
  Tr_vec <- if (is.null(seed)) {
    stats::rbinom(n = nrow(X), size = 1, prob = pi_vec)
  } else {
    withr::with_seed(seed + 2, stats::rbinom(n = nrow(X), size = 1, prob = pi_vec))
  }

  return(list(Tr = data.frame(Tr = Tr_vec), pi = pi_vec))
}


#' Convenience wrapper to generate a full simulated dataset
#'
#' Generates covariates, sample inclusion, treatment assignments, and observed outcomes
#' for a specified sample size. This wraps `gen_XY()`, `gen_S()`, and `gen_T()` in sequence.
#'
#' @param n Integer or numeric. Sample size (number of observations to generate).
#' @param seed Optional base seed for reproducibility. If provided, internal generators use
#'   offsets of this seed to ensure independent randomness. Default `NULL` means no explicit seeding.
#' @return A list with two components:
#' \item{data}{A data frame of length `n` containing covariates `X0,...`, sample indicator `S`, treatment indicator `Tr`, and observed outcome `Yobs`.}
#' \item{Y}{A data frame of length `n` containing the potential outcomes `Y0` and `Y1` for each observation.}
#' @details This function first generates covariates and potential outcomes with `gen_XY`.
#'   It then generates `S` (sample inclusion) and `Tr` (treatment assignment). The observed outcome `Yobs` is computed as \eqn{Y_{\mathrm{obs}} = Tr * Y1 + (1 - Tr) * Y0} for each observation.
#' @examples
#' sim <- get_data(n = 100, seed = 599)
#' dim(sim$data)    # should be 100 x (p + 3) columns (p features + S + Tr + Yobs)
#' head(sim$data$Yobs)  # observed outcomes
#' head(sim$Y)     # potential outcomes corresponding to those observations
#' @export
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

  # Generate sample inclusion and treatment assignment using offset seeds
  S_df  <- gen_S(X, seed = seed)
  T_lst <- gen_T(X, S_df, seed = seed)

  # Construct the observed outcome
  Tr_vec <- as.numeric(T_lst$Tr$Tr)  # extract treatment vector
  X$Yobs <- Tr_vec * Y$Y1 + (1 - Tr_vec) * Y$Y0

  # Combine X, S, Tr into one data frame
  full_data <- cbind(X, S = S_df$S, Tr = Tr_vec)

  return(list(data = full_data, Y = Y))
}


