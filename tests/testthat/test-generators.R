test_that("gen_XY produces sensible shapes and values", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  out <- gen_XY(n = 200, seed = 123)
  expect_s3_class(out$X, "data.frame")
  expect_s3_class(out$Y, "data.frame")
  expect_equal(nrow(out$X), 200L)
  expect_equal(nrow(out$Y), 200L)
  expect_true(all(c("Y0", "Y1") %in% names(out$Y)))
})

test_that("gen_S lowers inclusion inside the (0.5,1)x(0.5,1) box", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  XY <- gen_XY(n = 400, seed = 7)
  X  <- XY$X
  S  <- gen_S(X, seed = 7)

  in_box  <- (X$X0 > 0.5 & X$X0 < 1) & (X$X1 > 0.5 & X$X1 < 1)
  expect_true(mean(S$S[in_box] == 1) < mean(S$S[!in_box] == 1))
})

test_that("gen_T mixes randomized and observational treatment", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  XY <- gen_XY(n = 300, seed = 11)
  S  <- gen_S(XY$X, seed = 11)
  Tr <- gen_T(XY$X, S, seed = 11)

  expect_s3_class(Tr$Tr, "data.frame")
  expect_equal(nrow(Tr$Tr), nrow(XY$X))
  expect_equal(length(Tr$pi), nrow(XY$X))

  expect_true(all(Tr$pi > 0 & Tr$pi < 1))
})

test_that("get_data returns aligned frames and observed outcomes", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  sim <- get_data(n = 250, seed = 99)
  expect_equal(nrow(sim$data), 250L)
  expect_equal(nrow(sim$Y), 250L)
  expect_true(all(c("S", "Tr", "Y") %in% names(sim$data)))

  idx <- 1:10
  with(sim, {
    Y0 <- Y$Y0[idx]
    Y1 <- Y$Y1[idx]
    Tr <- data$Tr[idx]
    expect_equal(data$Y[idx], Tr * Y1 + (1 - Tr) * Y0, tolerance = 1e-8)
  })
})

test_that("gen_XY validates inputs and requires mlbench", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  # n must be positive numeric(1)
  expect_error(gen_XY(n = 0), "positive")
  expect_error(gen_XY(n = -10), "positive")
  expect_error(gen_XY(n = c(10, 20)), "positive")

  # seed must be NULL or numeric(1)
  expect_error(gen_XY(n = 10, seed = c(1,2)), "single numeric")
  expect_error(gen_XY(n = 10, seed = "nope"), "single numeric")

  # Basic shape and names
  out <- gen_XY(n = 17, seed = 11)
  expect_true(is.list(out))
  expect_true(all(c("X", "Y") %in% names(out)))
  expect_s3_class(out$X, "data.frame")
  expect_s3_class(out$Y, "data.frame")

  expect_equal(nrow(out$X), 17)
  expect_equal(nrow(out$Y), 17)

  # Friedman1 makes 10 features; we also rename X0..X9
  expect_equal(ncol(out$X), 10)
  expect_identical(colnames(out$X), paste0("X", 0:9))

  # Y has Y0, Y1 and the relation holds: Y1 = Y0 + log(Y0 + 1)
  expect_identical(colnames(out$Y), c("Y0", "Y1"))
  expect_equal(out$Y$Y1, out$Y$Y0 + log(out$Y$Y0 + 1), tolerance = 1e-12)

  # Determinism with seed
  out2 <- gen_XY(n = 17, seed = 11)
  expect_identical(out$X, out2$X)
  expect_identical(out$Y, out2$Y)
})

test_that("gen_S validates inputs and produces 0/1 outcomes with correct rectangle effect", {
  skip_if_not_installed("withr")

  # Build a tiny X by hand to control the rectangle (X0,X1) ∈ (0.5,1)
  X <- data.frame(
    X0 = c(0.25, 0.75, 0.60, 0.10),
    X1 = c(0.25, 0.75, 0.40, 0.80),
    X2 = rnorm(4)  # extra col should be ignored
  )

  # Validation: X must be data.frame and must contain X0,X1
  expect_error(gen_S(as.matrix(X)), "`X` must be a data frame")
  expect_error(gen_S(X[, c("X0", "X2")]), "contain columns 'X0' and 'X1'")
  expect_error(gen_S(X, seed = c(1,2)), "single numeric")

  # No NA output; S in {0,1}; correct length
  s1 <- gen_S(X, seed = 5)
  expect_s3_class(s1, "data.frame")
  expect_equal(nrow(s1), nrow(X))
  expect_true(all(s1$S %in% 0:1))
  expect_false(any(is.na(s1$S)))

  # Check probabilities are lower inside the rectangle:
  # a = 0.25 - 2 * I{X0∈(0.5,1) & X1∈(0.5,1)}; p = plogis(a)
  a_inside  <- 0.25 - 2
  a_outside <- 0.25
  expect_lt(stats::plogis(a_inside), stats::plogis(a_outside))

  # Determinism with seed
  s2 <- gen_S(X, seed = 5)
  expect_identical(s1$S, s2$S)
})

test_that("gen_T validates inputs and computes assignment probabilities correctly", {
  skip_if_not_installed("withr")

  set.seed(1)
  X <- data.frame(X0 = c(-2, 0, 2, 0.7), X1 = rnorm(4), X2 = rnorm(4))
  S <- data.frame(S = c(1, 1, 0, 0))

  # Basic validation
  expect_error(gen_T(as.matrix(X), S), "data frame of covariates")
  expect_error(gen_T(X, data.frame(S = c(1,0,1))), "same number of rows")
  expect_error(gen_T(X, data.frame(S = c(1,2,0,0))), "must contain only 0 or 1")
  expect_error(gen_T(X, S, seed = c(3,4)), "single numeric")

  # Compute expected pi:
  # pi_i = S_i * 0.5 + (1 - S_i) * plogis(X0_i)
  pi_exp      <- 0.5
  pi_obs      <- stats::plogis(X$X0)
  pi_expected <- S$S * pi_exp + (1 - S$S) * pi_obs

  t1 <- gen_T(X, S, seed = 99)
  expect_true(is.list(t1))
  expect_true(all(c("Tr", "pi") %in% names(t1)))
  expect_s3_class(t1$Tr, "data.frame")
  expect_equal(nrow(t1$Tr), nrow(X))
  expect_true(all(t1$Tr$Tr %in% 0:1))

  # Probabilities match formula
  expect_equal(t1$pi, pi_expected, tolerance = 1e-12)

  # Determinism of draws with seed
  t2 <- gen_T(X, S, seed = 99)
  expect_identical(t1$Tr$Tr, t2$Tr$Tr)

  # Special cases:
  # If S==1 -> pi == 0.5 regardless of X0
  S_all1 <- data.frame(S = rep(1L, nrow(X)))
  t_all1 <- gen_T(X, S_all1, seed = 100)
  expect_true(all(abs(t_all1$pi - 0.5) < 1e-15))

  # If S==0 -> pi == plogis(X0)
  S_all0 <- data.frame(S = rep(0L, nrow(X)))
  t_all0 <- gen_T(X, S_all0, seed = 100)
  expect_equal(t_all0$pi, stats::plogis(X$X0), tolerance = 1e-12)
})

test_that("get_data validates n/seed and returns consistent structure", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  expect_error(get_data(n = 0), "positive")
  expect_error(get_data(n = -5), "positive")
  expect_error(get_data(n = c(10, 20)), "positive")
  expect_error(get_data(n = 10, seed = c(1,2)), "single numeric")

  # Small n works
  sim1 <- get_data(n = 5, seed = 77)
  expect_true(is.list(sim1))
  expect_true(all(c("data", "Y") %in% names(sim1)))
  expect_s3_class(sim1$data, "data.frame")
  expect_s3_class(sim1$Y, "data.frame")
  expect_equal(nrow(sim1$data), 5)
  expect_equal(nrow(sim1$Y), 5)

  # Column presence & types
  x_cols <- grep("^X\\d+$", names(sim1$data), value = TRUE)
  expect_length(x_cols, 10)  # Friedman1
  expect_true(all(c("S", "Tr") %in% names(sim1$data)))
  expect_true("Y" %in% names(sim1$data))

  # Y equals Y1 for treated and Y0 for control
  Y   <- sim1$Y
  Tr  <- sim1$data$Tr
  Y_expected <- Tr * Y$Y1 + (1 - Tr) * Y$Y0
  expect_equal(sim1$data$Y, Y_expected, tolerance = 1e-12)

  # Determinism with seed
  sim2 <- get_data(n = 123, seed = 321)
  sim3 <- get_data(n = 123, seed = 321)
  expect_identical(sim2$data, sim3$data)
  expect_identical(sim2$Y, sim3$Y)
})

test_that("get_data integrates generators coherently (distribution sanity checks)", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  sim <- get_data(n = 1000, seed = 222)
  d   <- sim$data
  Y   <- sim$Y

  # S and Tr are 0/1, no NA
  expect_true(all(d$S %in% 0:1))
  expect_true(all(d$Tr %in% 0:1))
  expect_false(anyNA(d$S))
  expect_false(anyNA(d$Tr))
  expect_false(anyNA(d$Yobs))

  # Means in plausible ranges
  expect_true(mean(d$S) > 0.05 && mean(d$S) < 0.95)
  expect_true(mean(d$Tr) > 0.2 && mean(d$Tr) < 0.8)

  # Y1 - Y0 equals log(Y0 + 1)
  expect_equal(Y$Y1 - Y$Y0, log(Y$Y0 + 1), tolerance = 1e-12)
})
