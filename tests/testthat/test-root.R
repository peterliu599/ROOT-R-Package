test_that("ROOT runs in generalizability_path mode and returns structured outputs", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("rpart")
  skip_if_not_installed("withr")
  skip_if_not_installed("mlbench")

  sim <- get_data(n = 600, seed = 20)
  df  <- sim$data

  set.seed(42)
  out <- ROOT(
    data                  = df,
    generalizability_path = TRUE,
    seed                  = 99,
    num_trees             = 6,
    top_k_trees           = TRUE,
    k                     = 3,
    feature_est           = "Ridge",
    verbose               = TRUE
  )

  expect_s3_class(out, "ROOT")
  expect_true(all(c("D_rash", "D_forest", "w_forest", "rashomon_set", "f", "testing_data", "estimate") %in% names(out)))
  expect_equal(length(out$w_forest), 6L)
  expect_true(length(out$rashomon_set) >= 1L)
  expect_true(all(out$D_rash$w_opt %in% 0:1))

  est <- out$estimate
  expect_true(all(c("estimand_unweighted", "value_unweighted", "se_unweighted",
                    "estimand_weighted", "value_weighted", "se_weighted",
                    "n_analysis", "sum_w") %in% names(est)))
  expect_true(is.numeric(est$value_unweighted))
  expect_true(is.numeric(est$value_weighted) || is.na(est$value_weighted))
})

test_that("ROOT general optimization mode works when generalizability_path=FALSE", {
  skip_if_not_installed("mlbench")

  sim <- get_data(n = 300, seed = 21)
  df <- sim$data

  # Add v and vsq columns for general optimization
  df$v <- rnorm(nrow(df))
  df$vsq <- df$v^2

  out <- ROOT(
    data                  = df,
    generalizability_path = FALSE,
    num_trees             = 4,
    vote_threshold        = 0.6
  )
  expect_s3_class(out, "ROOT")
  expect_false(out$generalizability_path)
})

test_that("ROOT input validation catches errors", {
  skip_if_not_installed("mlbench")

  sim <- get_data(n = 100, seed = 5)
  df  <- sim$data

  expect_error(ROOT(as.matrix(df)), "`data` must be a data frame")
  expect_error(ROOT(df, generalizability_path = TRUE, num_trees = 0), "`num_trees` must be positive")
  expect_error(ROOT(df, generalizability_path = TRUE, vote_threshold = 1.1), "in \\(0, 1\\]")
})

mk_data_two_sample <- function(n = 200, p = 5, seed = 11) {
  set.seed(seed)
  X <- as.data.frame(matrix(rnorm(n * p), n, p))
  names(X) <- paste0("X", seq_len(p))
  Tr <- rbinom(n, 1, 0.5)
  S  <- rbinom(n, 1, 0.5)
  Y  <- rnorm(n)
  data.frame(X, Tr = Tr, S = S, Y = Y)
}

test_that("ROOT() validates arguments", {
  d <- mk_data_two_sample()

  expect_error(ROOT(as.matrix(d), generalizability_path = TRUE), "`data` must be a data frame")
  expect_error(ROOT(d, generalizability_path = TRUE, leaf_proba = -0.1), "between 0 and 1")
  expect_error(ROOT(d, generalizability_path = TRUE, num_trees = 0), "must be positive")
  expect_error(ROOT(d, generalizability_path = TRUE, vote_threshold = 0), "in \\(0, 1\\]")
  expect_error(ROOT(d, generalizability_path = TRUE, explore_proba = 2), "between 0 and 1")
  expect_error(ROOT(d, generalizability_path = TRUE, feature_est = 123), "must be \"Ridge\", \"GBM\", or a function")
  expect_error(ROOT(d, generalizability_path = TRUE, top_k_trees = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(ROOT(d, generalizability_path = TRUE, k = 0), "positive integer")
  expect_error(ROOT(d, generalizability_path = TRUE, cutoff = "nope"), "must be \"baseline\" or numeric")
  expect_error(ROOT(d, generalizability_path = TRUE, verbose = 1), "TRUE or FALSE")
  expect_error(ROOT(d, generalizability_path = TRUE, global_objective_fn = 1), "must be a function")
})

test_that("ROOT generalizability_path requires Y, Tr, S columns", {
  d <- mk_data_two_sample()
  d_no_Y <- d[, setdiff(names(d), "Y")]
  d_no_Tr <- d[, setdiff(names(d), "Tr")]
  d_no_S <- d[, setdiff(names(d), "S")]

  expect_error(ROOT(d_no_Y, generalizability_path = TRUE), "Missing: Y")
  expect_error(ROOT(d_no_Tr, generalizability_path = TRUE), "Missing: Tr")
  expect_error(ROOT(d_no_S, generalizability_path = TRUE), "Missing: S")
})

test_that("ROOT feature_est = Ridge / GBM / custom; bad custom errors", {
  d <- mk_data_two_sample(n = 220, p = 6)

  skip_if_not_installed("MASS")
  skip_if_not_installed("gbm")
  skip_if_not_installed("mlbench")

  set.seed(1)
  rR <- ROOT(d, generalizability_path = TRUE, num_trees = 2, feature_est = "Ridge")
  expect_s3_class(rR, "ROOT")

  set.seed(2)
  rG <- ROOT(d, generalizability_path = TRUE, num_trees = 2, feature_est = "GBM")
  expect_s3_class(rG, "ROOT")

  ok_imp <- function(X, y, ...) {
    setNames(abs(colSums(X^2)) + 1e-6, colnames(X))
  }
  set.seed(3)
  rC <- ROOT(d, generalizability_path = TRUE, num_trees = 2, feature_est = ok_imp)
  expect_s3_class(rC, "ROOT")

  bad1 <- function(X, y, ...) {
    as.numeric(abs(colSums(X)))
  }
  expect_error(
    ROOT(d, generalizability_path = TRUE, num_trees = 1, feature_est = bad1),
    "named numeric vector"
  )

  bad2 <- function(X, y, ...) {
    nm <- colnames(X)
    nm[length(nm)] <- "WRONG"
    setNames(abs(colSums(X)), nm)
  }
  expect_error(
    ROOT(d, generalizability_path = TRUE, num_trees = 1, feature_est = bad2),
    "Importance missing for some X_df columns"
  )

  bad3 <- function(X, y, ...) {
    setNames(c(-1, rep(1, ncol(X) - 1)), colnames(X))
  }
  expect_error(
    ROOT(d, generalizability_path = TRUE, num_trees = 1, feature_est = bad3),
    "must be non-negative"
  )
})

test_that("Rashomon selection: top-k, baseline cutoff, and empty set warning", {
  skip_if_not_installed("mlbench")

  d <- mk_data_two_sample(n = 260, p = 6)

  set.seed(11)
  r_topk <- ROOT(d, generalizability_path = TRUE, num_trees = 5, top_k_trees = TRUE, k = 3)
  expect_length(r_topk$rashomon_set, 3)

  set.seed(12)
  r_base <- ROOT(d, generalizability_path = TRUE, num_trees = 5, top_k_trees = FALSE, cutoff = "baseline")
  expect_true(length(r_base$rashomon_set) >= 0)

  set.seed(13)
  expect_warning(
    r_empty <- ROOT(d, generalizability_path = TRUE, num_trees = 3, top_k_trees = FALSE, cutoff = -1e9),
    "No trees selected"
  )
  expect_length(r_empty$rashomon_set, 0)
  expect_true("w_opt" %in% names(r_empty$D_rash))
})

test_that("ROOT estimate fields populated in generalizability_path mode", {
  skip_if_not_installed("mlbench")

  d <- mk_data_two_sample(n = 280, p = 6)

  set.seed(21)
  r <- ROOT(d, generalizability_path = TRUE, num_trees = 4)
  e <- r$estimate
  expect_true(all(c("estimand_unweighted", "value_unweighted", "se_unweighted",
                    "estimand_weighted", "value_weighted", "se_weighted",
                    "se_weighted_note", "n_analysis", "sum_w") %in% names(e)))
  expect_true(is.numeric(e$se_weighted) || is.na(e$se_weighted))
})

test_that("ROOT seed yields reproducible w_forest & estimates", {
  skip_if_not_installed("mlbench")

  d <- mk_data_two_sample(n = 200, p = 5)

  r1 <- ROOT(d, generalizability_path = TRUE, num_trees = 3, seed = 123, verbose = FALSE)
  r2 <- ROOT(d, generalizability_path = TRUE, num_trees = 3, seed = 123, verbose = FALSE)

  expect_equal(
    lapply(r1$w_forest, `[[`, "local objective"),
    lapply(r2$w_forest, `[[`, "local objective")
  )
  expect_equal(r1$estimate$value_unweighted, r2$estimate$value_unweighted, tolerance = 1e-12)
  expect_equal(r1$estimate$value_weighted, r2$estimate$value_weighted, tolerance = 1e-12)
})

test_that("ROOT runs with custom objective and zero-sum feature importances", {
  set.seed(10)
  n  <- 180
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  S  <- rbinom(n, 1, plogis(0.4 * X1 - 0.2 * X2))
  Tr <- rbinom(n, 1, plogis(0.3 + 0.7 * X1 - 0.2 * X2))
  Y  <- 0.5 + 1.1 * Tr + 0.5 * X1 - 0.3 * X2 + rnorm(n, 0.7)

  dat <- data.frame(Y = Y, Tr = Tr, S = S, X1 = X1, X2 = X2)

  zero_imp <- function(X, y, ...) stats::setNames(rep(0, ncol(X)), colnames(X))
  my_obj   <- function(D) mean(D$vsq, na.rm = TRUE) + 1e-8

  fit <- suppressWarnings(ROOT(
    data                  = dat,
    generalizability_path = TRUE,
    seed                  = 123,
    num_trees             = 5,
    top_k_trees           = TRUE,
    k                     = 5,
    vote_threshold        = 2/3,
    feature_est           = zero_imp,
    verbose               = FALSE,
    global_objective_fn   = my_obj
  ))
  expect_s3_class(fit, "ROOT")
})

test_that("ROOT argument guards", {
  set.seed(1)
  df <- data.frame(
    Y  = rnorm(40),
    Tr = sample(0:1, 40, TRUE),
    S  = sample(0:1, 40, TRUE),
    X  = rnorm(40)
  )

  expect_error(
    ROOT(as.matrix(df), generalizability_path = TRUE),
    "`data` must be a data frame"
  )

  expect_error(
    ROOT(df, generalizability_path = TRUE, leaf_proba = -0.1),
    "`leaf_proba` must be between 0 and 1"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, seed = c(1, 2)),
    "`seed` must be NULL or a single numeric value"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, num_trees = 0),
    "`num_trees` must be positive"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, vote_threshold = 0),
    "`vote_threshold` must be in \\(0, 1\\]"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, explore_proba = 2),
    "`explore_proba` must be between 0 and 1"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, feature_est = 123),
    "`feature_est` must be \"Ridge\", \"GBM\", or a function"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, k = 0),
    "`k` must be a positive integer"
  )
  expect_error(
    ROOT(df, generalizability_path = TRUE, cutoff = "not-baseline"),
    "`cutoff` must be \"baseline\" or numeric"
  )
})

test_that("Rashomon top-k warning when k > num_trees", {
  set.seed(11)
  n  <- 60
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  S  <- rbinom(n, 1, 0.5)
  Tr <- rbinom(n, 1, plogis(0.6 * X1))
  Y  <- 3 + 0 * Tr + 0 * X1 + rnorm(n, sd = 1e-8)

  dat <- data.frame(Y = Y, Tr = Tr, S = S, X1 = X1, X2 = X2)

  expect_warning(
    fit1 <- ROOT(dat, generalizability_path = TRUE, seed = 321, num_trees = 3,
                 feature_est = "Ridge", top_k_trees = TRUE, k = 999),
    "k > num_trees; using k = num_trees"
  )
  expect_s3_class(fit1, "ROOT")
})

test_that("ROOT informs when no summary tree is available (single-class w_opt)", {
  set.seed(13)
  n  <- 70
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  S  <- rbinom(n, 1, 0.5)
  Tr <- rbinom(n, 1, plogis(X1 - 0.2 * X2))
  Y  <- 0.2 + Tr + X1 + rnorm(n)

  dat <- data.frame(Y = Y, Tr = Tr, S = S, X1 = X1, X2 = X2)

  expect_message(
    fit <- ROOT(dat, generalizability_path = TRUE, seed = 99, num_trees = 4,
                vote_threshold = 1, verbose = FALSE),
    "No summary tree available to plot"
  )
  expect_s3_class(fit, "ROOT")
})

test_that("ROOT coerce01 handles bad sample values", {
  base <- data.frame(
    Y  = rnorm(8),
    Tr = sample(0:1, 8, TRUE),
    S  = c("yes", "no", "maybe", "no", "yes", "unknown", "0", "1"),
    X  = rnorm(8)
  )

  expect_error(ROOT(base, generalizability_path = TRUE), "Non 0/1 values found")
})
