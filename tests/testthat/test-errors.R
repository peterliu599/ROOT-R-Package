test_that("ROOT enforces NA rules in generalizability_path mode", {
  # NA in Y among S == 1
  d <- data.frame(
    Y  = c(NA, 2, 3, 4, 5, 6),
    Tr = c(1, 0, 1, 0, 1, 0),
    S  = c(1, 1, 1, 1, 0, 0),
    X1 = rnorm(6),
    X2 = runif(6)
  )

  # This should error because Y has NA in the analysis sample
  expect_error(
    ROOT(d, generalizability_path = TRUE),
    regexp = "NA|missing|non-finite",
    ignore.case = TRUE
  )
})

test_that("ROOT enforces Tr must be 0/1 in generalizability_path mode", {
  d <- data.frame(
    Y  = rnorm(6),
    Tr = c(1, 0, 2, 0, 1, 0),  # Invalid: contains 2
    S  = c(1, 1, 1, 1, 0, 0),
    X1 = rnorm(6)
  )

  expect_error(
    ROOT(d, generalizability_path = TRUE),
    "Non 0/1 values found"
  )
})

test_that("ROOT enforces S must be 0/1 in generalizability_path mode", {
  d <- data.frame(
    Y  = rnorm(6),
    Tr = c(1, 0, 1, 0, 1, 0),
    S  = c(1, 1, 2, 1, 0, 0),  # Invalid: contains 2
    X1 = rnorm(6)
  )

  expect_error(
    ROOT(d, generalizability_path = TRUE),
    "Non 0/1 values found"
  )
})

test_that("characterize_tree() fails on non-binary w and length mismatch", {
  skip_if_not_installed("rpart")

  X <- data.frame(X1 = rnorm(5))

  # Non-binary w (more than two distinct values)
  expect_error(
    characterize_tree(X, w = c(0, 1, 2, 0, 1)),
    "exactly two classes"
  )

  # Length mismatch between X and w
  expect_error(
    characterize_tree(X, w = c(0, 1, 1)),
    "Length of `w` must equal the number of rows in `X`"
  )
})

test_that("ROOT requires at least one covariate column", {
  d <- data.frame(
    Y  = rnorm(10),
    Tr = sample(0:1, 10, TRUE),
    S  = sample(0:1, 10, TRUE)
  )

  # No covariate columns - only Y, Tr, S
  # The function should handle this gracefully (leaf-only trees)
  # or error if it requires covariates
  result <- tryCatch(

    ROOT(d, generalizability_path = TRUE, num_trees = 2),
    error = function(e) e
  )

  # Either succeeds or gives informative error

  expect_true(
    inherits(result, "ROOT") ||
      (inherits(result, "error") && grepl("covariate|feature", result$message, ignore.case = TRUE))
  )
})

test_that("compute_transport_scores requires covariates", {
  cts <- getFromNamespace("compute_transport_scores", "ROOT")

  # Data with only outcome, treatment, sample - no covariates
  d <- data.frame(
    Y  = rnorm(10),
    Tr = sample(0:1, 10, TRUE),
    S  = sample(0:1, 10, TRUE)
  )

  expect_error(
    cts(d, outcome = "Y", treatment = "Tr", sample = "S"),
    "need at least one covariate"
  )
})

test_that("compute_transport_scores requires S==1 rows", {
  cts <- getFromNamespace("compute_transport_scores", "ROOT")

  # Data with no trial rows (all S == 0)
  d <- data.frame(
    Y  = rnorm(10),
    Tr = sample(0:1, 10, TRUE),
    S  = rep(0, 10),
    X1 = rnorm(10)
  )

  expect_error(
    cts(d, outcome = "Y", treatment = "Tr", sample = "S"),
    "no S == 1"
  )
})

test_that(".check_no_na properly detects missing values", {
  check_fn <- getFromNamespace(".check_no_na", "ROOT")

  df_clean <- data.frame(a = 1:5, b = 6:10)
  expect_invisible(check_fn(df_clean, c("a", "b")))

  df_na <- data.frame(a = c(1, NA, 3), b = 4:6)
  expect_error(
    check_fn(df_na, c("a", "b")),
    "contains missing values"
  )
})
