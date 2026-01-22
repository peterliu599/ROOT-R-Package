test_that("characterizing_underrep integrates ROOT and returns leaf summaries", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("mlbench")
  skip_if_not_installed("MASS")
  skip_if_not_installed("withr")

  sim  <- get_data(n = 500, seed = 1234)
  full <- sim$data

  out <- characterizing_underrep(
    data                  = full,
    generalizability_path = TRUE,
    seed                  = 99,
    num_trees             = 5,
    top_k_trees           = TRUE,
    k                     = 3,
    feature_est           = "Ridge"
  )

  expect_s3_class(out, "characterizing_underrep")
  expect_true(all(c("root", "combined", "leaf_summary") %in% names(out)))
  expect_equal(nrow(out$combined), nrow(full))

  expect_silent(capture.output(summary(out)))

  if (!is.null(out$leaf_summary)) {
    expect_true(all(c("rule", "predicted_w", "n", "pct", "label") %in% names(out$leaf_summary)))
  } else {
    expect_true(is.null(out$root$f))
  }
})

test_that("characterizing_underrep validation errors on missing columns", {
  skip_if_not_installed("mlbench")
  skip_if_not_installed("withr")

  sim  <- get_data(n = 200, seed = 66)
  full <- sim$data

  # Remove Y column
  bad_data <- full[, setdiff(names(full), "Y")]

  expect_error(
    characterizing_underrep(
      data                  = bad_data,
      generalizability_path = TRUE
    ),
    "Missing: Y"
  )
})

test_that("characterizing_underrep works in general optimization mode", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("MASS")

  set.seed(123)
  n <- 100
  data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    v  = rnorm(n)
  )
  data$vsq <- data$v^2

  out <- characterizing_underrep(
    data                  = data,
    generalizability_path = FALSE,
    seed                  = 42,
    num_trees             = 3
  )

  expect_s3_class(out, "characterizing_underrep")
  expect_s3_class(out$root, "ROOT")
})

test_that("characterizing_underrep validates data is a data.frame", {
  expect_error(
    characterizing_underrep(data = matrix(1:10, nrow = 2)),
    "`data` must be a data.frame"
  )
})

test_that("characterizing_underrep validates generalizability_path argument", {
  df <- data.frame(X1 = 1:10, Y = rnorm(10), Tr = rep(0:1, 5), S = rep(1, 10))

  expect_error(
    characterizing_underrep(data = df, generalizability_path = "yes"),
    "`generalizability_path` must be TRUE or FALSE"
  )
})
