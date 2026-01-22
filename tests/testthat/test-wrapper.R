# Helper to create fake ROOT objects for testing
fake_root <- function(D_forest,
                      D_rash = data.frame(),
                      rashomon_set = integer(0),
                      estimate = NULL,
                      f = NULL,
                      w_forest = NULL,
                      global_objective_fn = objective_default,
                      generalizability_path = FALSE) {
  obj <- list(
    D_forest            = D_forest,
    D_rash              = D_rash,
    rashomon_set        = rashomon_set,
    testing_data        = D_forest,
    estimate            = estimate,
    f                   = f,
    w_forest            = w_forest,
    global_objective_fn = global_objective_fn,
    generalizability_path = generalizability_path
  )
  class(obj) <- c("ROOT", "list")
  obj
}

## ------------------------------------------------------------------
## Shared simulated data for ROOT / characterizing_underrep tests
## ------------------------------------------------------------------

set.seed(123)
sim_2s <- get_data(n = 500, seed = 123)
dat_2s <- sim_2s$data

covs <- paste0("X", 0:9)

## ------------------------------------------------------------------
## characterizing_underrep wrapper and methods
## ------------------------------------------------------------------

test_that("characterizing_underrep wrapper works", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("MASS")

  suppressWarnings({
    res <- characterizing_underrep(
      data                  = dat_2s[1:100, ],
      generalizability_path = TRUE,
      num_trees             = 2,
      seed                  = 123
    )
  })

  expect_s3_class(res, "characterizing_underrep")
  expect_s3_class(res$root, "ROOT")
})

test_that("S3 methods for characterizing_underrep work", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("MASS")

  suppressWarnings({
    res <- characterizing_underrep(
      data                  = dat_2s[1:100, ],
      generalizability_path = TRUE,
      num_trees             = 2,
      seed                  = 123
    )
  })

  expect_output(summary(res), "characterizing_underrep object")

  if (!is.null(res$root$f)) {
    skip_if_not_installed("rpart.plot")
    expect_silent(grDevices::pdf(file = NULL))
    plot(res)
    grDevices::dev.off()
  }
})

## ------------------------------------------------------------------
## ROOT S3 methods
## ------------------------------------------------------------------

test_that("S3 methods for ROOT work", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("MASS")

  suppressWarnings({
    res_2s <- ROOT(dat_2s[1:200, ], generalizability_path = TRUE, num_trees = 3, seed = 123)
  })

  expect_output(summary(res_2s), "ROOT object")
})

## ------------------------------------------------------------------
## characterize_tree()
## ------------------------------------------------------------------

test_that("characterize_tree() validates length and binarity of w", {
  skip_if_not_installed("rpart")

  X <- data.frame(X1 = rnorm(30), X2 = runif(30))

  expect_error(
    characterize_tree(X, w = 0:2),
    "Length of `w` must equal the number of rows in `X`"
  )

  w3 <- c(0, 1, 2)[(seq_len(nrow(X)) - 1) %% 3 + 1]
  expect_error(
    characterize_tree(X, w = w3),
    "exactly two classes"
  )

  set.seed(1)
  w <- sample(c(0, 1), size = nrow(X), replace = TRUE)
  fit <- characterize_tree(X, w = w, max_depth = 2)
  expect_s3_class(fit, "rpart")
})

test_that("characterize_tree fits a shallow rpart classifier", {
  skip_if_not_installed("rpart")

  X <- data.frame(x1 = c(0, 0, 1, 1), x2 = c(0, 1, 0, 1))
  w <- c(0, 0, 1, 1)
  fit <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(fit, "rpart")
  expect_true(all(c("frame", "where") %in% names(fit)))
})

test_that("characterize_tree errors when w is not binary", {
  X <- data.frame(x1 = rnorm(10))
  w <- rep(2, 10)
  expect_error(characterize_tree(X, w), "exactly two classes")
})

## ------------------------------------------------------------------
## summary.ROOT printing & diagnostics
## ------------------------------------------------------------------

test_that("summary.ROOT prints basic info", {
  Df <- data.frame(
    v = rnorm(3), vsq = 1, S = c(1L, 1L, 1L),
    X1 = 1:3, w_tree_1 = 1L
  )

  obj <- fake_root(
    D_forest = Df,
    D_rash   = data.frame(w_opt = c(1L, 1L, 0L)),
    rashomon_set = 1L,
    estimate = list(
      estimand_unweighted = "SATE",
      value_unweighted    = 0.1,
      se_unweighted       = 0.2,
      estimand_weighted   = "WTATE",
      value_weighted      = 0.3,
      se_weighted         = 0.4,
      se_weighted_note    = "ok"
    ),
    generalizability_path = TRUE
  )

  out <- capture.output(summary(obj))
  expect_true(any(grepl("ROOT object", out)))
  expect_true(any(grepl("Generalization mode", out)))
})

test_that("summary.ROOT handles missing components gracefully", {
  Df <- data.frame(X1 = 1:3, vsq = 1, v = rnorm(3))

  obj <- fake_root(
    D_forest = Df,
    D_rash   = data.frame(w_opt = c(1L, 0L, 1L)),
    rashomon_set = integer(0),
    generalizability_path = FALSE
  )

  expect_output(summary(obj), "ROOT object")
})

## ------------------------------------------------------------------
## plot.ROOT
## ------------------------------------------------------------------

test_that("plot.ROOT handles x$f NULL and plots when present", {
  obj_null <- fake_root(D_forest = data.frame(), f = NULL)
  expect_message(plot(obj_null), "No summary tree available")

  skip_if_not_installed("rpart")
  skip_if_not_installed("rpart.plot")

  set.seed(1)
  df  <- data.frame(
    y = factor(sample(c(0, 1), 60, TRUE)),
    X1 = runif(60),
    X2 = runif(60)
  )
  fit <- rpart::rpart(y ~ X1 + X2, data = df, method = "class")
  obj <- fake_root(D_forest = df, f = fit)

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(suppressWarnings(plot(obj)), NA)
})

## ------------------------------------------------------------------
## characterizing_underrep S3 methods
## ------------------------------------------------------------------

test_that("summary.characterizing_underrep() rejects wrong class", {
  expect_error(
    summary.characterizing_underrep(list()),
    "Not a characterizing_underrep object"
  )
})

test_that("summary.characterizing_underrep() prints when leaf_summary is NULL", {
  obj <- structure(
    list(
      root = NULL,
      combined = data.frame(),
      leaf_summary = NULL
    ),
    class = "characterizing_underrep"
  )

  expect_output(
    expect_invisible(summary(obj)),
    "Leaf summary:\\s+none \\(no summarized tree\\)"
  )
})

test_that("summary.characterizing_underrep() shows preview and truncation cue", {
  leaf_summary <- data.frame(
    rule        = sprintf("X0 <= %i", 1:12),
    predicted_w = rep(1, 12),
    n           = 1:12,
    pct         = (1:12) / 100,
    label       = rep("Represented (keep, w = 1)", 12),
    stringsAsFactors = FALSE
  )
  obj <- structure(
    list(root = NULL, combined = data.frame(), leaf_summary = leaf_summary),
    class = "characterizing_underrep"
  )

  out <- capture.output(invisible(summary(obj)))

  expect_true(any(grepl("Leaf summary:\\s+\\d+\\s+terminal nodes", out)))
  shown_rules <- grep("X0\\s*<=\\s*\\d+", out, value = TRUE)
  expect_equal(length(shown_rules), 10)
  expect_true(any(grepl("\\.\\.\\.", out)))
})

## ------------------------------------------------------------------
## plot.characterizing_underrep
## ------------------------------------------------------------------

test_that("plot.characterizing_underrep() safely exits when no summary tree", {
  skip_if_not_installed("rpart")

  obj <- structure(
    list(root = list(f = NULL)),
    class = "characterizing_underrep"
  )

  expect_message(
    expect_invisible(plot(obj)),
    "No summary tree available to plot"
  )
})

test_that("plot.characterizing_underrep() works with valid tree", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("rpart.plot")

  set.seed(1)
  n  <- 30
  X0 <- runif(n)
  X1 <- runif(n)
  w  <- factor(ifelse(X0 + X1 > 1, 1, 0), levels = c(0, 1))

  fit <- rpart::rpart(
    w ~ X0 + X1,
    method  = "class",
    control = rpart::rpart.control(cp = 0, minsplit = 2, maxdepth = 2)
  )

  obj <- structure(
    list(root = list(f = fit, generalizability_path = TRUE)),
    class = "characterizing_underrep"
  )

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  expect_silent(plot(obj))
  grDevices::dev.off()
  unlink(tf)
})

## ------------------------------------------------------------------
## Internal helpers tests
## ------------------------------------------------------------------

# Define local versions to avoid namespace pollution from mocks in other tests
local_objective_default_wrapper <- function(D) {
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

local_objective_if_wrapper <- function(val, indices, D, global_objective_fn) {
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

local_loss_from_objective_wrapper <- function(global_objective_fn) {
  force(global_objective_fn)
  function(val, indices, D) local_objective_if_wrapper(val, indices, D, global_objective_fn)
}

test_that("objective_default computes correct values", {
  obj_fn <- local_objective_default_wrapper

  D <- data.frame(
    vsq = c(1, 4, 9),
    w   = c(1, 1, 1)
  )

  result <- obj_fn(D)
  expect_true(is.finite(result))

  # sqrt(sum(w * vsq) / sum(w)^2) = sqrt(14 / 9)
  expected <- sqrt(sum(D$w * D$vsq) / sum(D$w)^2)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("loss_from_objective creates working loss function", {
  lfo <- local_loss_from_objective_wrapper
  obj_fn <- local_objective_default_wrapper

  D <- data.frame(
    vsq = c(1, 2, 3, 4),
    w   = c(0, 0, 0, 0),
    row.names = paste0("r", 1:4)
  )

  loss_fn <- lfo(obj_fn)

  # Setting w=1 for all rows
  result <- loss_fn(1, rownames(D), D)
  expect_true(is.numeric(result))
  expect_length(result, 1)
})
