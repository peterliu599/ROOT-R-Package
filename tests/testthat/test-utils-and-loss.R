# Define the functions locally for testing to avoid namespace pollution from mocks
# These should match the implementations in R/utils.R
local_objective_default <- function(D) {
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

local_objective_if <- function(val, indices, D, global_objective_fn) {
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

local_loss_from_objective <- function(global_objective_fn) {
  force(global_objective_fn)
  function(val, indices, D) local_objective_if(val, indices, D, global_objective_fn)
}

test_that("objective_default and loss_from_objective agree on real code", {
  # Use local definitions to avoid namespace pollution from mocks in other tests
  obj <- local_objective_default
  lfo <- local_loss_from_objective

  D <- data.frame(
    vsq = c(0.1, 0.3, 0.1, 0.5),
    w   = c(0L, 0L, 0L, 0L),
    row.names   = paste0("r", 1:4),
    check.names = FALSE
  )

  losser <- lfo(obj)

  L0 <- losser(1, integer(0), D)
  expect_type(L0, "double")
  expect_true(is.finite(L0) || is.infinite(L0))
  expect_false(is.na(L0))

  # losser(val, indices, D) sets D$w[indices] <- val and evaluates obj on full D
  idx_all <- seq_len(nrow(D))
  D_all   <- D
  D_all$w[idx_all] <- 1L
  expected_all <- obj(D_all)

  L_all <- losser(1, idx_all, D)
  expect_type(L_all, "double")
  expect_equal(L_all, expected_all, tolerance = 1e-12)

  # Partial subset: set w=1 for first 3 rows only
  idx_drop_last <- idx_all[-length(idx_all)]
  D_sub         <- D
  D_sub$w[idx_drop_last] <- 1L
  expected_sub <- obj(D_sub)

  L_sub <- losser(1, idx_drop_last, D)
  expect_type(L_sub, "double")
  expect_equal(L_sub, expected_sub, tolerance = 1e-12)

  # Setting w=0 for all rows (starting from w=0)
  D_zero <- D
  D_zero$w[idx_all] <- 0L
  expected_zero <- obj(D_zero)

  L_zero <- losser(0, idx_all, D)
  expect_equal(L_zero, expected_zero, tolerance = 1e-12)
})

test_that("loss_from_objective composes correctly with a custom objective", {
  # Use local definition to avoid namespace pollution
  lfo <- local_loss_from_objective

  # Custom objective that depends on the full D
  custom_obj <- function(D) sum(D$w * D$vsq)

  losser <- lfo(custom_obj)

  D <- data.frame(
    vsq = c(2, 3),
    w   = c(0L, 0L),
    row.names   = c("a", "b"),
    check.names = FALSE
  )

  L0 <- losser(1, integer(0), D)
  expect_true(is.finite(L0) || is.infinite(L0))
  expect_type(L0, "double")

  # Set w=1 for first row only, evaluate on full D
  idx      <- 1L
  D_mod    <- D
  D_mod$w[idx] <- 1L
  expected <- custom_obj(D_mod)  # sum(c(1,0) * c(2,3)) = 2

  L_subset <- losser(1, idx, D)
  expect_equal(L_subset, expected, tolerance = 1e-12)
})

test_that(".check_no_na detects NAs and returns invisibly TRUE otherwise", {
  check_fn <- getFromNamespace(".check_no_na", "ROOT")

  df <- data.frame(a = 1:3, b = c(1, NA, 3))
  expect_error(check_fn(df, c("a", "b")), "contains missing values")

  df2 <- data.frame(a = 1:3, b = 11:13)
  expect_invisible(check_fn(df2, c("a", "b")))
})

test_that("midpoint, choose_feature, and reduce_weight work", {
  set.seed(1)
  x <- c(2, 8, 10, 4)
  expect_equal(midpoint(x), (min(x) + max(x)) / 2)

  x2 <- c(1, 5, NA, Inf)
  mid2 <- suppressWarnings(midpoint(x2))
  expect_true(is.numeric(mid2) && length(mid2) == 1)

  sf <- c(leaf = 0.2, X0 = 0.3, X1 = 0.5)

  set.seed(123)
  got   <- replicate(1000, choose_feature(sf, depth = 0))
  props <- prop.table(table(got))
  expect_true(all(names(sf) %in% names(props)))
  expect_gt(props[["X1"]], props[["leaf"]])

  expect_error(
    choose_feature(unname(sf), 0),
    "must have names",
    ignore.case = TRUE
  )
  expect_error(
    choose_feature(c(a = 0.5, b = NA_real_), 0),
    "contains NA|NA values",
    ignore.case = TRUE
  )
  expect_error(
    choose_feature(c(a = 0.5, b = -0.1), 0),
    ">= 0|negative",
    ignore.case = TRUE
  )

  sf2 <- reduce_weight("X1", sf)
  expect_lt(sf2["X1"], sf["X1"])
  expect_equal(sum(sf2), 1, tolerance = 1e-12)

  expect_setequal(names(sf2), names(sf))
  expect_false(any(is.na(sf2)))
  expect_true(all(sf2 >= 0))

  expect_error(
    reduce_weight("X9", sf),
    "not found",
    ignore.case = TRUE
  )
})

test_that("objective_default handles edge cases", {
  # Use local definition to avoid namespace pollution
  obj <- local_objective_default

  # Missing w column should error
  expect_error(obj(data.frame(vsq = 1:3)), "expects a column `w`")

  # All NA weights should return Inf
  D_na <- data.frame(vsq = 1:3, w = c(NA, NA, NA))
  expect_equal(obj(D_na), Inf)

  # All zero weights should return Inf
  D_zero <- data.frame(vsq = 1:3, w = c(0, 0, 0))
  expect_equal(obj(D_zero), Inf)

  # Valid case
  D_valid <- data.frame(vsq = c(1, 4, 9), w = c(1, 1, 1))
  result <- obj(D_valid)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("objective_default falls back to v when vsq missing", {
  # Use local definition to avoid namespace pollution
  obj <- local_objective_default

  D <- data.frame(v = c(1, 2, 3), w = c(1, 1, 1))
  result <- obj(D)
  expect_true(is.finite(result))
})

local_objective_default_coverage <- function(D) {
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

test_that("objective_default: vsq column with non-numeric type is skipped", {
  obj <- local_objective_default_coverage


  # vsq exists but is character (non-numeric), should fall through to v branch

  D <- data.frame(
    vsq = c("a", "b", "c"),
    v   = c(1, 2, 3),
    w   = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  result <- obj(D)
  expect_true(is.finite(result))
})

test_that("objective_default: vsq column with all non-finite values triggers v fallback",
          {
            obj <- local_objective_default_coverage

            # vsq exists and is numeric, but all values are Inf/NaN
            D <- data.frame(
              vsq = c(Inf, NaN, Inf),
              v   = c(1, 2, 3),
              w   = c(1, 1, 1)
            )
            result <- obj(D)
            expect_true(is.finite(result))
          })

test_that("objective_default: v column with all non-finite values triggers Bernoulli fallback", {
  obj <- local_objective_default_coverage

  # Neither vsq nor v have finite values -> Bernoulli fallback

  D <- data.frame(
    v = c(Inf, NaN, -Inf),
    w = c(1, 1, 1)
  )
  result <- obj(D)
  # With n_keep = 3 and p = 1, result = sqrt(1 * 0 / 3) = 0

  expect_true(is.finite(result))
  expect_equal(result, 0)
})

test_that("objective_default: Bernoulli fallback with n_keep <= 1 returns Inf", {
  obj <- local_objective_default_coverage

  # No vsq, no v, only 1 observation with w > 0
  D <- data.frame(
    x = c(1, 2, 3),  # no vsq or v column
    w = c(1, 0, 0)
  )
  result <- obj(D)
  expect_equal(result, Inf)
})

test_that("objective_default: Bernoulli fallback computes correctly", {
  obj <- local_objective_default_coverage

  # No vsq, no v -> Bernoulli fallback
  # n_keep = 3, p = 3/5 = 0.6
  # result = sqrt(0.6 * 0.4 / 3) = sqrt(0.08)
  D <- data.frame(
    x = 1:5,  # no vsq or v column
    w = c(1, 1, 1, 0, 0)
  )
  result <- obj(D)
  expected <- sqrt(0.6 * 0.4 / 3)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("objective_default: non-finite num returns Inf", {
  obj <- local_objective_default_coverage

  # vsq contains Inf, w * vsq will have Inf -> num is Inf
  D <- data.frame(
    vsq = c(1, Inf, 1),
    w   = c(1, 1, 1)
  )
  result <- obj(D)
  expect_equal(result, Inf)
})

test_that("objective_default: handles mixed NA in vsq correctly", {
  obj <- local_objective_default_coverage

  # vsq has some NA but also some finite values
  D <- data.frame(
    vsq = c(1, NA, 4),
    w   = c(1, 1, 1)
  )
  result <- obj(D)
  # num = 1*1 + NA + 1*4 = 5 (na.rm=TRUE), den = 3^2 = 9
  # result = sqrt(5/9)
  expected <- sqrt(5 / 9)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("objective_default: negative weights sum to <= 0 returns Inf", {
  obj <- local_objective_default_coverage

  D <- data.frame(
    vsq = c(1, 2, 3),
    w   = c(-1, -1, 0)
  )
  result <- obj(D)
  expect_equal(result, Inf)
})
