test_that("ROOT + Characterization run end-to-end on Box DGP (100-D)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("rpart")
  skip_if_not_installed("withr")

  set.seed(20251028)

  expit <- function(z) 1 / (1 + exp(-z))

  simulate_box_dgp <- function(n = 600L, p = 100L,
                               s_shift      = 0.10,
                               s_box        = -1.0,
                               s_noise_sd   = 0.25,
                               s_label_flip = 0.00) {
    stopifnot(n > 50, p >= 5)

    X <- as.data.frame(matrix(runif(n * p), nrow = n, ncol = p))
    names(X) <- paste0("X", 0:(p - 1))

    in_box  <- (X$X0 > 0.5 & X$X0 < 1) & (X$X1 > 0.5 & X$X1 < 1)
    linpred <- s_shift + s_box * as.integer(in_box) + rnorm(n, 0, s_noise_sd)
    p_s     <- expit(linpred)
    S       <- rbinom(n, 1, p_s)

    if (s_label_flip > 0) {
      flip <- rbinom(n, 1, s_label_flip)
      S[flip == 1] <- 1L - S[flip == 1]
    }

    Tr <- rep(NA_integer_, n)
    Tr[S == 1] <- rbinom(sum(S == 1), 1, 0.5)

    eps <- rnorm(n, 0, 1)
    Y0  <- 10 * sin(pi * X$X0 * X$X1) +
      20 * (X$X2 - 0.5)^2 +
      10 * X$X3 +
      5 * X$X4 + eps

    tau <- log(pmax(Y0, -0.999) + 1)
    Y1  <- Y0 + tau

    Y <- rep(NA_real_, n)
    idx  <- which(S == 1)
    Y[idx] <- Tr[idx] * Y1[idx] + (1 - Tr[idx]) * Y0[idx]

    list(
      data_full = cbind(X, S = S, Tr = Tr, Y = Y),
      X        = X,
      S        = S,
      Tr       = Tr,
      Y0       = Y0,
      Y1       = Y1,
      Y        = Y,
      in_box   = in_box
    )
  }

  sim <- simulate_box_dgp(n = 600, p = 100)

  inside_rate  <- mean(sim$S[sim$in_box] == 1)
  outside_rate <- mean(sim$S[!sim$in_box] == 1)
  expect_lt(inside_rate, outside_rate)

  root_out <- ROOT(
    data                  = sim$data_full,
    generalizability_path = TRUE,
    seed                  = 99,
    num_trees             = 5,
    top_k_trees           = TRUE,
    k                     = 3,
    feature_est           = "Ridge",
    verbose               = FALSE
  )

  expect_s3_class(root_out, "ROOT")
  expect_true(isTRUE(root_out$generalizability_path))

  expect_type(root_out, "list")
  expect_true(all(c(
    "D_rash", "D_forest", "w_forest",
    "rashomon_set", "f", "testing_data", "estimate"
  ) %in% names(root_out)))

  expect_equal(length(root_out$w_forest), 5L)
  expect_true(length(root_out$rashomon_set) >= 1L &&
                length(root_out$rashomon_set) <= 5L)

  expect_equal(nrow(root_out$D_rash), nrow(root_out$testing_data))
  expect_true("w_opt" %in% names(root_out$D_rash))
  expect_true(all(root_out$D_rash$w_opt %in% 0:1))

  if (!is.null(root_out$f)) {
    expect_s3_class(root_out$f, "rpart")
  }

  char_out <- characterizing_underrep(
    data                  = sim$data_full,
    generalizability_path = TRUE,
    seed                  = 99,
    num_trees             = 5,
    top_k_trees           = TRUE,
    k                     = 3,
    feature_est           = "Ridge",
    verbose               = FALSE
  )

  expect_type(char_out, "list")
  expect_s3_class(char_out, "characterizing_underrep")
  expect_true(all(c("root", "combined", "leaf_summary") %in% names(char_out)))
  expect_equal(nrow(char_out$combined), nrow(sim$data_full))

  if (!is.null(char_out$leaf_summary)) {
    expect_true(all(c("rule", "predicted_w", "n", "pct", "label")
                    %in% names(char_out$leaf_summary)))
  } else {
    expect_true(is.null(char_out$root$f))
  }
})

test_that("ROOT works with small dataset", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("rpart")

  set.seed(42)
  n <- 100
  df <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n),
    Tr = sample(0:1, n, TRUE),
    S  = sample(0:1, n, TRUE),
    Y  = rnorm(n)
  )

  result <- ROOT(
    data                  = df,
    generalizability_path = TRUE,
    num_trees             = 3,
    seed                  = 123
  )

  expect_s3_class(result, "ROOT")
  expect_true("estimate" %in% names(result))
  expect_true("D_rash" %in% names(result))
})
