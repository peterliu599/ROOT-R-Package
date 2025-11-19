test_that("ROOT + Characterization run end-to-end on Box DGP (100-D)", {
  skip_on_cran()
  skip_if_not_installed("MASS")   # lm.ridge for feature_est = "Ridge"
  skip_if_not_installed("rpart")  # characterize_tree()
  skip_if_not_installed("withr")  # ROOT uses withr::with_seed when seed is set

  set.seed(20251028)

  expit <- function(z) 1 / (1 + exp(-z))

  simulate_box_dgp <- function(n = 600L, p = 100L,
                               s_shift = 0.10,    # baseline logit intercept
                               s_box   = -1.0,    # box penalty (weaker than -2)
                               s_noise_sd = 0.25, # Gaussian noise in logit for S
                               s_label_flip = 0.00  # optional tiny label noise
  ) {
    stopifnot(n > 50, p >= 5)
    X <- as.data.frame(matrix(runif(n * p), nrow = n, ncol = p))
    names(X) <- paste0("X", 0:(p - 1))

    in_box  <- (X$X0 > 0.5 & X$X0 < 1) & (X$X1 > 0.5 & X$X1 < 1)
    linpred <- s_shift + s_box * as.integer(in_box) + rnorm(n, 0, s_noise_sd)
    p_s     <- expit(linpred)
    S       <- rbinom(n, 1, p_s)

    # optional: tiny symmetric label noise for S to further break separation
    if (s_label_flip > 0) {
      flip <- rbinom(n, 1, s_label_flip)
      S[flip == 1] <- 1L - S[flip == 1]
    }

    Tr <- rep(NA_integer_, n)
    Tr[S == 1] <- rbinom(sum(S == 1), 1, 0.5)

    eps <- rnorm(n, 0, 1)
    Y0  <- 10 * sin(pi * X$X0 * X$X1) + 20 * (X$X2 - 0.5)^2 + 10 * X$X3 + 5 * X$X4 + eps
    tau <- log(pmax(Y0, -0.999) + 1)
    Y1  <- Y0 + tau

    Yobs <- rep(NA_real_, n)
    idx  <- which(S == 1)
    Yobs[idx] <- Tr[idx] * Y1[idx] + (1 - Tr[idx]) * Y0[idx]

    list(
      data_full = cbind(X, S = S, Tr = Tr, Yobs = Yobs),
      X = X, S = S, Tr = Tr, Y0 = Y0, Y1 = Y1, Yobs = Yobs, in_box = in_box
    )
  }

  # 1) Simulate
  sim <- simulate_box_dgp(n = 600, p = 100)
  inside_rate  <- mean(sim$S[sim$in_box] == 1)
  outside_rate <- mean(sim$S[!sim$in_box] == 1)
  expect_lt(inside_rate, outside_rate)

  # 2) ROOT
  root_out <- ROOT(
    data           = sim$data_full,
    outcome        = "Yobs",
    treatment      = "Tr",
    sample         = "S",
    seed           = 99,
    num_trees      = 5,
    top_k_trees    = TRUE,
    k              = 3,
    feature_est    = "Ridge",
    verbose        = FALSE
  )

  expect_type(root_out, "list")
  expect_true(all(c("D_rash", "D_forest", "w_forest", "rashomon_set", "f", "testing_data") %in% names(root_out)))
  expect_equal(length(root_out$w_forest), 5L)
  expect_true(length(root_out$rashomon_set) >= 1L && length(root_out$rashomon_set) <= 5L)
  expect_equal(nrow(root_out$D_rash), nrow(root_out$testing_data))
  expect_true(all(root_out$D_rash$w_opt %in% 0:1))
  expect_s3_class(root_out$f, "rpart")
  expect_true(any(root_out$D_rash$w_opt == 1) && any(root_out$D_rash$w_opt == 0))

  # 3) Characterization wrapper
  DataRCT    <- subset(sim$data_full, S == 1, select = c(grep("^X", names(sim$data_full), value = TRUE), "Tr", "Yobs"))
  DataTarget <- subset(sim$data_full, S == 0, select = grep("^X", names(sim$data_full), value = TRUE))
  covs <- grep("^X", names(DataRCT), value = TRUE)

  char_out <- characterizing_underrep(
    DataRCT               = DataRCT,
    covariateColName_RCT     = covs,
    trtColName_RCT     = "Tr",
    outcomeColName_RCT       = "Yobs",
    DataTarget            = DataTarget,
    covariateColName_TargetData  = covs,
    seed                  = 99,
    num_trees             = 5,
    top_k_trees           = TRUE,
    k                     = 3,
    feature_est           = "Ridge",
    verbose               = FALSE
  )

  expect_type(char_out, "list")
  expect_true(all(c("root","combined","leaf_summary") %in% names(char_out)))
  expect_equal(nrow(char_out$combined), nrow(sim$data_full))
  if (!is.null(char_out$leaf_summary)) {
    expect_true(all(c("rule","predicted_w","n","pct","label") %in% names(char_out$leaf_summary)))
    expect_true(all(char_out$leaf_summary$label %in% c("Under-represented (drop, w=0)", "Represented (keep, w=1)")))
  } else {
    expect_null(char_out$root$f)
  }
})
