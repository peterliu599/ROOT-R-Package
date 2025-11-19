test_that("ROOT runs in two-sample mode and returns structured outputs", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("rpart")
  skip_if_not_installed("withr")
  skip_if_not_installed("gbm")   # used if feature_est = "GBM" in user code paths
  skip_if_not_installed("mlbench")

  sim <- get_data(n = 600, seed = 20)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  set.seed(42)
  out <- ROOT(
    data        = df,
    outcome     = "Yobs",
    treatment   = "Tr",
    sample      = "S",
    seed        = 99,
    num_trees   = 6,
    top_k_trees = TRUE,
    k           = 3,
    feature_est = "Ridge",
    verbose     = TRUE
  )

  expect_s3_class(out, "ROOT")
  expect_true(all(c("D_rash","D_forest","w_forest","rashomon_set","f","testing_data","estimate") %in% names(out)))
  expect_equal(length(out$w_forest), 6L)
  expect_true(length(out$rashomon_set) >= 1L)
  expect_true(all(out$D_rash$w_opt %in% 0:1))

  # Estimates present with labels and SDs (keys updated from se_* to sd_*)
  est <- out$estimate
  expect_true(all(c("estimand_unweighted","value_unweighted","se_unweighted",
                    "estimand_weighted","value_weighted","se_weighted",
                    "n_analysis","sum_w") %in% names(est)))
  expect_true(is.numeric(est$value_unweighted))
  expect_true(is.numeric(est$value_weighted) || is.na(est$value_weighted))

  # summary.ROOT prints without error
  expect_invisible(summary(out))
})

test_that("ROOT single-sample mode works when sample=NULL", {
  skip_if_not_installed("mlbench")
  sim <- get_data(n = 300, seed = 21)
  # Single-sample: keep only S==1 and drop sample
  dfS <- subset(sim$data, S == 1)
  dfS$Yobs <- dfS$Yobs

  out <- ROOT(
    data        = dfS[, c(grep("^X", names(dfS), value = TRUE), "Tr", "Yobs")],
    outcome     = "Yobs",
    treatment   = "Tr",
    sample      = NULL,
    num_trees   = 4,
    vote_threshold = 0.6
  )
  expect_s3_class(out, "ROOT")
  expect_true(out$single_sample_mode)
  expect_equal(out$estimate$estimand_unweighted, "SATE")
  expect_equal(out$estimate$estimand_weighted,   "WATE")
})

test_that("ROOT input validation catches errors", {
  sim <- get_data(n = 100, seed = 5)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  expect_error(ROOT(df, "Yobs", "Tr", sample = "S_missing"), "not found")
  expect_error(ROOT(df, "Yobs", "Tr", sample = "S", num_trees = 0), "must be positive")
  expect_error(ROOT(df, "Yobs", "Tr", sample = "S", vote_threshold = 1.1), "in \\(0, 1\\]")
})
