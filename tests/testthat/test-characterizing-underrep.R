test_that("characterizing_underrep integrates ROOT and returns leaf summaries", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("mlbench")

  sim <- get_data(n = 500, seed = 1234)
  full <- data.frame(sim$data, Yobs = sim$data$Yobs)
  covs <- grep("^X", names(full), value = TRUE)

  DataRCT    <- subset(full, S == 1, select = c(covs, "Tr", "Yobs"))
  DataTarget <- subset(full, S == 0, select = covs)

  out <- characterizing_underrep(
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
  )

  expect_s3_class(out, "characterizing_underrep")
  expect_true(all(c("root","combined","leaf_summary") %in% names(out)))
  expect_equal(nrow(out$combined), nrow(full))

  # summary.characterizing_underrep prints without error
  expect_invisible(summary(out))

  if (!is.null(out$leaf_summary)) {
    expect_true(all(c("rule","predicted_w","n","pct","label") %in% names(out$leaf_summary)))
    expect_true(all(out$leaf_summary$label %in% c("Under-represented (drop, w=0)", "Represented (keep, w=1)")))
  }
})

test_that("characterizing_underrep validation errors on mismatched covariates", {
  sim <- get_data(n = 200, seed = 66)
  full <- data.frame(sim$data, Yobs = sim$data$Yobs)
  covs <- grep("^X", names(full), value = TRUE)

  DataRCT    <- subset(full, S == 1, select = c(covs, "Tr", "Yobs"))
  DataTarget <- subset(full, S == 0, select = covs[-1])  # drop one to force warning/stop

  expect_error(
    characterizing_underrep(
      DataRCT, covariateColName_RCT = covs, trtColName_RCT = "Tr", outcomeColName_RCT = "Yobs",
      DataTarget, covariateColName_TargetData = covs
    ),
    "Missing target covariates"
  )
})
