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
  #expect_visible(summary(out))
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

mk_data_two_sample <- function(n=200, p=5, seed=11) {
  set.seed(seed)
  X <- as.data.frame(matrix(rnorm(n*p), n, p)); names(X) <- paste0("X", seq_len(p))
  Tr <- rbinom(n, 1, 0.5)
  S  <- rbinom(n, 1, 0.5)
  Y  <- rnorm(n)
  data.frame(X, Tr=Tr, S=S, Yobs=Y)
}

test_that("ROOT() validates arguments and 0/1 coercion", {
  d <- mk_data_two_sample()

  expect_error(ROOT(as.matrix(d), "Yobs", "Tr", "S"), "`data` must be a data frame")
  expect_error(ROOT(d, "Nope", "Tr", "S"), "`outcome` must be a single column name")
  expect_error(ROOT(d, "Yobs", "ZZ", "S"), "`treatment` must be a single column name")
  expect_error(ROOT(d, "Yobs", "Tr", "SS"), "`sample` column not found")

  expect_error(ROOT(d, "Yobs", "Tr", "S", leaf_proba = -0.1), "between 0 and 1")
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees = 0), "must be positive")
  expect_error(ROOT(d, "Yobs", "Tr", "S", vote_threshold = 0), "in \\(0, 1\\]")
  expect_error(ROOT(d, "Yobs", "Tr", "S", explore_proba = 2), "between 0 and 1")
  expect_error(ROOT(d, "Yobs", "Tr", "S", feature_est = 123), "must be \"Ridge\", \"GBM\", or a function")
  expect_error(ROOT(d, "Yobs", "Tr", "S", top_k_trees = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(ROOT(d, "Yobs", "Tr", "S", k = 0), "positive integer")
  expect_error(ROOT(d, "Yobs", "Tr", "S", cutoff = "nope"), "must be \"baseline\" or numeric")
  expect_error(ROOT(d, "Yobs", "Tr", "S", verbose = 1), "TRUE or FALSE")
  expect_error(ROOT(d, "Yobs", "Tr", "S", global_objective_fn = 1), "must be a function")

  # OK: treatment/sample coercion of strings
  d2 <- d
  d2$Tr <- sample(c("treated","control","yes","no","True","False"), nrow(d2), TRUE)
  d2$S  <- sample(c("1","0","yes","no","true","false"), nrow(d2), TRUE)
  # Won't error, but will stop later if estimate_dml needs numeric; just ensure coercion branch runs:
  expect_error(
    ROOT(d2, "Yobs", "Tr", "S", num_trees=1),
    regexp = NA
  )
})

test_that("ROOT single-sample mode via sample=NULL and via constant S", {
  d <- mk_data_two_sample()
  dS1 <- d; dS1$S <- 1L

  skip_if_not_installed("mlbench")

  # (a) sample=NULL path (SATE/WATE)
  set.seed(7)
  r1 <- ROOT(dS1[, c(grep("^X", names(dS1), value=TRUE), "Tr", "Yobs")],
             outcome="Yobs", treatment="Tr", sample=NULL, num_trees=2)
  expect_true(r1$single_sample_mode)
  expect_identical(r1$estimate$estimand_unweighted, "SATE")
  expect_identical(r1$estimate$estimand_weighted, "WATE")

  # (b) sample column present but constant -> also single_sample_mode TRUE
  set.seed(8)
  r2 <- ROOT(dS1, outcome="Yobs", treatment="Tr", sample="S", num_trees=2)
  expect_true(r2$single_sample_mode)
})

test_that("ROOT feature_est = Ridge / GBM / custom; bad custom errors", {
  d <- mk_data_two_sample(n=220, p=6)

  skip_if_not_installed("MASS")
  skip_if_not_installed("gbm")
  skip_if_not_installed("mlbench")

  # Ridge
  set.seed(1)
  rR <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, feature_est="Ridge")
  expect_s3_class(rR, "ROOT")

  # GBM
  set.seed(2)
  rG <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, feature_est="GBM")
  expect_s3_class(rG, "ROOT")

  # custom ok: returns named nonnegative importances for all X columns
  ok_imp <- function(X, y, ...) { setNames(abs(colSums(X^2)) + 1e-6, colnames(X)) }
  set.seed(3)
  rC <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, feature_est=ok_imp)
  expect_s3_class(rC, "ROOT")

  # custom bad: unnamed -> error from .norm_feat_prob
  bad1 <- function(X, y, ...) { as.numeric(abs(colSums(X))) }
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees=1, feature_est=bad1),
               "named numeric vector")

  # custom bad: missing a name -> error
  bad2 <- function(X, y, ...) {
    nm <- colnames(X); nm[length(nm)] <- "WRONG"
    setNames(abs(colSums(X)), nm)
  }
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees=1, feature_est=bad2),
               "Importance missing for some X_df columns")

  # custom bad: negative -> error
  bad3 <- function(X, y, ...) { setNames(c(-1, rep(1, ncol(X)-1)), colnames(X)) }
  expect_error(ROOT(d, "Yobs", "Tr", "S", num_trees=1, feature_est=bad3),
               "must be non-negative")
})

test_that("Rashomon selection: top-k, baseline cutoff, and empty set warning", {
  skip_if_not_installed("mlbench")
  d <- mk_data_two_sample(n=260, p=6)

  set.seed(11)
  r_topk <- ROOT(d, "Yobs", "Tr", "S", num_trees=5, top_k_trees=TRUE, k=3)
  expect_length(r_topk$rashomon_set, 3)

  set.seed(12)
  r_base <- ROOT(d, "Yobs", "Tr", "S", num_trees=5, top_k_trees=FALSE, cutoff="baseline")
  expect_true(length(r_base$rashomon_set) >= 1)

  set.seed(13)
  # absurdly small cutoff forces empty set → warning; w_opt should exist but be all 0
  expect_warning({
    r_empty <- ROOT(d, "Yobs", "Tr", "S", num_trees=3, top_k_trees=FALSE, cutoff = -1e9)
    expect_length(r_empty$rashomon_set, 0)
    expect_true("w_opt" %in% names(r_empty$D_rash))
    expect_true(all(r_empty$D_rash$w_opt == 0))
  }, "No trees selected")
})

test_that("ROOT estimate fields populated; binary vs non-binary weighted SE behavior", {
  skip_if_not_installed("mlbench")
  d <- mk_data_two_sample(n=280, p=6)

  set.seed(21)
  r <- ROOT(d, "Yobs", "Tr", "S", num_trees=4)
  e <- r$estimate
  expect_true(all(c("estimand_unweighted","value_unweighted","se_unweighted",
                    "estimand_weighted","value_weighted","se_weighted",
                    "se_weighted_note","n_analysis","sum_w") %in% names(e)))
  # Binary case: either valid SE or NA when n_A=1
  expect_true(is.numeric(e$se_weighted) || is.na(e$se_weighted))

  # Force non-binary to hit omit-SE path in summary/print
  r2 <- r
  r2$D_rash$w_opt <- as.numeric(r2$D_rash$w_opt) + runif(nrow(r2$D_rash), 0, 0.1)
  # Call summary/print to ensure they don't error
  expect_silent(capture.output(summary(r2)))
  expect_silent(capture.output(print(r2)))
})

test_that("ROOT seed yields reproducible w_forest & estimates; verbose prints", {
  skip_if_not_installed("mlbench")
  d <- mk_data_two_sample(n=200, p=5)

  r1 <- ROOT(d, "Yobs", "Tr", "S", num_trees=3, seed=123, verbose=TRUE)
  r2 <- ROOT(d, "Yobs", "Tr", "S", num_trees=3, seed=123, verbose=TRUE)

  expect_equal(lapply(r1$w_forest, `[[`, "local objective"),
               lapply(r2$w_forest, `[[`, "local objective"))
  expect_equal(r1$estimate$value_unweighted, r2$estimate$value_unweighted, tolerance=1e-12)
  expect_equal(r1$estimate$value_weighted,  r2$estimate$value_weighted,  tolerance=1e-12)

  # verbose -> messages (do not assert exact strings, just that something prints)
  expect_output({
    r3 <- ROOT(d, "Yobs", "Tr", "S", num_trees=2, seed=124, verbose=TRUE)
    print(r3)  # ensure print path runs too
  })
})

test_that("ROOT runs (two-sample) with custom objective and zero-sum feature importances", {
  set.seed(10)
  n  <- 180
  X1 <- rnorm(n); X2 <- rnorm(n)
  S  <- rbinom(n, 1, plogis(0.4*X1 - 0.2*X2))
  Tr <- rbinom(n, 1, plogis(0.3 + 0.7*X1 - 0.2*X2))
  Y  <- 0.5 + 1.1*Tr + 0.5*X1 - 0.3*X2 + rnorm(n, 0.7)

  dat <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)

  # Zero-sum feature importance → code falls back to uniform
  zero_imp <- function(X, y, ...) stats::setNames(rep(0, ncol(X)), colnames(X))
  my_obj   <- function(D) mean(D$vsq, na.rm = TRUE) + 1e-8

  # Use top_k_trees to guarantee a non-empty Rashomon set (and thus nontrivial votes)
  fit <- ROOT(
    data = dat, outcome = "Y", treatment = "Tr", sample = "S",
    seed = 123, num_trees = 5, top_k_trees = TRUE, k = 5,
    vote_threshold = 2/3, feature_est = zero_imp, verbose = FALSE,
    global_objective_fn = my_obj
  )
  expect_s3_class(fit, "ROOT")

  est <- fit$estimate
  # Weighted-SE note must contain the standard WTATE formula OR the "no kept observations" text
  expect_true(
    grepl("Calculation of SE for WTATE uses sqrt", est$se_weighted_note, fixed = TRUE) ||
      grepl("no kept observations", est$se_weighted_note, ignore.case = TRUE)
  )
})

test_that("ROOT argument guards and coercions", {
  set.seed(1)
  df <- data.frame(
    Y  = rnorm(40),
    Tr = sample(0:1, 40, TRUE),
    S  = sample(0:1, 40, TRUE),
    X  = rnorm(40)
  )

  # Type/name checks
  expect_error(ROOT(as.matrix(df), "Y", "Tr", "S"),
               "`data` must be a data frame.", fixed = TRUE)
  expect_error(ROOT(df, "ZZ", "Tr", "S"),
               "`outcome` must be a single column name present in `data`.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "TT", "S"),
               "`treatment` must be a single column name present in `data`.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "SSS"),
               "`sample` column not found; pass `sample = NULL` to run single-sample mode.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", sample = c("S","S")),
               "`sample` must be NULL or a single column name string.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", sample = 1L),
               "`sample` must be NULL or a single column name string.", fixed = TRUE)

  # No covariates
  df_nocov <- df[, c("Y","Tr","S")]
  expect_error(
    ROOT(df_nocov, "Y", "Tr", "S"),
    "ROOT\\(\\): no covariate columns found \\(need at least one feature besides outcome/treatment/sample\\)\\."
  )

  # Scalar-argument guards (now assert exact messages rather than snapshots)
  expect_error(ROOT(df, "Y", "Tr", "S", leaf_proba = -0.1),
               "`leaf_proba` must be between 0 and 1.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", seed = c(1,2)),
               "`seed` must be NULL or a single numeric value.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", num_trees = 0),
               "`num_trees` must be positive.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", vote_threshold = 0),
               "`vote_threshold` must be in (0, 1].", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", explore_proba = 2),
               "`explore_proba` must be between 0 and 1.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", feature_est = 123),
               "`feature_est` must be \"Ridge\", \"GBM\", or a function.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", k = 0),
               "`k` must be a positive integer.", fixed = TRUE)
  expect_error(ROOT(df, "Y", "Tr", "S", cutoff = "not-baseline"),
               "`cutoff` must be \"baseline\" or numeric\\.", perl = TRUE)

  # Bad sample values (exercises coerce01 NA-forbidden path)
  df_badS <- df
  df_badS$S <- c("yes","maybe", rep("no", nrow(df)-2))
  expect_error(ROOT(df_badS, "Y", "Tr", "S"),
               "Non 0/1 values found\\.", perl = TRUE)

  # Treatment with no variation among S==1 → DML training error from train()
  df_badT <- df
  df_badT$S  <- ifelse(runif(nrow(df)) < 0.7, 1L, 0L)
  df_badT$Tr[df_badT$S == 1L] <- 1L  # no variation in treated among S==1
  expect_error(
    ROOT(df_badT, "Y", "Tr", "S"),
    "Cannot train model: treatment has no variation among S==1 observations\\.", perl = TRUE
  )

  # All S == 1 should quietly switch to single-sample mode; avoid ridge path by supplying a custom FE
  df_all1 <- transform(
    df,
    S = 1L,
    X2 = rnorm(nrow(df))
  )
  uni_imp <- function(X, y, ...) stats::setNames(rep(1, ncol(X)), colnames(X))

  # Suppress incidental GLM warnings from tiny synthetic data
  fit <- suppressWarnings(
    ROOT(df_all1, "Y", "Tr", "S",
         feature_est = uni_imp, seed = 99, num_trees = 3, verbose = FALSE)
  )
  expect_true(isTRUE(fit$single_sample_mode))
})
test_that("Rashomon top-k warning, cutoff non-finite, and ridge-uniform path", {
  set.seed(11)
  n  <- 60
  X1 <- rnorm(n); X2 <- rnorm(n)
  S  <- rbinom(n, 1, 0.5)
  Tr <- rbinom(n, 1, plogis(0.6*X1))
  # Make Y so that v-squared is (nearly) zero after centering -> triggers ridge uniform branch
  Y  <- 3 + 0*Tr + 0*X1 + rnorm(n, sd = 1e-8)

  dat <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)

  # top_k_trees with k > num_trees emits a warning and clamps k
  expect_warning(
    fit1 <- ROOT(dat, "Y","Tr","S", seed = 321, num_trees = 3,
                 feature_est = "Ridge", top_k_trees = TRUE, k = 999),
    "k > num_trees; using k = num_trees."
  )
  expect_s3_class(fit1, "ROOT")

  # cutoff = NA_real_ -> non-finite cutoff handled (set to Inf)
  expect_silent(
    fit2 <- ROOT(dat, "Y","Tr","S", seed = 322, num_trees = 3,
                 feature_est = "Ridge", cutoff = NA_real_)
  )
  expect_s3_class(fit2, "ROOT")
})

test_that("ROOT single-sample mode uses SATE/WATE labels and returns estimates", {
  set.seed(12)
  n  <- 80
  X1 <- rnorm(n); X2 <- rnorm(n)
  Tr <- rbinom(n, 1, plogis(0.4*X1 - 0.1*X2))
  Y  <- 1 + Tr + 0.5*X1 + rnorm(n)

  dat <- data.frame(Y=Y, Tr=Tr, X1=X1, X2=X2)

  res <- ROOT(dat, outcome = "Y", treatment = "Tr", sample = NULL, seed = 42, num_trees = 4)
  expect_true(res$single_sample_mode)
  expect_identical(res$estimate$estimand_unweighted, "SATE")
  expect_identical(res$estimate$estimand_weighted,   "WATE")
  expect_true(is.numeric(res$estimate$value_unweighted))
  expect_true(is.numeric(res$estimate$value_weighted))
})

test_that("ROOT informs when no summary tree is available (single-class w_opt)", {
  set.seed(13)
  n  <- 70
  X1 <- rnorm(n); X2 <- rnorm(n)
  S  <- rbinom(n, 1, 0.5)
  Tr <- rbinom(n, 1, plogis(X1 - 0.2*X2))
  Y  <- 0.2 + Tr + X1 + rnorm(n)

  dat <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)

  # vote_threshold = 1 makes it very hard for any row to be kept by ALL trees -> usually all 0
  expect_message(
    fit <- ROOT(dat, "Y","Tr","S", seed=99, num_trees=4, vote_threshold = 1, verbose = FALSE),
    "No summary tree available to plot"
  )
  expect_s3_class(fit, "ROOT")
})

test_that("coerce01() NA-forbidden path is exercised via sample column", {
  base <- data.frame(
    Y=rnorm(8), Tr=sample(0:1, 8, TRUE),
    S=c("yes","no","maybe","no","yes","unknown","0","1"),
    X=rnorm(8)
  )
  # sample contains values not mapped to 0/1 -> coerce01(..., allow_na=FALSE) must error
  expect_snapshot_error(ROOT(base, "Y","Tr","S"))
})
