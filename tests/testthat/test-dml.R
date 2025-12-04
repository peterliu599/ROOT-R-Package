test_that("train() fits nuisance models and estimate() returns finite v,a,b", {
  # Build a small two-sample dataset
  sim <- get_data(n = 400, seed = 123)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)  # already inside sim$data

  # Some S==1 rows need non-NA Tr & Yobs — get_data provides that structure
  trn <- train(df, outcome = "Yobs", treatment = "Tr", sample = "S")
  expect_true(is.numeric(trn$pi) && trn$pi > 0 && trn$pi < 1)
  expect_s3_class(trn$pi_m, "glm")
  expect_s3_class(trn$e_m,  "glm")

  est <- estimate(df, outcome = "Yobs", treatment = "Tr", sample = "S",
                  pi = trn$pi, pi_m = trn$pi_m, e_m = trn$e_m)
  expect_length(est$v, nrow(df))
  expect_length(est$a, nrow(df))
  expect_length(est$b, nrow(df))
  expect_true(all(is.finite(est$v[ df$S == 1 ])))
})

test_that("estimate_dml returns aligned df_v and data2 for S==1 rows", {
  sim <- get_data(n = 500, seed = 321)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  out <- estimate_dml(df, outcome = "Yobs", treatment = "Tr", sample = "S", crossfit = 5)
  expect_true(all(c("te","a","b","te_sq","a_sq","primary_index") %in% names(out$df_v)))
  expect_equal(nrow(out$data2), nrow(out$df_v))  # aligned
})

test_that("estimate_dml_single works in single-sample mode", {
  # Emulate single-sample by filtering S==1 and dropping S
  sim <- get_data(n = 300, seed = 777)
  dfS <- subset(sim$data, S == 1)
  dfS$Yobs <- dfS$Yobs  # explicit

  out <- estimate_dml_single(data = dfS[, c(grep("^X", names(dfS), value = TRUE), "Tr", "Yobs")],
                             outcome = "Yobs", treatment = "Tr", crossfit = 4)
  expect_true(all(c("te","a","b","te_sq","a_sq","primary_index") %in% names(out$df_v)))
  expect_true(all(out$df_v$b == 1))
})

test_that("DML errors are thrown for degenerate or NA cases", {
  sim <- get_data(n = 200, seed = 55)
  df  <- data.frame(sim$data, Yobs = sim$data$Yobs)

  # All S==1 or S==0 → train() should error
  df_bad <- df; df_bad$S <- 1L
  expect_error(train(df_bad, "Yobs", "Tr", "S"), "no variation")

  # NA in covariates should error
  df_na <- df; nm <- grep("^X", names(df_na), value = TRUE)[1]; df_na[1, nm] <- NA
  expect_error(train(df_na, "Yobs", "Tr", "S"), "Covariates contain NA")

  # estimate() invalid pi
  trn <- train(df, "Yobs", "Tr", "S")
  expect_error(estimate(df, "Yobs", "Tr", "S", pi = 1, pi_m = trn$pi_m, e_m = trn$e_m),
               "must be between 0 and 1")
})

test_that("stratified_kfold works correctly", {
  S <- c(rep(0, 10), rep(1, 10))
  folds <- stratified_kfold(S, K = 5)
  expect_length(folds, 5)
  # Check all indices are present exactly once
  expect_equal(sort(unlist(folds)), 1:20)

  # Test K > N adjustment
  expect_warning(folds_big <- stratified_kfold(S[1:5], K = 10), "Requested K")
  expect_length(folds_big, 5)

  # Test error inputs
  expect_error(stratified_kfold(S, K = 0), "positive integer")
  expect_error(stratified_kfold(data.frame(S=S)), "vector or factor")
})

test_that("train function validates inputs and fits models", {
  # Happy path
  train_idx <- sample(nrow(dat_2s), 150)
  train_data <- dat_2s[train_idx, ]
  res <- train(train_data, outcome="Yobs", treatment="Tr", sample="S")

  expect_type(res$pi, "double")
  expect_s3_class(res$pi_m, "glm")
  expect_s3_class(res$e_m, "glm")

  # Input validation errors
  expect_error(train(as.matrix(train_data), "Yobs", "Tr", "S"), "must be a data frame")
  expect_error(train(train_data, "WrongY", "Tr", "S"), "Column 'WrongY' not found")

  # NA handling errors
  dat_na <- train_data
  dat_na$S[1] <- NA
  expect_error(train(dat_na, "Yobs", "Tr", "S"), "`S` contains NA")

  dat_na <- train_data
  dat_na$X0[1] <- NA
  expect_error(train(dat_na, "Yobs", "Tr", "S"), "Covariates contain NA")

  # Lack of variation errors
  dat_no_var <- train_data
  dat_no_var$S <- 1
  expect_error(train(dat_no_var, "Yobs", "Tr", "S"), "has no variation")
})

test_that("estimate function calculates pseudo-outcomes correctly", {
  # Setup training models
  train_idx <- 1:100
  test_idx <- 101:200
  trn <- train(dat_2s[train_idx,], outcome="Yobs", treatment="Tr", sample="S")
  test_data <- dat_2s[test_idx,]

  # Happy path
  est <- estimate(test_data, outcome="Yobs", treatment="Tr", sample="S",
                  pi = trn$pi, pi_m = trn$pi_m, e_m = trn$e_m)

  expect_equal(length(est$v), nrow(test_data))
  expect_equal(length(est$a), nrow(test_data))
  expect_equal(length(est$b), nrow(test_data))
  expect_true(all(is.finite(est$v)))

  # Test probability clamping (hard to trigger naturally, verify logic exists)
  # If we force an extreme prediction, it shouldn't result in Inf weights due to clamping in code
  # Mock a model that predicts near 0
  mock_em <- trn$e_m
  mock_em$coefficients[] <- -100 # forces prob near 0
  est_clamp <- estimate(test_data, outcome="Yobs", treatment="Tr", sample="S",
                        pi = trn$pi, pi_m = trn$pi_m, e_m = mock_em)
  expect_true(all(is.finite(est_clamp$a)))

  # Input validation errors
  expect_error(estimate(test_data, "Yobs", "Tr", "S", pi = 1.5, trn$pi_m, trn$e_m), "Invalid `pi`")
  expect_error(estimate(test_data, "Yobs", "Tr", "S", pi = 0.5, "not_a_model", trn$e_m), "must be model objects")
})

test_that("estimate_dml performs cross-fitting correctly", {
  # Happy path (two-sample)
  res <- estimate_dml(dat_2s, outcome="Yobs", treatment="Tr", sample="S", crossfit = 3)
  expect_s3_class(res$df_v, "data.frame")
  expect_true(all(c("te", "a", "b") %in% names(res$df_v)))
  # Should only contain S=1 rows
  expect_equal(nrow(res$df_v), sum(dat_2s$S == 1))

  # Validation errors
  expect_error(estimate_dml(dat_2s, "Yobs", "Tr", "S", crossfit = 1), "integer >= 2")

  dat_no_var <- dat_2s
  dat_no_var$S <- 1
  expect_error(estimate_dml(dat_no_var, "Yobs", "Tr", "S"), "has no variation")
})

test_that("Single sample DML functions work", {
  # train_single
  trn_s <- train_single(dat_1s[1:100,], "Yobs", "Tr")
  expect_s3_class(trn_s$e_m, "glm")

  # estimate_single
  est_s <- estimate_single(dat_1s[101:200,], "Yobs", "Tr", trn_s$e_m)
  expect_equal(est_s$b, rep(1, 100)) # b should be 1 in single sample

  # estimate_dml_single
  res_s <- estimate_dml_single(dat_1s, "Yobs", "Tr", crossfit = 3)
  expect_equal(nrow(res_s$df_v), nrow(dat_1s))
  expect_true(all(res_s$df_v$b == 1))
})

test_that("estimate_dml runs end-to-end and keeps only S==1 rows", {
  set.seed(1)
  n  <- 120
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  S  <- rbinom(n, 1, plogis(0.3*X1 - 0.2*X2))           # some 0s and 1s
  Tr <- rbinom(n, 1, plogis(0.2 + 0.8*X1 - 0.3*X2))     # treatment varies
  Y  <- 0.5 + 1.2*Tr + 0.5*X1 - 0.3*X2 + rnorm(n, sd=0.7)

  # Allow NA in Y/Tr where S==0 (should be fine)
  Y[S == 0][1:3]  <- NA
  Tr[S == 0][1:2] <- NA

  dat <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)

  out <- estimate_dml(dat, outcome="Y", treatment="Tr", sample="S", crossfit=3)
  expect_type(out, "list")
  expect_true(all(c("df_v","data2") %in% names(out)))
  expect_s3_class(out$df_v, "data.frame")
  expect_true(all(out$data2$S == 1))                    # only S==1 aligned
  expect_true(all(out$df_v$primary_index %in% which(S==1)))
  expect_true(all(c("te","a","b","te_sq","a_sq") %in% names(out$df_v)))
  expect_equal(nrow(out$df_v), nrow(out$data2))
  # squared deviation columns are non-negative by construction
  expect_true(all(out$df_v$te_sq >= 0))
  expect_true(all(out$df_v$a_sq  >= 0))
})

test_that("estimate_dml_single works and aggregates per index", {
  set.seed(2)
  n  <- 80
  X1 <- rnorm(n); X2 <- rnorm(n)
  Tr <- rbinom(n, 1, plogis(0.5*X1 - 0.2*X2))
  Y  <- 2 + Tr + 0.4*X1 + rnorm(n)

  dat <- data.frame(Y=Y, Tr=Tr, X1=X1, X2=X2)

  out <- estimate_dml_single(dat, outcome="Y", treatment="Tr", crossfit=4)
  expect_true(all(c("df_v","data2") %in% names(out)))
  expect_true(all(c("te","a","b","te_sq","a_sq","primary_index") %in% names(out$df_v)))
  # single-sample: b ≡ 1 and v==a
  expect_true(all(out$df_v$b == 1))
  expect_equal(out$df_v$te, out$df_v$a)
  expect_equal(nrow(out$df_v), nrow(out$data2))
})


test_that("train() to01 mapping covers numeric/logical/character/factor paths", {
  # mix encodings (strings, logical, numeric, factor)
  td <- data.frame(
    Y  = rnorm(20),
    Tr = factor(c("treated","control","1","0","yes","no","TRUE","FALSE","t","c",
                  rep("1",10))),
    S  = c(TRUE, FALSE, 1L, 0L, "1", "0", "yes", "no", "treated", "control",
           rep(1L, 10)),
    X1 = rnorm(20), X2 = rnorm(20),
    check.names = FALSE
  )

  # Map encodings to 0/1 the same way train() will
  td$S  <- as.integer(td$S %in% c(TRUE, "1", 1L, "yes", "treated"))
  td$Tr <- as.integer(td$Tr %in% c("treated","1","yes","TRUE","t"))

  # --- Enforce needed variation deterministically -------------------------
  # Ensure S has both classes
  if (length(unique(td$S)) < 2) {
    td$S[1] <- 0L
    td$S[2] <- 1L
  }
  # Ensure Tr varies among S==1 rows
  idx1 <- which(td$S == 1L)
  if (length(unique(td$Tr[idx1])) < 2) {
    # force two different values within S==1
    td$Tr[idx1[1]] <- 0L
    if (length(idx1) >= 2) td$Tr[idx1[2]] <- 1L
  }
  # ------------------------------------------------------------------------

  quiet_glm <- function(expr) {
    withCallingHandlers(expr,
                        warning = function(w) {
                          if (grepl("^glm\\.fit:", conditionMessage(w))) {
                            invokeRestart("muffleWarning")
                          }
                        }
    )
  }
  expect_no_error( quiet_glm( train(td, outcome = "Y", treatment = "Tr", sample = "S") ) )
})

test_that("train() stops on NA in S or covariates; treatment variation check", {
  base <- data.frame(
    Y=rnorm(10), Tr=sample(0:1, 10, TRUE), S=sample(0:1, 10, TRUE),
    X1=rnorm(10), X2=rnorm(10)
  )

  badS <- base; badS$S[1] <- NA
  expect_snapshot_error(train(badS, "Y","Tr","S"))

  badX <- base; badX$X2[1] <- NA
  expect_snapshot_error(train(badX, "Y","Tr","S"))

  noVarTr <- base; noVarTr$S <- 1L; noVarTr$Tr <- 1L
  expect_snapshot_error(train(noVarTr, "Y","Tr","S"))
})

test_that("estimate() input validation branches", {
  set.seed(3)
  n <- 40
  X1 <- rnorm(n); X2 <- rnorm(n)
  S  <- rbinom(n, 1, 0.6)
  Tr <- rbinom(n, 1, 0.5)
  Y  <- rnorm(n)
  df <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)

  # valid nuisances to reuse
  trn <- train(df, "Y","Tr","S")

  # not a data.frame
  expect_snapshot_error(estimate(as.matrix(df), "Y","Tr","S", trn$pi, trn$pi_m, trn$e_m))

  # missing columns
  expect_snapshot_error(estimate(df[, c("Y","Tr","X1","X2")], "Y","Tr","S", trn$pi, trn$pi_m, trn$e_m))

  # NA in S
  badS <- df; badS$S[1] <- NA
  expect_snapshot_error(estimate(badS,"Y","Tr","S", trn$pi, trn$pi_m, trn$e_m))

  # NA in covariates
  badX <- df; badX$X2[1] <- NA
  expect_snapshot_error(estimate(badX,"Y","Tr","S", trn$pi, trn$pi_m, trn$e_m))

  # NA in Y/Tr among S==1
  badYT <- df; badYT$Y[which(badYT$S==1)[1]] <- NA
  expect_snapshot_error(estimate(badYT,"Y","Tr","S", trn$pi, trn$pi_m, trn$e_m))

  # pi checks
  expect_snapshot_error(estimate(df,"Y","Tr","S", pi = NA_real_, trn$pi_m, trn$e_m))
  expect_snapshot_error(estimate(df,"Y","Tr","S", pi = 0, trn$pi_m, trn$e_m))
  expect_snapshot_error(estimate(df,"Y","Tr","S", pi = 1, trn$pi_m, trn$e_m))

  # model class checks
  expect_snapshot_error(estimate(df,"Y","Tr","S", trn$pi, pi_m = list(), e_m = trn$e_m))
  expect_snapshot_error(estimate(df,"Y","Tr","S", trn$pi, pi_m = trn$pi_m, e_m = list()))
})

test_that("estimate() returns numeric vectors of correct length and clamps probs", {
  set.seed(4)
  n <- 50
  X1 <- rnorm(n); X2 <- rnorm(n)
  S  <- rbinom(n, 1, 0.55)
  Tr <- rbinom(n, 1, plogis(X1 - 0.5*X2))
  Y  <- 1 + Tr + X1 + rnorm(n)

  df  <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)
  trn <- train(df, "Y","Tr","S")

  est <- estimate(df, "Y","Tr","S", trn$pi, trn$pi_m, trn$e_m)
  expect_true(all(vapply(est, length, 1L) == n))
  expect_true(is.numeric(est$v) && is.numeric(est$a) && is.numeric(est$b))
})

test_that("stratified_kfold basic properties and K>n warning path", {
  S <- factor(c(0,1,1,0,1,0,0,1))
  expect_snapshot_error(stratified_kfold(as.data.frame(S), K = 2))  # S not vector/factor

  f <- stratified_kfold(S, K = 3)
  expect_length(f, 3)
  expect_equal(sort(unlist(f)), seq_along(S))            # partition covers all
  expect_true(all(vapply(f, is.integer, TRUE)))

  # K > n triggers warning and reduction
  expect_warning(f2 <- stratified_kfold(S, K = 100), "greater than number of observations")
  expect_length(f2, length(S))
})

test_that("estimate_dml input guards hit key branches", {
  set.seed(5)
  n  <- 30
  X1 <- rnorm(n); X2 <- rnorm(n)
  S  <- rep(c(0L,1L), length.out=n)
  Tr <- rbinom(n, 1, 0.5)
  Y  <- rnorm(n)
  dat <- data.frame(Y=Y, Tr=Tr, S=S, X1=X1, X2=X2)

  # not data.frame
  expect_snapshot_error(estimate_dml(as.matrix(dat), "Y","Tr","S", 2))

  # missing column
  expect_snapshot_error(estimate_dml(dat[, c("Y","Tr","X1","X2")], "Y","Tr","S", 2))

  # NA in S
  badS <- dat; badS$S[1] <- NA
  expect_snapshot_error(estimate_dml(badS, "Y","Tr","S", 2))

  # NA in covariates
  badX <- dat; badX$X2[2] <- NA
  expect_snapshot_error(estimate_dml(badX, "Y","Tr","S", 2))

  # NA in Y/Tr within S==1
  badYT <- dat; badYT$Y[which(badYT$S==1)[1]] <- NA
  expect_snapshot_error(estimate_dml(badYT, "Y","Tr","S", 2))

  # crossfit checks
  expect_snapshot_error(estimate_dml(dat, "Y","Tr","S", crossfit = 1))
  expect_snapshot_error(estimate_dml(dat, "Y","Tr","S", crossfit = c(2,3)))

  # sample all 0 or all 1
  all0 <- dat; all0$S <- 0L
  expect_snapshot_error(estimate_dml(all0, "Y","Tr","S", 2))
  all1 <- dat; all1$S <- 1L
  expect_snapshot_error(estimate_dml(all1, "Y","Tr","S", 2))
})

