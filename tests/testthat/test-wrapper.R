# helper in tests/test-helper.R (or at top of the test file)
fake_root <- function(D_forest,
                      D_rash = data.frame(),
                      rashomon_set = integer(0),
                      estimate = NULL,
                      f = NULL,
                      w_forest = NULL,
                      global_objective_fn = objective_default,
                      single_sample_mode = NULL) {
  obj <- list(
    D_forest = D_forest,
    D_rash = D_rash,
    rashomon_set = rashomon_set,
    testing_data = D_forest,           # used by diagnostics
    estimate = estimate,
    f = f,
    w_forest = w_forest,
    global_objective_fn = global_objective_fn,
    single_sample_mode = single_sample_mode
  )
  class(obj) <- "ROOT"
  obj
}

dat_rct <- dat_2s[dat_2s$S == 1, ]
dat_tgt <- dat_2s[dat_2s$S == 0, setdiff(names(dat_2s), c("Yobs", "Tr", "S"))]
covs <- paste0("X", 0:9)

test_that("characterizing_underrep wrapper works", {
  suppressWarnings({
    res <- characterizing_underrep(
      DataRCT = dat_rct[1:50,],
      covariateColName_RCT = covs,
      trtColName_RCT = "Tr",
      outcomeColName_RCT = "Yobs",
      DataTarget = dat_tgt[1:50,],
      covariateColName_TargetData = covs,
      num_trees = 2, seed = 123
      # removed max_depth argument which was causing error
    )
  })

  expect_s3_class(res, "characterizing_underrep")
  expect_s3_class(res$root, "ROOT")
})

test_that("S3 methods for characterizing_underrep work", {
  suppressWarnings({
    res <- characterizing_underrep(
      DataRCT = dat_rct[1:50,], covs, "Tr", "Yobs",
      DataTarget = dat_tgt[1:50,], covs,
      num_trees = 2, seed = 123
    )
  })

  expect_output(print(res), "characterizing_underrep object")
  expect_output(summary(res), "ROOT summary")

  if (!is.null(res$root$f)) {
    expect_silent(grDevices::pdf(file = NULL))
    plot(res)
    grDevices::dev.off()
  }
})

test_that("S3 methods for ROOT work", {
  # Use dat_2s to ensure stable fit
  suppressWarnings({
    res_2s <- ROOT(dat_2s, "Yobs", "Tr", "S", num_trees = 3, seed = 123)
  })

  # Internals check
  expect_false(.root_is_single_sample(res_2s))

  # Methods check
  expect_output(print(res_2s), "TATE.+unweighted")
  expect_output(summary(res_2s), "Rashomon size")
})

test_that("characterize_tree() validates length and binarity of w", {
  skip_if_not_installed("rpart")

  X <- data.frame(X1 = rnorm(30), X2 = runif(30))

  # 1) Length mismatch should error with the length message
  expect_error(
    characterize_tree(X, w = 0:2),
    "Length of `w` must equal the number of rows in `X`"
  )

  # 2) Same length but 3 classes should trip the 'exactly two classes' check
  w3 <- c(0,1,2)[(seq_len(nrow(X)) - 1) %% 3 + 1]   # length 30, classes {0,1,2}
  expect_error(
    characterize_tree(X, w = w3),
    "exactly two classes"
  )

  # 3) Happy path
  set.seed(1)
  w <- sample(c(0,1), size = nrow(X), replace = TRUE)
  fit <- characterize_tree(X, w = w, max_depth = 2)
  expect_s3_class(fit, "rpart")
})

test_that("characterize_tree fits a shallow rpart classifier", {
  skip_if_not_installed("rpart")
  # Simple separable data
  X <- data.frame(x1 = c(0,0,1,1), x2 = c(0,1,0,1))
  w <- c(0,0,1,1)
  fit <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(fit, "rpart")
  expect_true(all(c("frame","where") %in% names(fit)))
})

test_that("characterize_tree errors when w is not binary", {
  X <- data.frame(x1 = rnorm(10))
  w <- rep(2, 10)  # single level
  expect_error(characterize_tree(X, w), "exactly two classes")
})

test_that(".root_is_single_sample respects flag and lX fallback", {
  obj_flag_true  <- list(single_sample_mode = TRUE)
  obj_flag_false <- list(single_sample_mode = FALSE)
  expect_true(ROOT:::.root_is_single_sample(obj_flag_true))
  expect_false(ROOT:::.root_is_single_sample(obj_flag_false))

  # Fallback: all lX NA -> TRUE; any non-NA -> FALSE
  Dna  <- data.frame(lX = c(NA_real_, NA_real_))
  Dmix <- data.frame(lX = c(NA_real_,  0.5))
  obj_na  <- list(D_forest = Dna)
  obj_mix <- list(D_forest = Dmix)
  expect_true(ROOT:::.root_is_single_sample(obj_na))
  expect_false(ROOT:::.root_is_single_sample(obj_mix))
})

test_that(".root_covariate_names excludes v,vsq,S,lX and w_tree_*", {
  Df <- data.frame(
    X1 = 1, X2 = 2, v = 3, vsq = 4, S = 1, lX = NA_real_,
    w_tree_1 = 0L, w_tree_2 = 1L
  )
  obj <- list(D_forest = Df)
  covs <- ROOT:::.root_covariate_names(obj)
  expect_equal(sort(covs), c("X1", "X2"))
})

test_that(".root_baseline_loss uses sqrt(sum(vsq)/n^2) with na.rm", {
  Df <- data.frame(vsq = c(1, 4, NA_real_), X1 = 1:3)
  obj <- list(D_forest = Df)
  got <- ROOT:::.root_baseline_loss(obj)
  expect_equal(got, sqrt(sum(c(1,4), na.rm = TRUE) / (nrow(Df)^2)))
})

test_that(".root_selected_objectives handles empty and names outputs", {
  # Empty forest
  obj_empty <- list(w_forest = list(), rashomon_set = integer(0))
  expect_length(ROOT:::.root_selected_objectives(obj_empty), 0)

  # Some trees missing the field -> NA, but rashomon_set empty -> numeric(0)
  obj_none <- list(
    w_forest = list(list(), list(`local objective` = 0.25)),
    rashomon_set = integer(0)
  )
  expect_length(ROOT:::.root_selected_objectives(obj_none), 0)

  # Select 2nd only; ensure name = "w_tree_2"
  obj_sel <- list(
    w_forest = list(list(`local objective` = 0.5),
                    list(`local objective` = 0.25),
                    list(`local objective` = 0.75)),
    rashomon_set = 2L
  )
  out <- ROOT:::.root_selected_objectives(obj_sel)
  expect_equal(unname(out), 0.25)
  expect_equal(names(out), "w_tree_2")
})

test_that("summary.ROOT prints precomputed estimates (with/without SE)", {
  Df <- data.frame(v = rnorm(3), vsq = 1, S = c(1L,1L,1L), lX = NA_real_,
                   X1 = 1:3, w_tree_1 = 1L)

  obj_with_se <- fake_root(
    D_forest = Df,
    D_rash   = data.frame(w_opt = c(1L,1L,0L)),
    rashomon_set = 1L,
    estimate = list(
      estimand_unweighted = "ATE in RCT",
      value_unweighted    = 0.1,
      se_unweighted       = 0.2,
      estimand_weighted   = "WATE",
      value_weighted      = 0.3,
      se_weighted         = 0.4,
      se_weighted_note    = "ok"
    )
  )
  expect_output(summary(obj_with_se), "ATE in RCT \\(unweighted\\).*SE = 0.200000")
  expect_output(summary(obj_with_se), "WATE \\(weighted\\).*SE = 0.400000")
  expect_output(summary(obj_with_se), "Note: ok")

  obj_no_se <- obj_with_se
  obj_no_se$estimate$se_weighted <- NA_real_
  expect_output(summary(obj_no_se), "WATE \\(weighted\\)   = 0.300000")
})

test_that("summary.ROOT fallback paths: binary vs nonbinary and no-kept", {
  Df <- data.frame(
    v = c(1, 2, 3, 4), vsq = 1, S = c(1L,1L,1L,0L), lX = 0,
    X1 = 0, w_tree_1 = 1L
  )

  # (A) Binary weights with some kept
  objA <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(1L,0L,1L,NA)))
  expect_output(summary(objA), "WTATE \\(weighted\\)")

  # (B) Binary weights but no kept
  objB <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(0L,0L,0L,NA)))
  expect_output(summary(objB), "NA \\(no kept observations\\)")
  expect_output(summary(objB), "SE omitted because no kept observations")

  # (C) Non-binary weights
  objC <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(0.2, 0.8, 0.6, NA)))
  expect_output(summary(objC), "WTATE \\(weighted\\)   = ")
  expect_output(summary(objC), "SE omitted: non-binary w_opt detected")

  # (D) Custom objective -> advisory note
  custom_obj <- function(D) mean(D$vsq)
  objD <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(1L,0L,1L,NA)),
                    global_objective_fn = custom_obj)
  expect_output(summary(objD), "custom global_objective_fn")
})

test_that("summary.ROOT diagnostics block prints baseline/selected loss and kept %", {
  Df <- data.frame(v = rnorm(5), vsq = c(1,2,3,4,5), S = c(1,1,1,1,0), lX = 0,
                   X1 = 0, w_tree_1 = 1L, w_tree_2 = 0L)
  obj <- fake_root(
    D_forest = Df,
    D_rash   = data.frame(w_opt = c(1,0,1,0,NA)),
    rashomon_set = c(1L, 2L),
    w_forest = list(
      list(`local objective` = 0.33),
      list(`local objective` = 0.10)
    )
  )
  out <- capture.output(summary(obj))
  expect_true(any(grepl("^  Baseline loss:\\s+\\d", out)))
  expect_true(any(grepl("^  Selected loss:\\s+min/median = ", out)))
  expect_true(any(grepl("^  Kept \\(w_opt=1\\):", out)))
})

test_that("print.ROOT covers weighted/unweighted branches and kept%", {
  Df <- data.frame(v = c(5,6,7), vsq = 1, S = c(1L,1L,0L), lX = 0,
                   X1 = 0, w_tree_1 = 1L)
  obj <- fake_root(D_forest = Df, D_rash = data.frame(w_opt = c(1L, 0L, NA)),
                   rashomon_set = 1L)
  expect_output(print(obj), "TATE \\(unweighted\\)|ATE in RCT \\(unweighted\\)")
  expect_output(print(obj), "Kept \\(w_opt=1\\):")
})

test_that("plot.ROOT handles x$f NULL and plots when present", {
  # NULL branch still returns invisibly with a message
  obj_null <- fake_root(D_forest = data.frame(), f = NULL)
  expect_invisible(plot(obj_null))  # message is fine

  skip_if_not_installed("rpart")
  skip_if_not_installed("rpart.plot")

  # Use a classification tree so default extra=109 is valid
  set.seed(1)
  df  <- data.frame(y = factor(sample(c(0, 1), 60, TRUE)),
                    X1 = runif(60), X2 = runif(60))
  fit <- rpart::rpart(y ~ X1 + X2, data = df, method = "class")
  obj <- fake_root(D_forest = df, f = fit)

  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off(), add = TRUE)

  # We only care that it doesn't error; suppress incidental layout warnings.
  expect_error(suppressWarnings(plot(obj)), NA)
})

test_that("summary.characterizing_underrep() rejects wrong class", {
  # Call method directly to exercise the guard
  expect_error(
    summary.characterizing_underrep(list()),
    "Not a characterizing_underrep object\\."
  )
})

test_that("summary.characterizing_underrep() prints when leaf_summary is NULL", {
  # Minimal object: root=NULL is OK (summary(NULL) is harmless)
  obj <- structure(
    list(
      root = NULL,
      combined = data.frame(),   # not used here, but present in real objects
      leaf_summary = NULL
    ),
    class = "characterizing_underrep"
  )

  # Should print a 'none (no summarized tree)' line and return invisibly
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
    pct         = (1:12)/100,
    label       = rep("Represented (keep, w=1)", 12),
    stringsAsFactors = FALSE
  )
  obj <- structure(
    list(root = NULL, combined = data.frame(), leaf_summary = leaf_summary),
    class = "characterizing_underrep"
  )

  out <- capture.output(invisible(summary(obj)))

  # 1) Header line appears somewhere (don’t anchor to line start)
  expect_true(any(grepl("Leaf summary:\\s+\\d+\\s+terminal nodes", out)))

  # 2) Exactly 10 preview lines printed (look for the 'rule' values anywhere in a row)
  #    We search for rows that contain 'X0 <=' followed by a number.
  shown_rules <- grep("X0\\s*<=\\s*\\d+", out, value = TRUE)
  expect_equal(length(shown_rules), 10)

  # 3) Ellipsis cue printed because there were >10 rows
  expect_true(any(grepl("\\.\\.\\.", out)))
})

test_that("print.characterizing_underrep() rejects wrong class", {
  expect_error(
    print.characterizing_underrep(list()),
    "Not a characterizing_underrep object\\."
  )
})

test_that("print.characterizing_underrep() prints header and delegates to print(root)", {
  obj <- structure(
    list(root = NULL, combined = data.frame(), leaf_summary = NULL),
    class = "characterizing_underrep"
  )
  # With root=NULL, base::print(NULL) is fine; just ensure the header appears
  expect_output(
    invisible(print(obj)),
    "characterizing_underrep object\\n\\s+--- ROOT summary ---"
  )
})

test_that("plot.characterizing_underrep() safely exits when no summary tree", {
  skip_if_not_installed("rpart")    # keep parity with the rest of the suite

  obj <- structure(
    list(root = list(f = NULL)),    # triggers the early guard
    class = "characterizing_underrep")
  # Expect a message and invisible(NULL)
  expect_message(
    expect_invisible(plot(obj)),
    "No summary tree available to plot"
  )
})

test_that("plot.characterizing_underrep() uses f$ylevels when attr is NULL (classification)", {
  skip_if_not_installed("rpart"); skip_if_not_installed("rpart.plot")

  set.seed(1)
  n  <- 30
  X0 <- runif(n); X1 <- runif(n)
  w  <- factor(ifelse(X0 + X1 > 1, 1, 0), levels = c(0, 1))  # levels "0","1"

  fit <- rpart::rpart(w ~ X0 + X1,
                      method  = "class",
                      control = rpart::rpart.control(cp = 0, minsplit = 2, maxdepth = 2))

  # Normalize levels layout to hit the intended branch:
  attr(fit, "ylevels") <- NULL
  fit$ylevels <- levels(w)

  obj <- structure(list(root = list(f = fit)), class = "characterizing_underrep")

  tf <- tempfile(fileext = ".pdf"); grDevices::pdf(tf)
  expect_silent(plot(obj))  # exercises the f$ylevels fallback
  grDevices::dev.off()
})

test_that("plot.characterizing_underrep() falls back to regression threshold when no levels", {
  skip_if_not_installed("rpart"); skip_if_not_installed("rpart.plot")

  set.seed(2)
  n  <- 30
  X0 <- runif(n); X1 <- runif(n)
  y  <- X0 - X1 + rnorm(n, sd = 0.1)

  fit <- rpart::rpart(y ~ X0 + X1,
                      method  = "anova",
                      control = rpart::rpart.control(cp = 0, minsplit = 2, maxdepth = 2))

  # Ensure there are truly no class levels:
  attr(fit, "ylevels") <- NULL
  fit$ylevels <- NULL

  obj <- structure(list(root = list(f = fit)), class = "characterizing_underrep")

  tf <- tempfile(fileext = ".pdf"); grDevices::pdf(tf)
  expect_silent(plot(obj))  # hits the regression threshold path
  grDevices::dev.off()
})

test_that("internal helpers behave as expected", {
  # Minimal ROOT-like scaffold
  mk_root <- function(D_forest, D_rash = data.frame(), w_forest = NULL,
                      rashomon_set = integer(), testing_data = NULL,
                      global_objective_fn = objective_default,
                      single_sample_mode = NULL, f = NULL, estimate = NULL) {
    structure(list(
      D_forest = D_forest,
      D_rash   = D_rash,
      w_forest = w_forest,
      rashomon_set = rashomon_set,
      testing_data = if (is.null(testing_data)) D_forest else testing_data,
      global_objective_fn = global_objective_fn,
      single_sample_mode = single_sample_mode,
      f = f,
      estimate = estimate
    ), class = c("ROOT", "list"))
  }

  set.seed(1)
  n <- 8
  Df <- data.frame(
    x1  = rnorm(n), x2 = rbinom(n, 1, 0.5),
    v   = rnorm(n), vsq = abs(rnorm(n)),
    S   = sample(c(0L, 1L), n, TRUE),
    lX  = NA_real_,
    w_tree_1 = rbinom(n, 1, 0.5),
    stringsAsFactors = FALSE
  )

  # .root_is_single_sample: explicit flag first, then inferred via lX
  r1 <- mk_root(Df, single_sample_mode = TRUE)
  r2 <- mk_root(Df, single_sample_mode = FALSE)
  r3 <- mk_root(transform(Df, lX = NA_real_))  # all-NA lX -> TRUE
  r4 <- mk_root(transform(Df, lX = c(NA, 1, rep(NA, n - 2))))  # not all NA -> FALSE

  expect_true( ROOT:::.root_is_single_sample(r1) )
  expect_false(ROOT:::.root_is_single_sample(r2))
  expect_true( ROOT:::.root_is_single_sample(r3) )
  expect_false(ROOT:::.root_is_single_sample(r4))

  # .root_covariate_names excludes internals and w_tree_*
  covs <- ROOT:::.root_covariate_names(r1)
  expect_setequal(covs, c("x1", "x2"))

  # .root_baseline_loss uses sqrt(sum(vsq)/n^2)
  expect_equal(
    ROOT:::.root_baseline_loss(r1),
    sqrt(sum(Df$vsq) / (n^2)),
    tolerance = 1e-12
  )

  # .root_selected_objectives: empty, and with names for selected
  expect_equal(ROOT:::.root_selected_objectives(r1), numeric(0))
  w_forest <- list(
    list("local objective" = 0.2),
    list("local objective" = 0.1),
    list("local objective" = NA_real_)
  )
  r5 <- mk_root(Df, w_forest = w_forest, rashomon_set = c(2L, 1L))
  so <- ROOT:::.root_selected_objectives(r5)
  expect_equal(unname(so), c(0.1, 0.2))
  expect_identical(names(so), c("w_tree_2", "w_tree_1"))
})

test_that("summary.ROOT: guards and precomputed-estimate printing", {
  # Not a ROOT object -> error
  expect_error(summary.ROOT(list()), "Not a ROOT object.")

  n <- 6
  Df <- data.frame(
    x1 = rnorm(n), v = rnorm(n), vsq = abs(rnorm(n)),
    S = rep(1L, n), lX = NA_real_, w_tree_1 = rbinom(n, 1, 0.5)
  )
  # Precomputed estimate with SE present
  est1 <- list(
    estimand_unweighted = "TATE",
    value_unweighted    = 0.12,
    se_unweighted       = 0.34,
    estimand_weighted   = "WTATE",
    value_weighted      = 0.56,
    se_weighted         = 0.78,
    se_weighted_note    = "custom-note-here"
  )
  obj1 <- structure(list(
    D_forest = Df,
    D_rash   = transform(Df[, 0], w_opt = rbinom(n, 1, 0.5)),
    rashomon_set = integer(), w_forest = NULL,
    testing_data = Df, f = NULL,
    estimate = est1, single_sample_mode = FALSE,
    global_objective_fn = objective_default
  ), class = c("ROOT", "list"))

  expect_output(
    summary(obj1),
    paste0(
      "TATE \\(unweighted\\) = 0\\.120000, SE = 0\\.340000\\n",
      "WTATE \\(weighted\\)   = 0\\.560000, SE = 0\\.780000\\n",
      "\\s+Note: custom-note-here"
    )
  )

  # Precomputed estimate with weighted SE omitted (NA) but note printed
  est2 <- est1; est2$se_weighted <- NA_real_; est2$se_weighted_note <- "note-NA-path"
  obj2 <- obj1; obj2$estimate <- est2
  expect_output(
    summary(obj2),
    paste0("WTATE \\(weighted\\)\\s+= 0\\.560000\\n\\s+Note: note-NA-path")
  )
})

test_that("summary.ROOT: fallback recompute branches incl. custom-objective notes", {
  set.seed(2)
  n <- 10
  # Two-sample (S==1 analysis), binary w with kept > 0
  Df <- data.frame(
    x1 = rnorm(n), v = rnorm(n), vsq = abs(rnorm(n)),
    S = c(rep(1L, n-2), 0L, 0L), lX = 1,  # lX non-NA => two-sample path
    w_tree_1 = rbinom(n, 1, 0.5)
  )
  w_opt <- rbinom(n, 1, 1/2)
  # Make sure at least some kept among S==1
  if (!any(w_opt[Df$S == 1L] == 1L)) w_opt[which(Df$S == 1L)[1]] <- 1L

  # Custom objective -> extra paragraph appended
  g_custom <- function(D) 999

  obj_bin <- structure(list(
    D_forest = Df,
    D_rash   = transform(Df[, 0], w_opt = w_opt),
    rashomon_set = 1L, w_forest = list(list("local objective" = 0.1)),
    testing_data = Df,
    f = NULL, estimate = NULL, single_sample_mode = NULL,
    global_objective_fn = g_custom
  ), class = c("ROOT", "list"))

  out <- testthat::capture_output(summary(obj_bin))
  expect_match(out, "TATE \\(unweighted\\) = ", perl = TRUE)
  expect_match(out, "WTATE \\(weighted\\)   = .* SE = ", perl = TRUE)
  expect_match(out, "Calculation of SE for WTATE uses sqrt\\(", perl = TRUE)
  expect_match(out, "You supplied a custom global_objective_fn; please verify this SE matches your", fixed = TRUE)

  # Binary w but NO kept among S==1 -> NA line and explicit empty-note
  w_none <- replace(w_opt, Df$S == 1L, 0L)
  obj_none <- obj_bin; obj_none$D_rash$w_opt <- w_none
  out2 <- testthat::capture_output(summary(obj_none))
  expect_match(out2, "WTATE \\(weighted\\)\\s+= NA \\(no kept observations\\)")
  expect_match(out2, "SE omitted because no kept observations", perl = TRUE)

  # Non-binary w path with custom-objective addendum
  w_nonbin <- ifelse(Df$S == 1L, runif(n), NA_real_)
  obj_nb <- obj_bin; obj_nb$D_rash$w_opt <- w_nonbin
  out3 <- testthat::capture_output(summary(obj_nb))
  expect_match(out3, "WTATE \\(weighted\\)\\s+= ", perl = TRUE)
  expect_match(out3, "SE omitted: non-binary w_opt detected", fixed = TRUE)
  expect_match(out3, "please ensure your variance method matches that estimand", fixed = TRUE)
})

test_that("print.ROOT mirrors summary logic and shows notes", {
  n <- 7
  Df <- data.frame(
    x1 = rnorm(n), v = rnorm(n), vsq = abs(rnorm(n)),
    S = rep(1L, n), lX = NA_real_, w_tree_1 = rbinom(n, 1, 0.5)
  )
  # Precomputed estimate with weighted SE omitted and a note
  est <- list(
    estimand_unweighted = "ATE in RCT",
    value_unweighted    = 1.23,
    se_unweighted       = 0.11,
    estimand_weighted   = "WATE",
    value_weighted      = -0.5,
    se_weighted         = NA_real_,
    se_weighted_note    = "non-binary weights -> omit SE"
  )
  obj <- structure(list(
    D_forest = Df,
    D_rash   = transform(Df[, 0], w_opt = rbinom(n, 1, 0.5)),
    rashomon_set = integer(), w_forest = NULL,
    testing_data = Df, f = NULL,
    estimate = est, single_sample_mode = TRUE,
    global_objective_fn = objective_default
  ), class = c("ROOT", "list"))

  expect_output(
    print(obj),
    paste0(
      "ATE in RCT \\(unweighted\\) = 1\\.230000, SE = 0\\.110000\\n",
      "WATE \\(weighted\\)\\s+= -0\\.500000\\n\\s+Note: non-binary weights -> omit SE"
    )
  )

  # Fallback path with binary kept > 0 and custom objective -> addendum is printed
  obj2 <- obj
  obj2$estimate <- NULL
  obj2$single_sample_mode <- FALSE
  obj2$D_forest$lX <- 1  # force two-sample
  obj2$global_objective_fn <- function(D) 123
  obj2$D_rash$w_opt <- as.integer(runif(n) > 0.4)
  if (!any(obj2$D_rash$w_opt == 1L)) obj2$D_rash$w_opt[1] <- 1L

  out <- testthat::capture_output(print(obj2))
  expect_match(out, "WTATE \\(weighted\\)   = .* SE = ", perl = TRUE)
  expect_match(out, "You supplied a custom global_objective_fn; please verify this SE matches your", fixed = TRUE)
})

test_that("plot.ROOT works with multiclass rpart when box.palette is valid", {
  skip_if_not_installed("rpart.plot")
  skip_if_not_installed("rpart")

  f <- rpart::rpart(
    Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
    data = iris
  )

  fake_root <- list(
    D_forest     = data.frame(v = 0, vsq = 0, S = 1L, lX = NA_real_),
    D_rash       = data.frame(),
    rashomon_set = integer(),
    f            = f,
    testing_data = data.frame()
  )
  class(fake_root) <- "ROOT"

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp); on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)

  # Either auto palette…
  expect_warning(plot(fake_root, box.palette = "auto"), regexp = NA)
  # …or explicit list with 3 entries (for 3 classes)
  expect_warning(plot(fake_root, box.palette = list("pink", "lightblue", "lightgray")), regexp = NA)
  # (If available) expect_no_warning(...)
})
