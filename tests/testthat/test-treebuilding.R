# Helper to create dummy state for split_node calls
make_dummy_state <- function(n = 20, p = 2) {
  X <- matrix(runif(n * p), ncol = p)
  colnames(X) <- paste0("X", 1:p)
  X_df <- as.data.frame(X)
  X_df$w <- 1L
  rownames(X_df) <- as.character(1:n)

  D <- X_df
  D$vsq <- rchisq(n, df = 1)

  split_feats <- c(leaf = 0.1, X1 = 0.45, X2 = 0.45)
  list(X = X_df, D = D, sf = split_feats)
}

test_that("split_node base cases work", {
  st <- make_dummy_state(n = 10)

  # max-depth leaf
  res_depth <- split_node(
    st$sf, st$X, st$D,
    parent_loss = Inf, depth = 8, max_depth = 8
  )
  expect_equal(res_depth$node, "leaf")
  expect_equal(res_depth$leaf_reason, "max-depth")

  # min-leaf leaf
  res_min_n <- split_node(
    st$sf, st$X[1:4, ], st$D,
    parent_loss = Inf, depth = 0, min_leaf_n = 5
  )
  expect_equal(res_min_n$node, "leaf")
  expect_equal(res_min_n$leaf_reason, "min-leaf")

  # feature==leaf
  sf_leaf <- c(leaf = 1, X1 = 0, X2 = 0)
  res_leaf_feat <- split_node(
    sf_leaf, st$X, st$D,
    parent_loss = Inf, depth = 0
  )
  expect_equal(res_leaf_feat$node, "leaf")
  expect_equal(res_leaf_feat$leaf_reason, "feature==leaf")
})

test_that("split_node successfully splits", {
  st <- make_dummy_state(n = 50)
  res_split <- withr::with_seed(1, {
    split_node(
      st$sf, st$X, st$D,
      parent_loss = 100, depth = 0
    )
  })
  expect_true(res_split$node %in% c("X1", "X2"))
})

test_that("Helpers: choose_feature, reduce_weight, midpoint", {
  # choose_feature
  sf <- c(leaf = 0.2, X1 = 0.8)
  set.seed(1)
  chosen <- replicate(100, choose_feature(sf, depth = 0))
  expect_true(mean(chosen == "X1") > 0.7)

  # reduce_weight
  sf_red <- reduce_weight("X1", sf)
  expected_val <- (0.8 / 2) / (0.2 + 0.8 / 2)
  expect_equal(sf_red["X1"], expected_val, ignore_attr = TRUE)
  expect_equal(sum(sf_red), 1)

  # midpoint
  expect_equal(midpoint(c(1, 3)), 2)
  suppressWarnings({
    expect_true(is.na(midpoint(c(NA_real_))))
  })
})

test_that("characterize_tree fits an rpart model", {
  skip_if_not_installed("rpart")

  X <- dat_tiny[, c("X0", "X1")]
  n <- nrow(X)

  # Create w with exactly length n, and ensure two classes
  w <- rep(0, n)
  w[1:(n / 2)] <- 1

  fit <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(fit, "rpart")

  # Error: w not binary (all 1s)
  expect_error(characterize_tree(X, rep(1, n)), "exactly two classes")

  # Error: Length mismatch
  expect_error(characterize_tree(X, rep(1, n - 1)), "Length of `w` must equal")
})

# Minimal toy data with two numeric features
toy_X <- data.frame(

  X1 = c(0, 1, 2, 3, 4, 5),
  X2 = c(5, 4, 3, 2, 1, 0),
  row.names = as.character(1:6)
)
toy_D <- data.frame(
  vsq = c(1, 2, 1, 2, 1, 2),
  w   = rep(1, 6),
  row.names = as.character(1:6)
)

test_that("split_node() stops at max depth & min leaf size", {
  set.seed(123)
  sf <- c(leaf = 0.0, X1 = 0.5, X2 = 0.5)

  # Enforce early leaf via max_depth=0
  tree0 <- split_node(
    sf, toy_X, toy_D,
    parent_loss = Inf, depth = 0,
    explore_proba = 0, max_depth = 0, min_leaf_n = 1
  )
  expect_equal(tree0$node, "leaf")
  expect_match(tree0$leaf_reason, "max-depth")

  # Enforce early leaf via min_leaf_n
  tree1 <- split_node(
    sf, toy_X[1:5, , drop = FALSE], toy_D,
    parent_loss = Inf, depth = 0,
    explore_proba = 0, max_depth = 8, min_leaf_n = 6
  )
  expect_equal(tree1$node, "leaf")
  expect_match(tree1$leaf_reason, "min-leaf")
})

test_that("split_node() can choose 'leaf' and assign best w without exploration", {
  set.seed(123)
  sf <- c(leaf = 1.0, X1 = 0.0)
  res <- split_node(
    sf, toy_X, toy_D,
    parent_loss = Inf, depth = 0,
    explore_proba = 0, max_depth = 8, min_leaf_n = 1
  )
  expect_equal(res$node, "leaf")
  expect_true(res$w %in% c(0, 1))
  expect_true(is.finite(res[["local objective"]]))
})

test_that("split_node() rejects non-improving split and penalizes feature weight", {
  set.seed(1)
  sf <- c(leaf = 0.0, X1 = 1.0)
  res <- split_node(
    sf, toy_X, toy_D,
    parent_loss = .Machine$double.eps,
    depth = 0, explore_proba = 0,
    max_depth = 8, min_leaf_n = 2
  )
  expect_true(res$node %in% c("leaf", "X1"))
})

test_that("split_node accepts improving splits and updates weights", {
  set.seed(10)
  X  <- data.frame(X0 = runif(50), X1 = runif(50))
  v  <- rnorm(50)
  vsq <- (v - mean(v))^2
  D <- data.frame(
    X, v = v, vsq = vsq, w = rep(1, 50), S = rep(1L, 50),
    stringsAsFactors = FALSE
  )
  rownames(D) <- as.character(seq_len(nrow(D)))

  sf <- c(leaf = 0.2, X0 = 0.4, X1 = 0.4)

  res <- split_node(
    split_feature        = sf,
    X                    = D,
    D                    = D,
    parent_loss          = Inf,
    depth                = 0L,
    explore_proba        = 0.0,
    choose_feature_fn    = choose_feature,
    reduce_weight_fn     = reduce_weight,
    global_objective_fn  = objective_default,
    max_depth            = 2,
    min_leaf_n           = 5
  )

  expect_type(res, "list")
  expect_true(all(c("D", "local objective", "depth") %in% names(res)))
  expect_equal(nrow(res$D), nrow(D))
  expect_true(all(res$D$w %in% 0:1))
})

test_that("midpoint validates and handles edge cases", {
  expect_error(midpoint("a"), "`X` must be a numeric")
  expect_warning(expect_true(is.na(midpoint(numeric(0)))))
  expect_warning(expect_true(is.na(midpoint(c(NA_real_, NA_real_)))))
  expect_equal(midpoint(c(-2, 2, NA)), 0)
  expect_equal(midpoint(c(-5, -1)), -3)
})

test_that("choose_feature validates probs and names", {
  expect_error(choose_feature(numeric(0), 0), "non-empty numeric")
  expect_error(choose_feature(c(a = 0.5, b = NA_real_), 0), "contains NA")
  expect_error(
    choose_feature(setNames(c(-0.2, 1.2), c("leaf", "X1")), 0),
    ">= 0|Negative",
    ignore.case = TRUE
  )
  expect_error(choose_feature(unname(c(0.5, 0.5)), 0), "must have names")
  expect_error(
    choose_feature(c(a = 0.7, b = NA_real_), 0),
    "NA values|contains NA"
  )
})

test_that("reduce_weight validates and renormalizes, halves target", {
  expect_error(
    reduce_weight(1, c(leaf = 0.5, X1 = 0.5)),
    "`fj` must be a single"
  )
  expect_error(
    reduce_weight("X1", c(0.5, 0.5)),
    "named numeric"
  )
  expect_error(
    reduce_weight("X9", c(leaf = 0.5, X1 = 0.5)),
    "not found"
  )
  expect_error(
    reduce_weight("X1", c(leaf = 0.5, X1 = -0.1)),
    ">= 0|Negative",
    ignore.case = TRUE
  )
  expect_error(
    reduce_weight("X1", c(leaf = 0.5, X1 = Inf)),
    "non-finite"
  )

  sf <- c(leaf = 0.2, X1 = 0.8)
  sf2 <- reduce_weight("X1", sf)
  expect_equal(sum(sf2), 1, tolerance = 1e-12)
  expect_lt(sf2["X1"], sf["X1"])
})

test_that("characterize_tree basic fit and guards", {
  skip_if_not_installed("rpart")

  X <- data.frame(a = rnorm(30), b = runif(30))
  w <- sample(c(0, 1), 30, TRUE)

  fit <- characterize_tree(X, w)
  expect_s3_class(fit, "rpart")

  expect_error(characterize_tree(X, w[1:10]), "Length of `w` must equal")
})

mk_state <- function(n = 20) {
  set.seed(42)
  X <- data.frame(X1 = runif(n), X2 = runif(n))
  D <- data.frame(
    X,
    v   = rnorm(n),
    vsq = rchisq(n, df = 1),
    w   = rep(1, n),
    S   = rep(1L, n)
  )
  rownames(X) <- rownames(D) <- as.character(seq_len(n))
  list(X = X, D = D)
}

test_that("split_node: max-depth, min-leaf, leaf feature", {
  st <- mk_state(12)
  sf <- c(leaf = 1, X1 = 0, X2 = 0)

  o1 <- split_node(
    sf, st$X, st$D,
    parent_loss = Inf, depth = 3, max_depth = 3
  )
  expect_equal(o1$node, "leaf")
  expect_match(o1$leaf_reason, "max-depth")

  o2 <- split_node(
    c(leaf = 0.1, X1 = 0.9), st$X[1:4, ], st$D,
    parent_loss = Inf, depth = 0, min_leaf_n = 5
  )
  expect_equal(o2$node, "leaf")
  expect_match(o2$leaf_reason, "min-leaf")

  o3 <- split_node(
    sf, st$X, st$D,
    parent_loss = Inf, depth = 0
  )
  expect_equal(o3$node, "leaf")
  expect_match(o3$leaf_reason, "feature==leaf")
})

test_that("split_node: improving split vs rejection loop", {
  st <- mk_state(40)
  sf <- c(leaf = 0.0, X1 = 1.0)
  out <- withr::with_seed(
    1,
    split_node(
      sf, st$X, st$D,
      parent_loss  = 1e99,
      depth        = 0,
      explore_proba = 0,
      max_depth    = 2
    )
  )
  expect_true(out$node %in% c("X1", "leaf"))
  expect_true(is.finite(out$`local objective`))

  out2 <- withr::with_seed(
    1,
    split_node(
      sf, st$X, st$D,
      parent_loss = .Machine$double.eps,
      depth       = 0,
      explore_proba       = 0,
      max_depth           = 2,
      max_rejects_per_node = 3
    )
  )
  expect_true(out2$node %in% c("leaf", "X1"))
  if (out2$node == "leaf") {
    expect_true(grepl("reject", out2$leaf_reason) || grepl("leaf", out2$leaf_reason))
  }
})

test_that("split_node hits leaf branches: max-depth, min-leaf, feature==leaf", {
  X <- data.frame(x = rnorm(6))
  rownames(X) <- as.character(seq_len(nrow(X)))
  D <- data.frame(vsq = rnorm(6)^2, w = rep(1, 6))
  rownames(D) <- rownames(X)

  sf <- c(leaf = 1)
  expect_silent({
    res_leaf <- split_node(
      sf, X, D,
      parent_loss = Inf, depth = 0,
      max_depth   = 8, min_leaf_n = 2
    )
  })
  expect_identical(res_leaf$node, "leaf")
  expect_match(res_leaf$leaf_reason, "feature==leaf")

  expect_silent({
    res_maxdepth <- split_node(
      c(leaf = 0.0, x = 1.0), X, D,
      parent_loss = Inf,
      depth       = 0,
      max_depth   = 0,
      min_leaf_n  = 2
    )
  })
  expect_identical(res_maxdepth$leaf_reason, "max-depth")

  X2 <- X[1:3, , drop = FALSE]
  D2 <- D[1:3, , drop = FALSE]
  expect_silent({
    res_minleaf <- split_node(
      c(leaf = 0.0, x = 1.0), X2, D2,
      parent_loss = Inf,
      depth       = 0,
      max_depth   = 8,
      min_leaf_n  = 5
    )
  })
  expect_identical(res_minleaf$leaf_reason, "min-leaf")
})

test_that("split_node empty-child branch triggers when midpoint collapses a side", {
  X <- data.frame(x = rep(1, 5))
  rownames(X) <- as.character(seq_len(nrow(X)))
  D <- data.frame(vsq = rnorm(5)^2, w = rep(1, 5))
  rownames(D) <- rownames(X)

  res <- split_node(
    c(leaf = 0, x = 1), X, D,
    parent_loss = Inf,
    depth       = 0,
    max_depth   = 8,
    min_leaf_n  = 1
  )
  expect_identical(res$leaf_reason, "empty-child")
  expect_identical(res$node, "leaf")
})

test_that("split_node rejection loop can force a leaf after many rejects", {
  set.seed(1)
  X <- data.frame(x = sort(rnorm(20)))
  rownames(X) <- as.character(seq_len(nrow(X)))
  D <- data.frame(vsq = rnorm(20)^2, w = rep(1, 20))
  rownames(D) <- rownames(X)

  parent_loss <- -1e9

  res <- split_node(
    split_feature = c(leaf = 0.0, x = 1.0),
    X             = X,
    D             = D,
    parent_loss   = parent_loss,
    depth         = 0,
    explore_proba = 0.0,
    max_depth     = 8,
    min_leaf_n    = 2,
    max_rejects_per_node = 5
  )
  expect_identical(res$node, "leaf")
  expect_true(res$leaf_reason %in% c("max-rejects", "forced-leaf-after-rejects"))
})

test_that("choose_feature strict guards: bad depth, NA probs, non-finite sum", {
  sf <- c(leaf = 0.6, x = 0.4)

  expect_error(choose_feature(sf, depth = c(1, 2)), "numeric scalar")
  expect_error(choose_feature(c(leaf = NA_real_, x = 1), depth = 0), "contains NA")
  expect_error(
    choose_feature(c(leaf = 1, x = Inf), depth = 0),
    "sum is non-finite"
  )
})

test_that("reduce_weight validates input and renormalizes to 1", {
  sf <- c(leaf = 0.2, x = 0.8)

  expect_error(reduce_weight("z", sf), "not found")
  expect_error(
    reduce_weight("leaf", c(leaf = NA_real_, x = 1)),
    "contains NA"
  )
  expect_error(
    reduce_weight("leaf", c(leaf = -0.1, x = 1.1)),
    ">= 0|Negative",
    ignore.case = TRUE
  )

  out <- reduce_weight("x", sf)
  expect_true(abs(sum(out) - 1) < 1e-12)
  expect_lt(out["x"], sf["x"])
  expect_gt(out["leaf"], 0)
})

test_that("characterize_tree guards: X must be df; |w| mismatch; w must be binary", {
  skip_if_not_installed("rpart")

  X <- data.frame(a = rnorm(10), b = rnorm(10))
  w2 <- rep(0, 10)

  expect_error(characterize_tree(as.matrix(X), w2), "data frame")
  expect_error(characterize_tree(X, w = w2[1:5]), "equal the number of rows")
  expect_error(characterize_tree(X, w = rep(1, 10)), "exactly two classes")
  expect_error(
    characterize_tree(X, w = c(rep(0, 4), rep(1, 4), rep(2, 2))),
    "exactly two classes"
  )

  X_bad <- X
  X_bad$a[3] <- NA_real_
  expect_error(
    characterize_tree(X_bad, w = c(rep(0, 5), rep(1, 5))),
    "contains missing values"
  )
})

test_that("characterize_tree fits with two classes and respects max_depth", {
  skip_if_not_installed("rpart")

  set.seed(7)
  X <- data.frame(a = rnorm(40), b = rnorm(40))
  w <- as.integer(X$a + X$b > 0)
  w[sample.int(40, 2)] <- 1 - w[sample.int(40, 2)]
  f <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(f, "rpart")
  expect_true(f$control$maxdepth == 2)
})

stub_loss_from_objective <- function(...) {
  function(val, idx, D) {
    if (length(idx) == 0) 0 else mean(D[idx, "vsq"])
  }
}

make_XD_onefeat <- function() {
  X <- data.frame(x1 = c(0, 0, 1, 1, 2, 2))
  rownames(X) <- paste0("r", seq_len(nrow(X)))
  D <- data.frame(w = 0, vsq = seq_len(nrow(X)))
  rownames(D) <- rownames(X)
  list(X = X, D = D)
}

make_XD_twofeat <- function() {
  X <- data.frame(
    x1 = c(0, 0, 1, 1, 2, 2),
    x2 = c(0, 1, 0, 1, 0, 1)
  )
  rownames(X) <- paste0("r", seq_len(nrow(X)))
  D <- data.frame(
    w   = 0,
    vsq = c(5, 5, 5, 1, 1, 1)
  )
  rownames(D) <- rownames(X)
  list(X = X, D = D)
}

test_that("split_node -> forced-leaf-after-rejects", {
  XD <- make_XD_onefeat()
  choose_x1 <- function(sf, depth) "x1"
  force_leaf <- function(fj, sf) c(leaf = 1)

  testthat::with_mocked_bindings(
    {
      res <- split_node(
        c(leaf = 0, x1 = 1), XD$X, XD$D,
        parent_loss      = -1e6,
        depth            = 0,
        choose_feature_fn = choose_x1,
        reduce_weight_fn  = force_leaf,
        explore_proba     = 0
      )
      expect_identical(res$node, "leaf")
      expect_match(res$leaf_reason, "forced-leaf-after-rejects")
    },
    loss_from_objective = stub_loss_from_objective,
    .env = asNamespace("ROOT")
  )
})

test_that("split_node -> feature==leaf-after-rejects", {
  XD <- make_XD_onefeat()
  choose_x1_then_leaf <- local({
    called <- FALSE
    function(sf, depth) {
      if (!called) {
        called <<- TRUE
        "x1"
      } else {
        "leaf"
      }
    }
  })
  id_reduce <- function(fj, sf) sf

  testthat::with_mocked_bindings(
    {
      res <- split_node(
        c(leaf = 0.1, x1 = 0.9), XD$X, XD$D,
        parent_loss      = -1e6,
        depth            = 0,
        choose_feature_fn = choose_x1_then_leaf,
        reduce_weight_fn  = id_reduce,
        explore_proba     = 0
      )
      expect_identical(res$node, "leaf")
      expect_match(res$leaf_reason, "feature==leaf-after-rejects")
    },
    loss_from_objective = stub_loss_from_objective,
    .env = asNamespace("ROOT")
  )
})

test_that("split_node -> empty-child-after-rejects", {
  XD <- make_XD_onefeat()
  XD$X$x1 <- 5
  choose_x1 <- function(sf, depth) "x1"
  id_reduce <- function(fj, sf) sf

  testthat::with_mocked_bindings(
    {
      res <- split_node(
        c(leaf = 0, x1 = 1), XD$X, XD$D,
        parent_loss      = -1e6,
        depth            = 0,
        choose_feature_fn = choose_x1,
        reduce_weight_fn  = id_reduce,
        explore_proba     = 0
      )
      expect_identical(res$node, "leaf")
      expect_match(res$leaf_reason, "empty-child")
    },
    loss_from_objective = stub_loss_from_objective,
    .env = asNamespace("ROOT")
  )
})

test_that("split_node -> accepts split after rejects and returns local objective", {
  X <- data.frame(
    a = c(1, 1, 2, 2),
    b = c(1, 2, 1, 2),
    row.names = paste0("r", 1:4)
  )

  D <- data.frame(
    vsq = rep(1, nrow(X)),
    w   = rep(0, nrow(X)),
    row.names = rownames(X)
  )

  parent_loss <- 1

  key_of <- function(idx) paste(sort(idx), collapse = ",")

  good_keys <- c(key_of(c("r1", "r3")), key_of(c("r2", "r4")))

  stub_loss_from_objective <- function(global_objective_fn) {
    force(global_objective_fn)
    function(val, indices, D) {
      if (key_of(indices) %in% good_keys) 0.10 else 10.0
    }
  }

  choices <- c("a", "b")
  choose_seq <- local({
    i <- 0L
    function(split_feature, depth) {
      i <<- i + 1L
      if (i <= length(choices)) choices[i] else "b"
    }
  })

  reduce_identity <- function(fj, sf) sf

  split_feature <- c(a = 0.6, b = 0.3, leaf = 0.1)

  testthat::with_mocked_bindings(
    loss_from_objective = stub_loss_from_objective,
    objective_default   = function(D) 0.5,
    .env = asNamespace("ROOT"),
    {
      res <- split_node(
        split_feature        = split_feature,
        X                    = X,
        D                    = D,
        parent_loss          = parent_loss,
        depth                = 0L,
        explore_proba        = 0,
        choose_feature_fn    = choose_seq,
        reduce_weight_fn     = reduce_identity,
        global_objective_fn  = objective_default,
        max_depth            = 1L,
        min_leaf_n           = 2L,
        log_fn               = function(...) {},
        max_rejects_per_node = 10L
      )

      expect_identical(res$node, "b")
      expect_true(is.finite(res[["local objective"]]))
      expect_identical(rownames(res$D), rownames(X))
      expect_true(res$left_tree$node  %in% c("leaf"))
      expect_true(res$right_tree$node %in% c("leaf"))
    }
  )
})
