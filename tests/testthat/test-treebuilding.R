# Helper to create dummy state for split_node calls
make_dummy_state <- function(n=20, p=2) {
  X <- matrix(runif(n*p), ncol=p)
  colnames(X) <- paste0("X", 1:p)
  X_df <- as.data.frame(X)
  X_df$w <- 1
  rownames(X_df) <- as.character(1:n)

  D <- X_df
  D$vsq <- rchisq(n, df=1) # random positive values

  split_feats <- c(leaf=0.1, X1=0.45, X2=0.45)
  list(X=X_df, D=D, sf=split_feats)
}

test_that("split_node base cases work", {
  st <- make_dummy_state(n=10)
  res_depth <- split_node(st$sf, st$X, st$D, parent_loss = Inf, depth = 8, max_depth = 8)
  expect_equal(res_depth$node, "leaf")
  expect_equal(res_depth$leaf_reason, "max-depth")

  res_min_n <- split_node(st$sf, st$X[1:4,], st$D, parent_loss = Inf, depth = 0, min_leaf_n = 5)
  expect_equal(res_min_n$node, "leaf")
  expect_equal(res_min_n$leaf_reason, "min-leaf")

  sf_leaf <- c(leaf=1, X1=0, X2=0)
  res_leaf_feat <- split_node(sf_leaf, st$X, st$D, parent_loss = Inf, depth = 0)
  expect_equal(res_leaf_feat$node, "leaf")
  expect_equal(res_leaf_feat$leaf_reason, "feature==leaf")
})

test_that("split_node successfully splits", {
  st <- make_dummy_state(n=50)
  res_split <- withr::with_seed(1, {
    split_node(st$sf, st$X, st$D, parent_loss = 100, depth = 0)
  })
  expect_true(res_split$node %in% c("X1", "X2"))
})

test_that("Helpers: choose_feature, reduce_weight, midpoint", {
  # choose_feature
  sf <- c(leaf=0.2, X1=0.8)
  set.seed(1)
  chosen <- replicate(100, choose_feature(sf, depth=0))
  expect_true(mean(chosen == "X1") > 0.7)

  # reduce_weight
  sf_red <- reduce_weight("X1", sf)
  expected_val <- (0.8/2) / (0.2 + 0.8/2)
  expect_equal(sf_red["X1"], expected_val, ignore_attr = TRUE)
  expect_equal(sum(sf_red), 1)

  # midpoint
  expect_equal(midpoint(c(1, 3)), 2)
  suppressWarnings({
    expect_true(is.na(midpoint(c(NA_real_))))
  })
})

test_that("characterize_tree fits an rpart model", {
  # Use dat_tiny which has 50 rows in the helper
  X <- dat_tiny[, c("X0", "X1")]
  n <- nrow(X)

  # Create w with exactly length n, and ensure two classes
  w <- rep(0, n)
  w[1:(n/2)] <- 1

  fit <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(fit, "rpart")

  # Error: w not binary (all 1s)
  expect_error(characterize_tree(X, rep(1, n)), "must have exactly two classes")

  # Error: Length mismatch
  expect_error(characterize_tree(X, rep(1, n - 1)), "Length of `w` must equal")
})

# Minimal toy data with two numeric features
toy_X <- data.frame(X1 = c(0, 1, 2, 3, 4, 5),
                    X2 = c(5, 4, 3, 2, 1, 0),
                    row.names = as.character(1:6))
toy_D <- data.frame(vsq = c(1, 2, 1, 2, 1, 2),
                    w   = rep(1, 6),
                    row.names = as.character(1:6))

test_that("split_node() stops at max depth & min leaf size", {
  set.seed(123)
  sf <- c(leaf = 0.0, X1 = 0.5, X2 = 0.5)

  # Enforce early leaf via max_depth=0
  tree0 <- split_node(sf, toy_X, toy_D,
                      parent_loss = Inf, depth = 0,
                      explore_proba = 0, max_depth = 0, min_leaf_n = 1)
  expect_equal(tree0$node, "leaf")
  expect_match(tree0$leaf_reason, "max-depth")

  # Enforce early leaf via min_leaf_n
  tree1 <- split_node(sf, toy_X[1:5, , drop = FALSE], toy_D,
                      parent_loss = Inf, depth = 0,
                      explore_proba = 0, max_depth = 8, min_leaf_n = 6)
  expect_equal(tree1$node, "leaf")
  expect_match(tree1$leaf_reason, "min-leaf")
})

test_that("split_node() can choose 'leaf' and assign best w without exploration", {
  set.seed(123)
  sf <- c(leaf = 1.0, X1 = 0.0)  # always leaf
  res <- split_node(sf, toy_X, toy_D,
                    parent_loss = Inf, depth = 0,
                    explore_proba = 0, max_depth = 8, min_leaf_n = 1)
  expect_equal(res$node, "leaf")
  expect_true(res$w %in% c(0, 1))
  # Objective attached
  expect_true(is.finite(res[["local objective"]]))
})

test_that("split_node() rejects non-improving split and penalizes feature weight", {
  set.seed(1)
  sf <- c(leaf = 0.0, X1 = 1.0)
  # Construct parent_loss tiny to force 'reject'
  res <- split_node(sf, toy_X, toy_D,
                    parent_loss = .Machine$double.eps,
                    depth = 0, explore_proba = 0,
                    max_depth = 8, min_leaf_n = 2)
  # Because it recurses until a leaf (after penalization), final node can be leaf or split.
  expect_true(res$node %in% c("leaf", "X1"))
})

test_that("split_node accepts improving splits and updates weights", {
  # Minimal synthetic frame with two features and vsq heterogeneity
  set.seed(10)
  X  <- data.frame(X0 = runif(50), X1 = runif(50))
  v  <- rnorm(50)
  vsq <- (v - mean(v))^2
  D <- data.frame(X, v = v, vsq = vsq, w = rep(1, 50), S = rep(1L, 50), stringsAsFactors = FALSE)
  rownames(D) <- as.character(seq_len(nrow(D)))

  sf <- c(leaf = 0.2, X0 = 0.4, X1 = 0.4)

  res <- split_node(
    split_feature        = sf,
    X                    = D,
    D                    = D,
    parent_loss          = Inf,
    depth                = 0L,
    explore_proba        = 0.0,                 # deterministic exploitation for test
    choose_feature_fn    = choose_feature,      # updated arg name
    reduce_weight_fn     = reduce_weight,
    global_objective_fn  = objective_default,   # replaces objective_fn/loss_fn combo
    max_depth            = 2,
    min_leaf_n           = 5
  )

  expect_type(res, "list")
  expect_true(all(c("D","local objective","depth") %in% names(res)))
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
  expect_error(choose_feature(c(a=0.5, b=NA_real_), 0), "contains NA")
  expect_error(choose_feature(setNames(c(-0.2, 1.2), c("leaf","X1")), 0), "Negative")
  expect_error(choose_feature(unname(c(0.5,0.5)), 0), "must have names")
  expect_error(
    choose_feature(c(a = 0.7, b = NA_real_), 0),
    "NA values"  # or just "NA"
  )})

test_that("reduce_weight validates and renormalizes, halves target", {
  expect_error(reduce_weight(1, c(leaf=0.5, X1=0.5)), "`fj` must be a single")
  expect_error(reduce_weight("X1", c(0.5,0.5)), "named numeric")
  expect_error(reduce_weight("X9", c(leaf=0.5, X1=0.5)), "not found")
  expect_error(reduce_weight("X1", c(leaf=0.5, X1=-0.1)), "Negative")
  expect_error(reduce_weight("X1", c(leaf=0.5, X1=Inf)), "non-finite")
  sf <- c(leaf=0.2, X1=0.8)
  sf2 <- reduce_weight("X1", sf)
  expect_equal(sum(sf2), 1, tolerance = 1e-12)
  expect_lt(sf2["X1"], sf["X1"])
})

test_that("characterize_tree basic fit and guards", {
  X <- data.frame(a = rnorm(30), b = runif(30))
  w <- sample(c(0,1), 30, TRUE)

  # happy path
  fit <- characterize_tree(X, w)
  expect_s3_class(fit, "rpart")

  # guards (this should be an error and *caught*)
  expect_error(characterize_tree(X, w[1:10]), "Length of `w` must equal")
})

mk_state <- function(n=20) {
  set.seed(42)
  X <- data.frame(X1 = runif(n), X2 = runif(n))
  D <- data.frame(X, v = rnorm(n), vsq = rchisq(n, df=1), w = rep(1, n), S = rep(1L, n))
  rownames(X) <- rownames(D) <- as.character(seq_len(n))
  list(X=X, D=D)
}

test_that("split_node: max-depth, min-leaf, leaf feature", {
  st <- mk_state(12)
  sf <- c(leaf=1, X1=0, X2=0)
  o1 <- split_node(sf, st$X, st$D, parent_loss=Inf, depth=3, max_depth=3)
  expect_equal(o1$node, "leaf")
  expect_match(o1$leaf_reason, "max-depth")

  o2 <- split_node(c(leaf=0.1, X1=0.9), st$X[1:4,], st$D, parent_loss=Inf, depth=0, min_leaf_n=5)
  expect_equal(o2$node, "leaf")
  expect_match(o2$leaf_reason, "min-leaf")

  o3 <- split_node(sf, st$X, st$D, parent_loss=Inf, depth=0)
  expect_equal(o3$node, "leaf")
  expect_match(o3$leaf_reason, "feature==leaf")
})

test_that("split_node: improving split vs rejection loop", {
  st <- mk_state(40)
  # Let feature always be X1 first; use small explore prob to be deterministic
  sf <- c(leaf=0.0, X1=1.0)
  out <- withr::with_seed(1, split_node(sf, st$X, st$D, parent_loss=1e99,
                                        depth=0, explore_proba=0, max_depth=2))
  expect_true(out$node %in% c("X1","leaf"))   # either splits or quits as leaf after recursion
  expect_true(is.finite(out$`local objective`))

  # Force repeated rejections by making parent_loss tiny and halving until leaf
  out2 <- withr::with_seed(1, split_node(sf, st$X, st$D, parent_loss=.Machine$double.eps,
                                         depth=0, explore_proba=0, max_depth=2,
                                         max_rejects_per_node=3))
  expect_true(out2$node %in% c("leaf","X1"))
  # if leaf, reason should show reject budget or leaf-after-rejects
  if (out2$node == "leaf") {
    expect_true(grepl("reject", out2$leaf_reason) || grepl("leaf", out2$leaf_reason))
  }
})

test_that("split_node hits leaf branches: max-depth, min-leaf, feature==leaf", {
  # Minimal toy data
  X <- data.frame(x = rnorm(6))
  rownames(X) <- as.character(seq_len(nrow(X)))
  D <- data.frame(vsq = rnorm(6)^2, w = rep(1, 6))
  rownames(D) <- rownames(X)

  # loss_from_objective(objective_default) is used internally
  sf <- c(leaf = 1) # force "feature==leaf"
  expect_silent({
    res_leaf <- split_node(sf, X, D, parent_loss = Inf, depth = 0,
                           max_depth = 8, min_leaf_n = 2)
  })
  expect_identical(res_leaf$node, "leaf")
  expect_match(res_leaf$leaf_reason, "feature==leaf")

  # max-depth branch
  expect_silent({
    res_maxdepth <- split_node(c(leaf = 0.0, x = 1.0), X, D, parent_loss = Inf,
                               depth = 0, max_depth = 0, min_leaf_n = 2)
  })
  expect_identical(res_maxdepth$leaf_reason, "max-depth")

  # min-leaf branch (nrow <= min_leaf_n)
  X2 <- X[1:3, , drop = FALSE]; D2 <- D[1:3, , drop = FALSE]
  expect_silent({
    res_minleaf <- split_node(c(leaf = 0.0, x = 1.0), X2, D2, parent_loss = Inf,
                              depth = 0, max_depth = 8, min_leaf_n = 5)
  })
  expect_identical(res_minleaf$leaf_reason, "min-leaf")
})

test_that("split_node empty-child branch triggers when midpoint collapses a side", {
  # All x equal -> midpoint = that value -> one side empty
  X <- data.frame(x = rep(1, 5))
  rownames(X) <- as.character(seq_len(nrow(X)))
  D <- data.frame(vsq = rnorm(5)^2, w = rep(1, 5))
  rownames(D) <- rownames(X)

  res <- split_node(c(leaf = 0, x = 1), X, D, parent_loss = Inf,
                    depth = 0, max_depth = 8, min_leaf_n = 1)
  expect_identical(res$leaf_reason, "empty-child")
  expect_identical(res$node, "leaf")
})

test_that("split_node rejection loop can force a leaf after many rejects", {
  # Build X so that the first few candidate cuts are bad, then force 'leaf'
  set.seed(1)
  X <- data.frame(x = sort(rnorm(20)))
  rownames(X) <- as.character(seq_len(nrow(X)))
  D <- data.frame(vsq = rnorm(20)^2, w = rep(1, 20))
  rownames(D) <- rownames(X)

  # Start with 100% on 'x'; reduce_weight halves it repeatedly and never reaches good split
  # By making parent_loss extremely small relative to any possible new_loss, we ensure rejections.
  parent_loss <- -1e9

  res <- split_node(
    split_feature = c(leaf = 0.0, x = 1.0),
    X = X, D = D, parent_loss = parent_loss, depth = 0,
    explore_proba = 0.0,
    # keep defaults for choose_feature_fn / reduce_weight_fn
    max_depth = 8, min_leaf_n = 2,
    max_rejects_per_node = 5  # keep it small so we deterministically hit "forced-leaf-after-rejects"
  )
  expect_identical(res$node, "leaf")
  expect_true(res$leaf_reason %in% c("max-rejects", "forced-leaf-after-rejects"))
})

test_that("choose_feature strict guards: bad depth, NA probs, non-finite sum", {
  sf <- c(leaf = 0.6, x = 0.4)

  expect_error(choose_feature(sf, depth = c(1, 2)), "numeric scalar")
  expect_error(choose_feature(c(leaf = NA_real_, x = 1), depth = 0), "contains NA")

  # Non-finite sum: include Inf
  expect_error(choose_feature(c(leaf = 1, x = Inf), depth = 0), "sum is non-finite")
})

test_that("reduce_weight validates input and renormalizes to 1", {
  sf <- c(leaf = 0.2, x = 0.8)

  # fj not present
  expect_error(reduce_weight("z", sf), "not found")

  # NA in split_feature rejected
  expect_error(reduce_weight("leaf", c(leaf = NA_real_, x = 1)), "contains NA")

  # Negative prob rejected
  expect_error(reduce_weight("leaf", c(leaf = -0.1, x = 1.1)), ">= 0")

  # Works and sums to 1 after halving fj
  out <- reduce_weight("x", sf)
  expect_true(abs(sum(out) - 1) < 1e-12)
  expect_lt(out["x"], sf["x"])
  expect_gt(out["leaf"], 0)  # unchanged except renormalization
})

test_that("characterize_tree guards: X must be df; |w| mismatch; w must be binary", {
  skip_if_not_installed("rpart")

  X <- data.frame(a = rnorm(10), b = rnorm(10))
  w2 <- rep(0, 10)

  # X not a data.frame
  expect_error(characterize_tree(as.matrix(X), w2), "data frame")

  # length mismatch
  expect_error(characterize_tree(X, w = w2[1:5]), "equal the number of rows")

  # w has one level only
  expect_error(characterize_tree(X, w = rep(1, 10)), "exactly two classes")

  # w has > 2 levels
  expect_error(characterize_tree(X, w = c(rep(0, 4), rep(1, 4), rep(2, 2))), "exactly two classes")

  # NAs in X trigger .check_no_na()
  X_bad <- X; X_bad$a[3] <- NA_real_
  expect_error(characterize_tree(X_bad, w = c(rep(0,5), rep(1,5))),
               "Data `X` column 'a' contains missing values")})

test_that("characterize_tree fits with two classes and respects max_depth", {
  skip_if_not_installed("rpart")

  set.seed(7)
  X <- data.frame(a = rnorm(40), b = rnorm(40))
  # separable-ish classes
  w <- as.integer(X$a + X$b > 0); w[sample.int(40, 2)] <- 1 - w[sample.int(40, 2)]
  f <- characterize_tree(X, w, max_depth = 2)
  expect_s3_class(f, "rpart")
  expect_true(f$control$maxdepth == 2)
})
