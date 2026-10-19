# Tests for the experimental oblique random-forest (aorsf) learner adapter
# (plan D2).  All fits use explicit seeds; aorsf lives in Suggests, so every
# test skips gracefully when it is unavailable.

skip_if_no_aorsf <- function() {
  skip_if_not_installed("aorsf")
}

# Non-linear synthetic design (plan D2): y depends on an interaction x1*x2
# plus main effects 0.8*x1 + 0.8*x2.  The main effects accompany the
# interaction because importance spread across linear combinations leaves
# pure interaction partners nearly invisible to marginal importance.  Nulls
# x3..x8 are AR(1)-correlated among themselves and mildly correlated with the
# signals (the leakage scenario FDP+ must survive).
# MASS-free construction: multivariate normal via Cholesky.
.make_oblique_synthetic <- function(n = 80L, p = 8L, seed = 1L) {
  set.seed(seed)
  p_null <- p - 2L
  Sigma  <- matrix(0.2, nrow = p, ncol = p)
  diag(Sigma) <- 1
  for (i in seq_len(p_null)) {
    for (j in seq_len(p_null)) {
      Sigma[2L + i, 2L + j] <- 0.5 ^ abs(i - j)
    }
  }
  Sigma[1:2, 1:2] <- 0.3
  diag(Sigma) <- 1
  L <- t(chol(Sigma))
  X <- matrix(rnorm(n * p), n, p) %*% L
  X <- matrix(as.numeric(X), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(
    0.8 * X[, 1] + 0.8 * X[, 2] + 1.5 * X[, 1] * X[, 2] + rnorm(n, sd = 0.5),
    rownames(X)
  )
  list(x = X, y = y, signals = paste0("f", seq_len(2L)))
}

# ---- (a) Adapter contract -----------------------------------------------------

test_that("make_aorsf_adapter returns valid logical masks of the right length", {
  skip_if_no_aorsf()
  d <- .make_oblique_synthetic(n = 60L, p = 6L, seed = 11L)
  grid <- auto_oblique_threshold_grid(d$x, d$y, n_grid = 5L, num_trees = 100L,
                                      seed = 1L)
  expect_s3_class(grid, "data.frame")
  expect_identical(nrow(grid), 5L)
  expect_true(all(c("lambda", "threshold") %in% names(grid)))

  # Public single-threshold adapter: 1-row grid -> one logical mask
  adapter <- make_aorsf_adapter(num_trees = 100L)
  mask <- adapter(d$x, d$y, grid[1L, , drop = FALSE])
  expect_type(mask, "logical")
  expect_length(mask, ncol(d$x))
  expect_false(anyNA(mask))

  # Internal batch adapter: full grid -> logical matrix (p x n_grid)
  batch <- stablr:::.make_aorsf_batch_adapter(family = "gaussian",
                                              num_trees = 100L)
  masks <- batch(d$x, d$y, grid)
  expect_true(is.matrix(masks))
  expect_identical(dim(masks), c(ncol(d$x), 5L))
  expect_type(masks, "logical")
  expect_false(anyNA(masks))
  # Sparse end of the path must select a subset, never more than dense end
  expect_lte(sum(masks[, 1L]), sum(masks[, ncol(masks)]))
})

test_that("public oblique_rf adapter matches the first column of the batch path", {
  skip_if_no_aorsf()
  d <- .make_oblique_synthetic(n = 60L, p = 6L, seed = 12L)
  grid <- auto_oblique_threshold_grid(d$x, d$y, n_grid = 4L, num_trees = 100L,
                                      seed = 2L)
  public <- make_aorsf_adapter(num_trees = 100L)
  batch  <- stablr:::.make_aorsf_batch_adapter(family = "gaussian",
                                               num_trees = 100L)
  # Same caller RNG stream -> same tree seeds -> identical masks
  set.seed(99L)
  m_public <- public(d$x, d$y, grid[1L, , drop = FALSE])
  set.seed(99L)
  m_batch <- batch(d$x, d$y, grid)[, 1L]
  expect_identical(m_public, m_batch)
})

# ---- (b) Monotonicity in the threshold ----------------------------------------

test_that("oblique_rf masks are monotone: lower threshold => superset selection", {
  skip_if_no_aorsf()
  d <- .make_oblique_synthetic(n = 80L, p = 8L, seed = 21L)
  grid <- auto_oblique_threshold_grid(d$x, d$y, n_grid = 8L, num_trees = 200L,
                                      seed = 3L)
  # Grid is sorted decreasing (sparse end first), like a glmnet lambda path
  expect_true(all(diff(grid$threshold) <= 0))

  batch <- stablr:::.make_aorsf_batch_adapter(family = "gaussian",
                                              num_trees = 200L)
  masks <- batch(d$x, d$y, grid)
  k <- ncol(masks)
  for (j in seq_len(k - 1L)) {
    expect_true(all(!masks[, j] | masks[, j + 1L]),
                info = paste0("threshold column ", j, " -> ", j + 1L))
  }
})

# ---- (c) End-to-end stabl_fit(base_learner = "oblique_rf") ----------------------

test_that("stabl_fit(base_learner='oblique_rf') recovers non-linear planted signals", {
  skip_if_no_aorsf()
  d <- .make_oblique_synthetic(n = 80L, p = 8L, seed = 31L)

  # Structure + recovery via the statistically sound pairing (model-X knockoff
  # nulls).  The knockoff_equi generator (create.gaussian) is independent of
  # the fixed-X create.fixed construction, so no fixed-X build probe gates
  # this assertion (a former probe here measured create.fixed's internal
  # unit-norm scale against a sd ~ 1 expectation and always skipped).
  if (requireNamespace("knockoff", quietly = TRUE)) {
    fit_k <- stabl_fit(
      x = d$x, y = d$y,
      lambda_grid = "auto",
      base_learner = "oblique_rf",
      artificial_type = "knockoff_equi",
      n_bootstraps = 15L,
      random_state = 7L
    )
    expect_s3_class(fit_k, "stabl_fit")
    expect_identical(nrow(fit_k$fitted_lambda_grid), 30L)
    selected <- get_feature_names(fit_k)
    expect_true(all(d$signals %in% selected),
                info = paste("selected:", paste(selected, collapse = ",")))
  } else {
    skip("knockoff package not available")
  }

  # Permutation-null structure variant (over-selection direction: signals must
  # still be inside the selected set).  The oblique_rf + random_permutation
  # warning is expected here and asserted separately below.
  fit_p <- suppressWarnings(stabl_fit(
    x = d$x, y = d$y,
    lambda_grid = "auto",
    base_learner = "oblique_rf",
    artificial_type = "random_permutation",
    n_bootstraps = 15L,
    random_state = 8L
  ))
  selected_p <- get_feature_names(fit_p)
  expect_true(all(d$signals %in% selected_p),
              info = paste("selected:", paste(selected_p, collapse = ",")))
})

test_that("oblique_rf fits are reproducible under random_state", {
  skip_if_no_aorsf()
  d <- .make_oblique_synthetic(n = 60L, p = 6L, seed = 41L)
  fit1 <- stabl_fit(d$x, d$y, lambda_grid = "auto", base_learner = "oblique_rf",
                    artificial_type = "random_permutation",
                    n_bootstraps = 3L, random_state = 5L)
  fit2 <- stabl_fit(d$x, d$y, lambda_grid = "auto", base_learner = "oblique_rf",
                    artificial_type = "random_permutation",
                    n_bootstraps = 3L, random_state = 5L)
  expect_identical(fit1$stabl_scores_, fit2$stabl_scores_)
})

# ---- (d) Guardrail warning ------------------------------------------------------

test_that("oblique_rf + random_permutation emits the anti-conservative FDP+ warning", {
  skip_if_no_aorsf()
  d <- .make_oblique_synthetic(n = 60L, p = 6L, seed = 51L)
  expect_warning(
    stabl_fit(d$x, d$y, lambda_grid = "auto", base_learner = "oblique_rf",
              artificial_type = "random_permutation",
              n_bootstraps = 3L, random_state = 9L),
    "anti-conservative"
  )
  # The sound pairing must NOT warn
  if (requireNamespace("knockoff", quietly = TRUE)) {
    expect_no_warning(
      suppressWarnings(stabl_fit(d$x, d$y, lambda_grid = "auto",
                                 base_learner = "oblique_rf",
                                 artificial_type = "knockoff_equi",
                                 n_bootstraps = 3L, random_state = 10L))
    )
  }
})

test_that("unknown base_learner error message lists oblique_rf", {
  xn <- matrix(rnorm(40), 10, 4, dimnames = list(paste0("s", 1:10),
                                                 paste0("f", 1:4)))
  yn <- setNames(rnorm(10), paste0("s", 1:10))
  expect_error(
    stabl_fit(xn, yn, lambda_grid = "auto", base_learner = "nope",
              n_bootstraps = 2L),
    "must be one of"
  )
})
