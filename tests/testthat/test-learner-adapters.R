# Tests for the three public single-lambda adapter factories (Item 15).
# stabl_fit() uses the batch variants internally; these factories are the
# user-facing single-lambda equivalents exported for custom loop usage.

.make_adapter_data <- function(n = 20L, p = 6L, seed = 1L) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", seq_len(p))))
  y <- rnorm(n)
  list(x = x, y = y, lambda_val = data.frame(lambda = 0.1))
}

# ---- make_glmnet_adapter -----------------------------------------------------

test_that("make_glmnet_adapter returns a selection mask of correct type/length", {
  skip_if_not_installed("glmnet")
  d <- .make_adapter_data()
  adapter <- make_glmnet_adapter(family = "gaussian")
  mask <- adapter(d$x, d$y, d$lambda_val)
  expect_type(mask, "logical")
  expect_length(mask, ncol(d$x))
})

test_that("make_glmnet_adapter alpha_fixed overrides lambda_val alpha column", {
  skip_if_not_installed("glmnet")
  d <- .make_adapter_data()
  d$lambda_val$alpha <- 0.5
  adapter_fixed <- make_glmnet_adapter(family = "gaussian", alpha_fixed = 0.0)
  adapter_col   <- make_glmnet_adapter(family = "gaussian")
  mask_fixed <- adapter_fixed(d$x, d$y, d$lambda_val)
  mask_col   <- adapter_col(d$x, d$y, d$lambda_val)
  expect_false(identical(mask_fixed, mask_col))
})

test_that("make_glmnet_adapter bootstrap_threshold controls sparsity", {
  skip_if_not_installed("glmnet")
  d <- .make_adapter_data()
  lax    <- make_glmnet_adapter(bootstrap_threshold = 1e-10)
  strict <- make_glmnet_adapter(bootstrap_threshold = 1e3)
  expect_gte(sum(lax(d$x, d$y, d$lambda_val)), sum(strict(d$x, d$y, d$lambda_val)))
})

test_that("public glmnet adapter matches one-row batch path exactly", {
  skip_if_not_installed("glmnet")
  d <- .make_adapter_data(seed = 11L)
  public <- make_glmnet_adapter(family = "gaussian")
  batch <- stablr:::.make_glmnet_batch_adapter(
    family = "gaussian", alpha_fixed = NULL,
    bootstrap_threshold = stablr:::.BOOTSTRAP_COEF_THRESHOLD
  )
  expect_identical(
    public(d$x, d$y, d$lambda_val),
    batch(d$x, d$y, d$lambda_val)[, 1L]
  )
})

# ---- make_adaptive_lasso_adapter ---------------------------------------------

test_that("make_adaptive_lasso_adapter returns a selection mask", {
  skip_if_not_installed("glmnet")
  d <- .make_adapter_data()
  adapter <- make_adaptive_lasso_adapter(family = "gaussian")
  mask <- adapter(d$x, d$y, d$lambda_val)
  expect_type(mask, "logical")
  expect_length(mask, ncol(d$x))
})

test_that("make_adaptive_lasso_adapter rejects non-positive gamma", {
  expect_error(make_adaptive_lasso_adapter(gamma = -1), "gamma")
  expect_error(make_adaptive_lasso_adapter(gamma = 0),  "gamma")
})

test_that("make_adaptive_lasso_adapter rejects non-positive epsilon", {
  expect_error(make_adaptive_lasso_adapter(epsilon = 0),  "epsilon")
  expect_error(make_adaptive_lasso_adapter(epsilon = -1), "epsilon")
})

test_that("public adaptive adapter matches one-row batch path exactly", {
  skip_if_not_installed("glmnet")
  d <- .make_adapter_data(seed = 12L)
  public <- make_adaptive_lasso_adapter(family = "gaussian")
  batch <- stablr:::.make_adaptive_lasso_batch_adapter(
    family = "gaussian", gamma = 1, epsilon = 1e-6,
    bootstrap_threshold = stablr:::.BOOTSTRAP_COEF_THRESHOLD
  )
  expect_identical(
    public(d$x, d$y, d$lambda_val),
    batch(d$x, d$y, d$lambda_val)[, 1L]
  )
})

# ---- make_sgl_adapter --------------------------------------------------------

test_that("make_sgl_adapter returns a selection mask", {
  skip_if_not_installed("sparsegl")
  d <- .make_adapter_data()
  groups  <- rep(1:3, each = 2L)
  adapter <- make_sgl_adapter(family = "gaussian", feature_groups = groups)
  mask <- adapter(d$x, d$y, d$lambda_val)
  expect_type(mask, "logical")
  expect_length(mask, ncol(d$x))
})

test_that("make_sgl_adapter rejects cox family", {
  skip_if_not_installed("sparsegl")
  expect_error(
    make_sgl_adapter(family = "cox", feature_groups = 1:3),
    "Cox family is not supported"
  )
})

test_that("make_sgl_adapter requires feature_groups", {
  skip_if_not_installed("sparsegl")
  expect_error(make_sgl_adapter(), "feature_groups")
})

# ---- C7: adaptive-lasso hyperparameters are validated on the batch path -------
# stabl_fit() reaches .make_adaptive_lasso_batch_adapter directly, bypassing the
# public single-lambda factory, so the guards must live in the batch factory.

test_that("C7 the adaptive batch adapter rejects non-positive gamma and epsilon", {
  skip_if_not_installed("glmnet")
  batch <- function(gamma, epsilon) {
    stablr:::.make_adaptive_lasso_batch_adapter(
      family = "gaussian", gamma = gamma, epsilon = epsilon,
      bootstrap_threshold = stablr:::.BOOTSTRAP_COEF_THRESHOLD
    )
  }
  for (bad in list(-1, 0, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(batch(gamma = bad, epsilon = 1e-6), "gamma")
  }
  for (bad in list(0, -1, NA_real_, Inf, c(1e-6, 1e-6), "1e-6")) {
    expect_error(batch(gamma = 1, epsilon = bad), "epsilon")
  }
  expect_no_error(batch(gamma = 1, epsilon = 1e-6))
})

test_that("C7 stabl_fit rejects invalid adaptive hyperparameters before bootstrapping", {
  d <- .make_adapter_data(n = 20L, p = 5L)
  x <- d$x
  rownames(x) <- paste0("s", seq_len(nrow(x)))
  y <- stats::setNames(d$y, rownames(x))
  expect_error(
    stabl_fit(x, y, lambda_grid = d$lambda_val, base_learner = "adaptive_lasso",
              adaptive_gamma = -1, n_bootstraps = 1L, artificial_type = NULL,
              hard_threshold = 0.5, sample_fraction = 1, random_state = 1L),
    "gamma"
  )
  expect_error(
    stabl_fit(x, y, lambda_grid = d$lambda_val, base_learner = "adaptive_lasso",
              adaptive_epsilon = 0, n_bootstraps = 1L, artificial_type = NULL,
              hard_threshold = 0.5, sample_fraction = 1, random_state = 1L),
    "epsilon"
  )
})

# ---- C8: the epsilon floor is on the scale of the initial scores ---------------

test_that("C8 the adaptive penalty floor scales with the response, not absolutely", {
  # Ridge scores on a small-scale response.  An absolute 1e-6 floor dominates
  # every score, so all penalty factors collapse to one value and adaptive lasso
  # degenerates to plain lasso without any signal that it has done so.
  scores <- c(3e-8, 1e-8, 5e-9, 2e-8, 8e-9, 1.2e-8, 4e-8, 6e-9)
  absolute <- 1 / ((scores + 1e-6)^1)
  expect_lt(max(absolute) / min(absolute), 1.05)

  relative <- stablr:::.adaptive_lasso_penalty_factor(scores, 1, 1e-6)
  expect_gt(length(unique(relative)), 1L)
  # The relative floor keeps a real spread and keeps the ordering: a larger
  # ridge score still earns a smaller penalty factor.
  expect_gt(max(relative) / min(relative), 2)
  expect_identical(order(relative), order(-scores))
})

test_that("C8 the adaptive penalty factor degenerates only when there is no signal", {
  # A non-positive median score means the ridge fit carries no scale
  # information; the documented fallback is the historical absolute floor, and
  # equal scores must still give equal weights.
  expect_identical(
    stablr:::.adaptive_lasso_penalty_factor(c(0, 0, 0), 1, 1e-6),
    rep(1e6, 3)
  )
  # A feature with no estimable ridge coefficient takes the maximum penalty
  # rather than producing a NaN weight that glmnet would reject.
  weights <- stablr:::.adaptive_lasso_penalty_factor(c(1, 2, NA_real_), 1, 1e-6)
  expect_true(all(is.finite(weights)))
  expect_true(all(weights > 0))
  expect_identical(weights[[3L]], max(weights))
})

test_that("C8 adaptivity actually changes the support relative to plain lasso", {
  skip_if_not_installed("glmnet")
  # The oracle-property claim has no test without this: identical supports would
  # mean the adaptive weights have no effect at all.
  set.seed(7)
  n <- 30L; p <- 20L
  x <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("f", seq_len(p))))
  y <- 2 * x[, 1L] + stats::rnorm(n, sd = 0.5)

  adaptive <- make_adaptive_lasso_adapter(family = "gaussian")
  plain <- make_glmnet_adapter(family = "gaussian")
  lambda_val <- data.frame(lambda = 0.08)

  adaptive_mask <- adaptive(x, y, lambda_val)
  plain_mask <- plain(x, y, lambda_val)

  expect_true(any(adaptive_mask != plain_mask))
  # Both retain the true signal feature, so the difference is in the
  # feature-specific weighting rather than in a wholesale failure.
  expect_true(adaptive_mask[[1L]])
  expect_true(plain_mask[[1L]])
})
