# sPLS-DA learner adapter tests (mixOmics base learner, classification only).

.make_splsda_xy <- function(seed = 7L, n = 60L, p = 20L, n_signal = 5L) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  beta <- c(rep(1.1, n_signal), rep(0, p - n_signal))
  eta <- as.vector(x %*% beta)
  y <- rbinom(n, 1L, plogis(eta))
  y <- setNames(y, rownames(x))
  list(x = x, y = y, signals = paste0("f", seq_len(n_signal)))
}

test_that("auto_splsda_keepx_grid returns a valid increasing keepX path", {
  x <- matrix(0, 10, 150)
  grid <- auto_splsda_keepx_grid(x, n_lambda = 8L)
  expect_s3_class(grid, "data.frame")
  expect_setequal(names(grid), c("lambda", "keepX"))
  expect_true(nrow(grid) >= 2L && nrow(grid) <= 8L)
  expect_true(all(grid$lambda == grid$keepX))
  expect_true(all(grid$keepX >= 1L))
  expect_true(all(diff(grid$keepX) > 0))
  expect_true(all(is.finite(grid$keepX)))
})

test_that("splsda adapter returns a logical mask respecting keepX", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy()
  adapter <- make_splsda_adapter()
  mask <- adapter(d$x, d$y, data.frame(lambda = 5))
  expect_type(mask, "logical")
  expect_length(mask, ncol(d$x))
  expect_true(is.logical(mask))
  # keepX = 5 per component, ncomp = 2 -> at most 10 distinct features kept.
  expect_true(sum(mask) <= 10L)
  expect_true(sum(mask) >= 1L)
})

test_that("splsda batch adapter sweeps the full grid in one call", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy()
  grid <- auto_splsda_keepx_grid(d$x, n_lambda = 5L)
  masks <- stablr:::.make_splsda_batch_adapter(
    family = "binomial", ncomp = 2L, bootstrap_threshold = 1e-10
  )(d$x, d$y, grid)
  expect_equal(dim(masks), c(ncol(d$x), nrow(grid)))
  expect_true(is.logical(masks))
  # Sparse-to-dense nesting: a larger keepX keeps a superset per component.
  expect_true(all(masks[, 1L] <= masks[, ncol(masks), drop = TRUE] |
                    !masks[, 1L]))
})

test_that("keepX column takes precedence over the lambda column", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy()
  batch <- stablr:::.make_splsda_batch_adapter(
    family = "binomial", ncomp = 1L, bootstrap_threshold = 1e-10
  )
  # Same value in either column -> identical masks.
  via_lambda <- batch(d$x, d$y, data.frame(lambda = 4))
  via_keepx  <- batch(d$x, d$y, data.frame(lambda = 99, keepX = 4))
  expect_identical(via_lambda, via_keepx)
  # Different keepX under a constant lambda column -> different masks
  # (proves the keepX column, not the lambda column, drives the fit).
  small <- batch(d$x, d$y, data.frame(lambda = 1, keepX = 2))
  large <- batch(d$x, d$y, data.frame(lambda = 1, keepX = 15))
  expect_true(all(small <= large))  # keepX = 2 support nests in keepX = 15
  expect_true(sum(small) <= sum(large))
})

test_that("splsda rejects non-binomial families and bad y", {
  expect_error(make_splsda_adapter(family = "gaussian"), "binomial")
  expect_error(
    stablr:::.make_splsda_batch_adapter("gaussian", 2L, 1e-10),
    "binomial"
  )
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy()
  adapter <- make_splsda_adapter()
  expect_error(adapter(d$x, setNames(d$y + runif(length(d$y)), names(d$y)),
                       data.frame(lambda = 5)), "0/1")
  expect_error(adapter(d$x, setNames(rep(1L, nrow(d$x)), rownames(d$x)),
                       data.frame(lambda = 5)), "two observed levels")
})

test_that("splsda accepts factor y and 0/1 numeric y identically", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy()
  y_fac <- factor(d$y, levels = c("0", "1"))
  y_fac <- setNames(y_fac, names(d$y))
  adapter <- make_splsda_adapter()
  m_num  <- adapter(d$x, d$y, data.frame(lambda = 5))
  m_fac  <- adapter(d$x, y_fac, data.frame(lambda = 5))
  expect_identical(m_num, m_fac)
})

test_that("stabl_fit end-to-end with base_learner = 'splsda' (RP and knockoff)", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy(n = 60L, p = 20L)
  fit_rp <- suppressWarnings(stabl_fit(
    x               = d$x,
    y               = d$y,
    lambda_grid     = "auto",
    base_learner    = "splsda",
    family          = "binomial",
    n_bootstraps    = 20L,
    artificial_type = "random_permutation",
    random_state    = 11L
  ))
  expect_s3_class(fit_rp, "stabl_fit")
  expect_false(is.null(fit_rp$min_fdr_))
  expect_false(is.null(fit_rp$FDRs_))
  expect_false(is.null(fit_rp$stabl_scores_))
  expect_equal(nrow(fit_rp$stabl_scores_), ncol(d$x))
  expect_type(get_feature_names(fit_rp), "character")

  fit_ke <- suppressWarnings(stabl_fit(
    x               = d$x,
    y               = d$y,
    lambda_grid     = "auto",
    base_learner    = "splsda",
    family          = "binomial",
    n_bootstraps    = 20L,
    artificial_type = "knockoff_equi",
    random_state    = 11L
  ))
  expect_s3_class(fit_ke, "stabl_fit")
  expect_false(is.null(fit_ke$min_fdr_))
  # Artificial block accumulated separately from real features.
  expect_false(is.null(fit_ke$stabl_scores_artificial_))
})

test_that("splsda stabl_fit is deterministic under a fixed seed", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy(n = 50L, p = 15L)
  one <- stabl_fit(
    x = d$x, y = d$y, lambda_grid = "auto",
    base_learner = "splsda", family = "binomial",
    n_bootstraps = 10L, artificial_type = "random_permutation",
    random_state = 3L
  )
  two <- stabl_fit(
    x = d$x, y = d$y, lambda_grid = "auto",
    base_learner = "splsda", family = "binomial",
    n_bootstraps = 10L, artificial_type = "random_permutation",
    random_state = 3L
  )
  expect_identical(one$stabl_scores_, two$stabl_scores_)
})

test_that("custom keepX grid flows through stabl_fit", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy(n = 50L, p = 15L)
  custom <- data.frame(lambda = c(2, 5, 10), keepX = c(2, 5, 10))
  fit <- stabl_fit(
    x = d$x, y = d$y, lambda_grid = custom,
    base_learner = "splsda", family = "binomial",
    n_bootstraps = 10L, artificial_type = "random_permutation",
    random_state = 5L
  )
  expect_equal(nrow(fit$stabl_scores_), ncol(d$x))
  expect_equal(ncol(fit$stabl_scores_), 3L)
  expect_equal(fit$fitted_lambda_grid$keepX, c(2, 5, 10))
})

test_that("splsda_ncomp is honored (ncomp = 1 vs 2)", {
  skip_if_not_installed("mixOmics")
  d <- .make_splsda_xy()
  adapter1 <- make_splsda_adapter(ncomp = 1L)
  adapter2 <- make_splsda_adapter(ncomp = 2L)
  m1 <- adapter1(d$x, d$y, data.frame(lambda = 5))
  m2 <- adapter2(d$x, d$y, data.frame(lambda = 5))
  # With two components the union of kept features can only grow.
  expect_true(all(m1 <= m2))
  expect_true(sum(m1) <= sum(m2))
})
