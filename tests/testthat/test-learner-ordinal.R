# WS-9 tests: ordinal (cumulative-logit) learner adapter.

.make_ordinal_xy <- function(seed = 5L, n = 120L, p = 8L, n_bands = 4L) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  # f1 and f2 drive the latent score; the bands are cuts of the latent score.
  eta <- 1.2 * x[, "f1"] - 0.9 * x[, "f2"] + rnorm(n, sd = 0.6)
  y <- cut(eta, breaks = n_bands, labels = paste0("L", seq_len(n_bands)),
           ordered_result = TRUE)
  y <- setNames(y, rownames(x))
  list(x = x, y = y, signals = c("f1", "f2"))
}

test_that("auto_ordinal_lambda_grid returns a decreasing finite path", {
  skip_if_not_installed("ordinalNet")
  d <- .make_ordinal_xy()
  grid <- auto_ordinal_lambda_grid(d$x, d$y, n_lambda = 8L)
  expect_s3_class(grid, "data.frame")
  expect_equal(names(grid), "lambda")
  expect_true(nrow(grid) >= 2L)
  expect_true(all(is.finite(grid$lambda)))
  expect_true(all(diff(grid$lambda) < 0))
})

test_that("ordinal adapter masks are nested supersets along the path", {
  skip_if_not_installed("ordinalNet")
  d <- .make_ordinal_xy()
  grid <- auto_ordinal_lambda_grid(d$x, d$y, n_lambda = 6L)
  adapter <- make_ordinal_adapter()
  # Public adapter contract: one grid row per call, logical vector returned.
  one <- adapter(d$x, d$y, grid[1L, , drop = FALSE])
  expect_type(one, "logical")
  expect_length(one, ncol(d$x))
  # Internal batch adapter: full grid in one fit, logical matrix returned.
  masks <- stablr:::.make_ordinal_batch_adapter(alpha = 1)(d$x, d$y, grid)
  expect_equal(dim(masks), c(ncol(d$x), nrow(grid)))
  # Sparse-to-dense nesting: selections only grow as the penalty relaxes.
  cum <- t(apply(masks, 1L, cummax))
  expect_identical(masks, cum == 1L)
  expect_identical(one, masks[, 1L])
})

test_that("stabl_fit with base_learner = 'ordinal' recovers planted signals", {
  skip_if_not_installed("ordinalNet")
  d <- .make_ordinal_xy()
  fit <- suppressWarnings(stabl_fit(
    x               = d$x,
    y               = d$y,
    lambda_grid     = "auto",
    base_learner    = "ordinal",
    family          = "ordinal",
    n_bootstraps    = 40L,
    artificial_type = "random_permutation",
    random_state    = 42L
  ))
  expect_s3_class(fit, "stabl_fit")
  selected <- get_feature_names(fit)
  expect_true(all(d$signals %in% selected))
  # FDP+ machinery is outcome-agnostic and must be computable.
  expect_true(!is.null(fit$fdr_min_threshold_))
  # The fitted grid is decreasing (canonical plot order).
  expect_true(all(diff(fit$fitted_lambda_grid$lambda) < 0))
})

test_that("ordinal learner works with knockoff decoys end to end", {
  skip_if_not_installed("ordinalNet")
  skip_if_not_installed("knockoff")
  d <- .make_ordinal_xy(n = 150L, p = 6L)
  fit <- suppressWarnings(stabl_fit(
    x               = d$x,
    y               = d$y,
    lambda_grid     = "auto",
    base_learner    = "ordinal",
    family          = "ordinal",
    n_bootstraps    = 30L,
    artificial_type = "knockoff_equi",
    random_state    = 7L
  ))
  expect_s3_class(fit, "stabl_fit")
  expect_equal(nrow(fit$stabl_scores_artificial_), ncol(d$x))
})

test_that("ordinal learner rejects non-factor outcomes and honours seeds", {
  skip_if_not_installed("ordinalNet")
  d <- .make_ordinal_xy()
  y_num <- setNames(as.numeric(d$y), rownames(d$x))
  expect_error(
    stabl_fit(d$x, y_num, lambda_grid = "auto",
              base_learner = "ordinal", family = "ordinal",
              n_bootstraps = 5L, random_state = 1L),
    regexp = "factor"
  )
  f1 <- suppressWarnings(stabl_fit(
    d$x, d$y, lambda_grid = "auto", base_learner = "ordinal",
    family = "ordinal", n_bootstraps = 10L, random_state = 99L))
  f2 <- suppressWarnings(stabl_fit(
    d$x, d$y, lambda_grid = "auto", base_learner = "ordinal",
    family = "ordinal", n_bootstraps = 10L, random_state = 99L))
  expect_identical(f1$stabl_scores_, f2$stabl_scores_)
})

test_that("ordinal learner is exposed through the public interface audit", {
  # Prefer the shipped registry (available in R CMD check trees); fall
  # back to the skip-guarded source-tree lookup for test_local runs.
  audit_path <- system.file(
    "release", "public-interface-audit.csv", package = "stablr"
  )
  if (!nzchar(audit_path)) {
    audit_path <- file.path(.interface_source_root(), "inst", "release",
                            "public-interface-audit.csv")
  }
  audit <- utils::read.csv(
    audit_path,
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = character()
  )
  expect_true("make_ordinal_adapter" %in% audit$interface[audit$kind == "export"])
  expect_true("auto_ordinal_lambda_grid" %in% audit$interface[audit$kind == "export"])
})
