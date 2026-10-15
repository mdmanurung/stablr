# Tests for fold-safe preprocessing (PR2).
# Core property: a transform fitted on training data must be completely
# insensitive to assessment values.

.make_pp_data <- function(seed = 1L, n_tr = 50L, n_as = 12L, p = 5L) {
  set.seed(seed)
  x_tr <- matrix(rnorm(n_tr * p), n_tr, p,
                 dimnames = list(paste0("s", seq_len(n_tr)),
                                 paste0("f", seq_len(p))))
  x_as <- matrix(rnorm(n_as * p), n_as, p,
                 dimnames = list(paste0("v", seq_len(n_as)),
                                 paste0("f", seq_len(p))))
  list(x_tr = x_tr, x_as = x_as)
}

test_that("center/scale standardise training data", {
  d <- .make_pp_data()
  pp <- stabl_preprocess_fit(d$x_tr, transforms = c("center", "scale"))
  xt <- stabl_preprocess_apply(pp, d$x_tr)
  expect_equal(colMeans(xt), setNames(rep(0, 5), colnames(xt)),
               tolerance = 1e-10)
  expect_equal(apply(xt, 2L, sd), setNames(rep(1, 5), colnames(xt)),
               tolerance = 1e-10)
})

test_that("assessment shift does not alter training-fitted parameters", {
  d <- .make_pp_data()
  pp <- stabl_preprocess_fit(d$x_tr, transforms = c("center", "scale"))
  params_before <- pp$params

  # Extreme assessment values must not change the fitted transform.
  x_as_extreme <- d$x_as
  x_as_extreme[, 1L] <- 1e6
  invisible(stabl_preprocess_apply(pp, x_as_extreme))
  expect_identical(pp$params, params_before)

  # Applying to assessment uses the training parameters, so the assessment
  # output is NOT standardised around its own mean.
  xa <- stabl_preprocess_apply(pp, x_as_extreme)
  expect_false(isTRUE(all.equal(unname(colMeans(xa)[1L]), 0,
                                tolerance = 1e-6)))
})

test_that("median imputation uses training medians on assessment data", {
  d <- .make_pp_data()
  d$x_tr[1L, 2L] <- NA
  d$x_tr[2L, 2L] <- NA
  pp <- stabl_preprocess_fit(d$x_tr, transforms = "impute_median")
  med <- stats::median(d$x_tr[, 2L], na.rm = TRUE)

  d$x_as[1L, 2L] <- NA
  xa <- stabl_preprocess_apply(pp, d$x_as)
  expect_equal(unname(xa[1L, 2L]), unname(med))
})

test_that("zero-variance filter drops constant training features", {
  d <- .make_pp_data()
  d$x_tr[, 3L] <- 7  # constant in training
  pp <- stabl_preprocess_fit(d$x_tr, transforms = "zv")
  expect_equal(pp$retained_features, c("f1", "f2", "f4", "f5"))
  xt <- stabl_preprocess_apply(pp, d$x_tr)
  expect_equal(ncol(xt), 4L)
  # The dropped feature is absent even if the assessment set varies in it.
  d$x_as[, 3L] <- rnorm(nrow(d$x_as))
  xa <- stabl_preprocess_apply(pp, d$x_as)
  expect_equal(ncol(xa), 4L)
})

test_that("NULL / empty transforms leave input unchanged", {
  d <- .make_pp_data()
  pp <- stabl_preprocess_fit(d$x_tr, transforms = character())
  expect_identical(stabl_preprocess_apply(pp, d$x_tr), d$x_tr)
})

test_that("apply rejects missing required features", {
  d <- .make_pp_data()
  pp <- stabl_preprocess_fit(d$x_tr, transforms = c("center", "scale"))
  expect_error(
    stabl_preprocess_apply(pp, d$x_as[, -1L]),
    "missing"
  )
})

test_that("fit rejects unnamed or nonnumeric input", {
  d <- .make_pp_data()
  x_nonames <- d$x_tr
  colnames(x_nonames) <- NULL
  expect_error(stabl_preprocess_fit(x_nonames), "column")
  expect_error(stabl_preprocess_fit(matrix(letters[1:10], 2, 5)), "numeric")
})

test_that("unsupported transforms are rejected", {
  d <- .make_pp_data()
  expect_error(
    stabl_preprocess_fit(d$x_tr, transforms = "pca"),
    "Unsupported"
  )
})

test_that("all-missing-in-training feature fails closed on apply", {
  d <- .make_pp_data()
  d$x_tr[, 4L] <- NA_real_  # entirely missing -> median NA
  pp <- stabl_preprocess_fit(
    d$x_tr, transforms = c("impute_median", "center", "scale")
  )
  # No training information exists for this feature; imputation cannot
  # invent one, so apply must fail closed rather than emit silent garbage
  # (regression guard for the Windows CI failure of 2026-10-15: the old
  # check only fired when the NEW data also had missing values).
  expect_error(
    suppressWarnings(stabl_preprocess_apply(pp, d$x_tr)),
    "entirely missing in training"
  )
  # Same fail-closed behavior on fully observed new data: the fitted
  # center/scale parameters for the feature are NA, so the output would be
  # an all-NA column.
  expect_error(
    stabl_preprocess_apply(pp, d$x_as),
    "entirely missing in training"
  )
})

test_that("constant-after-imputation feature scales safely", {
  d <- .make_pp_data()
  # f4: observed values all equal with some missing -> constant after
  # median imputation; the SD guard must keep scaling finite (no NaN/Inf).
  d$x_tr[, 4L] <- 2.5
  d$x_tr[sample(seq_len(nrow(d$x_tr)), 10L), 4L] <- NA_real_
  pp <- stabl_preprocess_fit(
    d$x_tr, transforms = c("impute_median", "center", "scale")
  )
  xt <- suppressWarnings(stabl_preprocess_apply(pp, d$x_tr))
  expect_true(all(is.finite(xt[, "f4"])))
  # (c - c) / 1: centered constant column is exactly zero.
  expect_true(all(xt[, "f4"] == 0))
})

test_that("feature order in assessment does not matter", {
  d <- .make_pp_data()
  pp <- stabl_preprocess_fit(d$x_tr, transforms = c("center", "scale"))
  xa1 <- stabl_preprocess_apply(pp, d$x_as)
  xa2 <- stabl_preprocess_apply(pp, d$x_as[, rev(seq_len(ncol(d$x_as)))])
  expect_identical(xa1, xa2)
})

test_that("refit on the same training data is deterministic", {
  d <- .make_pp_data()
  pp1 <- stabl_preprocess_fit(d$x_tr, transforms = c("impute_median",
                                                     "center", "scale"))
  pp2 <- stabl_preprocess_fit(d$x_tr, transforms = c("impute_median",
                                                     "center", "scale"))
  expect_identical(pp1$params, pp2$params)
  expect_identical(pp1$retained_features, pp2$retained_features)
})

test_that("recipes path requires the recipes package and preps on training", {
  skip_if_not_installed("recipes")
  d <- .make_pp_data()
  rec <- recipes::recipe(~ ., data = as.data.frame(d$x_tr)) |>
    recipes::step_center(recipes::all_predictors()) |>
    recipes::step_scale(recipes::all_predictors())
  pp <- stabl_preprocess_fit_recipes(rec, d$x_tr)
  expect_s3_class(pp, "stabl_preprocessing_recipes")
  xt <- stabl_preprocess_apply(pp, d$x_tr)
  expect_equal(colMeans(xt), setNames(rep(0, 5), colnames(xt)),
               tolerance = 1e-8)
})
