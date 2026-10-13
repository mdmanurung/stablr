# Workstream D tests: C3 accessors (get_fdrs, get_lambda_grid), C4 S3 methods
# (plot.stabl_fit, summary.stabl_fit, print.summary.stabl_fit), and the C5
# `...`-collision fix in the multiomic late-fusion path.

.ws_d_design <- function() {
  set.seed(20261012L)
  n <- 80L; p <- 12L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  list(x = x, y = y)
}

# ---- C3: get_fdrs -----------------------------------------------------------

test_that("C3: get_fdrs returns named FDP+ estimates matching the fit fields", {
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1, 0.05)),
                   n_bootstraps    = 20L,
                   artificial_type = "random_permutation",
                   random_state    = 11L)
  fdrs <- get_fdrs(fit)
  expect_true(is.numeric(fdrs))
  expect_false(is.null(names(fdrs)))
  expect_identical(unname(fdrs), fit$FDRs_)
  expect_identical(length(fdrs), length(fit$fdr_threshold_range))
  expect_identical(names(fdrs), as.character(fit$fdr_threshold_range))
  expect_identical(min(fdrs), fit$min_fdr_)
})

test_that("C3: get_fdrs errors informatively when no artificial features were used", {
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps    = 20L,
                   artificial_type = NULL,
                   hard_threshold  = 0.5,
                   random_state    = 11L)
  expect_error(get_fdrs(fit), "FDRs_|artificial")
})

# ---- C3: get_lambda_grid ----------------------------------------------------

test_that("C3: get_lambda_grid returns the fitted lambda data.frame", {
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps    = 20L,
                   artificial_type = NULL,
                   hard_threshold  = 0.5,
                   random_state    = 11L)
  expect_s3_class(get_lambda_grid(fit), "data.frame")
  expect_identical(get_lambda_grid(fit), fit$fitted_lambda_grid)
  expect_identical(nrow(get_lambda_grid(fit)), ncol(fit$stabl_scores_))
})

# ---- C4: plot.stabl_fit -----------------------------------------------------

test_that("C4: plot.stabl_fit delegates to plot_stabl_path and returns a ggplot invisibly", {
  skip_if_not_installed("ggplot2")
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps    = 20L,
                   artificial_type = NULL,
                   hard_threshold  = 0.5,
                   random_state    = 11L)
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
  expect_invisible(plot(fit))
  # ... forwarding: a custom threshold override reaches plot_stabl_path().
  p_strict <- plot(fit, new_hard_threshold = 0.9)
  expect_s3_class(p_strict, "ggplot")
})

# ---- C4: summary.stabl_fit / print.summary.stabl_fit ------------------------

test_that("C4: summary.stabl_fit returns the documented fields with correct values", {
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps    = 20L,
                   artificial_type = NULL,
                   hard_threshold  = 0.5,
                   random_state    = 11L)
  s <- summary(fit)
  expect_s3_class(s, "summary.stabl_fit")
  expect_type(s$selected_features, "character")
  expect_setequal(s$selected_features, get_feature_names(fit))
  expect_identical(s$n_selected, length(s$selected_features))
  expect_identical(s$stabl_threshold, fit$hard_threshold)
  expect_null(s$min_fdr)
  expect_identical(s$lambda_grid_used, dim(fit$fitted_lambda_grid))
  expect_identical(s$n_features_in, fit$n_features_in_)
  expect_identical(s$n_bootstraps, 20L)
})

test_that("C4: summary.stabl_fit resolves the FDP+ threshold when no hard threshold", {
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps    = 20L,
                   artificial_type = "random_permutation",
                   random_state    = 11L)
  s <- summary(fit)
  expect_identical(s$stabl_threshold, fit$fdr_min_threshold_)
  expect_identical(s$min_fdr, fit$min_fdr_)
})

test_that("C4: print.summary.stabl_fit runs invisibly and returns the summary", {
  d <- .ws_d_design()
  fit <- stabl_fit(d$x, d$y,
                   lambda_grid     = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps    = 20L,
                   artificial_type = NULL,
                   hard_threshold  = 0.5,
                   random_state    = 11L)
  s <- summary(fit)
  expect_output(out <- print(s), "stabl_fit summary")
  expect_s3_class(out, "summary.stabl_fit")
  expect_invisible(print(s))
})

# ---- C5: `...` collision fix in the late-fusion path ------------------------

.ws_d_multiomic_inputs <- function(n = 30L) {
  set.seed(20261013L)
  ids <- paste0("s", seq_len(n))
  x1 <- matrix(rnorm(n * 4L), n, 4L, dimnames = list(ids, paste0("g", 1:4)))
  x2 <- matrix(rnorm(n * 3L), n, 3L, dimnames = list(ids, paste0("p", 1:3)))
  y  <- setNames(rnorm(n), ids)
  list(x_list = list(omic1 = x1, omic2 = x2), y = y,
       lam = data.frame(lambda = c(0.2, 0.1)))
}

test_that("C5: former internal late-fusion formal names no longer collide through ...", {
  d <- .ws_d_multiomic_inputs()

  # `nfolds` used to be a formal of the internal .late_fusion_oof_fit(); a
  # user passing it through ... produced a confusing duplicate-formal error
  # ("formal argument matched by multiple actual arguments").  It must now
  # simply reach stabl_fit(), which rejects it as an unused argument.
  expect_error(
    stabl_multiomic_train_validate(
      x_train_list = d$x_list, y_train = d$y,
      lambda_grid = d$lam, late_fusion = TRUE, n_iter_lf = 50L,
      n_bootstraps = 5L, artificial_type = NULL, hard_threshold = 0.5,
      random_state = 1L, nfolds = 3L
    ),
    "unused argument"
  )
})

test_that("C5: valid stabl_fit args still forward through ... with late fusion", {
  d <- .ws_d_multiomic_inputs()
  fit_lf <- stabl_multiomic_train_validate(
    x_train_list = d$x_list, y_train = d$y,
    lambda_grid = d$lam, late_fusion = TRUE, n_iter_lf = 50L,
    n_bootstraps = 5L, artificial_type = NULL, hard_threshold = 0.5,
    random_state = 1L, sample_fraction = 0.75
  )
  expect_s3_class(fit_lf, "stabl_multiomic_fit")
  expect_false(is.null(fit_lf$late_fusion))
  expect_identical(ncol(fit_lf$fits$omic1$stabl_scores_), 2L)
  expect_identical(ncol(fit_lf$fits$omic2$stabl_scores_), 2L)
})

test_that("C5: stabl_multiomic_cv forwards ... to the fixed late-fusion path", {
  d <- .ws_d_multiomic_inputs(n = 24L)
  cv <- stabl_multiomic_cv(
    x_list = d$x_list, y = d$y, lambda_grid = d$lam, v = 2L,
    late_fusion = TRUE, n_iter_lf = 20L,
    n_bootstraps = 3L, artificial_type = NULL, hard_threshold = 0.5,
    random_state = 2L
  )
  expect_s3_class(cv, "stabl_multiomic_cv")
  expect_identical(length(cv$fold_results), 2L)
  expect_true(all(vapply(cv$fold_results, function(f) !is.null(f$late_fusion),
                         logical(1L))))
})
