# Tests for the gaussian validated path (1.0.0 spec A5; Phase 13 PR-F):
# stabl_multiomic_nested_cv(family = "gaussian") and
# stabl_analyze(mode = "validated", family = "gaussian").
#
# Coverage: synthetic linear benchmark, fold-isolated selection, negative
# R2 retention, empty-support deterministic fallback, study-disjoint outer
# folds, missing-outcome error, repeats and seed invariance, metric/family
# mismatch errors, event_class / return_probabilities rejection, and the
# stabl_analyze validated-gaussian surface (performance, held-out
# importance, final_refit, print).

# ---- shared fixtures --------------------------------------------------------

.make_gauss_data <- function(n = 120L, p = 30L, seed = 11L, sd = 2) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(ids, paste0("f", seq_len(p))))
  beta <- c(2, -1.5, 1, rep(0, p - 3L))
  y <- setNames(as.numeric(x %*% beta + rnorm(n, sd = sd)), ids)
  list(x = x, y = y, ids = ids)
}

grid3 <- data.frame(lambda = c(0.3, 0.2, 0.1))

# ---- stabl_multiomic_nested_cv gaussian ------------------------------------

test_that("gaussian nested CV recovers signal with honest held-out R2", {
  d <- .make_gauss_data()
  fit <- stabl_multiomic_nested_cv(
    x_list = list(x = d$x),
    y = d$y,
    candidates = list(list(name = "x", blocks = "x")),
    family = "gaussian",
    lambda_grid = grid3,
    outer_v = 4L, inner_v = 3L,
    n_bootstraps = 30L,
    random_state = 42L
  )
  expect_identical(fit$family, "gaussian")
  expect_true(is.numeric(fit$outer_predictions$truth))
  expect_true(is.numeric(fit$outer_predictions$predicted))
  expect_false("event_probability" %in% names(fit$outer_predictions))
  m <- .regression_metrics(fit$outer_predictions$truth,
                           fit$outer_predictions$predicted)
  # Loose bound: the 3-signal linear problem must beat the intercept-only
  # model by a clear margin, but we do not pin a point estimate.
  expect_gt(m$r2, 0.3)
  expect_lt(m$rmse, stats::sd(d$y))
  # Selection happened inside folds: per-fold selected panels exist and
  # the true signals dominate the selection counts.
  sel <- fit$selected_features
  expect_true(all(c("repeat_id", "fold_id", "feature") %in% names(sel)))
  top <- names(sort(table(sel$feature), decreasing = TRUE))
  expect_true(all(c("f1", "f2", "f3") %in% top[seq_len(3L)]))
})

test_that("gaussian nested CV is replay-identical under the same seed", {
  d <- .make_gauss_data(n = 90L, p = 20L)
  run <- function() {
    stabl_multiomic_nested_cv(
      x_list = list(x = d$x), y = d$y,
      candidates = list(list(name = "x", blocks = "x")),
      family = "gaussian",
      lambda_grid = grid3,
      outer_v = 3L, inner_v = 3L, n_bootstraps = 20L,
      random_state = 7L
    )
  }
  a <- run()
  b <- run()
  expect_identical(a$outer_predictions, b$outer_predictions)
  expect_identical(a$selected_features, b$selected_features)
})

test_that("gaussian nested CV supports repeated folds", {
  d <- .make_gauss_data(n = 80L, p = 15L)
  fit <- stabl_multiomic_nested_cv(
    x_list = list(x = d$x), y = d$y,
    candidates = list(list(name = "x", blocks = "x")),
    family = "gaussian",
    lambda_grid = grid3,
    outer_v = 3L, outer_repeats = 2L, inner_v = 3L,
    n_bootstraps = 15L, random_state = 9L
  )
  expect_equal(sort(unique(fit$outer_predictions$repeat_id)), c(1L, 2L))
  expect_equal(nrow(fit$outer_predictions), 2L * 80L)
})

test_that("empty support falls back to deterministic constant predictions", {
  d <- .make_gauss_data(n = 80L, p = 15L)
  fit <- stabl_multiomic_nested_cv(
    x_list = list(x = d$x), y = d$y,
    candidates = list(list(name = "x", blocks = "x")),
    family = "gaussian",
    lambda_grid = grid3,
    outer_v = 3L, inner_v = 3L, n_bootstraps = 15L,
    hard_threshold = 1,  # strict > comparator: nothing is ever selected
    random_state = 5L
  )
  preds <- fit$outer_predictions
  # Every fold predicts its own training mean: within-fold predictions are
  # constant, and the pooled R2 is <= 0 (kept, not clipped or NA-hidden).
  by_fold <- split(preds$predicted, preds$fold_id)
  expect_true(all(vapply(by_fold, function(v) length(unique(v)) == 1L,
                         logical(1L))))
  m <- .regression_metrics(preds$truth, preds$predicted)
  expect_true(is.finite(m$rmse))
  expect_lte(m$r2, 0)
})

test_that("negative R2 is reported as-is for a useless panel", {
  # Intercept-only predictions on held-out data: R2 <= 0 by construction
  # and must not be clipped to 0 or converted to NA.
  m <- .regression_metrics(truth = c(1, 2, 3, 10), predicted = c(5, 5, 5, 5))
  expect_lt(m$r2, 0)
  # Constant held-out truths: R2 is NA with an explicit reason.
  m_const <- .regression_metrics(truth = c(2, 2, 2), predicted = c(2, 2, 2))
  expect_true(is.na(m_const$r2))
  expect_equal(m_const$rmse, 0)
})

test_that("gaussian study-disjoint outer folds keep studies apart", {
  d <- .make_gauss_data(n = 80L, p = 15L)
  studies <- setNames(rep(paste0("study", 1:4), each = 20L), d$ids)
  fit <- stabl_multiomic_nested_cv(
    x_list = list(x = d$x), y = d$y,
    candidates = list(list(name = "x", blocks = "x")),
    family = "gaussian",
    lambda_grid = grid3,
    outer_v = 4L, inner_v = 2L, n_bootstraps = 15L,
    study_labels = studies,
    random_state = 13L
  )
  expect_identical(fit$fold_scheme, "study_disjoint")
  for (fr in fit$fold_results) {
    tr <- unique(as.character(studies[fr$outer_fold$train_ids]))
    va <- unique(as.character(studies[fr$outer_fold$valid_ids]))
    expect_length(intersect(tr, va), 0L)
  }
  expect_true(!is.null(fit$performance_by_study))
  expect_true(all(c("study", "mae", "rmse", "r2") %in%
                    names(fit$performance_by_study)))
})

test_that("gaussian nested CV rejects bad inputs fail-closed", {
  d <- .make_gauss_data(n = 60L, p = 10L)
  y_na <- d$y
  y_na[3L] <- NA
  expect_error(
    stabl_multiomic_nested_cv(
      x_list = list(x = d$x), y = y_na,
      candidates = list(list(name = "x", blocks = "x")),
      family = "gaussian",
      lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
      n_bootstraps = 10L, random_state = 1L
    ),
    "missing|NA|finite", ignore.case = TRUE
  )
  # metric/family mismatch
  expect_error(
    stabl_multiomic_nested_cv(
      x_list = list(x = d$x), y = d$y,
      candidates = list(list(name = "x", blocks = "x")),
      family = "gaussian",
      lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
      n_bootstraps = 10L, metric = "ber", random_state = 1L
    ),
    "should be one of"
  )
  # event_class is meaningless for regression
  expect_error(
    stabl_multiomic_nested_cv(
      x_list = list(x = d$x), y = d$y,
      candidates = list(list(name = "x", blocks = "x")),
      family = "gaussian",
      lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
      n_bootstraps = 10L, event_class = "case", random_state = 1L
    ),
    "event_class"
  )
  # return_probabilities is a classification-only request
  expect_error(
    stabl_multiomic_nested_cv(
      x_list = list(x = d$x), y = d$y,
      candidates = list(list(name = "x", blocks = "x")),
      family = "gaussian",
      lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
      n_bootstraps = 10L, return_probabilities = TRUE, random_state = 1L
    ),
    "return_probabilities"
  )
  # factor y under gaussian
  expect_error(
    stabl_multiomic_nested_cv(
      x_list = list(x = d$x),
      y = setNames(factor(rep(c("a", "b"), each = 30L)), d$ids),
      candidates = list(list(name = "x", blocks = "x")),
      lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
      n_bootstraps = 10L, family = "gaussian", random_state = 1L
    ),
    "numeric"
  )
})

test_that("gaussian metric selection supports rmse/mae/r2 and defaults", {
  d <- .make_gauss_data(n = 60L, p = 10L)
  fit_default <- stabl_multiomic_nested_cv(
    x_list = list(x = d$x), y = d$y,
    candidates = list(list(name = "x", blocks = "x")),
    family = "gaussian",
    lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
    n_bootstraps = 10L, random_state = 3L
  )
  expect_identical(fit_default$metric, "rmse")
  fit_r2 <- stabl_multiomic_nested_cv(
    x_list = list(x = d$x), y = d$y,
    candidates = list(list(name = "x", blocks = "x")),
    family = "gaussian",
    lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
    n_bootstraps = 10L, metric = "r2", random_state = 3L
  )
  expect_identical(fit_r2$metric, "r2")
})

# ---- stabl_analyze validated gaussian ---------------------------------------

test_that("stabl_analyze validated gaussian returns the full evidence surface", {
  d <- .make_gauss_data()
  val <- stabl_analyze(
    d$x, d$y, family = "gaussian", mode = "validated",
    preprocessing = c("center", "scale"),
    lambda_grid = grid3,
    outer_v = 4L, inner_v = 3L, n_bootstraps = 20L,
    random_state = 42L, n_repeats = 5L, final_refit = TRUE
  )
  expect_s3_class(val, "stabl_analysis")
  expect_identical(val$mode, "validated")
  expect_identical(val$family, "gaussian")
  expect_null(val$event_class)
  expect_identical(val$evidence_type, "cross_validated")
  expect_identical(val$evaluation_scope, "held_out_assessment")

  pooled <- val$performance$pooled
  expect_equal(sort(names(pooled)), sort(c("n", "mae", "rmse", "r2")))
  expect_equal(pooled$n, length(d$y))
  expect_gt(pooled$r2, 0.3)
  expect_equal(nrow(val$performance$by_repeat), 1L)

  imp <- val$heldout_importance$summary
  expect_true(all(c("feature", "mean_importance", "n_folds_selected") %in%
                    names(imp)))
  expect_true(all(is.na(imp$event_class)))
  selected_imp <- imp[!is.na(imp$mean_importance), ]
  expect_true(all(selected_imp$n_folds_selected >= 1L))
  # The three true signals dominate held-out importance on average
  # (per-fold importance is noisy; compare signal vs noise means).
  is_signal <- selected_imp$feature %in% c("f1", "f2", "f3")
  expect_true(all(c("f1", "f2", "f3") %in% selected_imp$feature))
  if (any(!is_signal)) {
    expect_gt(mean(selected_imp$mean_importance[is_signal]),
              mean(selected_imp$mean_importance[!is_signal]))
  }

  expect_identical(val$final_refit$evidence_type, "final_refit")
  expect_identical(val$final_refit$fit$family, "gaussian")
})

test_that("stabl_analyze validated gaussian is replay-identical", {
  d <- .make_gauss_data(n = 90L, p = 20L)
  run <- function() {
    stabl_analyze(
      d$x, d$y, family = "gaussian", mode = "validated",
      lambda_grid = grid3,
      outer_v = 3L, inner_v = 3L, n_bootstraps = 15L,
      random_state = 17L, n_repeats = 3L
    )
  }
  a <- run()
  b <- run()
  expect_identical(a$predictions, b$predictions)
  expect_identical(a$performance, b$performance)
  expect_identical(a$heldout_importance$summary, b$heldout_importance$summary)
})

test_that("stabl_analyze metric resolution and mismatch guards", {
  d <- .make_gauss_data(n = 60L, p = 10L)
  val <- stabl_analyze(
    d$x, d$y, family = "gaussian", mode = "validated",
    lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
    n_bootstraps = 10L, random_state = 1L, n_repeats = 2L
  )
  expect_identical(val$nested$metric, "rmse")
  expect_error(
    stabl_analyze(d$x, d$y, family = "gaussian", mode = "validated",
                  metric = "ber", lambda_grid = grid3,
                  outer_v = 3L, inner_v = 3L, n_bootstraps = 5L,
                  random_state = 1L),
    "should be one of"
  )
  expect_error(
    stabl_analyze(d$x, d$y, family = "gaussian", mode = "validated",
                  event_class = "case", lambda_grid = grid3,
                  outer_v = 3L, inner_v = 3L, n_bootstraps = 5L,
                  random_state = 1L),
    "`event_class` is only defined for family"
  )
  # Binomial default is unchanged (ber).
  set.seed(3L)
  yb <- setNames(factor(rep(c("A", "B"), each = 30L)), d$ids)
  xb <- d$x
  xb[, 1L] <- xb[, 1L] + 1.5 * as.numeric(yb == "B")
  valb <- suppressWarnings(stabl_analyze(
    xb, yb, family = "binomial", mode = "validated",
    lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
    n_bootstraps = 10L, random_state = 5L, n_repeats = 2L
  ))
  expect_identical(valb$nested$metric, "ber")
  expect_true(all(c("accuracy", "balanced_error_rate", "brier", "auroc") %in%
                    names(valb$performance$pooled)))
})

test_that("print and summary handle the gaussian validated surface", {
  d <- .make_gauss_data(n = 60L, p = 10L)
  val <- stabl_analyze(
    d$x, d$y, family = "gaussian", mode = "validated",
    lambda_grid = grid3, outer_v = 3L, inner_v = 3L,
    n_bootstraps = 10L, random_state = 2L, n_repeats = 2L
  )
  out <- capture.output(print(val))
  expect_true(any(grepl("family:      gaussian", out, fixed = TRUE)))
  expect_true(any(grepl("RMSE", out, fixed = TRUE)))
  expect_false(any(grepl("AUROC", out, fixed = TRUE)))
  s <- capture.output(invisible(summary(val)))
  expect_true(any(grepl("Held-out performance", s, fixed = TRUE)))
})
