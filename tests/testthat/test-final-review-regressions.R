# Final pre-publication review (2026-10-18) regression tests.
# Each block references the finding ID from
# docs/audit/2026-10-17-adversarial/final-review.md.

test_that("B-01: auto_rf_threshold_grid works with default seed = NULL", {
  # Regression: .rf_fit_permutation_importance computed
  # `seed + (attempt - 1L) * 104729L` on a NULL seed, producing numeric(0)
  # and a ranger error ("Expecting a single value: [extent=0]") on every
  # default-path call, including stabl_fit(base_learner = "rf",
  # lambda_grid = "auto") without random_state.
  skip_if_not_installed("ranger")
  set.seed(1L)
  n <- 20L
  x <- matrix(rnorm(n * 10L), n, 10L,
              dimnames = list(NULL, paste0("f", 1:10)))
  y <- x[, 1L] * x[, 2L] + rnorm(n)

  imp <- stablr:::.rf_fit_permutation_importance(x, y, 100L, NULL)
  expect_length(imp, 10L)
  expect_false(anyNA(imp))

  grid <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L)
  expect_s3_class(grid, "data.frame")
  expect_identical(nrow(grid), 5L)
  expect_true(all(is.finite(grid$lambda)))
})

test_that("B-02: threshold grid builder never returns a 0-row grid", {
  # The grid builder must stop with an informative message when the
  # importance helper returns NULL (ranger degenerate on every retry),
  # instead of silently returning a 0-row data.frame. The NULL path is not
  # reachable without mocking ranger; here we assert the 0-row symptom
  # cannot recur for valid importances.
  set.seed(2L)
  imp <- setNames(rnorm(10L), paste0("f", 1:10))
  grid <- .rf_threshold_grid_from_importance(imp, 5L)
  expect_gte(nrow(grid), 1L)
})

test_that("N-02: get_support validates new_hard_threshold in (0, 1]", {
  withr::local_seed(0)
  n <- 40L; p <- 6L
  x <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid = data.frame(lambda = c(0.2, 0.1, 0.05)),
                   n_bootstraps = 8L,
                   artificial_type = "random_permutation",
                   random_state = 1L)

  expect_error(get_support(fit, new_hard_threshold = 1.5), "new_hard_threshold")
  expect_error(get_support(fit, new_hard_threshold = -1), "new_hard_threshold")
  expect_error(get_support(fit, new_hard_threshold = 0), "new_hard_threshold")
  expect_error(
    get_feature_names(fit, new_hard_threshold = 5), "new_hard_threshold"
  )
  # Boundary values remain valid.
  expect_type(get_support(fit, new_hard_threshold = 1), "logical")
  expect_type(get_support(fit, new_hard_threshold = 0.5), "logical")
})

test_that("D-1: feature plots enforce length(y) == nrow(x)", {
  skip_if_not_installed("ggplot2")
  set.seed(3L)
  n <- 40L
  x <- matrix(rnorm(n * 4L), n, 4L,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", 1:4)))
  y_short <- setNames(rnorm(n / 2L), paste0("s", seq_len(n / 2L)))
  expect_error(
    boxplot_features(c("f1", "f2"), x, y_short),
    "same number of elements"
  )
  expect_error(
    scatterplot_features(c("f1", "f2"), x, y_short),
    "same number of elements"
  )
})

test_that("C-3: y_valid without x_valid_list fails closed", {
  set.seed(4L)
  n <- 40L
  x_list <- list(
    omic1 = matrix(rnorm(n * 5L), n, 5L,
                   dimnames = list(paste0("s", 1:n), paste0("g", 1:5))),
    omic2 = matrix(rnorm(n * 4L), n, 4L,
                   dimnames = list(paste0("s", 1:n), paste0("p", 1:4)))
  )
  y <- setNames(rnorm(n), paste0("s", 1:n))
  train <- 1:30
  lam <- data.frame(lambda = c(0.3, 0.1))
  expect_error(
    stabl_multiomic_train_validate(
      lapply(x_list, function(m) m[train, , drop = FALSE]),
      y[train],
      lambda_grid = lam,
      y_valid = y[31:40],
      family = "gaussian",
      n_bootstraps = 4L,
      hard_threshold = 0.2,
      random_state = 1L
    ),
    "y_valid` was supplied without `x_valid_list`"
  )
})

test_that("C-5: cooperative prediction matrix fails closed on row mismatch", {
  pred <- matrix(1:6, nrow = 3L)
  expect_error(
    stablr:::.cooperative_prediction_matrix(pred, n_rows = 4L),
    "expected 4"
  )
  out <- stablr:::.cooperative_prediction_matrix(pred, n_rows = 3L)
  expect_identical(nrow(out), 3L)
})
