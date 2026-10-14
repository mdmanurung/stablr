test_that("repeated stratified folds cover all samples once per repeat", {
  y <- setNames(
    factor(rep(c("A", "B", "C"), each = 6L)),
    paste0("s", seq_len(18L))
  )

  folds <- .make_repeated_cv_folds(
    y = y,
    v = 3L,
    repeats = 2L,
    stratified = TRUE,
    random_state = 42L
  )

  expect_length(folds, 6L)
  for (rep_i in 1:2) {
    rep_folds <- folds[vapply(folds, `[[`, integer(1L), "repeat") == rep_i]
    valid_ids <- unlist(lapply(rep_folds, `[[`, "valid_ids"))
    expect_equal(sort(valid_ids), sort(names(y)))
    expect_equal(max(table(valid_ids)), 1L)
    for (fold in rep_folds) {
      expect_equal(as.integer(table(y[fold$valid_ids])), c(2L, 2L, 2L))
    }
  }
})

test_that("stabl_multiomic_nested_cv returns deterministic diagnostics and predictions", {
  set.seed(101)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  signal <- model.matrix(~ y - 1)

  x_a <- matrix(rnorm(n * 6L, sd = 0.3), nrow = n,
                dimnames = list(ids, paste0("a", seq_len(6L))))
  x_b <- matrix(rnorm(n * 5L, sd = 0.3), nrow = n,
                dimnames = list(ids, paste0("b", seq_len(5L))))
  x_a[, 1:3] <- x_a[, 1:3] + signal
  x_b[, 1:3] <- x_b[, 1:3] + signal

  fit1 <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x_a, mirna = x_b),
    y = y,
    candidates = list(
      list(name = "mrna", blocks = "mrna"),
      list(name = "early_fusion", blocks = c("mrna", "mirna"))
    ),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L,
    outer_repeats = 1L,
    inner_v = 3L,
    metric = "ber",
    family = "multinomial",
    artificial_type = NULL,
    hard_threshold = 1e-9,
    n_bootstraps = 3L,
    sample_fraction = 1.0,
    random_state = 11L
  ))

  fit2 <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x_a, mirna = x_b),
    y = y,
    candidates = list(
      list(name = "mrna", blocks = "mrna"),
      list(name = "early_fusion", blocks = c("mrna", "mirna"))
    ),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L,
    outer_repeats = 1L,
    inner_v = 3L,
    metric = "ber",
    family = "multinomial",
    artificial_type = NULL,
    hard_threshold = 1e-9,
    n_bootstraps = 3L,
    sample_fraction = 1.0,
    random_state = 11L
  ))

  expect_s3_class(fit1, "stabl_multiomic_nested_cv")
  expect_length(fit1$outer_folds, 3L)
  expect_equal(nrow(fit1$outer_predictions), n)
  expect_equal(sort(fit1$outer_predictions$sample_id), sort(ids))
  expect_true(all(fit1$diagnostics$candidate %in% c("mrna", "early_fusion")))
  expect_named(fit1$performance,
               c("accuracy", "balanced_error_rate", "per_class_recall",
                 "macro_f1", "mcc", "confusion"))

  expect_equal(fit1$outer_predictions, fit2$outer_predictions)
  expect_equal(fit1$diagnostics, fit2$diagnostics)
})

test_that("stabl_multiomic_nested_cv supports custom categorical and numeric strata", {
  set.seed(202)
  n <- 36L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(rep(c("A", "B", "C"), each = 3L), 4L)), ids)
  batch <- setNames(rep(paste0("batch", 1:4), each = 9L), ids)
  age <- setNames(seq_len(n), ids)
  x <- matrix(rnorm(n * 5L), nrow = n,
              dimnames = list(ids, paste0("f", seq_len(5L))))

  fit_cat <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x),
    y = y,
    strata = batch,
    outer_v = 3L,
    outer_repeats = 1L,
    inner_v = 3L,
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    artificial_type = NULL,
    hard_threshold = 0.2,
    n_bootstraps = 2L,
    sample_fraction = 1.0,
    random_state = 3L
  ))

  expect_true(fit_cat$stratified)
  expect_equal(levels(fit_cat$strata), paste0("batch", 1:4))
  for (fold in fit_cat$outer_folds) {
    expect_true(all(table(fit_cat$strata[fold$valid_ids]) >= 1L))
  }

  fit_num <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x),
    y = y,
    strata = age,
    strata_bins = 3L,
    outer_v = 3L,
    outer_repeats = 1L,
    inner_v = 3L,
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    artificial_type = NULL,
    hard_threshold = 0.2,
    n_bootstraps = 2L,
    sample_fraction = 1.0,
    random_state = 4L,
    cv_workers = 1L
  ))

  expect_equal(nlevels(fit_num$strata), 3L)
  expect_equal(nrow(fit_num$outer_predictions), n)
})

test_that("stabl_multiomic_nested_cv rejects folds with too few class samples", {
  ids <- paste0("s", seq_len(8L))
  x <- matrix(rnorm(8L * 3L), nrow = 8L,
              dimnames = list(ids, paste0("f", seq_len(3L))))
  y <- setNames(factor(c(rep("A", 6L), rep("B", 2L))), ids)

  expect_error(
    stabl_multiomic_nested_cv(
      x_list = list(mrna = x),
      y = y,
      outer_v = 3L,
      inner_v = 3L,
      n_bootstraps = 2L,
      artificial_type = NULL,
      hard_threshold = 0.1
    ),
    "Each stratum"
  )
})

test_that("stabl_multiomic_nested_cv forwards l1_ratio to auto lambda grids", {
  set.seed(303)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  x <- matrix(rnorm(n * 5L), nrow = n,
              dimnames = list(ids, paste0("f", seq_len(5L))))

  fit <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x),
    y = y,
    candidates = list(list(name = "mrna", blocks = "mrna")),
    lambda_grid = "auto",
    outer_v = 3L,
    outer_repeats = 1L,
    inner_v = 3L,
    family = "multinomial",
    base_learner = "elastic_net",
    l1_ratio = 0.5,
    n_lambda = 2L,
    artificial_type = NULL,
    hard_threshold = 1,
    n_bootstraps = 2L,
    sample_fraction = 1.0,
    random_state = 303L
  ))

  grid <- fit$fold_results[[1L]]$fit$stabl_fit$fitted_lambda_grid
  expect_true("alpha" %in% names(grid))
  expect_equal(unique(grid$alpha), 0.5)
})

test_that("selected multinomial fitting rethrows programming errors", {
  x_train <- matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  x_valid <- matrix(rnorm(6L), nrow = 3L, ncol = 2L)
  y_train <- factor(c(rep("A", 6L), rep("B", 4L)))
  testthat::local_mocked_bindings(
    .with_glmnet_numerical_infeasibility = function(...) {
      stop("deliberate selected-model programming defect")
    },
    .package = "stablr"
  )

  expect_error(
    stablr:::.predict_selected_multinomial(
      x_train, y_train, x_valid, levels(y_train)
    ),
    "deliberate selected-model programming defect"
  )
})

test_that("selected multinomial fitting catches typed numerical failures", {
  x_train <- matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  x_valid <- matrix(rnorm(6L), nrow = 3L, ncol = 2L)
  y_train <- factor(c(rep("A", 6L), rep("B", 4L)))
  testthat::local_mocked_bindings(
    .with_glmnet_numerical_infeasibility = function(...) {
      stablr:::.abort_numerical_infeasibility(
        "stablr_learner_numerical_infeasibility",
        "deliberate selected-model numerical infeasibility"
      )
    },
    .package = "stablr"
  )

  expect_identical(
    stablr:::.predict_selected_multinomial(
      x_train, y_train, x_valid, levels(y_train)
    ),
    rep("A", nrow(x_valid))
  )
})

test_that("one selected feature uses a genuine one-feature binary fit", {
  # Strongly informative single feature: the one-feature glm must recover the
  # signal instead of falling back to the majority class.
  set.seed(7L)
  x_train <- matrix(rnorm(40L), nrow = 40L, ncol = 1L,
                    dimnames = list(paste0("s", 1:40), "f1"))
  y_train <- factor(ifelse(x_train[, 1L] > 0, "event", "control"),
                    levels = c("control", "event"))
  x_valid <- matrix(c(-2, -1, 1, 2), nrow = 4L, ncol = 1L,
                    dimnames = list(paste0("v", 1:4), "f1"))

  pred <- stablr:::.predict_selected_multinomial(
    x_train, y_train, x_valid, levels(y_train)
  )
  expect_identical(pred, c("control", "control", "event", "event"))
})

test_that("one selected feature falls back to majority on separation", {
  # Complete separation with a tiny training set: the glm cannot produce a
  # finite fit, so the training-derived majority class must be returned.
  x_train <- matrix(c(-3, -2, -1, 1, 2, 3), nrow = 6L, ncol = 1L,
                    dimnames = list(paste0("s", 1:6), "f1"))
  y_train <- factor(c(rep("control", 3L), rep("event", 3L)),
                    levels = c("control", "event"))
  x_valid <- matrix(0, nrow = 2L, ncol = 1L,
                    dimnames = list(c("v1", "v2"), "f1"))

  pred <- suppressWarnings(stablr:::.predict_selected_multinomial(
    x_train, y_train, x_valid, levels(y_train)
  ))
  expect_true(all(pred %in% levels(y_train)))
  expect_length(pred, 2L)
})

test_that("zero selected features keep the majority-class contract", {
  x_train <- matrix(numeric(0), nrow = 10L, ncol = 0L)
  x_valid <- matrix(numeric(0), nrow = 3L, ncol = 0L)
  y_train <- factor(c(rep("B", 7L), rep("A", 3L)), levels = c("A", "B"))

  expect_identical(
    stablr:::.predict_selected_multinomial(
      x_train, y_train, x_valid, levels(y_train)
    ),
    rep("B", 3L)
  )
})

test_that("one selected feature works for multinomial outcomes", {
  set.seed(11L)
  x_train <- matrix(rnorm(60L), nrow = 60L, ncol = 1L,
                    dimnames = list(paste0("s", 1:60), "f1"))
  y_train <- factor(
    ifelse(x_train[, 1L] < -0.5, "low",
           ifelse(x_train[, 1L] > 0.5, "high", "mid")),
    levels = c("low", "mid", "high")
  )
  x_valid <- matrix(c(-2, 0, 2), nrow = 3L, ncol = 1L,
                    dimnames = list(paste0("v", 1:3), "f1"))

  pred <- stablr:::.predict_selected_multinomial(
    x_train, y_train, x_valid, levels(y_train)
  )
  expect_true(all(pred %in% levels(y_train)))
  expect_length(pred, 3L)
  # The extreme validation points should map to the extreme classes.
  expect_identical(pred[1L], "low")
  expect_identical(pred[3L], "high")
})

test_that("one-feature prediction is deterministic across repeated calls", {
  set.seed(13L)
  x_train <- matrix(rnorm(30L), nrow = 30L, ncol = 1L,
                    dimnames = list(paste0("s", 1:30), "f1"))
  y_train <- factor(ifelse(x_train[, 1L] > 0, "event", "control"),
                    levels = c("control", "event"))
  x_valid <- matrix(rnorm(6L), nrow = 6L, ncol = 1L,
                    dimnames = list(paste0("v", 1:6), "f1"))

  p1 <- stablr:::.predict_selected_multinomial(
    x_train, y_train, x_valid, levels(y_train)
  )
  p2 <- stablr:::.predict_selected_multinomial(
    x_train, y_train, x_valid, levels(y_train)
  )
  expect_identical(p1, p2)
})

# ── PR3: preprocessing threading and held-out probabilities ─────────────────

test_that("preprocessing is fitted inside fold training (zv drops constants)", {
  set.seed(202)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B"), each = 9L)), ids)
  x <- matrix(rnorm(n * 5L, sd = 0.3), n, 5L,
              dimnames = list(ids, paste0("g", 1:5)))
  x[, 1] <- x[, 1] + as.numeric(y == "B")
  x[, 5] <- 7  # constant feature: must be filtered inside every fold

  fit <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(rna = x), y = y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
    random_state = 5L,
    preprocessing = c("zv", "center", "scale"),
    keep_fold_predictors = TRUE
  ))
  expect_false("g5" %in% fit$selected_features$feature)
  for (fr in fit$fold_results) {
    expect_s3_class(fr$fit$preprocessing, "stabl_preprocessing")
    expect_false("rna::g5" %in% fr$fit$preprocessing$retained_features)
  }
})

test_that("return_probabilities adds an aligned binary event probability", {
  set.seed(303)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B"), each = 9L)), ids)
  x <- matrix(rnorm(n * 5L, sd = 0.3), n, 5L,
              dimnames = list(ids, paste0("g", 1:5)))
  x[, 1] <- x[, 1] + as.numeric(y == "B")

  fit <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(rna = x), y = y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
    random_state = 6L,
    return_probabilities = TRUE
  ))
  preds <- fit$outer_predictions
  expect_true("event_probability" %in% names(preds))
  expect_true(all(preds$event_probability >= 0 & preds$event_probability <= 1))
  expect_equal(nrow(preds), n)
  expect_identical(fit$event_class, "B")
  # Legacy hard-class columns are untouched.
  expect_true(all(c("truth", "predicted", "selected_candidate") %in%
                    names(preds)))

  # Deterministic for a fixed seed.
  fit2 <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(rna = x), y = y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
    random_state = 6L,
    return_probabilities = TRUE
  ))
  expect_identical(fit$outer_predictions, fit2$outer_predictions)
})

test_that("legacy nested-CV output is unchanged when new options are omitted", {
  set.seed(404)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B"), each = 9L)), ids)
  x <- matrix(rnorm(n * 5L, sd = 0.3), n, 5L,
              dimnames = list(ids, paste0("g", 1:5)))
  x[, 1] <- x[, 1] + as.numeric(y == "B")

  fit <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(rna = x), y = y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
    random_state = 8L
  ))
  expect_false("event_probability" %in% names(fit$outer_predictions))
  expect_null(fit$event_class)
  expect_null(fit$preprocessing)
  for (fr in fit$fold_results) {
    expect_null(fr$fit$predictor)
    expect_null(fr$fit$preprocessing)
  }
})

test_that("return_probabilities fails clearly for multiclass outcomes", {
  set.seed(505)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  x <- matrix(rnorm(n * 4L, sd = 0.3), n, 4L,
              dimnames = list(ids, paste0("g", 1:4)))
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = list(rna = x), y = y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 9L,
      return_probabilities = TRUE
    )),
    "binary"
  )
})
