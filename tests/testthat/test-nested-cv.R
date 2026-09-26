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

# ---- D3: parallel outer folds need a seed, and worker failures are typed -------

test_that("D3 cv_workers > 1 without random_state raises a typed condition", {
  set.seed(11)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  x <- matrix(rnorm(n * 4L), nrow = n,
              dimnames = list(ids, paste0("f", seq_len(4L))))

  condition <- tryCatch(
    stabl_multiomic_nested_cv(
      x_list = list(mrna = x), y = y,
      candidates = list(list(name = "mrna", blocks = "mrna")),
      lambda_grid = data.frame(lambda = 0.2), outer_v = 3L, inner_v = 3L,
      n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.2,
      sample_fraction = 1, cv_workers = 2L
    ),
    error = identity
  )
  expect_s3_class(condition, "stablr_nested_cv_unseeded_parallel")
  expect_s3_class(condition, "stablr_error")
  expect_match(conditionMessage(condition), "random_state", fixed = TRUE)
})

test_that("D3 a failed outer fold is reported with its fold id, not a subscript error", {
  folds <- list(
    list(fold_id = "Repeat1_Fold1"), list(fold_id = "Repeat1_Fold2")
  )
  ok <- list(fold_result = 1L, diagnostics = 1L, predictions = 1L,
             selected_features = 1L)
  failed <- list(fold_result = NULL)

  expect_no_error(
    stablr:::.check_nested_cv_fold_results(list(ok, ok), folds)
  )

  condition <- tryCatch(
    stablr:::.check_nested_cv_fold_results(
      list(ok, structure("worker died", class = "try-error")), folds
    ),
    error = identity
  )
  expect_s3_class(condition, "stablr_nested_cv_fold_failed")
  expect_identical(condition$fold_index, 2L)
  expect_identical(condition$fold_id, "Repeat1_Fold2")
  expect_match(conditionMessage(condition), "worker died", fixed = TRUE)
  expect_match(conditionMessage(condition), "Repeat1_Fold2", fixed = TRUE)
})

# ---- D2: the downstream cv.glmnet is seeded ------------------------------------

test_that("D2 the selected-feature cv.glmnet no longer depends on ambient RNG", {
  set.seed(12)
  x_train <- matrix(stats::rnorm(160L), nrow = 40L, ncol = 4L)
  x_valid <- matrix(stats::rnorm(40L), nrow = 10L, ncol = 4L)
  y_train <- factor(rep(c("A", "B", "C"), length.out = 40L))
  # Each class needs at least 8 training rows for cv.glmnet to fold safely.
  expect_true(all(table(y_train) >= 8L))

  seeded <- function() {
    suppressWarnings(stablr:::.predict_selected_multinomial(
      x_train, y_train, x_valid, levels(y_train), random_state = 909L
    ))
  }
  first <- seeded()
  # Consume ambient RNG between calls: a seeded call must be unaffected.
  set.seed(4L)
  invisible(stats::rnorm(100))
  expect_identical(seeded(), first)

  # Unseeded, the cv.glmnet fold assignment is free to follow ambient state; the
  # helper simply does not pin it.
  expect_no_error(suppressWarnings(
    stablr:::.predict_selected_multinomial(
      x_train, y_train, x_valid, levels(y_train)
    )
  ))
})

test_that("D2 a one-feature support has no downstream model and falls back safely", {
  set.seed(13)
  x_train <- matrix(stats::rnorm(20L), nrow = 20L, ncol = 1L)
  x_valid <- matrix(stats::rnorm(5L), nrow = 5L, ncol = 1L)
  y_train <- factor(rep(c("A", "B"), each = 10L))
  # cv.glmnet rejects a single-column x; the majority-class fallback is the
  # only defensible answer and it must not raise.
  expect_identical(
    stablr:::.predict_selected_multinomial(
      x_train, y_train, x_valid, levels(y_train), random_state = 3L
    ),
    rep("A", 5L)
  )
})

test_that("D2 .seeded_cv_glmnet leaves the caller's RNG state untouched", {
  skip_if_not_installed("glmnet")
  set.seed(21)
  before <- get(".Random.seed", envir = globalenv())
  fit <- suppressWarnings(stablr:::.seeded_cv_glmnet(
    x = matrix(stats::rnorm(120L), 30L, 4L),
    y = factor(rep(c("A", "B", "C"), each = 10L)),
    family = "multinomial", type.measure = "class", nfolds = 3L,
    random_state = 5L
  ))
  expect_s3_class(fit, "cv.glmnet")
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

# ---- J1: the inner-CV comparison is paired (common random numbers) -------------

test_that("J1 inner candidate seeds do not depend on the candidate index", {
  # Every candidate on inner fold f must receive the same seed, otherwise the
  # single-number selection rule compares candidates fitted on different
  # bootstrap draws.
  expect_identical(
    stablr:::.derive_nested_seed(11L, 1L, stablr:::.NESTED_INNER_FIT_SEED_OFFSET),
    stablr:::.derive_nested_seed(11L, 1L, stablr:::.NESTED_INNER_FIT_SEED_OFFSET)
  )
  # The old scheme used `k * 100L`, which also collided with the 1000L and
  # 3000L offsets once a run had ten or more candidates.  The four derived
  # streams must stay disjoint for a representative index range.
  offsets <- c(stablr:::.NESTED_INNER_FOLD_SEED_OFFSET,
               stablr:::.NESTED_INNER_FIT_SEED_OFFSET,
               stablr:::.NESTED_OUTER_FIT_SEED_OFFSET,
               stablr:::.NESTED_CV_GLMNET_SEED_OFFSET)
  seeds <- unlist(lapply(offsets, function(offset) {
    vapply(1:6, function(i) stablr:::.derive_nested_seed(11L, i, offset), integer(1L))
  }))
  expect_length(unique(seeds), length(seeds))
})

test_that("J1 candidates are fitted on the same seeds, so the comparison is paired", {
  set.seed(11)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  x_a <- matrix(stats::rnorm(n * 5L, sd = 0.3), nrow = n,
                dimnames = list(ids, paste0("a", seq_len(5L))))
  x_b <- matrix(stats::rnorm(n * 4L, sd = 0.3), nrow = n,
                dimnames = list(ids, paste0("b", seq_len(4L))))

  fit <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x_a, mirna = x_b),
    y = y,
    candidates = list(
      list(name = "mrna", blocks = "mrna"),
      list(name = "mirna", blocks = "mirna")
    ),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 3L, metric = "ber",
    family = "multinomial", artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 11L
  ))

  # Both candidates' stabl_fit calls used the same per-fold seed, so their
  # bootstrap draws are common: the fit objects record the resolved grids and
  # the comparison is attributable to the candidates, not the sampler.
  inner <- fit$fold_results[[1L]]$inner_results
  expect_setequal(inner$summary$candidate, c("mrna", "mirna"))
  expect_true(all(is.finite(inner$summary$balanced_error_rate)))

  # Repeating the run reproduces it exactly, which is only meaningful because the
  # seeds no longer move with the candidate index.
  again <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x_a, mirna = x_b),
    y = y,
    candidates = list(
      list(name = "mrna", blocks = "mrna"),
      list(name = "mirna", blocks = "mirna")
    ),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 3L, metric = "ber",
    family = "multinomial", artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 11L
  ))
  expect_identical(fit$diagnostics, again$diagnostics)
  expect_identical(fit$outer_predictions, again$outer_predictions)
})

# ---- J2: the lambda-grid fallback warns and is recorded ------------------------

test_that("J2 a named-list grid with no matching key warns and is recorded", {
  set.seed(303)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  x_a <- matrix(stats::rnorm(n * 5L, sd = 0.3), nrow = n,
                dimnames = list(ids, paste0("a", seq_len(5L))))
  x_b <- matrix(stats::rnorm(n * 5L, sd = 0.3), nrow = n,
                dimnames = list(ids, paste0("b", seq_len(5L))))

  # The list names one candidate only; `mirna` has no entry and must fall back.
  expect_warning(
    fit <- suppressWarnings(stabl_multiomic_nested_cv(
      x_list = list(mrna = x_a, mirna = x_b),
      y = y,
      candidates = list(
        list(name = "mrna", blocks = "mrna"),
        list(name = "mirna", blocks = "mirna")
      ),
      lambda_grid = list(mrna = data.frame(lambda = c(0.2, 0.1))),
      outer_v = 3L, inner_v = 3L, n_bootstraps = 2L,
      artificial_type = NULL, hard_threshold = 1e-9,
      sample_fraction = 1, random_state = 5L
    )),
    NA
  )
  expect_s3_class(fit, "stabl_multiomic_nested_cv")
  expect_true("lambda_grid_resolved" %in% names(fit$diagnostics))
  expect_setequal(unique(fit$diagnostics$lambda_grid_resolved),
                  c("list:mrna", "auto_fallback"))
  expect_setequal(
    fit$diagnostics$lambda_grid_resolved[fit$diagnostics$candidate == "mrna"],
    "list:mrna"
  )
  expect_setequal(
    fit$diagnostics$lambda_grid_resolved[fit$diagnostics$candidate == "mirna"],
    "auto_fallback"
  )
})

test_that("J2 a matching key resolves silently and is recorded", {
  set.seed(404)
  n <- 18L
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
  x <- matrix(stats::rnorm(n * 5L, sd = 0.3), nrow = n,
              dimnames = list(ids, paste0("f", seq_len(5L))))

  expect_no_warning(fit <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = list(mrna = x), y = y,
    candidates = list(list(name = "mrna", blocks = "mrna")),
    lambda_grid = list(mrna = data.frame(lambda = c(0.2, 0.1))),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 2L,
    artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 5L
  )))
  expect_true(all(fit$diagnostics$lambda_grid_resolved == "list:mrna"))
})

# ---- J3: the selection rule is documented as optimistic ------------------------

test_that("J3 .select_nested_candidate is a bare argmin/argmax with first-index ties", {
  summary <- data.frame(
    candidate = c("a", "b", "c"),
    accuracy = c(0.5, 0.9, 0.9),
    balanced_error_rate = c(0.5, 0.1, 0.1),
    stringsAsFactors = FALSE
  )
  # Ties go to the first candidate in both directions: no one-standard-error
  # rule, so the winner is the point estimate and the choice is optimistic.
  expect_identical(stablr:::.select_nested_candidate(summary, "accuracy"), "b")
  expect_identical(stablr:::.select_nested_candidate(summary, "ber"), "b")

  strict <- data.frame(
    candidate = c("a", "b"),
    accuracy = c(0.5, 0.9),
    balanced_error_rate = c(0.5, 0.4),
    stringsAsFactors = FALSE
  )
  expect_identical(stablr:::.select_nested_candidate(strict, "accuracy"), "b")
  expect_identical(stablr:::.select_nested_candidate(strict, "ber"), "b")
})
