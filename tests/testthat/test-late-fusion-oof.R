test_that("OOF late fusion is invariant to named outcome order", {
  set.seed(2207)
  ids <- paste0("s", seq_len(18L))
  x1 <- matrix(rnorm(18L * 4L), 18L,
               dimnames = list(ids, paste0("a", 1:4)))
  x2 <- matrix(rnorm(18L * 3L), 18L,
               dimnames = list(ids, paste0("b", 1:3)))
  y <- setNames(0.7 * x1[, 1L] + rnorm(18L), ids)
  groups <- setNames(rep(paste0("g", 1:6), each = 3L), ids)
  strata <- setNames(rep(rep(c("site1", "site2", "site3"), each = 3L), 2L), ids)
  args <- list(
    x_train_list = list(a = x1, b = x2),
    groups_train = groups, bootstrap_strata_train = strata,
    lambda_grid = data.frame(lambda = 0.1),
    artificial_type = NULL, hard_threshold = 1,
    n_bootstraps = 2L, sample_fraction = 1,
    late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 20L, random_state = 77L
  )
  ordered <- do.call(stabl_multiomic_train_validate, c(args, list(y_train = y)))
  shuffled_args <- args
  shuffled_args$y_train <- y[sample(names(y))]
  shuffled_args$groups_train <- groups[sample(names(groups))]
  shuffled_args$bootstrap_strata_train <- strata[sample(names(strata))]
  shuffled <- do.call(stabl_multiomic_train_validate, shuffled_args)

  expect_identical(ordered$late_fusion$weights, shuffled$late_fusion$weights)
  expect_identical(ordered$late_fusion$train_predictions,
                   shuffled$late_fusion$train_predictions)
  expect_identical(ordered$late_fusion$provenance$fold_id,
                   shuffled$late_fusion$provenance$fold_id)
  expect_identical(
    ordered$late_fusion$oof_predictions,
    ordered$late_fusion$train_predictions
  )
  expect_identical(ordered$late_fusion$oof_score, ordered$late_fusion$score)
  expect_false("in_sample_predictions" %in% names(ordered$late_fusion))
  expect_identical(
    ordered$late_fusion$provenance$result_semantics$compatibility_aliases,
    c(train_predictions = "oof_predictions", score = "oof_score")
  )
})

test_that("OOF late fusion covers samples once and separates groups", {
  set.seed(2208)
  ids <- paste0("s", seq_len(18L))
  x <- matrix(rnorm(18L * 3L), 18L,
              dimnames = list(ids, paste0("f", 1:3)))
  y <- setNames(rnorm(18L), ids)
  groups <- setNames(rep(paste0("g", 1:6), each = 3L), ids)
  fit <- stabl_multiomic_train_validate(
    list(a = x), y, data.frame(lambda = 0.1),
    groups_train = groups, artificial_type = NULL, hard_threshold = 1,
    n_bootstraps = 2L, sample_fraction = 1,
    late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 10L, random_state = 78L
  )
  provenance <- fit$late_fusion$provenance
  expect_identical(names(provenance$fold_id), ids)
  expect_true(all(provenance$fold_id %in% 1:3))
  oof_columns <- setdiff(names(fit$late_fusion$train_predictions),
                         "Stacked Gen. Predictions")
  expect_true(all(is.finite(as.matrix(
    fit$late_fusion$train_predictions[, oof_columns, drop = FALSE]
  ))))
  for (fold in provenance$folds) {
    expect_false(any(groups[fold$train_ids] %in% groups[fold$assessment_ids]))
    expect_setequal(intersect(fold$train_ids, fold$assessment_ids), character())
  }
})

test_that("legacy late fusion remains explicitly available", {
  ids <- paste0("s", 1:12)
  x <- matrix(seq_len(36), 12, 3,
              dimnames = list(ids, paste0("f", 1:3)))
  y <- setNames(seq_len(12), ids)
  fit <- stabl_multiomic_train_validate(
    list(a = x), y, data.frame(lambda = 0.1),
    artificial_type = NULL, hard_threshold = 1,
    n_bootstraps = 2L, sample_fraction = 1,
    late_fusion = TRUE, late_fusion_training = "python_legacy",
    n_iter_lf = 10L, random_state = 79L
  )
  expect_identical(fit$late_fusion$provenance$training_mode, "python_legacy")
  expect_identical(
    fit$late_fusion$in_sample_predictions,
    fit$late_fusion$train_predictions
  )
  expect_identical(
    fit$late_fusion$in_sample_score,
    fit$late_fusion$score
  )
  expect_false("oof_predictions" %in% names(fit$late_fusion))
  expect_identical(
    fit$late_fusion$provenance$result_semantics$compatibility_aliases,
    c(
      train_predictions = "in_sample_predictions",
      score = "in_sample_score"
    )
  )
})

test_that("stacking rejects malformed scalar outcomes and infinite predictions", {
  p <- cbind(a = 1:4, b = 4:1)
  expect_error(stacked_multi_omic(p, 1:3, "regression"), "one value")
  expect_error(stacked_multi_omic(p, c(0, 0, 0, 0), "binary"), "both event")
  p[1, 1] <- Inf
  expect_error(stacked_multi_omic(p, 1:4, "regression"), "finite")
  p <- cbind(a = c(NA, 2:4), b = c(NA, 3:1))
  stacked <- stacked_multi_omic(p, 1:4, "regression", n_iter = 2L,
                                random_state = 1L)
  expect_true(is.na(stacked$predictions[["Stacked Gen. Predictions"]][[1L]]))
})

test_that("late fusion validates training outcomes without a validation split", {
  ids <- paste0("s", 1:8)
  x <- matrix(rnorm(16), 8, 2, dimnames = list(ids, c("x1", "x2")))
  common <- list(
    x_train_list = list(a = x), lambda_grid = data.frame(lambda = 0.1),
    artificial_type = NULL, hard_threshold = 1, n_bootstraps = 2L,
    late_fusion = TRUE, late_fusion_nfolds = 2L, n_iter_lf = 2L
  )
  expect_error(do.call(stabl_multiomic_train_validate, c(
    common, list(y_train = setNames(letters[1:8], ids), family = "gaussian")
  )), "finite numeric")
  expect_error(do.call(stabl_multiomic_train_validate, c(
    common, list(y_train = setNames(factor(rep("A", 8)), ids),
                 family = "multinomial")
  )), "at least two")
})

test_that("OOF binomial fusion maps named factor events deterministically", {
  set.seed(2210)
  ids <- paste0("s", seq_len(18L))
  x <- matrix(rnorm(18L * 3L), 18L,
              dimnames = list(ids, paste0("f", 1:3)))
  y <- setNames(factor(rep(c("control", "event"), 9L),
                       levels = c("control", "event")), ids)
  fit <- suppressWarnings(stabl_multiomic_train_validate(
    list(a = x), y, data.frame(lambda = 0.1), family = "binomial",
    artificial_type = NULL, hard_threshold = 1,
    n_bootstraps = 2L, sample_fraction = 1,
    late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 10L, random_state = 80L
  ))
  expect_identical(
    fit$late_fusion$provenance$binary_event_mapping,
    c(control = 0L, event = 1L)
  )
})

test_that("multiclass stacking rejects incomplete or unnormalised probabilities", {
  ids <- paste0("s", 1:3)
  y <- factor(c("A", "B", "C"))
  probs <- matrix(
    c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
    3L, byrow = TRUE, dimnames = list(ids, levels(y))
  )
  bad_na <- probs
  bad_na[1, 1] <- NA_real_
  bad_sum <- probs
  bad_sum[1, ] <- c(0.2, 0.2, 0.2)
  expect_error(stacked_multi_omic(list(a = bad_na), y, "multiclass"),
               "complete and finite")
  expect_error(stacked_multi_omic(list(a = bad_sum), y, "multiclass"),
               "sum to one")
})

test_that("OOF assessment predictions do not use their fold outcomes", {
  set.seed(2211)
  ids <- paste0("s", seq_len(18L))
  x <- matrix(rnorm(18L * 3L), 18L,
              dimnames = list(ids, paste0("f", 1:3)))
  y <- setNames(rnorm(18L), ids)
  args <- list(
    x_train_list = list(a = x), lambda_grid = data.frame(lambda = 0.1),
    artificial_type = NULL, hard_threshold = 1, n_bootstraps = 2L,
    sample_fraction = 1, late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 10L, random_state = 81L
  )
  original <- do.call(stabl_multiomic_train_validate, c(args, list(y_train = y)))
  fold_one <- names(original$late_fusion$provenance$fold_id)[
    original$late_fusion$provenance$fold_id == 1L
  ]
  perturbed_y <- y
  perturbed_y[fold_one] <- perturbed_y[fold_one] + 1000
  perturbed <- do.call(
    stabl_multiomic_train_validate, c(args, list(y_train = perturbed_y))
  )
  expect_equal(
    original$late_fusion$train_predictions[fold_one, "a"],
    perturbed$late_fusion$train_predictions[fold_one, "a"]
  )
})

test_that("validation outcomes cannot influence OOF weights or predictions", {
  set.seed(2212)
  train_ids <- paste0("tr", seq_len(18L))
  valid_ids <- paste0("va", seq_len(8L))
  x_train <- matrix(rnorm(18L * 3L), 18L,
                    dimnames = list(train_ids, paste0("f", 1:3)))
  x_valid <- matrix(rnorm(8L * 3L), 8L,
                    dimnames = list(valid_ids, paste0("f", 1:3)))
  y_train <- setNames(rnorm(18L), train_ids)
  y_valid <- setNames(rnorm(8L), valid_ids)
  args <- list(
    x_train_list = list(a = x_train), y_train = y_train,
    lambda_grid = data.frame(lambda = 0.1), x_valid_list = list(a = x_valid),
    artificial_type = NULL, hard_threshold = 1, n_bootstraps = 2L,
    sample_fraction = 1, late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 10L, random_state = 82L
  )
  first <- do.call(stabl_multiomic_train_validate, c(args, list(y_valid = y_valid)))
  second <- do.call(
    stabl_multiomic_train_validate,
    c(args, list(y_valid = setNames(rev(unname(y_valid)), valid_ids)))
  )
  expect_equal(first$late_fusion$weights, second$late_fusion$weights)
  expect_equal(first$late_fusion$valid_predictions,
               second$late_fusion$valid_predictions)
})

test_that("validation omics are finite numeric matrices with one sample-ID set", {
  ids <- paste0("s", 1:8)
  valid_ids <- paste0("v", 1:4)
  x <- matrix(rnorm(16), 8, 2, dimnames = list(ids, c("x1", "x2")))
  y <- setNames(rnorm(8), ids)
  base <- list(
    x_train_list = list(a = x, b = x), y_train = y,
    lambda_grid = data.frame(lambda = 0.1),
    artificial_type = NULL, hard_threshold = 1, n_bootstraps = 2L
  )
  va <- matrix(rnorm(8), 4, 2, dimnames = list(valid_ids, c("x1", "x2")))
  vb <- va[rev(valid_ids), , drop = FALSE]
  expect_silent(do.call(
    stabl_multiomic_train_validate,
    c(base, list(x_valid_list = list(a = va, b = vb)))
  ))
  expect_error(do.call(
    stabl_multiomic_train_validate,
    c(base, list(x_valid_list = list(a = va, b = vb[-1, , drop = FALSE])))
  ), "same sample-ID set")
  bad <- va
  bad[1, 1] <- Inf
  expect_error(do.call(
    stabl_multiomic_train_validate,
    c(base, list(x_valid_list = list(a = va, b = bad)))
  ), "finite")
  bad_text <- as.data.frame(va)
  bad_text[[1]] <- letters[1:4]
  expect_error(do.call(
    stabl_multiomic_train_validate,
    c(base, list(x_valid_list = list(a = va, b = bad_text)))
  ), "numeric")
})

test_that("late-fusion validation labels must belong to training classes", {
  ids <- paste0("s", 1:12)
  valid_ids <- paste0("v", 1:4)
  x <- matrix(rnorm(24), 12, 2, dimnames = list(ids, c("x1", "x2")))
  xv <- matrix(rnorm(8), 4, 2, dimnames = list(valid_ids, c("x1", "x2")))
  common <- list(
    x_train_list = list(a = x), lambda_grid = data.frame(lambda = 0.1),
    x_valid_list = list(a = xv), artificial_type = NULL,
    n_bootstraps = 2L, late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 2L
  )
  y_binary <- setNames(factor(rep(c("control", "event"), 6)), ids)
  expect_error(do.call(
    stabl_multiomic_train_validate,
    c(common, list(y_train = y_binary, y_valid = setNames(
      factor(c("control", "event", "other", "control")), valid_ids
    ), family = "binomial"))
  ), "training binary event levels")
  y_multi <- setNames(factor(rep(c("A", "B", "C"), 4)), ids)
  expect_error(do.call(
    stabl_multiomic_train_validate,
    c(common, list(y_train = y_multi, y_valid = setNames(
      factor(c("A", "B", "D", "A")), valid_ids
    ), family = "multinomial"))
  ), "training multiclass levels")
})

test_that("classification metrics never drop rows and preserve declared levels", {
  truth <- factor(c("A", "A"), levels = c("A", "B", "C"))
  predicted <- factor(c("A", "C"), levels = c("A", "B", "C"))
  metrics <- stablr:::.classification_metrics(truth, predicted)
  expect_equal(sum(metrics$confusion), 2)
  expect_identical(dim(metrics$confusion), c(3L, 3L))
  expect_error(stablr:::.classification_metrics(c("A", NA), c("A", "A")),
               "complete")
  expect_error(stablr:::.classification_metrics(c("A", "B"), c("A")),
               "equal length")
})

test_that("OOF late fusion is identical for sequential and parallel bootstraps", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  skip_if_not(.can_open_test_server_socket(),
              "future multisession requires opening a local server socket")
  set.seed(2213)
  ids <- paste0("s", seq_len(18L))
  x <- matrix(rnorm(18L * 3L), 18L,
              dimnames = list(ids, paste0("f", 1:3)))
  y <- setNames(rnorm(18L), ids)
  args <- list(
    x_train_list = list(a = x), y_train = y,
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    artificial_type = "random_permutation", n_bootstraps = 4L,
    late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 10L, random_state = 83L
  )
  sequential <- do.call(stabl_multiomic_train_validate,
                        c(args, list(workers = 1L)))
  old_plan <- future::plan(future::multisession, workers = 2L)
  withr::defer(future::plan(old_plan))
  parallel <- suppressWarnings(do.call(
    stabl_multiomic_train_validate, c(args, list(workers = 2L))
  ))
  expect_identical(sequential$late_fusion$weights,
                   parallel$late_fusion$weights)
  expect_identical(sequential$late_fusion$train_predictions,
                   parallel$late_fusion$train_predictions)
  expect_identical(sequential$late_fusion$valid_predictions,
                   parallel$late_fusion$valid_predictions)
})

test_that("OOF selector configuration errors are rethrown", {
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(24), 12, 2, dimnames = list(ids, c("x1", "x2")))
  y <- setNames(rnorm(12), ids)
  expect_error(
    stabl_multiomic_train_validate(
      x_train_list = list(a = x), y_train = y,
      lambda_grid = data.frame(lambda = 0.1),
      base_learner = "not-a-learner", family = "gaussian",
      n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 1,
      late_fusion = TRUE, late_fusion_nfolds = 3L,
      n_iter_lf = 2L, random_state = 1L
    ),
    "base_learner"
  )
})

test_that("OOF selector programming errors escape the public workflow", {
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(24), 12, 2, dimnames = list(ids, c("x1", "x2")))
  y <- setNames(rnorm(12), ids)
  testthat::local_mocked_bindings(
    .fit_multiomic_per_omic = function(...) {
      stop("deliberate selector programming defect")
    },
    .package = "stablr"
  )

  expect_error(
    stabl_multiomic_train_validate(
      x_train_list = list(a = x), y_train = y,
      lambda_grid = data.frame(lambda = 0.1),
      n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 1,
      late_fusion = TRUE, late_fusion_nfolds = 3L,
      n_iter_lf = 2L, random_state = 2L
    ),
    "deliberate selector programming defect"
  )
})

test_that("OOF typed selector infeasibility uses a structured fallback", {
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(24), 12, 2, dimnames = list(ids, c("x1", "x2")))
  y <- setNames(rnorm(12), ids)
  x_valid <- matrix(
    rnorm(4), 2, 2,
    dimnames = list(c("v1", "v2"), colnames(x))
  )
  testthat::local_mocked_bindings(
    .fit_multiomic_per_omic = function(...) {
      stablr:::.abort_numerical_infeasibility(
        "stablr_learner_numerical_infeasibility",
        "deliberate selector numerical infeasibility"
      )
    },
    .package = "stablr"
  )

  fit <- stabl_multiomic_train_validate(
    x_train_list = list(a = x), y_train = y,
    lambda_grid = data.frame(lambda = 0.1),
    x_valid_list = list(a = x_valid),
    n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 1,
    late_fusion = TRUE, late_fusion_nfolds = 3L,
    n_iter_lf = 2L, random_state = 3L
  )

  captured <- fit$late_fusion$provenance$folds[[1L]]
  expect_s3_class(
    captured$selector_conditions$a,
    "stablr_numerical_infeasibility"
  )
  expect_match(captured$selector_errors$a, "deliberate selector")
  expect_identical(captured$selected_features$a, character())
  expect_match(captured$fallback_reasons[["a"]], "selector_fit_error")
  expect_equal(
    fit$late_fusion$valid_predictions,
    setNames(rep(mean(y), nrow(x_valid)), rownames(x_valid)),
    ignore_attr = TRUE
  )
  full <- fit$late_fusion$provenance$full_refit
  expect_s3_class(
    full$selector_conditions$a,
    "stablr_numerical_infeasibility"
  )
  expect_identical(full$fallback_conditions$a, full$selector_conditions$a)
})

test_that("OOF downstream programming errors are rethrown", {
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(24), 12, 2, dimnames = list(ids, c("x1", "x2")))
  y <- setNames(rnorm(12), ids)
  testthat::local_mocked_bindings(
    .late_fusion_fit_omic = function(...) {
      stop("deliberate downstream programming defect")
    },
    .package = "stablr"
  )

  expect_error(
    stablr:::.late_fusion_fit_omic_safe(x, y, NULL, "regression"),
    "deliberate downstream programming defect"
  )
})

test_that("OOF downstream numerical infeasibility uses a typed fallback", {
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(24), 12, 2, dimnames = list(ids, c("x1", "x2")))
  y <- setNames(rnorm(12), ids)
  original_fit <- stablr:::.late_fusion_fit_omic
  testthat::local_mocked_bindings(
    .late_fusion_fit_omic = function(x_train_sel, y_train,
                                     x_valid_sel = NULL, y_train_mean,
                                     task_type, levels = NULL) {
      if (ncol(x_train_sel) > 0L) {
        stablr:::.abort_numerical_infeasibility(
          "stablr_downstream_numerical_infeasibility",
          "deliberate downstream numerical infeasibility"
        )
      }
      original_fit(
        x_train_sel, y_train, x_valid_sel, y_train_mean, task_type, levels
      )
    },
    .package = "stablr"
  )

  pred <- stablr:::.late_fusion_fit_omic_safe(x, y, NULL, "regression")
  expect_equal(pred$train_preds, rep(mean(y), length(y)))
  expect_match(pred$fallback_reason, "downstream_fit_error")
  expect_s3_class(
    pred$fallback_condition,
    "stablr_downstream_numerical_infeasibility"
  )
})

# ---------------------------------------------------------------------------
# H8 / H9: the downstream multinomial fit must not read ambient RNG state, and
# both fusion branches must encode the same event from the same `y`.
# ---------------------------------------------------------------------------

.lf_multiclass_design <- function(n = 45L, seed = 8100L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  x <- matrix(rnorm(n * 4L), n, dimnames = list(ids, paste0("f", 1:4)))
  classes <- rep(c("A", "B", "C"), length.out = n)
  for (level in c("A", "B", "C")) {
    x[classes == level, 1:2] <- x[classes == level, 1:2] + 2
  }
  y <- factor(classes, levels = c("A", "B", "C"))
  names(y) <- ids
  list(x = x, y = y)
}

test_that("H8 a seeded downstream multinomial fit ignores ambient RNG state", {
  d <- .lf_multiclass_design()
  fit_with <- function(ambient_seed) {
    set.seed(ambient_seed)
    stablr:::.late_fusion_fit_omic_safe(
      x_train_sel = d$x, y_train = d$y, x_valid_sel = NULL,
      task_type = "multiclass", levels = levels(d$y), random_state = 4242L
    )
  }

  first <- fit_with(1L)
  second <- fit_with(987654L)
  expect_identical(first$train_preds, second$train_preds)
  expect_null(first$fallback_reason)
})

test_that("H8 an unseeded downstream multinomial fit still reads ambient state", {
  # The honest consequence of random_state = NULL: the caller asked for an
  # unseeded run, so the fit is free to depend on ambient RNG state.
  d <- .lf_multiclass_design()
  fit_with <- function(ambient_seed) {
    set.seed(ambient_seed)
    stablr:::.late_fusion_fit_omic_safe(
      x_train_sel = d$x, y_train = d$y, x_valid_sel = NULL,
      task_type = "multiclass", levels = levels(d$y)
    )
  }
  expect_false(identical(fit_with(1L)$train_preds, fit_with(999L)$train_preds))
})

test_that("H8 seeded python_legacy multinomial late fusion is RNG independent", {
  d <- .lf_multiclass_design(n = 30L)
  ids <- rownames(d$x)
  run <- function(ambient_seed) {
    set.seed(ambient_seed)
    suppressWarnings(stabl_multiomic_train_validate(
      x_train_list = list(a = d$x[, 1:3, drop = FALSE],
                          b = d$x[, c(1L, 2L, 4L), drop = FALSE]),
      y_train = d$y, lambda_grid = data.frame(lambda = 0.05),
      family = "multinomial", artificial_type = NULL, hard_threshold = 1e-9,
      n_bootstraps = 2L, sample_fraction = 1, late_fusion = TRUE,
      late_fusion_training = "python_legacy", n_iter_lf = 20L,
      random_state = 5150L
    ))
  }
  first <- run(2L)
  second <- run(7777L)
  expect_identical(first$late_fusion$train_predictions,
                   second$late_fusion$train_predictions)
  expect_identical(first$late_fusion$weights, second$late_fusion$weights)
  expect_identical(first$late_fusion$score, second$late_fusion$score)
  expect_identical(rownames(first$late_fusion$train_predictions), ids)
})

test_that("H8 seeded OOF multinomial late fusion is RNG independent", {
  d <- .lf_multiclass_design(n = 30L)
  run <- function(ambient_seed) {
    set.seed(ambient_seed)
    stabl_multiomic_train_validate(
      x_train_list = list(a = d$x[, 1:3, drop = FALSE],
                          b = d$x[, c(1L, 2L, 4L), drop = FALSE]),
      y_train = d$y, lambda_grid = data.frame(lambda = 0.05),
      family = "multinomial", artificial_type = NULL, hard_threshold = 1e-9,
      n_bootstraps = 2L, sample_fraction = 1, late_fusion = TRUE,
      late_fusion_nfolds = 3L, n_iter_lf = 10L, random_state = 5151L
    )
  }
  first <- run(3L)
  second <- run(8888L)
  expect_identical(first$late_fusion$train_predictions,
                   second$late_fusion$train_predictions)
  expect_identical(first$late_fusion$valid_predictions,
                   second$late_fusion$valid_predictions)
  expect_identical(first$late_fusion$score, second$late_fusion$score)
})

test_that("H9 one mapper serves both branches and returns the mapping", {
  expect_false(exists(".cooperative_coerce_binary_outcome",
                      envir = asNamespace("stablr"), inherits = FALSE))

  logical_y <- c(TRUE, FALSE, TRUE)
  numeric_y <- c(0, 1, 0, 1)
  factor_sorted <- factor(c("control", "event", "control", "event"))
  factor_declared <- factor(c("control", "event", "control", "event"),
                            levels = c("event", "control"))

  expect_identical(
    stablr:::.multiomic_binary_outcome(logical_y)$y,
    c(1L, 0L, 1L)
  )
  expect_identical(
    stablr:::.multiomic_binary_outcome(logical_y)$mapping,
    c("FALSE" = 0L, "TRUE" = 1L)
  )
  expect_identical(
    stablr:::.multiomic_binary_outcome(numeric_y)$y,
    c(0L, 1L, 0L, 1L)
  )
  expect_identical(
    stablr:::.multiomic_binary_outcome(factor_sorted)$mapping,
    c(control = 0L, event = 1L)
  )
  # A declared level order is the caller's event declaration, not a hint to be
  # re-sorted: the second declared level is the event.
  expect_identical(
    stablr:::.multiomic_binary_outcome(factor_declared)$mapping,
    c(event = 0L, control = 1L)
  )
  expect_error(
    stablr:::.multiomic_binary_outcome(c(0, 1, 2)),
    "exactly two complete event levels"
  )
  expect_error(
    stablr:::.multiomic_binary_outcome(c(0L, 0L, 0L)),
    "exactly two complete event levels"
  )
})

test_that("H9 the cooperative AUC metric maps a labelled outcome, not NA", {
  # A character validation outcome used to reach `sort(unique(as.numeric(y)))`,
  # which drops the NA and leaves an empty vector, so every candidate was
  # scored as all-NA and then as a no-information 0.5 AUC.
  truth <- c("event", "control", "event", "control")
  aligned <- c(0.8, 0.2, 0.7, 0.3)
  for (true_y in list(factor(truth, levels = c("control", "event")), truth)) {
    expect_equal(
      stablr:::.cooperative_validation_metric(
        true_y = true_y, pred_y = aligned, family = "binomial",
        type_measure = "auc"
      ),
      1
    )
  }
  # Reversing the ranking reverses the AUC, which is what proves the two
  # branches now encode the same event rather than opposite ones.
  expect_equal(
    stablr:::.cooperative_validation_metric(
      true_y = truth, pred_y = rev(aligned), family = "binomial",
      type_measure = "auc"
    ),
    0
  )
  # A single-class validation outcome cannot tune an AUC at all; scoring it at
  # the no-information 0.5 would silently select a candidate.
  expect_error(
    stablr:::.cooperative_validation_metric(
      true_y = rep("event", 4L), pred_y = aligned, family = "binomial",
      type_measure = "auc"
    ),
    "exactly two complete event levels"
  )
})

test_that("H9 a declared non-alphabetical event order is honoured end to end", {
  d <- .lf_multiclass_design(n = 24L, seed = 8101L)
  ids <- rownames(d$x)
  set.seed(8101L)
  events <- rbinom(length(ids), 1L, 0.5)
  y <- factor(ifelse(events == 1L, "event", "control"),
              levels = c("event", "control"))
  names(y) <- ids
  x <- matrix(rnorm(length(ids) * 4L), length(ids),
              dimnames = list(ids, paste0("f", 1:4)))

  fit <- suppressWarnings(stabl_multiomic_train_validate(
    x_train_list = list(a = x[, 1:2], b = x[, 3:4]), y_train = y,
    lambda_grid = data.frame(lambda = 0.1), family = "binomial",
    artificial_type = NULL, hard_threshold = 1, n_bootstraps = 2L,
    sample_fraction = 1, random_state = 8102L, late_fusion = TRUE,
    late_fusion_nfolds = 3L, n_iter_lf = 10L
  ))
  expect_identical(
    fit$late_fusion$provenance$binary_event_mapping,
    c(event = 0L, control = 1L)
  )
})
