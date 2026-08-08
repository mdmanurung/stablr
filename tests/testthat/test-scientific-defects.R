test_that("SCI-01 classification repairs guarantee two observations per class", {
  y <- factor(c(rep("major", 8L), rep("minor", 2L)))

  draw_sequence <- function() {
    sampler <- stablr:::.classification_bootstrap_sampler(
      y = y,
      groups = NULL,
      n_subsamples = 5L,
      replace = FALSE,
      strata_ids = NULL,
      family = "binomial",
      class_levels = levels(y)
    )
    set.seed(741L)
    lapply(seq_len(20L), sampler)
  }

  first <- draw_sequence()
  second <- draw_sequence()
  expect_identical(first, second)
  expect_true(all(vapply(first, function(indices) {
    all(table(factor(y[indices], levels = levels(y))) >= 2L)
  }, logical(1L))))
})

test_that("SCI-01 minimum stratified allocation is deterministic", {
  strata <- factor(c(rep("major", 8L), rep("minor", 2L)))
  counts <- stablr:::.stratified_counts(
    strata,
    n_subsamples = 6L,
    replace = FALSE,
    min_per_stratum = 2L
  )

  expect_identical(counts, c(major = 4L, minor = 2L))
})

test_that("SCI-01 impossible bootstrap size raises a typed condition", {
  set.seed(742L)
  x <- matrix(rnorm(18L), nrow = 6L, ncol = 3L)
  rownames(x) <- paste0("s", seq_len(nrow(x)))
  y <- factor(rep(c("A", "B", "C"), each = 2L))
  names(y) <- rownames(x)

  condition <- tryCatch(
    stabl_fit(
      x,
      y,
      lambda_grid = data.frame(lambda = 0.1),
      family = "multinomial",
      n_bootstraps = 1L,
      artificial_type = NULL,
      hard_threshold = 0.5,
      sample_fraction = 0.5,
      random_state = 742L
    ),
    error = identity
  )

  expect_s3_class(condition, "stablr_bootstrap_infeasible")
  expect_s3_class(condition, "stablr_error")
  expect_identical(condition$n_subsamples, 3L)
  expect_identical(condition$minimum_per_class, 2L)
  expect_match(conditionMessage(condition), "every class", fixed = TRUE)
})

test_that("SCI-01 sparse training classes fail before learner fitting", {
  set.seed(743L)
  x <- matrix(rnorm(28L), nrow = 7L, ncol = 4L)
  rownames(x) <- paste0("s", seq_len(nrow(x)))
  y <- factor(c(rep("A", 3L), rep("B", 3L), "C"))
  names(y) <- rownames(x)

  expect_error(
    stabl_fit(
      x,
      y,
      lambda_grid = data.frame(lambda = 0.1),
      family = "multinomial",
      n_bootstraps = 1L,
      artificial_type = NULL,
      hard_threshold = 0.5,
      sample_fraction = 1,
      random_state = 743L
    ),
    class = "stablr_bootstrap_infeasible"
  )
})

test_that("SCI-01 imbalanced binomial fitting is deterministic and feasible", {
  set.seed(744L)
  x <- matrix(rnorm(80L), nrow = 20L, ncol = 4L)
  rownames(x) <- paste0("s", seq_len(nrow(x)))
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  y <- factor(c(rep("major", 18L), rep("minor", 2L)))
  names(y) <- rownames(x)
  args <- list(
    x = x,
    y = y,
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    family = "binomial",
    n_bootstraps = 4L,
    artificial_type = NULL,
    hard_threshold = 0.5,
    sample_fraction = 0.25,
    replace = FALSE,
    random_state = 744L
  )

  first <- suppressWarnings(do.call(stabl_fit, args))
  second <- suppressWarnings(do.call(stabl_fit, args))
  expect_s3_class(first, "stabl_fit")
  expect_identical(first$stabl_scores_, second$stabl_scores_)
})

test_that("SCI-02 named stacking aligns outcomes to prediction rows", {
  ids <- paste0("s", seq_len(8L))
  y_ordered <- rep(c(0, 1), 4L)
  predictions <- matrix(
    c(y_ordered + 0.1, rev(y_ordered) + 0.2),
    nrow = length(ids),
    dimnames = list(ids, c("omic_a", "omic_b"))
  )
  y_named <- stats::setNames(y_ordered, ids)
  shuffled <- y_named[c(5L, 2L, 8L, 1L, 7L, 4L, 6L, 3L)]

  expected <- stacked_multi_omic(
    predictions,
    y_named,
    task_type = "binary",
    n_iter = 100L,
    random_state = 802L
  )
  actual <- stacked_multi_omic(
    predictions,
    shuffled,
    task_type = "binary",
    n_iter = 100L,
    random_state = 802L
  )

  expect_identical(actual, expected)
  expect_identical(rownames(actual$predictions), ids)
})

test_that("SCI-02 multiclass stacking aligns every named prediction matrix", {
  ids <- paste0("s", seq_len(6L))
  y <- factor(rep(c("A", "B", "C"), each = 2L))
  names(y) <- ids
  probabilities <- matrix(
    0.05,
    nrow = length(ids),
    ncol = 3L,
    dimnames = list(ids, levels(y))
  )
  probabilities[cbind(seq_along(y), as.integer(y))] <- 0.9
  reversed <- probabilities[rev(ids), , drop = FALSE]

  aligned <- stacked_multi_omic(
    list(omic_a = probabilities, omic_b = reversed),
    y[rev(ids)],
    task_type = "multiclass",
    n_iter = 30L,
    random_state = 803L
  )

  expect_identical(rownames(aligned$predictions), ids)
  expect_identical(aligned$predictions$predicted_class, as.character(y))
})

test_that("SCI-02 stacking rejects mixed, duplicate, and mismatched IDs", {
  ids <- paste0("s", seq_len(4L))
  predictions <- matrix(
    seq_len(8L) / 10,
    nrow = 4L,
    dimnames = list(ids, c("a", "b"))
  )
  y <- stats::setNames(c(0, 1, 0, 1), ids)

  expect_error(
    stacked_multi_omic(predictions, unname(y), n_iter = 2L),
    "present on both"
  )
  expect_error(
    stacked_multi_omic(unname(predictions), y, n_iter = 2L),
    "present on both"
  )

  duplicate_predictions <- predictions
  rownames(duplicate_predictions)[[4L]] <- ids[[3L]]
  expect_error(
    stacked_multi_omic(duplicate_predictions, y, n_iter = 2L),
    "unique sample IDs"
  )
  duplicate_y <- y
  names(duplicate_y)[[4L]] <- ids[[3L]]
  expect_error(
    stacked_multi_omic(predictions, duplicate_y, n_iter = 2L),
    "unique sample IDs"
  )

  mismatched_y <- y
  names(mismatched_y)[[4L]] <- "foreign"
  expect_error(
    stacked_multi_omic(predictions, mismatched_y, n_iter = 2L),
    "same sample IDs"
  )
})

test_that("SCI-03 stratified folds share loads and remain learner-feasible", {
  ids <- paste0("s", seq_len(6L))
  strata <- stats::setNames(rep(c("A", "B"), each = 3L), ids)

  first <- stablr:::.make_multiomic_cv_folds(
    sample_ids = ids,
    groups = NULL,
    v = 5L,
    random_state = 804L,
    strata = strata
  )
  second <- stablr:::.make_multiomic_cv_folds(
    sample_ids = ids,
    groups = NULL,
    v = 5L,
    random_state = 804L,
    strata = strata
  )

  expect_identical(first, second)
  expect_identical(
    sort(vapply(first, function(fold) length(fold$valid_ids), integer(1L))),
    c(1L, 1L, 1L, 1L, 2L)
  )
  expect_true(all(vapply(first, function(fold) {
    all(table(factor(strata[fold$train_ids], levels = c("A", "B"))) >= 2L)
  }, logical(1L))))
  expect_setequal(unlist(lapply(first, `[[`, "valid_ids")), ids)
})

test_that("SCI-03 infeasible classification folds raise a typed condition", {
  ids <- paste0("s", seq_len(6L))
  strata <- stats::setNames(c("rare", "rare", rep("common", 4L)), ids)

  condition <- tryCatch(
    stablr:::.make_multiomic_cv_folds(
      sample_ids = ids,
      groups = NULL,
      v = 3L,
      random_state = 805L,
      strata = strata
    ),
    error = identity
  )

  expect_s3_class(condition, "stablr_cv_infeasible")
  expect_s3_class(condition, "stablr_error")
  expect_identical(condition$minimum_per_class, 2L)
  expect_match(conditionMessage(condition), "modeled class", fixed = TRUE)
})

test_that("SCI-03 fold counts exceeding grouped units fail closed", {
  ids <- paste0("s", seq_len(6L))
  groups <- stats::setNames(rep(paste0("g", 1:2), each = 3L), ids)

  condition <- tryCatch(
    stablr:::.make_multiomic_cv_folds(
      sample_ids = ids,
      groups = groups,
      v = 3L,
      random_state = 806L
    ),
    error = identity
  )

  expect_s3_class(condition, "stablr_cv_infeasible")
  expect_identical(condition$requested_folds, 3L)
  expect_identical(condition$available_units, 2L)
})
