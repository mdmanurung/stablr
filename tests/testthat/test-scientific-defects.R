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
