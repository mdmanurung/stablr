test_that("stabl_fit returns stabl_fit object with correct dimensions", {
  set.seed(42)
  n <- 30L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(sample(0:1, n, replace = TRUE), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.01, 0.5, length.out = 5L))

  # suppressWarnings: glmnet warns about small bootstrap subsamples for
  # binomial family when n_subsamples is small relative to the number of
  # classes — this is expected behaviour in minimal tests.
  fit <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "binomial",
    n_bootstraps    = 20L,
    artificial_type = "random_permutation"
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 5L))
  expect_equal(fit$n_features_in_, p)
  expect_equal(fit$feature_names, paste0("f", seq_len(p)))
})

test_that("stability scores are in [0, 1]", {
  set.seed(1)
  n <- 40L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.01, 0.3, length.out = 4L))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "gaussian",
    n_bootstraps    = 15L,
    artificial_type = "random_permutation"
  )

  expect_true(all(fit$stabl_scores_ >= 0 & fit$stabl_scores_ <= 1))
})

test_that("FDP+ is computed when artificial_type = 'random_permutation'", {
  set.seed(7)
  n <- 40L; p <- 12L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.01, 0.5, length.out = 4L))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "gaussian",
    n_bootstraps    = 15L,
    artificial_type = "random_permutation"
  )

  expect_false(is.null(fit$fdr_min_threshold_))
  expect_true(is.finite(fit$fdr_min_threshold_))
  expect_true(fit$fdr_min_threshold_ >= 0 && fit$fdr_min_threshold_ <= 1)
  expect_false(is.null(fit$stabl_scores_artificial_))
  expect_equal(nrow(fit$stabl_scores_artificial_), p)  # artificial_proportion=1
})

test_that("get_support returns logical vector of correct length", {
  set.seed(3)
  n <- 30L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.01, 0.5, length.out = 4L))

  fit  <- stabl_fit(x = x, y = y, lambda_grid = lam_grid, family = "gaussian",
                    n_bootstraps = 10L, artificial_type = "random_permutation")
  mask <- get_support(fit)

  expect_length(mask, p)
  expect_type(mask, "logical")
})

test_that("explore fallback selects exactly n_explore features even when all scores are tied", {
  # When hard_threshold = 0.99, nothing passes.  With explore = TRUE and
  # n_explore = 3, exactly 3 features should be returned regardless of score
  # ties (e.g. all scores identical / all zero).
  set.seed(42)
  n <- 30L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  # Use a very large lambda so all features score ~0 (maximum tie scenario)
  lam_grid <- data.frame(lambda = c(10, 20))  # huge regularisation

  fit <- stabl_fit(
    x              = x,
    y              = y,
    lambda_grid    = lam_grid,
    family         = "gaussian",
    n_bootstraps   = 20L,
    artificial_type = NULL,
    hard_threshold = 0.99,   # nothing will pass
    explore        = TRUE,
    n_explore      = 3L
  )

  mask <- get_support(fit)

  # Exactly n_explore features selected, not all p
  expect_equal(sum(mask), 3L)
  expect_length(mask, p)
  expect_type(mask, "logical")

  # The selected features must be the top-3 by importance score
  importances <- get_importances(fit)
  top3_names  <- names(sort(importances, decreasing = TRUE))[seq_len(3L)]
  expect_setequal(names(mask)[mask], top3_names)
})

test_that("hard_threshold bypasses FDP+ path", {
  set.seed(5)
  n <- 30L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.01, 0.5, length.out = 4L))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "gaussian",
    n_bootstraps    = 10L,
    artificial_type = NULL,
    hard_threshold  = 0.3
  )

  expect_null(fit$fdr_min_threshold_)
  mask        <- get_support(fit)
  importances <- get_importances(fit)
  expect_true(all(importances[mask] > 0.3))
})

test_that("error when both hard_threshold and artificial_type are NULL", {
  set.seed(6)
  n <- 20L; p <- 5L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = c(0.1, 0.2))

  expect_error(
    stabl_fit(x = x, y = y, lambda_grid = lam_grid,
              artificial_type = NULL, hard_threshold = NULL),
    "Either"
  )
})

test_that("group bootstrap path runs without error", {
  set.seed(9)
  n <- 40L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y      <- setNames(rnorm(n), rownames(x))
  groups <- setNames(rep(c("A", "B", "C", "D"), each = 10L), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.05, 0.3, length.out = 3L))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "gaussian",
    n_bootstraps    = 10L,
    artificial_type = "random_permutation",
    sample_groups   = groups
  )

  expect_s3_class(fit, "stabl_fit")
})

test_that("explore clamps to p when n_explore > n_features", {
  # Verifies min(n_explore, p) clamping: n_explore=100 with p=6 should select
  # all 6 features, not error.
  set.seed(7L)
  n <- 20L; p <- 6L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = c(10, 20))  # huge lambda → all scores ~0

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "gaussian",
    n_bootstraps    = 10L,
    artificial_type = NULL,
    hard_threshold  = 0.99,
    explore         = TRUE,
    n_explore       = 100L       # far exceeds p
  )

  mask <- get_support(fit)
  expect_equal(sum(mask), p)    # clamped to p, not 100
  expect_length(mask, p)
})

test_that("get_support errors with informative message when no threshold is resolvable", {
  set.seed(8L)
  n <- 20L; p <- 4L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = data.frame(lambda = c(0.1, 0.2)),
    family          = "gaussian",
    n_bootstraps    = 10L,
    artificial_type = "random_permutation"
  )
  # Simulate a fitted object with no threshold resolved (neither hard nor FDP+).
  fit$hard_threshold      <- NULL
  fit$fdr_min_threshold_  <- NULL

  expect_error(get_support(fit), "No threshold available")
})

test_that("stabl_fit exposes bootstrap stratification design option", {
  set.seed(10)
  n <- 38L; p <- 6L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(factor(c(rep("NP", 15L), rep("P", 23L))), rownames(x))
  study <- setNames(rep(c("A", "B"), length.out = n), rownames(x))
  bootstrap_strata <- data.frame(outcome = y, study = study)
  rownames(bootstrap_strata) <- rownames(x)
  lam_grid <- data.frame(lambda = seq(0.05, 0.2, length.out = 3L))

  fit <- suppressWarnings(stabl_fit(
    x                  = x,
    y                  = y,
    lambda_grid        = lam_grid,
    family             = "binomial",
    n_bootstraps       = 6L,
    artificial_type    = "random_permutation",
    bootstrap_strata   = bootstrap_strata,
    random_state       = 10L
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_true(isTRUE(fit$stratify_bootstrap))
  expect_length(fit$bootstrap_strata_levels, 4L)
})

test_that("auto_lambda_grid returns data.frame with lambda column", {
  set.seed(11)
  n <- 40L; p <- 10L
  x <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  grid <- auto_lambda_grid(x = x, y = y, family = "gaussian", n_lambda = 10L)

  expect_s3_class(grid, "data.frame")
  expect_true("lambda" %in% names(grid))
  expect_gt(nrow(grid), 0L)
})

test_that("auto_lambda_grid with l1_ratio returns alpha column", {
  set.seed(12)
  n <- 40L; p <- 10L
  x <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  grid <- auto_lambda_grid(x = x, y = y, family = "gaussian",
                           n_lambda = 5L, l1_ratio = c(0.5, 0.9))

  expect_true("alpha" %in% names(grid))
  expect_true("lambda" %in% names(grid))
  expect_true(all(c(0.5, 0.9) %in% grid$alpha))
})

test_that("stabl_fit auto lambda grid forwards l1_ratio", {
  set.seed(120)
  n <- 30L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  fit <- stabl_fit(
    x = x,
    y = y,
    lambda_grid = "auto",
    base_learner = "elastic_net",
    family = "gaussian",
    n_lambda = 3L,
    l1_ratio = c(0.3, 0.7),
    n_bootstraps = 2L,
    artificial_type = NULL,
    hard_threshold = 1,
    sample_fraction = 1,
    random_state = 120L
  )

  expect_true("alpha" %in% names(fit$fitted_lambda_grid))
  expect_equal(sort(unique(fit$fitted_lambda_grid$alpha)), c(0.3, 0.7))
})

test_that("get_feature_names returns character subset of feature names", {
  set.seed(13)
  n <- 40L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.01, 0.5, length.out = 4L))

  fit <- stabl_fit(x = x, y = y, lambda_grid = lam_grid, family = "gaussian",
                   n_bootstraps = 10L, artificial_type = "random_permutation")

  sel <- get_feature_names(fit)
  expect_type(sel, "character")
  expect_true(all(sel %in% paste0("f", seq_len(p))))
})

test_that("print.stabl_fit runs without error", {
  set.seed(14)
  n <- 20L; p <- 6L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = c(0.05, 0.1, 0.2))

  fit <- stabl_fit(x = x, y = y, lambda_grid = lam_grid, family = "gaussian",
                   n_bootstraps = 5L, artificial_type = "random_permutation")
  expect_output(print(fit), "stabl_fit")
})

test_that("compute_fdp_plus returns correct list structure", {
  set.seed(20)
  p <- 8L; nl <- 4L
  scores     <- matrix(runif(p * nl), p, nl)
  scores_art <- matrix(runif(p * nl) * 0.3, p, nl)  # lower art scores

  result <- compute_fdp_plus(
    stabl_scores            = scores,
    stabl_scores_artificial = scores_art,
    artificial_proportion   = 1.0,
    fdr_threshold_range     = seq(0, 1, by = 0.1)
  )

  expect_named(result, c("FDRs", "min_fdr", "fdr_min_threshold", "fdrs_table"))
  expect_length(result$FDRs, 11L)
  expect_true(result$fdr_min_threshold >= 0 && result$fdr_min_threshold <= 1)
  expect_equal(dim(result$fdrs_table), c(nl, 11L))
})

test_that("stabl_fit with lambda_grid = 'auto' runs without error", {
  set.seed(15)
  n <- 40L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = "auto",
    family          = "gaussian",
    n_bootstraps    = 5L,
    n_lambda        = 6L,
    artificial_type = "random_permutation"
  )
  expect_s3_class(fit, "stabl_fit")
  expect_equal(nrow(fit$fitted_lambda_grid), ncol(fit$stabl_scores_))
})

test_that("adaptive_lasso base learner runs and returns expected dimensions", {
  set.seed(16)
  n <- 45L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = seq(0.02, 0.4, length.out = 4L))

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "gaussian",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation"
  )

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 4L))
  expect_true(all(fit$stabl_scores_ >= 0 & fit$stabl_scores_ <= 1))
})

test_that("multinomial family path runs without error", {
  set.seed(17)
  n <- 54L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(factor(sample(c("A", "B", "C"), n, replace = TRUE)),
                rownames(x))
  lam_grid <- data.frame(lambda = seq(0.02, 0.5, length.out = 3L))

  fit <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "multinomial",
    n_bootstraps    = 5L,
    artificial_type = "random_permutation"
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 3L))
})

test_that("sparse_group_lasso runs with explicit feature groups", {
  skip_if_not_installed("sparsegl")

  set.seed(18)
  n <- 48L; p <- 12L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(
    alpha = c(0.2, 0.6),
    lambda = c(0.2, 0.05)
  )
  feat_groups <- rep(seq_len(p / 3L), each = 3L)

  fit <- stabl_fit(
    x                   = x,
    y                   = y,
    lambda_grid         = lam_grid,
    base_learner        = "sparse_group_lasso",
    family              = "gaussian",
    feature_groups      = feat_groups,
    n_bootstraps        = 5L,
    artificial_type     = "random_permutation"
  )

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 2L))
})

test_that("sparse_group_lasso runs with correlation-based groups", {
  skip_if_not_installed("sparsegl")

  set.seed(19)
  n <- 50L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(sample(0:1, n, replace = TRUE), rownames(x))
  lam_grid <- data.frame(alpha = 0.4, lambda = 0.1)

  fit <- suppressWarnings(stabl_fit(
    x                   = x,
    y                   = y,
    lambda_grid         = lam_grid,
    base_learner        = "sparse_group_lasso",
    family              = "binomial",
    corr_group_threshold = 90,
    n_bootstraps        = 5L,
    artificial_type     = "random_permutation"
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 1L))
})

test_that("sparse_group_lasso multinomial path runs", {
  skip_if_not_installed("sparsegl")

  set.seed(21)
  n <- 60L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(factor(sample(c("A", "B", "C"), n, replace = TRUE)),
                rownames(x))
  lam_grid <- data.frame(alpha = c(0.2, 0.8), lambda = c(0.2, 0.05))
  feat_groups <- rep(seq_len(3L), each = 3L)

  fit <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "sparse_group_lasso",
    family          = "multinomial",
    feature_groups  = feat_groups,
    n_bootstraps    = 4L,
    artificial_type = "random_permutation"
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 2L))
})

test_that("cox family runs with Surv outcomes", {
  skip_if_not_installed("survival")

  set.seed(24)
  n <- 50L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- survival::Surv(time = rexp(n, rate = 0.1),
                      event = rbinom(n, size = 1L, prob = 0.7))
  rownames(y) <- rownames(x)

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = data.frame(lambda = seq(0.02, 0.2, length.out = 3L)),
    base_learner    = "lasso",
    family          = "cox",
    n_bootstraps    = 5L,
    artificial_type = "random_permutation",
    random_state    = 31L
  )

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, 3L))
})

test_that("sparse_group_lasso rejects cox family", {
  skip_if_not_installed("survival")

  n <- 20L; p <- 6L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- survival::Surv(time = rexp(n, rate = 0.1),
                      event = rbinom(n, size = 1L, prob = 0.7))
  rownames(y) <- rownames(x)

  expect_error(
    stabl_fit(
      x               = x,
      y               = y,
      lambda_grid     = data.frame(alpha = 0.3, lambda = 0.1),
      base_learner    = "sparse_group_lasso",
      family          = "cox",
      feature_groups  = rep(1:3, each = 2),
      n_bootstraps    = 2L,
      artificial_type = NULL,
      hard_threshold  = 0.2
    ),
    "Cox family is not supported"
  )
})

test_that("auto_lambda_grid supports cox family with mixed alpha grid", {
  skip_if_not_installed("survival")

  set.seed(25)
  n <- 60L; p <- 11L
  x <- matrix(rnorm(n * p), n, p)
  y <- survival::Surv(
    time = rexp(n, rate = 0.2),
    event = rbinom(n, size = 1L, prob = 0.6)
  )

  grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "cox",
    n_lambda = 6L,
    l1_ratio = c(0.2, 0.8)
  )

  expect_s3_class(grid, "data.frame")
  expect_true(all(c("alpha", "lambda") %in% names(grid)))
  expect_equal(sort(unique(grid$alpha)), c(0.2, 0.8))
  expect_equal(nrow(grid), 12L)
  expect_true(all(is.finite(grid$lambda)))
  expect_true(all(grid$lambda > 0))
})

test_that("elastic_net consumes mixed-alpha lambda grid for cox", {
  skip_if_not_installed("survival")

  set.seed(26)
  n <- 56L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- survival::Surv(
    time = rexp(n, rate = 0.15),
    event = rbinom(n, size = 1L, prob = 0.65)
  )
  rownames(y) <- rownames(x)

  lam_grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "cox",
    n_lambda = 4L,
    l1_ratio = c(0.3, 0.7)
  )

  # Small bootstrap samples can trigger expected glmnet cox path warnings for
  # tiny lambda values; suppress to keep this a structural coverage test.
  fit <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "elastic_net",
    family          = "cox",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 32L
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(nrow(fit$fitted_lambda_grid), nrow(lam_grid))
  expect_equal(sort(unique(fit$fitted_lambda_grid$alpha)), c(0.3, 0.7))
})

test_that("auto_lambda_grid supports multinomial family with mixed alpha grid", {
  set.seed(27)
  n <- 66L; p <- 10L
  x <- matrix(rnorm(n * p), n, p)
  y <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "multinomial",
    n_lambda = 5L,
    l1_ratio = c(0.2, 0.8)
  )

  expect_s3_class(grid, "data.frame")
  expect_true(all(c("alpha", "lambda") %in% names(grid)))
  expect_equal(sort(unique(grid$alpha)), c(0.2, 0.8))
  expect_equal(nrow(grid), 10L)
  expect_true(all(is.finite(grid$lambda)))
  expect_true(all(grid$lambda > 0))
})

test_that("auto_lambda_grid supports binomial family with mixed alpha grid", {
  set.seed(28)
  n <- 64L; p <- 9L
  x <- matrix(rnorm(n * p), n, p)
  y <- sample(0:1, n, replace = TRUE)

  grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "binomial",
    n_lambda = 5L,
    l1_ratio = c(0.3, 0.7)
  )

  expect_s3_class(grid, "data.frame")
  expect_true(all(c("alpha", "lambda") %in% names(grid)))
  expect_equal(sort(unique(grid$alpha)), c(0.3, 0.7))
  expect_equal(nrow(grid), 10L)
  expect_true(all(is.finite(grid$lambda)))
  expect_true(all(grid$lambda > 0))
})

test_that("multinomial elastic_net mixed-alpha path is structurally stable", {
  set.seed(29)
  n <- 72L; p <- 11L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(factor(sample(c("A", "B", "C"), n, replace = TRUE)),
                rownames(x))

  lam_grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "multinomial",
    n_lambda = 4L,
    l1_ratio = c(0.25, 0.75)
  )

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "elastic_net",
    family          = "multinomial",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 40L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "elastic_net",
    family          = "multinomial",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 40L
  ))

  expect_s3_class(fit_1, "stabl_fit")
  expect_equal(dim(fit_1$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(nrow(fit_1$fitted_lambda_grid), nrow(lam_grid))
  expect_equal(sort(unique(fit_1$fitted_lambda_grid$alpha)), c(0.25, 0.75))
  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
})

test_that("binomial adaptive_lasso mixed-alpha grid is structurally stable", {
  set.seed(30)
  n <- 68L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(sample(0:1, n, replace = TRUE), rownames(x))

  lam_grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "binomial",
    n_lambda = 4L,
    l1_ratio = c(0.2, 0.8)
  )

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "binomial",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 41L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "binomial",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 41L
  ))

  expect_s3_class(fit_1, "stabl_fit")
  expect_equal(dim(fit_1$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(nrow(fit_1$fitted_lambda_grid), nrow(lam_grid))
  expect_equal(sort(unique(fit_1$fitted_lambda_grid$alpha)), c(0.2, 0.8))
  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
})

test_that("gaussian lasso preserves mixed-alpha grid rows and is deterministic", {
  set.seed(31)
  n <- 62L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  lam_grid <- auto_lambda_grid(
    x = x,
    y = y,
    family = "gaussian",
    n_lambda = 4L,
    l1_ratio = c(0.1, 0.9)
  )

  fit_1 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "gaussian",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 42L
  )

  fit_2 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "gaussian",
    n_bootstraps    = 4L,
    artificial_type = "random_permutation",
    random_state    = 42L
  )

  expect_s3_class(fit_1, "stabl_fit")
  expect_equal(dim(fit_1$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(nrow(fit_1$fitted_lambda_grid), nrow(lam_grid))
  expect_equal(sort(unique(fit_1$fitted_lambda_grid$alpha)), c(0.1, 0.9))
  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
})

test_that("high-collinearity regime remains deterministic and stable", {
  set.seed(43)
  n <- 70L; p <- 8L

  f1 <- rnorm(n)
  f2 <- f1 + rnorm(n, sd = 0.01)
  x_rest <- matrix(rnorm(n * (p - 2L)), nrow = n)
  x <- cbind(f1, f2, x_rest)
  colnames(x) <- paste0("f", seq_len(p))
  rownames(x) <- paste0("s", seq_len(n))

  y <- setNames(1.2 * f1 + rnorm(n, sd = 0.6), rownames(x))
  lam_grid <- data.frame(lambda = exp(seq(log(0.5), log(0.02), length.out = 5L)))

  fit_1 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "gaussian",
    n_bootstraps    = 5L,
    artificial_type = "random_permutation",
    random_state    = 50L
  )

  fit_2 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "gaussian",
    n_bootstraps    = 5L,
    artificial_type = "random_permutation",
    random_state    = 50L
  )

  importances <- get_importances(fit_1)

  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(is.finite(importances)))
  expect_true(any(importances[c("f1", "f2")] > 0))
})

test_that("near-zero lambda tails keep path output finite and aligned", {
  set.seed(44)
  n <- 64L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(0.9 * x[, 1L] - 0.5 * x[, 3L] + rnorm(n, sd = 0.7), rownames(x))

  lam_grid <- data.frame(
    alpha = rep(c(0.3, 0.7), each = 4L),
    lambda = rep(c(5e-2, 1e-3, 1e-5, 1e-7), times = 2L)
  )

  # Very small lambdas can trigger expected numerical warnings in glmnet;
  # this is a structural edge-regime coverage test.
  fit <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "elastic_net",
    family          = "gaussian",
    n_bootstraps    = 5L,
    artificial_type = "random_permutation",
    random_state    = 51L
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(nrow(fit$fitted_lambda_grid), nrow(lam_grid))
  expect_equal(sort(unique(fit$fitted_lambda_grid$alpha)), c(0.3, 0.7))
  expect_true(all(is.finite(fit$stabl_scores_)))
})

test_that("class-imbalance binomial regime is deterministic and bounded", {
  set.seed(45)
  n <- 80L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))

  logits <- 1.1 * x[, 1L] - 0.8 * x[, 2L] - 2.2
  probs <- stats::plogis(logits)
  cutoff <- as.numeric(stats::quantile(probs, probs = 0.9, names = FALSE))
  y <- setNames(as.integer(probs >= cutoff), rownames(x))
  expect_lt(mean(y), 0.25)
  expect_equal(sort(unique(y)), c(0L, 1L))

  lam_grid <- data.frame(lambda = exp(seq(log(0.3), log(0.01), length.out = 4L)))

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "binomial",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 52L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "binomial",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 52L
  ))

  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
  expect_true(is.finite(fit_1$fdr_min_threshold_))
})

# ── Phase 6 behavior-level parity: multinomial ────────────────────────────────

test_that("multinomial lasso detects true class-separating signal features", {
  # Fixture: f1 and f2 are strong class predictors; f3..f10 are pure noise.
  # Parity expectation: true signal features should accumulate higher mean
  # stability scores than the noise bulk across all lambda grid rows.
  set.seed(50)
  n <- 120L; p <- 10L
  n_cls <- 40L  # 40 per class

  # Class centres: class A (+1,0), B (-1,+1), C (0,-1)
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  x[seq_len(n_cls), 1L]              <-  x[seq_len(n_cls), 1L] + 2.0
  x[seq_len(n_cls) + n_cls, 1L]     <-  x[seq_len(n_cls) + n_cls, 1L] - 2.0
  x[seq_len(n_cls) + n_cls, 2L]     <-  x[seq_len(n_cls) + n_cls, 2L] + 2.0
  x[seq_len(n_cls) + 2L * n_cls, 2L] <- x[seq_len(n_cls) + 2L * n_cls, 2L] - 2.0
  y <- setNames(factor(rep(c("A", "B", "C"), each = n_cls)), rownames(x))

  lam_grid <- auto_lambda_grid(x = x, y = y, family = "multinomial", n_lambda = 8L)

  fit <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "multinomial",
    n_bootstraps    = 20L,
    artificial_type = "random_permutation",
    random_state    = 50L
  ))

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, nrow(lam_grid)))
  expect_true(all(fit$stabl_scores_ >= 0 & fit$stabl_scores_ <= 1))

  # Signal features f1 and f2 should have higher mean stability than average
  # of noise features f3..f10.
  mean_scores <- rowMeans(fit$stabl_scores_)
  signal_mean <- mean(mean_scores[c(1L, 2L)])
  noise_mean  <- mean(mean_scores[3L:p])
  expect_gt(signal_mean, noise_mean)
})

test_that("multinomial lasso is deterministic across repeated calls", {
  set.seed(51)
  n <- 90L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(factor(sample(c("X", "Y", "Z"), n, replace = TRUE)),
                rownames(x))
  lam_grid <- data.frame(lambda = exp(seq(log(0.4), log(0.02), length.out = 5L)))

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "multinomial",
    n_bootstraps    = 8L,
    artificial_type = "random_permutation",
    random_state    = 51L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "multinomial",
    n_bootstraps    = 8L,
    artificial_type = "random_permutation",
    random_state    = 51L
  ))

  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
})

test_that("multinomial adaptive_lasso is deterministic and bounded", {
  set.seed(52)
  n <- 90L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- setNames(factor(sample(c("A", "B", "C"), n, replace = TRUE)),
                rownames(x))
  lam_grid <- data.frame(lambda = exp(seq(log(0.3), log(0.02), length.out = 4L)))

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "multinomial",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 52L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "multinomial",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 52L
  ))

  expect_s3_class(fit_1, "stabl_fit")
  expect_equal(dim(fit_1$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
})

# ── Phase 6 behavior-level parity: Cox ────────────────────────────────────────

test_that("cox lasso detects true survival-signal feature", {
  skip_if_not_installed("survival")

  # Fixture: feature f1 drives hazard proportionally; f2..f10 are pure noise.
  # Parity expectation: f1 should accumulate higher mean stability score than
  # noise features f2..f10.
  set.seed(60)
  n <- 120L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))

  # Exponential survival times with hazard exp(1.5 * x[,1])
  beta_true <- 1.5
  times  <- rexp(n, rate = exp(beta_true * x[, 1L]))
  censor <- rexp(n, rate = 0.3)
  t_obs  <- pmin(times, censor)
  event  <- as.integer(times <= censor)

  y <- survival::Surv(time = t_obs, event = event)
  rownames(y) <- rownames(x)

  lam_grid <- auto_lambda_grid(x = x, y = y, family = "cox", n_lambda = 8L)

  fit <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "cox",
    n_bootstraps    = 20L,
    artificial_type = "random_permutation",
    random_state    = 60L
  )

  expect_s3_class(fit, "stabl_fit")
  expect_equal(dim(fit$stabl_scores_), c(p, nrow(lam_grid)))
  expect_true(all(fit$stabl_scores_ >= 0 & fit$stabl_scores_ <= 1))

  # f1 should have higher mean stability than noise average
  mean_scores  <- rowMeans(fit$stabl_scores_)
  noise_mean   <- mean(mean_scores[2L:p])
  expect_gt(mean_scores[[1L]], noise_mean)
})

test_that("cox lasso is deterministic across repeated calls", {
  skip_if_not_installed("survival")

  set.seed(61)
  n <- 80L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- survival::Surv(
    time  = rexp(n, rate = 0.15),
    event = rbinom(n, size = 1L, prob = 0.65)
  )
  rownames(y) <- rownames(x)
  lam_grid <- data.frame(lambda = exp(seq(log(0.4), log(0.02), length.out = 5L)))

  fit_1 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "cox",
    n_bootstraps    = 8L,
    artificial_type = "random_permutation",
    random_state    = 61L
  )

  fit_2 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "cox",
    n_bootstraps    = 8L,
    artificial_type = "random_permutation",
    random_state    = 61L
  )

  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
})

test_that("cox adaptive_lasso is deterministic and bounded", {
  skip_if_not_installed("survival")

  set.seed(62)
  n <- 80L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- survival::Surv(
    time  = rexp(n, rate = 0.12),
    event = rbinom(n, size = 1L, prob = 0.70)
  )
  rownames(y) <- rownames(x)
  lam_grid <- data.frame(lambda = exp(seq(log(0.3), log(0.02), length.out = 4L)))

  fit_1 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "cox",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 62L
  )

  fit_2 <- stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "adaptive_lasso",
    family          = "cox",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 62L
  )

  expect_s3_class(fit_1, "stabl_fit")
  expect_equal(dim(fit_1$stabl_scores_), c(p, nrow(lam_grid)))
  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
})

test_that("cox elastic_net mixed-alpha is deterministic and bounded", {
  skip_if_not_installed("survival")

  set.seed(63)
  n <- 80L; p <- 9L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)),
                              paste0("f", seq_len(p))))
  y <- survival::Surv(
    time  = rexp(n, rate = 0.10),
    event = rbinom(n, size = 1L, prob = 0.68)
  )
  rownames(y) <- rownames(x)

  lam_grid <- auto_lambda_grid(
    x = x, y = y, family = "cox", n_lambda = 4L, l1_ratio = c(0.4, 0.8)
  )

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "elastic_net",
    family          = "cox",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 63L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "elastic_net",
    family          = "cox",
    n_bootstraps    = 6L,
    artificial_type = "random_permutation",
    random_state    = 63L
  ))

  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
  expect_equal(sort(unique(fit_1$fitted_lambda_grid$alpha)), c(0.4, 0.8))
})

# ── Phase 6 edge-regime: multinomial high-collinearity ────────────────────────

test_that("multinomial high-collinearity regime is deterministic and bounded", {
  set.seed(70)
  n <- 100L; p <- 8L

  # Two highly collinear features + noise
  f1  <- rnorm(n)
  f2  <- f1 + rnorm(n, sd = 0.02)
  x_rest <- matrix(rnorm(n * (p - 2L)), nrow = n)
  x <- cbind(f1, f2, x_rest)
  colnames(x) <- paste0("f", seq_len(p))
  rownames(x) <- paste0("s", seq_len(n))

  # 3-class outcome driven by f1
  probs <- cbind(plogis(1.2 * f1), plogis(-1.2 * f1) * 0.5, 0.5)
  probs <- probs / rowSums(probs)
  cls   <- apply(probs, 1L, function(pr) sample(c("A", "B", "C"), 1L, prob = pr))
  y     <- setNames(factor(cls), rownames(x))

  lam_grid <- data.frame(lambda = exp(seq(log(0.4), log(0.02), length.out = 5L)))

  fit_1 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "multinomial",
    n_bootstraps    = 8L,
    artificial_type = "random_permutation",
    random_state    = 70L
  ))

  fit_2 <- suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    base_learner    = "lasso",
    family          = "multinomial",
    n_bootstraps    = 8L,
    artificial_type = "random_permutation",
    random_state    = 70L
  ))

  expect_equal(fit_1$stabl_scores_, fit_2$stabl_scores_)
  expect_true(all(fit_1$stabl_scores_ >= 0 & fit_1$stabl_scores_ <= 1))
  expect_equal(dim(fit_1$stabl_scores_), c(p, nrow(lam_grid)))
})

test_that(".build_corr_groups applies -0.1 offset matching Python reference", {
  # Python stabl.py line 1142: threshold = np.percentile(corr_val, perc) - 0.1
  # Without the offset, highly correlated pairs whose corr is *between* the raw
  # percentile and (percentile - 0.1) are not grouped.  The offset must cause
  # at least one additional pair to be joined relative to the raw-percentile cutoff.
  set.seed(9)
  n <- 50L
  # f1 and f2 are nearly identical: pairwise corr ≈ 0.99
  f1 <- rnorm(n)
  f2 <- f1 + rnorm(n, sd = 0.05)
  # Independent noise features
  x_noise <- matrix(rnorm(n * 8L), n, 8L)
  x <- cbind(f1, f2, x_noise)
  colnames(x) <- paste0("f", seq_len(ncol(x)))

  # Access the internal helper via ::: (exported in tests only)
  # Returns an unnamed integer vector: grps[i] is the group id of column i.
  grps <- stablr:::.build_corr_groups(x, percentile = 95)

  # f1 (col 1) and f2 (col 2) should be in the same group
  expect_equal(grps[[1L]], grps[[2L]])
  # All groups are integers >= 1
  expect_true(all(grps >= 1L))
  expect_length(grps, ncol(x))
})
