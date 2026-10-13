# Audit M-7 / M-8 / V-8: edge-case input validation and degenerate fits.

test_that("stabl_fit errors when sample_fraction > 1 and replace = FALSE", {
  withr::local_seed(0)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  expect_error(
    stabl_fit(x, y,
              lambda_grid = data.frame(lambda = 0.1),
              sample_fraction = 1.5, replace = FALSE,
              n_bootstraps = 3L,
              artificial_type = "random_permutation",
              random_state = 1L),
    "`replace = FALSE`"
  )
})

test_that("stabl_fit handles n_subsamples == n_samples + replace = FALSE", {
  withr::local_seed(0)
  n <- 20L; p <- 4L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  # sample_fraction = 1 → n_subsamples = n_samples; with replace = FALSE every
  # bootstrap is the full dataset (degenerate but legal).  Should not error.
  expect_no_error({
    fit <- stabl_fit(x, y,
                     lambda_grid = data.frame(lambda = c(0.2, 0.1)),
                     sample_fraction = 1, replace = FALSE,
                     n_bootstraps = 3L,
                     artificial_type = "random_permutation",
                     random_state = 1L)
  })
  expect_s3_class(fit, "stabl_fit")
})

test_that("stabl_fit accepts a single-lambda grid", {
  withr::local_seed(0)
  n <- 30L; p <- 4L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))

  fit <- stabl_fit(x, y,
                   lambda_grid = data.frame(lambda = 0.1),
                   n_bootstraps = 5L,
                   artificial_type = "random_permutation",
                   random_state = 1L)

  expect_s3_class(fit, "stabl_fit")
  expect_equal(ncol(fit$stabl_scores_), 1L)
})

test_that("stabl_fit tolerates a zero-variance column", {
  withr::local_seed(0)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  x[, 3L] <- 0  # zero-variance feature
  y <- setNames(rnorm(n), rownames(x))

  expect_no_error({
    fit <- stabl_fit(x, y,
                     lambda_grid = data.frame(lambda = c(0.2, 0.1)),
                     n_bootstraps = 5L,
                     artificial_type = "random_permutation",
                     random_state = 1L)
  })
  # The zero-variance column should not be selected (no variability to score).
  expect_false("f3" %in% get_feature_names(fit))
})
