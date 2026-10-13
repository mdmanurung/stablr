# Audit V-10: assert that the public accessors round-trip consistently.
# Specifically, get_feature_names(fit) must equal
# names(get_support(fit))[get_support(fit)] for any threshold-resolution path.

test_that("get_feature_names equals names(get_support())[get_support()]", {
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

  mask <- get_support(fit)
  expect_identical(get_feature_names(fit), names(mask)[mask])
})

test_that("accessor round-trip holds with hard_threshold override", {
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

  mask_at_zero <- get_support(fit, new_hard_threshold = 0.001)
  expect_identical(get_feature_names(fit, new_hard_threshold = 0.001),
                   names(mask_at_zero)[mask_at_zero])
})

# E1: get_stabl_scores ---------------------------------------------------------

test_that("E1: get_stabl_scores returns a numeric matrix with correct dimensions", {
  withr::local_seed(42L)
  n <- 30L; p <- 5L; L <- 3L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid    = data.frame(lambda = c(0.3, 0.15, 0.05)),
                   n_bootstraps   = 4L,
                   hard_threshold = 0.5,
                   random_state   = 1L)

  sc <- get_stabl_scores(fit)

  expect_true(is.matrix(sc))
  expect_true(is.numeric(sc))
  expect_equal(nrow(sc), p)
  expect_equal(ncol(sc), L)
})

test_that("E1: get_stabl_scores row names match feature names", {
  withr::local_seed(42L)
  n <- 30L; p <- 4L
  feat_names <- paste0("gene", seq_len(p))
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), feat_names))
  y <- setNames(rnorm(n), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid    = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps   = 4L,
                   hard_threshold = 0.5,
                   random_state   = 1L)

  expect_equal(rownames(get_stabl_scores(fit)), feat_names)
})

test_that("E1: get_stabl_scores errors on an unfitted object", {
  not_fit <- list(stabl_scores_ = NULL)
  class(not_fit) <- "stabl_fit"
  expect_error(get_stabl_scores(not_fit), regexp = "stabl_scores_")
})

test_that("transform_stabl subsets matrix input in fitted feature order and preserves row names", {
  fit <- structure(
    list(
      stabl_scores_ = matrix(c(0.1, 0.9, 0.7), ncol = 1L,
                             dimnames = list(c("b", "a", "c"), NULL)),
      hard_threshold = 0.5,
      fdr_min_threshold_ = NULL,
      explore = FALSE,
      n_explore = 1L,
      feature_names = c("b", "a", "c"),
      n_features_in_ = 3L
    ),
    class = "stabl_fit"
  )
  x <- matrix(
    1:12,
    nrow = 3L,
    dimnames = list(c("s1", "s2", "s3"), c("a", "b", "c", "d"))
  )

  out <- transform_stabl(fit, x)

  expect_true(is.matrix(out))
  expect_identical(rownames(out), rownames(x))
  expect_identical(colnames(out), c("a", "c"))
  expect_equal(out, x[, c("a", "c"), drop = FALSE])
})

test_that("transform_stabl handles data.frame input, zero selections, and threshold override", {
  fit <- structure(
    list(
      stabl_scores_ = matrix(c(0.2, 0.4), ncol = 1L,
                             dimnames = list(c("a", "b"), NULL)),
      hard_threshold = 0.9,
      fdr_min_threshold_ = NULL,
      explore = FALSE,
      n_explore = 1L,
      feature_names = c("a", "b"),
      n_features_in_ = 2L
    ),
    class = "stabl_fit"
  )
  x <- data.frame(b = 1:3, a = 4:6, row.names = paste0("s", 1:3))

  none <- transform_stabl(fit, x)
  expect_s3_class(none, "data.frame")
  expect_equal(dim(none), c(3L, 0L))
  expect_identical(rownames(none), rownames(x))

  selected <- transform_stabl(fit, x, new_hard_threshold = 0.1)
  expect_s3_class(selected, "data.frame")
  expect_identical(names(selected), c("a", "b"))
})

test_that("transform_stabl rejects missing or duplicate input feature names", {
  fit <- structure(
    list(
      stabl_scores_ = matrix(c(0.8, 0.1), ncol = 1L,
                             dimnames = list(c("a", "b"), NULL)),
      hard_threshold = 0.5,
      fdr_min_threshold_ = NULL,
      explore = FALSE,
      n_explore = 1L,
      feature_names = c("a", "b"),
      n_features_in_ = 2L
    ),
    class = "stabl_fit"
  )

  expect_error(transform_stabl(fit, data.frame(b = 1:2)), "missing.*a", ignore.case = TRUE)

  x_dup <- matrix(1:4, nrow = 2L, dimnames = list(c("s1", "s2"), c("a", "a")))
  expect_error(transform_stabl(fit, x_dup), "duplicate", ignore.case = TRUE)
})
