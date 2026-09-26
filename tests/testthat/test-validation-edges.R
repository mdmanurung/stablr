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
  expect_false("f3" %in% get_feature_names_out(fit))
})

# ---- G9: permuted validation omics are repaired by name, not trusted --------

test_that("G9 permuted validation omics and outcome cannot misalign the fit", {
  set.seed(2201L)
  n_train <- 18L
  n_valid <- 8L
  train_ids <- paste0("t", seq_len(n_train))
  valid_ids <- paste0("v", seq_len(n_valid))

  x_a <- rbind(
    matrix(rnorm(n_train * 3L), n_train, 3L, dimnames = list(train_ids, paste0("a", 1:3))),
    matrix(rnorm(n_valid * 3L), n_valid, 3L, dimnames = list(valid_ids, paste0("a", 1:3)))
  )
  x_b <- rbind(
    matrix(rnorm(n_train * 2L), n_train, 2L, dimnames = list(train_ids, paste0("b", 1:2))),
    matrix(rnorm(n_valid * 2L), n_valid, 2L, dimnames = list(valid_ids, paste0("b", 1:2)))
  )
  y_all <- setNames(0.9 * x_a[, 1L] - 0.5 * x_b[, 2L] + rnorm(nrow(x_a), sd = 0.3),
                    rownames(x_a))
  grid <- data.frame(lambda = c(0.2, 0.1))

  aligned <- stabl_multiomic_train_validate(
    x_train_list = list(omic_a = x_a[train_ids, ], omic_b = x_b[train_ids, ]),
    y_train = y_all[train_ids],
    lambda_grid = grid,
    x_valid_list = list(omic_a = x_a[valid_ids, ], omic_b = x_b[valid_ids, ]),
    y_valid = y_all[valid_ids],
    artificial_type = NULL, hard_threshold = 0.3, n_bootstraps = 2L,
    early_fusion = TRUE, random_state = 11L
  )

  order <- c(5L, 2L, 8L, 1L, 7L, 3L, 6L, 4L)
  permuted <- stabl_multiomic_train_validate(
    x_train_list = list(omic_a = x_a[train_ids, ], omic_b = x_b[train_ids, ]),
    y_train = y_all[train_ids],
    lambda_grid = grid,
    x_valid_list = list(omic_a = x_a[valid_ids, ][order, ], omic_b = x_b[valid_ids, ][order, ]),
    y_valid = y_all[valid_ids][order],
    artificial_type = NULL, hard_threshold = 0.3, n_bootstraps = 2L,
    early_fusion = TRUE, random_state = 11L
  )

  # `.validate_multiomic_validation_inputs()` is order-tolerant, and its sole
  # caller normalises every validation table and the validation outcome to the
  # FIRST validation omic's row order.  So the returned validation matrices follow
  # the arrival order rather than a canonical one; what matters is that each row
  # is still matched to its own sample, which is what re-indexing by ID shows.
  by_id <- function(m) m[valid_ids, , drop = FALSE]

  expect_equal(permuted$early_fusion$selected_features, aligned$early_fusion$selected_features)
  expect_setequal(rownames(permuted$early_fusion$selected_valid), valid_ids)
  expect_equal(
    by_id(permuted$early_fusion$selected_valid),
    by_id(aligned$early_fusion$selected_valid)
  )
  expect_setequal(rownames(permuted$selected_valid$omic_a), valid_ids)
  expect_equal(
    by_id(permuted$selected_valid$omic_a),
    by_id(aligned$selected_valid$omic_a)
  )
  expect_equal(permuted$selected_features, aligned$selected_features)
})

test_that("G9 training omics in a different order are rejected before any fitting", {
  set.seed(2202L)
  n <- 12L
  ids <- paste0("s", seq_len(n))
  x_a <- matrix(rnorm(n * 3L), n, 3L, dimnames = list(ids, paste0("a", 1:3)))
  x_b <- matrix(rnorm(n * 2L), n, 2L, dimnames = list(rev(ids), paste0("b", 1:2)))
  y <- setNames(rnorm(n), ids)

  expect_error(
    stabl_multiomic_train_validate(
      x_train_list = list(omic_a = x_a, omic_b = x_b),
      y_train = y,
      lambda_grid = data.frame(lambda = 0.1),
      artificial_type = NULL, hard_threshold = 0.3, n_bootstraps = 2L,
      random_state = 1L
    ),
    "identical sample order"
  )
})

# ---- G2: the training side now matches the validation side -------------------

test_that("G2 non-finite training omics are rejected like non-finite validation omics", {
  set.seed(2203L)
  n <- 14L
  ids <- paste0("s", seq_len(n))
  x_a <- matrix(rnorm(n * 3L), n, 3L, dimnames = list(ids, paste0("a", 1:3)))
  x_b <- matrix(rnorm(n * 2L), n, 2L, dimnames = list(ids, paste0("b", 1:2)))
  y <- setNames(rnorm(n), ids)

  x_b_bad <- x_b
  x_b_bad[4L, 1L] <- Inf

  expect_error(
    stabl_multiomic_train_validate(
      x_train_list = list(omic_a = x_a, omic_b = x_b_bad),
      y_train = y,
      lambda_grid = data.frame(lambda = 0.1),
      artificial_type = NULL, hard_threshold = 0.3, n_bootstraps = 2L,
      random_state = 1L
    ),
    "finite"
  )
})
