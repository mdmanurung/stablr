test_that("validate_sample_alignment accepts aligned inputs", {
  x <- data.frame(a = c(1, 2), b = c(3, 4), row.names = c("s1", "s2"))
  y <- c(s1 = 0, s2 = 1)
  g <- c(s1 = "p1", s2 = "p2")

  expect_invisible(validate_sample_alignment(x, y, g))
})

test_that("validate_sample_alignment fails on id mismatch", {
  x <- data.frame(a = c(1, 2), row.names = c("s1", "s2"))
  y <- c(s1 = 0, s3 = 1)

  expect_error(validate_sample_alignment(x, y), "Sample mismatch")
})

test_that("validate_sample_alignment rejects duplicate sample ids before matching", {
  x <- matrix(c(1, 2), ncol = 1L, dimnames = list(c("s1", "s1"), "a"))
  y <- c(s1 = 0, s2 = 1)

  expect_error(validate_sample_alignment(x, y), "x.*row names.*duplicate", ignore.case = TRUE)
  expect_error(
    validate_sample_alignment(
      data.frame(a = c(1, 2), row.names = c("s1", "s2")),
      c(s1 = 0, s1 = 1)
    ),
    "y.*duplicate",
    ignore.case = TRUE
  )
  expect_error(
    validate_sample_alignment(
      data.frame(a = c(1, 2), row.names = c("s1", "s2")),
      c(s1 = 0, s2 = 1),
      c(s1 = "a", s1 = "b")
    ),
    "groups.*duplicate",
    ignore.case = TRUE
  )
})

test_that("validate_sample_alignment rejects duplicate feature names", {
  x <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(c("s1", "s2"), c("gene", "gene"))
  )
  y <- c(s1 = 0, s2 = 1)

  expect_error(validate_sample_alignment(x, y), "feature names.*duplicate", ignore.case = TRUE)
})

test_that("validate_multiomic_inputs enforces named list and sample order", {
  x1 <- data.frame(a = c(1, 2), row.names = c("s1", "s2"))
  x2 <- data.frame(a = c(2, 3), row.names = c("s2", "s1"))
  y <- c(s1 = 0, s2 = 1)

  expect_error(
    validate_multiomic_inputs(list(omic1 = x1, omic2 = x2), y),
    "identical sample order"
  )
})

test_that("validate_multiomic_inputs rejects duplicate omic names and duplicate omic columns", {
  x1 <- data.frame(a = c(1, 2), row.names = c("s1", "s2"))
  x2 <- matrix(1:4, nrow = 2L, dimnames = list(c("s1", "s2"), c("g", "g")))
  x_dup_rows <- matrix(1:4, nrow = 2L, dimnames = list(c("s1", "s1"), c("a", "b")))
  y <- c(s1 = 0, s2 = 1)

  expect_error(validate_multiomic_inputs(list(omic = x1, omic = x1), y), "omic names.*duplicate", ignore.case = TRUE)
  expect_error(validate_multiomic_inputs(list(omic1 = x1, omic2 = x2), y), "omic2.*feature names.*duplicate", ignore.case = TRUE)
  expect_error(validate_multiomic_inputs(list(omic1 = x1, omic2 = x_dup_rows), y), "x.*row names.*duplicate", ignore.case = TRUE)
})

test_that("multiomic validation inputs reject duplicate validation IDs and features", {
  x_valid <- matrix(1:4, nrow = 2L, dimnames = list(c("v1", "v1"), c("a", "b")))
  expect_error(
    stablr:::.validate_multiomic_validation_inputs(
      x_valid_list = list(omic1 = x_valid),
      y_valid = NULL,
      train_omic_names = "omic1"
    ),
    "validation omic.*row names.*duplicate",
    ignore.case = TRUE
  )

  x_valid_cols <- matrix(1:4, nrow = 2L, dimnames = list(c("v1", "v2"), c("a", "a")))
  expect_error(
    stablr:::.validate_multiomic_validation_inputs(
      x_valid_list = list(omic1 = x_valid_cols),
      y_valid = NULL,
      train_omic_names = "omic1"
    ),
    "validation omic.*feature names.*duplicate",
    ignore.case = TRUE
  )

  y_valid <- c(v1 = 0, v1 = 1)
  expect_error(
    stablr:::.validate_multiomic_validation_inputs(
      x_valid_list = list(omic1 = matrix(1:4, nrow = 2L, dimnames = list(c("v1", "v2"), c("a", "b")))),
      y_valid = y_valid,
      train_omic_names = "omic1"
    ),
    "y.*duplicate",
    ignore.case = TRUE
  )
})

test_that("stabl_fit scalar validators reject non-finite, fractional, vector, and string inputs", {
  x <- matrix(rnorm(20), nrow = 5L,
              dimnames = list(paste0("s", 1:5), paste0("f", 1:4)))
  y <- setNames(rnorm(5L), rownames(x))
  lam_grid <- data.frame(lambda = c(0.1, 0.2))

  base_call <- function(lambda_grid = lam_grid, ...) stabl_fit(
    x = x, y = y, lambda_grid = lam_grid, family = "gaussian",
    artificial_type = NULL, hard_threshold = 0.3, ...
  )

  expect_error(base_call(n_bootstraps = NA_integer_), "n_bootstraps")
  expect_error(base_call(n_bootstraps = Inf), "n_bootstraps")
  expect_error(base_call(n_bootstraps = c(2L, 3L)), "n_bootstraps")
  expect_error(base_call(n_bootstraps = 2.5), "n_bootstraps")
  expect_error(base_call(n_bootstraps = "2"), "n_bootstraps")
  expect_error(base_call(n_bootstraps = 0L), "n_bootstraps")

  expect_error(base_call(sample_fraction = NA_real_), "sample_fraction")
  expect_error(base_call(sample_fraction = Inf), "sample_fraction")
  expect_error(base_call(sample_fraction = c(0.5, 0.6)), "sample_fraction")
  expect_error(base_call(sample_fraction = "0.5"), "sample_fraction")

  expect_error(base_call(hard_threshold = NA_real_), "hard_threshold")
  expect_error(base_call(hard_threshold = Inf), "hard_threshold")
  expect_error(base_call(hard_threshold = c(0.2, 0.3)), "hard_threshold")
  expect_error(base_call(hard_threshold = "0.3"), "hard_threshold")

  expect_error(base_call(lambda_grid = "auto", n_lambda = "3"), "n_lambda")
  expect_error(base_call(workers = 0L), "workers")
  expect_error(base_call(random_state = 1.5), "random_state")
})

test_that("stabl_fit errors early when sample_fraction > 1 and replace = FALSE", {
  # Fix 7: validator must fire before the bootstrap loop to give a clear message.
  set.seed(1L)
  n <- 20L; p <- 5L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam_grid <- data.frame(lambda = c(0.1, 0.2))

  expect_error(
    stabl_fit(
      x               = x,
      y               = y,
      lambda_grid     = lam_grid,
      family          = "gaussian",
      n_bootstraps    = 5L,
      artificial_type = NULL,
      hard_threshold  = 0.3,
      sample_fraction = 1.5,   # > 1 with replace = FALSE -> n_subsamples > n_samples
      replace         = FALSE
    ),
    "cannot exceed 1.*replace = FALSE|reduce.*sample_fraction|replace = TRUE",
    ignore.case = TRUE
  )
})

# C3 characterization: .resolve_cooperation_type_measure default mapping
# and .supported_cooperation_type_measures behaviour for unknown families.
test_that(".resolve_cooperation_type_measure maps 'default' to the correct measure per family", {
  # Pin the four supported family defaults at tolerance = 0.
  expect_identical(stablr:::.resolve_cooperation_type_measure("gaussian", "default"), "mse")
  expect_identical(stablr:::.resolve_cooperation_type_measure("binomial", "default"), "deviance")
  expect_identical(stablr:::.resolve_cooperation_type_measure("poisson",  "default"), "deviance")
  expect_identical(stablr:::.resolve_cooperation_type_measure("cox",      "default"), "deviance")
})

test_that(".resolve_cooperation_type_measure errors on an unsupported family", {
  # An unknown family must fail loudly (via .supported_cooperation_type_measures
  # returning character(0), making "default" %in% character(0) FALSE).
  expect_error(
    stablr:::.resolve_cooperation_type_measure("multinomial", "default"),
    "cooperation_type_measure"
  )
})

test_that(".supported_cooperation_type_measures returns character(0) for unknown family", {
  expect_identical(stablr:::.supported_cooperation_type_measures("unknown_family"), character(0))
})


test_that("multi-omic workflow scalar validators reject invalid values", {
  set.seed(441)
  n <- 12L
  ids <- paste0("s", seq_len(n))
  x <- matrix(rnorm(n * 4L), nrow = n,
              dimnames = list(ids, paste0("f", seq_len(4L))))
  y <- setNames(rnorm(n), ids)
  x_list <- list(omic_a = x)
  lambda_grid <- data.frame(lambda = 0.2)

  expect_error(stabl_multiomic_cv(x_list, y, lambda_grid, v = NA_integer_), "v")
  expect_error(stabl_multiomic_cv(x_list, y, lambda_grid, v = Inf), "v")
  expect_error(stabl_multiomic_cv(x_list, y, lambda_grid, v = 1.5), "v")
  expect_error(stabl_multiomic_cv(x_list, y, lambda_grid, v = "2"), "v")
  expect_error(stabl_multiomic_cv(x_list, y, lambda_grid, v = 1L), "v")

  expect_error(
    stabl_multiomic_train_validate(
      x_train_list = x_list,
      y_train = y,
      lambda_grid = lambda_grid,
      late_fusion = TRUE,
      n_iter_lf = NA_integer_
    ),
    "n_iter_lf"
  )
  expect_error(
    stabl_multiomic_train_validate(x_list, y, lambda_grid, n_iter_lf = 0L),
    "n_iter_lf"
  )
  expect_error(
    stabl_multiomic_train_validate(x_list, y, lambda_grid, n_iter_lf = 1.5),
    "n_iter_lf"
  )
  expect_error(
    stabl_multiomic_train_validate(x_list, y, lambda_grid, n_iter_lf = "5"),
    "n_iter_lf"
  )
})

test_that("nested multi-omic CV scalar validators reject invalid values", {
  set.seed(442)
  n <- 12L
  ids <- paste0("s", seq_len(n))
  x <- matrix(rnorm(n * 4L), nrow = n,
              dimnames = list(ids, paste0("f", seq_len(4L))))
  y <- setNames(rep(c("A", "B"), each = n / 2L), ids)
  x_list <- list(omic_a = x)
  lambda_grid <- data.frame(lambda = 0.2)
  base_call <- function(...) stabl_multiomic_nested_cv(
    x_list = x_list,
    y = y,
    lambda_grid = lambda_grid,
    outer_v = 2L,
    inner_v = 2L,
    outer_repeats = 1L,
    strata_bins = 2L,
    cv_workers = 1L,
    workers = 1L,
    n_lambda = 2L,
    n_bootstraps = 2L,
    artificial_type = NULL,
    hard_threshold = 1,
    ...
  )

  expect_error(base_call(outer_v = NA_integer_), "outer_v")
  expect_error(base_call(outer_v = Inf), "outer_v")
  expect_error(base_call(outer_v = 1.5), "outer_v")
  expect_error(base_call(outer_v = "2"), "outer_v")
  expect_error(base_call(outer_v = 1L), "outer_v")
  expect_error(base_call(inner_v = 1L), "inner_v")
  expect_error(base_call(outer_repeats = 0L), "outer_repeats")
  expect_error(base_call(strata_bins = 0L), "strata_bins")
  expect_error(base_call(cv_workers = 0L), "cv_workers")
  expect_error(base_call(workers = 0L), "workers")
  expect_error(base_call(n_lambda = 0L), "n_lambda")
  expect_error(base_call(random_state = 1.5), "random_state")
})
