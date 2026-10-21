# Tests for the ten deprecated shims introduced by Phase 14 PR-B
# (mixOmics-style API streamlining):
#   similarity:    jaccard_similarity, jaccard_similarity_matrix,
#                  adjusted_similarity, adjusted_similarity_pairwise,
#                  adjusted_similarity_summary, pearson_similarity,
#                  pearson_similarity_pairwise, pearson_similarity_summary
#   exports:       export_stabl_to_csv
#   preprocessing: stabl_preprocess_fit_recipes
#
# Each shim must (a) fire exactly one .Deprecated warning naming the
# canonical replacement and (b) forward identically to that replacement.
# The legacy-argument mapping of the similarity family is pinned separately
# in test-api-deprecation.R (block C2).  All shim calls capture their
# warnings via .catch_warnings() so the suite-wide warning budget is
# unchanged.

# Helper: capture all warning messages from an expression (same per-file
# definition as test-api-deprecation.R / test-stabl-wrapper.R).
.catch_warnings <- function(expr) {
  w <- list()
  withCallingHandlers(
    force(expr),
    warning = function(cm) {
      w[[length(w) + 1L]] <<- conditionMessage(cm)
      invokeRestart("muffleWarning")
    }
  )
  w
}

# -- similarity shims ----------------------------------------------------------

test_that("jaccard_similarity shim warns once and forwards identically", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")
  w <- .catch_warnings(old <- jaccard_similarity(a, b))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'jaccard_similarity' is deprecated", fixed = TRUE)
  expect_match(w[[1L]], "stabl_similarity", fixed = TRUE)
  expect_identical(old, stabl_similarity(a, b, method = "jaccard"))
})

test_that("jaccard_similarity_matrix shim warns once and forwards identically", {
  sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
  w <- .catch_warnings(old <- jaccard_similarity_matrix(sets))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'jaccard_similarity_matrix' is deprecated", fixed = TRUE)
  expect_identical(
    old,
    stabl_similarity(sets, method = "jaccard", as_matrix = TRUE, remove_diag = TRUE)
  )
  w2 <- .catch_warnings(old2 <- jaccard_similarity_matrix(sets, remove_diag = FALSE))
  expect_length(w2, 1L)
  expect_identical(
    old2,
    stabl_similarity(sets, method = "jaccard", as_matrix = TRUE, remove_diag = FALSE)
  )
})

test_that("adjusted_similarity shim warns once and forwards identically", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")
  w <- .catch_warnings(old <- adjusted_similarity(a, b, n_universe = 10L))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'adjusted_similarity' is deprecated", fixed = TRUE)
  expect_identical(old, stabl_similarity(a, b, method = "adjusted", n_universe = 10L))
})

test_that("adjusted_similarity_pairwise shim warns once and forwards identically", {
  sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
  w <- .catch_warnings(old <- adjusted_similarity_pairwise(sets, 10L))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'adjusted_similarity_pairwise' is deprecated", fixed = TRUE)
  expect_identical(old, stabl_similarity(sets, method = "adjusted", n_universe = 10L))
})

test_that("adjusted_similarity_summary shim warns once and forwards identically", {
  sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
  w <- .catch_warnings(old <- adjusted_similarity_summary(sets, 10L))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'adjusted_similarity_summary' is deprecated", fixed = TRUE)
  expect_match(w[[1L]], "stabl_similarity_summary", fixed = TRUE)
  expect_identical(
    old,
    stabl_similarity_summary(sets, method = "adjusted", n_universe = 10L)
  )
  w2 <- .catch_warnings(old2 <- adjusted_similarity_summary(sets, 10L, stat = "mean"))
  expect_length(w2, 1L)
  expect_identical(
    old2,
    stabl_similarity_summary(sets, method = "adjusted", n_universe = 10L, stat = "mean")
  )
})

test_that("pearson_similarity shim warns once and forwards identically", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")
  w <- .catch_warnings(old <- pearson_similarity(a, b, n_universe = 10L))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'pearson_similarity' is deprecated", fixed = TRUE)
  expect_identical(old, stabl_similarity(a, b, method = "pearson", n_universe = 10L))
})

test_that("pearson_similarity_pairwise shim warns once and forwards identically", {
  sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
  w <- .catch_warnings(old <- pearson_similarity_pairwise(sets, 10L))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'pearson_similarity_pairwise' is deprecated", fixed = TRUE)
  expect_identical(old, stabl_similarity(sets, method = "pearson", n_universe = 10L))
})

test_that("pearson_similarity_summary shim warns once and forwards identically", {
  sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
  w <- .catch_warnings(old <- pearson_similarity_summary(sets, 10L))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'pearson_similarity_summary' is deprecated", fixed = TRUE)
  expect_match(w[[1L]], "stabl_similarity_summary", fixed = TRUE)
  expect_identical(
    old,
    stabl_similarity_summary(sets, method = "pearson", n_universe = 10L)
  )
})

# -- exports shim --------------------------------------------------------------

test_that("export_stabl_to_csv shim warns once and forwards identically", {
  set.seed(1L)
  x <- matrix(rnorm(40 * 6), 40, 6,
              dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
  y <- setNames(rnorm(40), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid  = data.frame(lambda = c(0.3, 0.1, 0.05)),
                   n_bootstraps = 6L, hard_threshold = 0.3, random_state = 1L)

  dir_old <- tempfile("shim_export_old")
  dir_new <- tempfile("shim_export_new")
  on.exit(unlink(c(dir_old, dir_new), recursive = TRUE), add = TRUE)

  w <- .catch_warnings(ret_old <- export_stabl_to_csv(fit, dir_old))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'export_stabl_to_csv' is deprecated", fixed = TRUE)
  expect_match(w[[1L]], "save_stabl_results", fixed = TRUE)

  ret_new <- save_stabl_results(fit, dir_new, what = "csv")
  expect_identical(ret_old, normalizePath(dir_old))
  expect_identical(ret_new, normalizePath(dir_new))

  files_old <- sort(list.files(dir_old, recursive = TRUE))
  files_new <- sort(list.files(dir_new, recursive = TRUE))
  expect_identical(files_old, files_new)
  expect_true(length(files_old) >= 2L)
  for (f in files_old) {
    expect_identical(
      utils::read.csv(file.path(dir_old, f)),
      utils::read.csv(file.path(dir_new, f)),
      info = f
    )
  }
})

# -- preprocessing shim --------------------------------------------------------

test_that("stabl_preprocess_fit_recipes shim warns once and forwards identically", {
  skip_if_not_installed("recipes")
  set.seed(2L)
  x <- matrix(rnorm(40 * 4), 40, 4,
              dimnames = list(paste0("s", 1:40), paste0("f", 1:4)))
  rec <- recipes::recipe(~ ., data = as.data.frame(x))
  rec <- recipes::step_center(rec, recipes::all_predictors())

  w <- .catch_warnings(old <- stabl_preprocess_fit_recipes(rec, x))
  expect_length(w, 1L)
  expect_match(w[[1L]], "'stabl_preprocess_fit_recipes' is deprecated", fixed = TRUE)
  expect_match(w[[1L]], "stabl_preprocess_fit", fixed = TRUE)

  new <- stabl_preprocess_fit(x, transforms = rec)
  expect_identical(class(old), class(new))
  expect_identical(old$retained_features, new$retained_features)
  expect_identical(stabl_preprocess_apply(old, x), stabl_preprocess_apply(new, x))
})
