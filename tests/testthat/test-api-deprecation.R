# Tests for the API deprecation shims introduced by the 0.2.0 renames:
#   * C1: sample_bootstrap_indices() / group_bootstrap_indices()
#         `seed` -> `random_state`, plus RNG hygiene (withr::local_seed())
#   * C2: adjusted_similarity() / pearson_similarity()
#         `list1`/`list2`/`nb_total_elements` and `list_i`/`list_j`/`d`
#         -> `set1`/`set2`/`n_universe`
#
# All seeds are explicit so the expectations are deterministic.

# Helper: capture all warning messages from an expression (robust across
# testthat editions; also lets us assert the exact warning count).
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

# -- C1: sample_bootstrap_indices() ------------------------------------------

test_that("sample_bootstrap_indices: deprecated `seed` warns once and matches `random_state`", {
  y <- c(rep(0, 15), rep(1, 15))

  w <- .catch_warnings(
    old_style <- sample_bootstrap_indices(y, n_subsamples = 20L, seed = 1L)
  )
  expect_length(w, 1L)
  expect_match(
    w[[1L]],
    "`seed` is deprecated; use `random_state` instead.",
    fixed = TRUE
  )

  new_style <- sample_bootstrap_indices(y, n_subsamples = 20L, random_state = 1L)
  expect_identical(old_style, new_style)
})

test_that("sample_bootstrap_indices: explicit seed = NULL still warns", {
  y <- c(rep(0, 15), rep(1, 15))
  w <- .catch_warnings(sample_bootstrap_indices(y, n_subsamples = 5L, seed = NULL))
  expect_length(w, 1L)
  expect_match(w[[1L]], "deprecated", fixed = TRUE)
})

test_that("sample_bootstrap_indices: new name does not warn", {
  y <- c(rep(0, 15), rep(1, 15))
  w <- .catch_warnings(sample_bootstrap_indices(y, n_subsamples = 5L, random_state = 1L))
  expect_length(w, 0L)
})

test_that("sample_bootstrap_indices with random_state leaves the global RNG state untouched", {
  # Ensure .Random.seed exists before capturing (it may not in a fresh session).
  if (!exists(".Random.seed", globalenv(), inherits = FALSE)) runif(1)
  before <- get(".Random.seed", globalenv())

  idx <- sample_bootstrap_indices(
    y = c(rep(0, 10), rep(1, 10)), n_subsamples = 8, random_state = 42
  )
  expect_identical(get(".Random.seed", globalenv()), before)
  expect_length(idx, 8L)

  # Same for the deprecated spelling (deprecation warning muffled; asserted
  # in the dedicated warning tests above).
  idx_old <- suppressWarnings(sample_bootstrap_indices(
    y = c(rep(0, 10), rep(1, 10)), n_subsamples = 8, seed = 42
  ))
  expect_identical(get(".Random.seed", globalenv()), before)
  expect_identical(idx_old, idx)
})

test_that("sample_bootstrap_indices: random_state reproduces the pre-rename set.seed() convention", {
  y <- c(rep(0, 10), rep(1, 10))
  # Pre-rename convention: set.seed(seed) followed by the plain sample.int draw.
  set.seed(99)
  expected <- sample.int(20, 8, replace = TRUE)
  idx <- sample_bootstrap_indices(y, n_subsamples = 8, random_state = 99)
  expect_identical(as.integer(idx), as.integer(expected))
})

test_that("sample_bootstrap_indices: ambient RNG stream is used when random_state = NULL (stabl_fit caller convention)", {
  y <- c(rep(0, 10), rep(1, 10))
  set.seed(5)
  expected <- sample.int(20, 8, replace = TRUE)
  set.seed(5)
  idx <- sample_bootstrap_indices(y, n_subsamples = 8)
  expect_identical(as.integer(idx), as.integer(expected))
})

# -- C1: group_bootstrap_indices() --------------------------------------------

test_that("group_bootstrap_indices: deprecated `seed` warns once and matches `random_state`", {
  y      <- c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
  groups <- rep(c("G1", "G2", "G3"), each = 4L)

  w <- .catch_warnings(
    old_style <- group_bootstrap_indices(y, groups, n_subsamples = 8L, seed = 1L)
  )
  expect_length(w, 1L)
  expect_match(
    w[[1L]],
    "`seed` is deprecated; use `random_state` instead.",
    fixed = TRUE
  )

  new_style <- group_bootstrap_indices(y, groups, n_subsamples = 8L, random_state = 1L)
  expect_identical(old_style, new_style)
})

test_that("group_bootstrap_indices with random_state leaves the global RNG state untouched", {
  if (!exists(".Random.seed", globalenv(), inherits = FALSE)) runif(1)
  before <- get(".Random.seed", globalenv())

  y      <- c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
  groups <- rep(c("G1", "G2", "G3"), each = 4L)
  idx <- group_bootstrap_indices(y, groups, n_subsamples = 8L, random_state = 7L)
  expect_identical(get(".Random.seed", globalenv()), before)
  expect_true(all(table(groups[idx]) > 0L))

  # Reproducible across calls with the same random_state.
  idx2 <- group_bootstrap_indices(y, groups, n_subsamples = 8L, random_state = 7L)
  expect_identical(idx, idx2)
})

# -- C2: adjusted_similarity() ------------------------------------------------

test_that("adjusted_similarity: deprecated names warn once and map to the new names", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")

  w <- .catch_warnings(
    old_style <- adjusted_similarity(list1 = a, list2 = b, nb_total_elements = 10L)
  )
  expect_length(w, 1L)
  expect_match(
    w[[1L]],
    "`list1`, `list2` and `nb_total_elements` are deprecated; use `set1`, `set2` and `n_universe` instead.",
    fixed = TRUE
  )

  new_style <- adjusted_similarity(set1 = a, set2 = b, n_universe = 10L)
  expect_identical(old_style, new_style)

  # Mixed positional/new + named-old also warns and maps.
  w2 <- .catch_warnings(
    mixed <- adjusted_similarity(a, list2 = b, nb_total_elements = 10L)
  )
  expect_length(w2, 1L)
  expect_identical(mixed, new_style)
})

test_that("adjusted_similarity: positional calls are unchanged and silent", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")
  w <- .catch_warnings(pos <- adjusted_similarity(a, b, 10L))
  expect_length(w, 0L)
  expect_identical(pos, adjusted_similarity(set1 = a, set2 = b, n_universe = 10L))
  expect_identical(pos, adjusted_similarity(a, b, n_universe = 10L))
})

test_that("adjusted_similarity: edge cases and errors identical under the new names", {
  expect_identical(
    adjusted_similarity(character(0), c("f1"), n_universe = 10L), 0
  )
  # Identical sets: chance-corrected value (0.8 here), same under old spelling.
  expect_equal(
    adjusted_similarity(c("f1", "f2"), c("f1", "f2"), n_universe = 10L), 0.8
  )
  expect_equal(
    suppressWarnings(
      adjusted_similarity(list1 = c("f1", "f2"), list2 = c("f1", "f2"),
                          nb_total_elements = 10L)
    ),
    0.8
  )
  expect_error(
    adjusted_similarity(letters, letters, n_universe = 5L),
    "exceeds n_universe",
    fixed = TRUE
  )
  expect_error(
    adjusted_similarity(set1 = letters, set2 = letters),
    "n_universe"
  )
})

# -- C2: pearson_similarity() -------------------------------------------------

test_that("pearson_similarity: deprecated names warn once and map to the new names", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")

  w <- .catch_warnings(
    old_style <- pearson_similarity(list_i = a, list_j = b, d = 10L)
  )
  expect_length(w, 1L)
  expect_match(
    w[[1L]],
    "`list_i`, `list_j` and `d` are deprecated; use `set1`, `set2` and `n_universe` instead.",
    fixed = TRUE
  )

  new_style <- pearson_similarity(set1 = a, set2 = b, n_universe = 10L)
  expect_identical(old_style, new_style)
})

test_that("pearson_similarity: positional calls are unchanged and silent", {
  a <- c("f1", "f2", "f3")
  b <- c("f2", "f3", "f4")
  w <- .catch_warnings(pos <- pearson_similarity(a, b, 10L))
  expect_length(w, 0L)
  expect_identical(pos, pearson_similarity(set1 = a, set2 = b, n_universe = 10L))
  expect_identical(pos, pearson_similarity(a, b, n_universe = 10L))
})

test_that("pearson_similarity: edge cases identical under the new names", {
  expect_identical(
    pearson_similarity(character(0), c("f1", "f2"), n_universe = 10L), 0
  )
  expect_identical(
    pearson_similarity(c("f1", "f2"), c("f1", "f2"), n_universe = 10L), 1
  )
  expect_error(
    pearson_similarity(set1 = letters, set2 = letters),
    "n_universe"
  )
})
