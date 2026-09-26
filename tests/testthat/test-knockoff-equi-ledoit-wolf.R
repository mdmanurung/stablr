.lw_fixture_dir <- function() {
  testthat::test_path("fixtures", "ledoit_wolf")
}

.read_named_fixture_matrix <- function(name) {
  frame <- utils::read.csv(
    file.path(.lw_fixture_dir(), name),
    row.names = 1L,
    check.names = FALSE
  )
  as.matrix(frame)
}

.fixture_sha256 <- function(path) {
  unname(tools::sha256sum(path))
}

test_that("Ledoit-Wolf fixture closure is intact", {
  fixture_dir <- .lw_fixture_dir()
  manifest <- readLines(file.path(fixture_dir, "manifest.sha256"), warn = FALSE)
  expect_true(all(grepl("^[0-9a-f]{64}  [A-Za-z0-9_.-]+$", manifest)))
  expected <- substr(manifest, 1L, 64L)
  files <- substring(manifest, 67L)
  expect_equal(anyDuplicated(files), 0L)
  expect_true(all(file.exists(file.path(fixture_dir, files))))
  observed <- vapply(
    file.path(fixture_dir, files),
    .fixture_sha256,
    character(1L)
  )
  expect_identical(unname(observed), expected)
})

test_that("Ledoit-Wolf covariance matches scikit-learn 1.7.1 fixture", {
  x <- .read_named_fixture_matrix("input.csv")
  expected <- .read_named_fixture_matrix("covariance.csv")
  metadata <- jsonlite::read_json(
    file.path(.lw_fixture_dir(), "metadata.json"),
    simplifyVector = TRUE
  )

  observed <- stablr:::.estimate_ledoit_wolf_sigma(x)
  expect_equal(observed$covariance, expected, tolerance = 2e-14)
  expect_equal(
    observed$shrinkage,
    metadata$ledoit_wolf$shrinkage,
    tolerance = 2e-14
  )
  expect_equal(observed$location, unname(colMeans(x)), tolerance = 0)
  expect_gt(observed$minimum_eigenvalue, 0)
})

test_that("equicorrelated S-matrix matches knockpy 1.3.5 fixture", {
  covariance <- .read_named_fixture_matrix("covariance.csv")
  expected <- .read_named_fixture_matrix("s_matrix.csv")
  metadata <- jsonlite::read_json(
    file.path(.lw_fixture_dir(), "metadata.json"),
    simplifyVector = TRUE
  )

  observed <- stablr:::.knockpy_equicorrelated_s(covariance)
  expect_equal(unname(diag(observed$diag_s)), unname(expected), tolerance = 2e-14)
  expect_equal(observed$tolerance, 1e-5, tolerance = 0)
  expect_identical(observed$iterations, 10L)
  expect_equal(
    observed$conditional_minimum_eigenvalue,
    metadata$equicorrelated_s$minimum_conditional_eigenvalue,
    tolerance = 2e-12
  )
  expect_gte(observed$conditional_minimum_eigenvalue, observed$tolerance)
})

test_that("Ledoit-Wolf covariance is positive definite when p greatly exceeds n", {
  set.seed(1101L)
  x <- matrix(rnorm(12L * 80L), nrow = 12L)
  observed <- stablr:::.estimate_ledoit_wolf_sigma(x)

  expect_equal(dim(observed$covariance), c(80L, 80L))
  expect_true(all(is.finite(observed$covariance)))
  expect_gt(min(eigen(observed$covariance, symmetric = TRUE)$values), 0)
  expect_true(observed$shrinkage >= 0 && observed$shrinkage <= 1)
})

test_that("invalid Ledoit-Wolf covariance fails with a typed condition", {
  x <- matrix(1, nrow = 8L, ncol = 5L)
  expect_error(
    stablr:::.estimate_ledoit_wolf_sigma(x),
    class = "stablr_knockoff_infeasible"
  )
})

test_that("artificial-feature local seeds preserve caller RNG state", {
  set.seed(1102L)
  x <- matrix(rnorm(60L), nrow = 12L)
  caller_state <- .Random.seed

  first <- make_artificial_features(
    x,
    n_injected = 3L,
    type = "random_permutation",
    random_state = 88L
  )
  expect_identical(.Random.seed, caller_state)
  second <- make_artificial_features(
    x,
    n_injected = 3L,
    type = "random_permutation",
    random_state = 88L
  )
  expect_identical(first, second)
  expect_identical(.Random.seed, caller_state)
})

test_that("valid p-greater-than-n equicorrelated generation has no fallback", {
  skip_if_not_installed("knockoff")
  set.seed(1103L)
  x <- matrix(
    rnorm(16L * 40L),
    nrow = 16L,
    dimnames = list(paste0("s", 1:16), paste0("f", 1:40))
  )
  caller_state <- .Random.seed

  first <- expect_warning(
    make_artificial_features(
      x,
      n_injected = 20L,
      type = "knockoff_equi",
      random_state = 1103L
    ),
    NA
  )
  second <- expect_warning(
    make_artificial_features(
      x,
      n_injected = 20L,
      type = "knockoff_equi",
      random_state = 1103L
    ),
    NA
  )

  expect_identical(first, second)
  expect_identical(.Random.seed, caller_state)
  expect_identical(
    first$artificial_provenance$actual_type,
    "knockoff_equi"
  )
  expect_equal(nrow(first$artificial_provenance$fallback_history), 0L)
  expect_identical(
    first$artificial_provenance$covariance_estimator$name,
    "sklearn.covariance.LedoitWolf"
  )
  expect_identical(
    first$artificial_provenance$covariance_estimator$compatibility_version,
    "1.7.1"
  )
  expect_identical(
    first$artificial_provenance$equicorrelated_s_matrix$name,
    "knockpy.smatrix.equicorrelated"
  )
})

# ── A11: overflow-scale designs must be rejected, not overflowed ─────────────

test_that("overflow-scale predictors are rejected before the crossproduct", {
  condition <- expect_error(
    stablr:::.estimate_ledoit_wolf_sigma(matrix(1e200, 8L, 5L)),
    class = "stablr_knockoff_infeasible"
  )
  expect_s3_class(condition, "stablr_numerical_infeasibility")
  expect_match(conditionMessage(condition), "overflow")
  expect_false(grepl("positive-definite", conditionMessage(condition),
                     fixed = TRUE))
})

test_that("overflow-scale rejection does not depend on error classification", {
  # The guard must be self-sufficient: with base::eigen's message no longer
  # classified, the overflow-scale input must still fail typed.
  testthat::local_mocked_bindings(
    .is_linear_algebra_numerical_error = function(condition) FALSE,
    .package = "stablr"
  )
  expect_error(
    stablr:::.estimate_ledoit_wolf_sigma(matrix(1e200, 8L, 5L)),
    class = "stablr_knockoff_infeasible"
  )
  expect_error(
    stablr:::.estimate_ledoit_wolf_sigma(matrix(1e300, 8L, 5L)),
    class = "stablr_knockoff_infeasible"
  )
})

test_that("ordinary-scale predictors are unaffected by the overflow guard", {
  set.seed(1104L)
  # A single predictor has zero Ledoit-Wolf shrinkage, so the estimate is the
  # raw maximum-likelihood covariance and the guard's inertness is observable.
  x <- matrix(rnorm(24L), ncol = 1L)
  observed <- stablr:::.estimate_ledoit_wolf_sigma(x)
  expect_identical(observed$shrinkage, 0)
  expect_equal(
    observed$covariance,
    matrix(sum((x - mean(x))^2) / nrow(x), 1L, 1L),
    tolerance = 1e-15
  )

  multi <- matrix(rnorm(24L * 6L), nrow = 24L)
  multi_fit <- stablr:::.estimate_ledoit_wolf_sigma(multi)
  expect_true(all(is.finite(multi_fit$covariance)))
  expect_equal(multi_fit$covariance, t(multi_fit$covariance), tolerance = 0)
  expect_gte(multi_fit$shrinkage, 0)
  expect_lte(multi_fit$shrinkage, 1)
  expect_gt(
    min(eigen(multi_fit$covariance, symmetric = TRUE, only.values = TRUE)$values),
    0
  )
})

test_that("overflow-scale equicorrelated design falls back with a reason", {
  skip_if_not_installed("knockoff")
  set.seed(1105L)
  x <- matrix(
    rnorm(20L * 4L) * 1e200,
    nrow = 20L,
    dimnames = list(paste0("s", seq_len(20L)), paste0("f", seq_len(4L)))
  )

  expect_warning(
    result <- make_artificial_features(
      x, n_injected = 2L, type = "knockoff_equi", random_state = 1105L
    ),
    "falling back to random permutation"
  )
  expect_identical(
    result$artificial_provenance$actual_type,
    "random_permutation"
  )
  expect_match(
    result$artificial_provenance$chunks$fallback_reason,
    "overflow"
  )
})
