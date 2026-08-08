test_that("MVR solver failure falls back to equi with recorded provenance", {
  skip_if_not_installed("knockoff")

  testthat::local_mocked_bindings(
    solve_mvr = function(...) {
      stablr:::.abort_numerical_infeasibility(
        "stablr_mvr_infeasible",
        "deliberate MVR solver failure"
      )
    },
    .package = "stablr"
  )
  set.seed(4201L)
  x <- matrix(
    rnorm(80L * 5L),
    nrow = 80L,
    dimnames = list(paste0("s", seq_len(80L)), paste0("f", seq_len(5L)))
  )

  result <- NULL
  expect_warning(
    result <- make_artificial_features(
      x = x,
      n_injected = 3L,
      type = "knockoff_mvr",
      random_state = 4201L
    ),
    "solve_mvr failed; using equi S"
  )

  provenance <- result$artificial_provenance
  expect_equal(provenance$requested_type, "knockoff_mvr")
  expect_equal(provenance$actual_type, "knockoff_equi")
  expect_equal(provenance$n_chunks, 1L)
  expect_equal(provenance$fallback_counts[["knockoff_equi"]], 1L)
  expect_equal(provenance$selected_type_counts[["knockoff_equi"]], 3L)
  expect_true(provenance$chunks$fallback)
  expect_match(
    provenance$chunks$fallback_reason,
    "deliberate MVR solver failure",
    fixed = TRUE
  )
  expect_identical(provenance$fallback_history$event, 1L)
  expect_identical(provenance$fallback_history$from_type, "knockoff_mvr")
  expect_identical(provenance$fallback_history$to_type, "knockoff_equi")
})

test_that("high-dimensional knockoff construction fails closed", {
  x <- matrix(
    seq_len(2L * 3001L),
    nrow = 2L,
    dimnames = list(c("s1", "s2"), paste0("f", seq_len(3001L)))
  )

  set.seed(4202L)
  before <- .Random.seed
  for (type in c("knockoff", "knockoff_equi", "knockoff_mvr")) {
    condition <- tryCatch(
      make_artificial_features(
        x = x,
        n_injected = 3L,
        type = type,
        random_state = 4202L
      ),
      error = identity
    )
    expect_s3_class(
      condition,
      "stablr_high_dimensional_knockoff_unsupported"
    )
    expect_identical(condition$requested_type, type)
    expect_identical(condition$n_features, 3001L)
    expect_identical(condition$maximum_supported_features, 3000L)
  }
  expect_identical(.Random.seed, before)
})
