test_that("make_artificial_features records fixed-X knockoff fallback provenance", {
  skip_if_not_installed("knockoff")

  set.seed(1001L)
  x <- matrix(
    rnorm(9L),
    nrow = 3L,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:3))
  )

  saw_fallback_warning <- FALSE
  result <- withCallingHandlers(
    make_artificial_features(
      x = x,
      n_injected = 2L,
      type = "knockoff",
      random_state = 1001L
    ),
    warning = function(w) {
      if (grepl("falling back to random permutation", conditionMessage(w),
                fixed = TRUE)) {
        saw_fallback_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(saw_fallback_warning)
  expect_equal(result$artificial_provenance$requested_type, "knockoff")
  expect_equal(
    result$artificial_provenance$actual_type,
    "random_permutation"
  )
  expect_identical(
    result$artificial_provenance$actual_types,
    "random_permutation"
  )
  expect_equal(result$artificial_provenance$n_chunks, 1L)
  expect_equal(result$artificial_provenance$fallback_counts[["random_permutation"]], 1L)
  expect_equal(result$artificial_provenance$selected_type_counts[["random_permutation"]], 2L)
  expect_true(all(result$artificial_provenance$chunks$fallback))
  expect_match(
    result$artificial_provenance$chunks$fallback_reason,
    "Input X must have dimensions n > p"
  )
  expect_identical(
    result$artificial_provenance$fallback_history$event,
    1L
  )
  expect_identical(
    result$artificial_provenance$fallback_history$from_type,
    "knockoff"
  )
  expect_identical(
    result$artificial_provenance$fallback_history$to_type,
    "random_permutation"
  )
})

test_that("fixed-X row augmentation falls back instead of truncating invalid knockoffs", {
  skip_if_not_installed("knockoff")

  # For p < n < 2p, knockoff::create.fixed() augments the design to 2p rows.
  # Truncating only Xk back to n rows destroys the fixed-X Gram identities, so
  # the current stablr interface must reject that result and use its documented
  # random-permutation fallback.
  set.seed(1003L)
  x <- matrix(
    rnorm(15L * 10L),
    nrow = 15L,
    dimnames = list(paste0("s", 1:15), paste0("f", 1:10))
  )

  saw_fallback_warning <- FALSE
  result <- withCallingHandlers(
    make_artificial_features(
      x = x,
      n_injected = 10L,
      type = "knockoff",
      random_state = 1003L
    ),
    warning = function(w) {
      if (grepl("falling back to random permutation", conditionMessage(w),
                fixed = TRUE)) {
        saw_fallback_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(saw_fallback_warning)
  expect_equal(dim(result$x_augmented), c(nrow(x), 2L * ncol(x)))
  expect_equal(
    result$artificial_provenance$selected_type_counts[["random_permutation"]],
    ncol(x)
  )
  expect_match(
    result$artificial_provenance$chunks$fallback_reason,
    "augmented.*rows"
  )
})

test_that("stabl_fit preserves artificial-feature fallback provenance", {
  skip_if_not_installed("knockoff")

  set.seed(1002L)
  x <- matrix(
    rnorm(9L),
    nrow = 3L,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:3))
  )
  y <- setNames(rnorm(3L), rownames(x))

  saw_fallback_warning <- FALSE
  fit <- withCallingHandlers(
    stabl_fit(
      x = x,
      y = y,
      lambda_grid = data.frame(lambda = 0.1),
      family = "gaussian",
      n_bootstraps = 1L,
      artificial_type = "knockoff",
      artificial_proportion = 2 / 3,
      sample_fraction = 1,
      random_state = 1002L
    ),
    warning = function(w) {
      if (grepl("falling back to random permutation", conditionMessage(w),
                fixed = TRUE)) {
        saw_fallback_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(saw_fallback_warning)
  expect_s3_class(fit, "stabl_fit")
  expect_equal(fit$artificial_type, "knockoff")
  expect_equal(fit$artificial_provenance$requested_type, "knockoff")
  expect_equal(fit$artificial_provenance$actual_type, "random_permutation")
  expect_equal(fit$artificial_provenance$n_generated, 2L)
  expect_equal(fit$artificial_provenance$fallback_counts[["random_permutation"]], 1L)
  expect_equal(fit$artificial_provenance$selected_type_counts[["random_permutation"]], 2L)
})

test_that("artificial fallback history retains nested event order", {
  chunks <- stablr:::.artificial_chunk_record(
    chunk = 1L,
    requested_type = "knockoff_mvr",
    actual_type = "random_permutation",
    n_columns = 4L,
    fallback_reason = "sampling failed"
  )
  first <- simpleError("solver failed")
  class(first) <- c("stablr_mvr_infeasible", class(first))
  second <- simpleError("sampling failed")
  history <- rbind(
    stablr:::.artificial_fallback_event(
      "knockoff_mvr",
      "knockoff_equi",
      first
    ),
    stablr:::.artificial_fallback_event(
      "knockoff_equi",
      "random_permutation",
      second
    )
  )

  provenance <- stablr:::.make_artificial_provenance(
    requested_type = "knockoff_mvr",
    selected_types = rep("random_permutation", 2L),
    chunks = chunks,
    fallback_histories = list(history)
  )

  expect_identical(provenance$actual_type, "random_permutation")
  expect_identical(provenance$fallback_history$event, 1:2)
  expect_identical(provenance$fallback_history$step, 1:2)
  expect_identical(
    provenance$fallback_history$from_type,
    c("knockoff_mvr", "knockoff_equi")
  )
  expect_identical(
    provenance$fallback_history$to_type,
    c("knockoff_equi", "random_permutation")
  )
  expect_match(
    provenance$fallback_history$condition_class[[1L]],
    "stablr_mvr_infeasible",
    fixed = TRUE
  )
})
