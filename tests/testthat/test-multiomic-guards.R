# Audit V-9: hard guards on multi-omic input contracts.
# (Existing test-input-validation.R covers happy-path basics; this file
# pins the failure paths the audit flagged as untested.)

test_that("validate_multiomic_inputs rejects mismatched sample order across omics", {
  s_ids <- paste0("s", 1:6)
  x1 <- matrix(rnorm(6 * 3), nrow = 6, dimnames = list(s_ids, paste0("a", 1:3)))
  x2 <- matrix(rnorm(6 * 3), nrow = 6, dimnames = list(rev(s_ids), paste0("b", 1:3)))
  y  <- setNames(rnorm(6), s_ids)

  expect_error(
    validate_multiomic_inputs(list(omicA = x1, omicB = x2), y = y),
    "identical sample order"
  )
})

test_that("cooperative_fusion + family = 'cox' is rejected in native v1", {
  skip_if_not_installed("survival")

  s_ids <- paste0("s", 1:30)
  x1 <- matrix(rnorm(30 * 4), nrow = 30, dimnames = list(s_ids, paste0("a", 1:4)))
  x2 <- matrix(rnorm(30 * 4), nrow = 30, dimnames = list(s_ids, paste0("b", 1:4)))
  y  <- survival::Surv(time = abs(rnorm(30)) + 0.1,
                        event = sample(0:1, 30, replace = TRUE))
  rownames(y) <- s_ids

  expect_error(
    stabl_multiomic_train_validate(
      x_train_list           = list(omicA = x1, omicB = x2),
      y_train                = y,
      lambda_grid            = data.frame(lambda = c(0.2, 0.1)),
      family                 = "cox",
      cooperative_fusion     = TRUE,
      n_bootstraps           = 5L,
      artificial_type        = "random_permutation",
      random_state           = 1L
    ),
    "gaussian.*binomial|binomial.*gaussian"
  )
})

test_that("cooperative_fusion errors when cooperative backend is unavailable", {
  s_ids <- paste0("s", 1:30)
  x1 <- matrix(rnorm(30 * 4), nrow = 30, dimnames = list(s_ids, paste0("a", 1:4)))
  x2 <- matrix(rnorm(30 * 4), nrow = 30, dimnames = list(s_ids, paste0("b", 1:4)))
  y  <- setNames(rnorm(30), s_ids)

  testthat::local_mocked_bindings(
    .has_cooperative_backend = function() FALSE,
    .package = "stablr"
  )

  expect_error(
    stabl_multiomic_train_validate(
      x_train_list      = list(omicA = x1, omicB = x2),
      y_train           = y,
      lambda_grid       = data.frame(lambda = c(0.2, 0.1)),
      family            = "gaussian",
      cooperative_fusion = TRUE,
      n_bootstraps      = 5L,
      artificial_type   = "random_permutation",
      random_state      = 1L
    ),
    "Cooperative fusion backend is unavailable"
  )
})


# ---------------------------------------------------------------------------
# Stacking scorability guards.  A candidate set that cannot be scored is a
# typed error, never a returned -Inf score with an all-NA stacked column and
# never a no-information metric value masquerading as evidence.
# ---------------------------------------------------------------------------

.stacking_guard_predictions <- function(complete_rows) {
  values <- rbind(
    c(0.2, 0.3), c(0.9, 0.8), c(0.4, 0.6), c(0.7, 0.1)
  )
  out <- data.frame(
    rna = rep(NA_real_, 4L), protein = rep(NA_real_, 4L),
    row.names = paste0("s", 1:4)
  )
  for (row in complete_rows) out[row, ] <- values[row, ]
  out
}

.stacking_guard_outcome <- function(values) stats::setNames(values, paste0("s", 1:4))

test_that("no scored row at all aborts and names the observed count", {
  predictions <- .stacking_guard_predictions(integer(0))
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, .stacking_guard_outcome(c(0, 1, 0, 1)),
      task_type = "binary", n_iter = 25L, random_state = 1L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
  expect_s3_class(condition, "stablr_error")
  expect_identical(condition$complete_rows, 0L)
  expect_identical(condition$n_samples, 4L)
  expect_identical(condition$task_type, "binary")
  expect_match(conditionMessage(condition), "0 of 4", fixed = TRUE)
})

test_that("one scored row aborts and names the observed count", {
  predictions <- .stacking_guard_predictions(2L)
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, .stacking_guard_outcome(c(0, 1, 0, 1)),
      task_type = "binary", n_iter = 25L, random_state = 1L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
  expect_identical(condition$complete_rows, 1L)
  expect_match(conditionMessage(condition), "1 of 4", fixed = TRUE)
})

test_that("a single-class scored subset aborts instead of reporting AUC 0.5", {
  predictions <- .stacking_guard_predictions(c(1L, 3L))
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, .stacking_guard_outcome(c(0, 1, 0, 1)),
      task_type = "binary", n_iter = 25L, random_state = 1L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
  expect_identical(condition$complete_rows, 2L)
  expect_match(conditionMessage(condition), "only 0", fixed = TRUE)
  expect_identical(as.integer(condition$scored_class_counts), c(2L, 0L))
})

test_that("a constant scored outcome aborts instead of reporting R2 = 0", {
  predictions <- .stacking_guard_predictions(c(1L, 3L))
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, .stacking_guard_outcome(c(0, 1, 0, 1)),
      task_type = "regression", n_iter = 25L, random_state = 1L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
  expect_identical(condition$task_type, "regression")
  expect_match(conditionMessage(condition), "constant")
})

test_that("the guard fires on scorability, not on the number of rows alone", {
  two_rows_two_classes <- .stacking_guard_predictions(c(1L, 3L))
  expect_no_error(
    stacked_multi_omic(
      two_rows_two_classes, .stacking_guard_outcome(c(0, 1, 1, 0)),
      task_type = "binary", n_iter = 25L, random_state = 1L
    )
  )
  expect_no_error(
    stacked_multi_omic(
      two_rows_two_classes, .stacking_guard_outcome(c(0.5, 1.5, 2.5, 3.5)),
      task_type = "regression", n_iter = 25L, random_state = 1L
    )
  )
  complete <- .stacking_guard_predictions(1:4)
  expect_no_error(
    stacked_multi_omic(
      complete, .stacking_guard_outcome(c(0, 1, 0, 1)),
      task_type = "binary", n_iter = 25L, random_state = 1L
    )
  )
})

test_that("the scorability guard leaves the caller's RNG stream untouched", {
  set.seed(4242L)
  before <- .Random.seed
  try(
    stacked_multi_omic(
      .stacking_guard_predictions(integer(0)),
      .stacking_guard_outcome(c(0, 1, 0, 1)),
      task_type = "binary", n_iter = 25L, random_state = 1L
    ),
    silent = TRUE
  )
  expect_identical(before, .Random.seed)
})
