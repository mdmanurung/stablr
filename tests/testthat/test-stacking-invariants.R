# Invariants of stacked_multi_omic() that the rest of the suite cannot see.
#
# The suite runs every stacking call with n_iter < 513, so all of it takes the
# single-chunk branch of the candidate loop.  The invariants below decide
# whether a returned `score` is evidence at all: an unscoreable search, a
# no-information metric value competing with genuine candidates, and
# chunk-boundary equivalence.

.stacking_sparse_fixture <- function(n_samples = 4L, complete_rows = 1L) {
  predictions <- matrix(
    NA_real_, nrow = n_samples, ncol = 2L,
    dimnames = list(paste0("s", seq_len(n_samples)), c("rna", "protein"))
  )
  if (complete_rows > 0L) {
    observed <- seq_len(complete_rows)
    predictions[observed, 1L] <- 0.2
    predictions[observed, 2L] <- 0.3
  }
  as.data.frame(predictions)
}

test_that("H1 stacked_multi_omic aborts instead of returning score = -Inf", {
  for (task_type in c("binary", "regression")) {
    for (complete_rows in 0:1) {
      predictions <- .stacking_sparse_fixture(complete_rows = complete_rows)
      outcome <- if (identical(task_type, "binary")) {
        c(0, 1, 0, 1)
      } else {
        c(0.5, 1.5, 0.5, 1.5)
      }
      names(outcome) <- rownames(predictions)

      condition <- tryCatch(
        stacked_multi_omic(
          predictions, outcome, task_type = task_type,
          n_iter = 20L, random_state = 1L
        ),
        condition = function(cnd) cnd
      )

      expect_s3_class(condition, "stablr_stacking_infeasible")
      expect_s3_class(condition, "stablr_error")
      expect_identical(condition$complete_rows, complete_rows)
      expect_identical(condition$n_samples, 4L)
      expect_identical(condition$task_type, task_type)
      expect_match(
        conditionMessage(condition),
        sprintf("%d of 4", complete_rows),
        fixed = TRUE
      )
    }
  }
})

test_that("H1 the guard is about scorability, not prediction validity", {
  predictions <- .stacking_sparse_fixture(complete_rows = 1L)
  scored <- cbind(predictions, lipid = c(NA, 0.1, 0.2, 0.3))
  expect_no_error(
    stacked_multi_omic(
      scored, setNames(c(0, 1, 0, 1), rownames(scored)),
      task_type = "binary", n_iter = 5L, random_state = 2L
    )
  )
})

test_that("H1 two scored rows of one class are infeasible for a binary task", {
  predictions <- .stacking_sparse_fixture(complete_rows = 1L)
  predictions[2L, ] <- c(0.4, 0.6)
  single_class_y <- setNames(c(0, 0, 1, 1), rownames(predictions))
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, single_class_y, task_type = "binary", n_iter = 20L,
      random_state = 3L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
  expect_identical(condition$complete_rows, 2L)
  expect_match(conditionMessage(condition), "class")
})

test_that("H2 .r_auc and .r_auc_batch return NA, not the no-information 0.5", {
  scores <- cbind(
    a = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    b = c(0.9, 0.7, 0.5, 0.3, 0.2, 0.1)
  )

  expect_true(is.na(stablr:::.r_auc(rep(1L, 6L), scores[, 1L])))
  expect_true(is.na(stablr:::.r_auc(rep(0L, 6L), scores[, 1L])))
  expect_identical(
    stablr:::.r_auc_batch(rep(1L, 6L), scores),
    c(NA_real_, NA_real_)
  )
  expect_equal(stablr:::.r_auc(c(0L, 1L), c(0.1, 0.9)), 1)
})

test_that("H2 a no-information AUC can no longer win the argmax", {
  # Both scored rows are y = 0, so no candidate carries event information.
  predictions <- data.frame(
    rna = c(0.1, 0.9, NA, NA),
    protein = c(0.8, 0.1, NA, NA)
  )
  rownames(predictions) <- paste0("s", 1:4)
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, setNames(c(0L, 0L, 1L, 1L), rownames(predictions)),
      task_type = "binary", n_iter = 50L, random_state = 2L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
})

test_that("H3 .r_squared and .r_squared_batch return NA when ss_tot is 0", {
  constant <- rep(2, 5)
  predictions <- cbind(a = rep(2, 5), b = rep(1, 5))

  expect_true(is.na(stablr:::.r_squared(constant, predictions[, 1L])))
  expect_identical(
    stablr:::.r_squared_batch(constant, predictions),
    c(NA_real_, NA_real_)
  )
  expect_equal(stablr:::.r_squared(1:5, rep(3, 5)), 0)
  expect_equal(stablr:::.r_squared(1:5, 1:5), 1)
})

test_that("H3 a constant scored outcome aborts instead of reporting R2 = 0", {
  # The two complete rows are both y = 1, so R^2 is undefined for every
  # candidate; the old floor reported 0 and won the argmax.
  predictions <- data.frame(
    rna = c(NA, 0.2, NA, 0.8),
    protein = c(NA, 0.3, NA, 0.7)
  )
  rownames(predictions) <- paste0("s", 1:4)
  condition <- tryCatch(
    stacked_multi_omic(
      predictions, setNames(c(0, 1, 0, 1), rownames(predictions)),
      task_type = "regression", n_iter = 50L, random_state = 1L
    ),
    condition = function(cnd) cnd
  )
  expect_s3_class(condition, "stablr_stacking_infeasible")
  expect_match(conditionMessage(condition), "constant")
})

test_that("H4 chunked and unchunked searches agree across a chunk boundary", {
  set.seed(4321L)
  n <- 24L
  y <- rbinom(n, 1L, 0.5)
  predictions <- data.frame(
    omic_a = y + rnorm(n, sd = 0.4),
    omic_b = y + rnorm(n, sd = 0.9)
  )
  predictions[c(3L, 11L), "omic_b"] <- NA_real_
  predictions[7L, ] <- NA_real_

  baseline <- stacked_multi_omic(
    predictions, y, task_type = "binary", n_iter = 40L, random_state = 61L
  )
  for (chunk_size in c(1L, 2L, 3L, 7L, 8L, 13L, 39L, 40L, 41L, 512L)) {
    chunked <- stablr:::.stacked_multi_omic_scalar(
      predictions = predictions,
      y = y,
      task_type = "binary",
      n_iter = 40L,
      random_state = 61L,
      chunk_size = chunk_size
    )
    expect_equal(chunked$score, baseline$score, tolerance = 0)
    expect_equal(chunked$weights, baseline$weights, tolerance = 0)
    expect_equal(
      chunked$predictions[["Stacked Gen. Predictions"]],
      baseline$predictions[["Stacked Gen. Predictions"]],
      tolerance = 0
    )
  }
})

test_that("H4 the 512/513 chunk boundary does not perturb the incumbent", {
  set.seed(4322L)
  n <- 20L
  y <- rbinom(n, 1L, 0.5)
  predictions <- data.frame(
    omic_a = y + rnorm(n, sd = 0.4),
    omic_b = y + rnorm(n, sd = 0.8)
  )

  unchunked_513 <- stablr:::.stacked_multi_omic_scalar(
    predictions, y, "binary", n_iter = 513L, random_state = 62L,
    chunk_size = 1024L
  )
  chunked_513 <- stablr:::.stacked_multi_omic_scalar(
    predictions, y, "binary", n_iter = 513L, random_state = 62L
  )
  expect_equal(chunked_513$score, unchunked_513$score, tolerance = 0)
  expect_equal(chunked_513$weights, unchunked_513$weights, tolerance = 0)
  expect_equal(chunked_513$predictions, unchunked_513$predictions, tolerance = 0)

  chunked_512 <- stablr:::.stacked_multi_omic_scalar(
    predictions, y, "binary", n_iter = 512L, random_state = 62L
  )
  # The 513th draw is a fresh candidate, so it can only help.
  expect_gte(chunked_513$score, chunked_512$score)
  expect_true(is.finite(chunked_513$score))
})

test_that("H4 chunk size is bounded by n_iter and must be positive", {
  predictions <- data.frame(a = c(0.1, 0.9, 0.2, 0.8), b = c(0.2, 0.8, 0.3, 0.7))
  y <- c(0L, 1L, 0L, 1L)
  baseline <- stablr:::.stacked_multi_omic_scalar(
    predictions, y, "binary", n_iter = 3L, random_state = 4L
  )
  for (chunk_size in c(0L, -1L, 3L, 100L)) {
    attempt <- stablr:::.stacked_multi_omic_scalar(
      predictions, y, "binary", n_iter = 3L, random_state = 4L,
      chunk_size = chunk_size
    )
    expect_equal(attempt$score, baseline$score, tolerance = 0)
    expect_equal(attempt$weights, baseline$weights, tolerance = 0)
  }
})

test_that("H5 weights report the raw draw and its normalised share", {
  set.seed(99L)
  n <- 30L
  y <- rbinom(n, 1L, 0.5)
  predictions <- data.frame(
    omic_a = y + rnorm(n, sd = 0.4),
    omic_b = y + rnorm(n, sd = 0.8)
  )
  res <- stacked_multi_omic(
    predictions, y, task_type = "binary", n_iter = 500L, random_state = 64L
  )

  expect_named(res$weights, c("Associated_weight", "Normalised_weight"))
  expect_equal(
    res$weights$Normalised_weight,
    res$weights$Associated_weight / sum(res$weights$Associated_weight)
  )
  expect_equal(sum(res$weights$Normalised_weight), 1)
  expect_true(all(res$weights$Associated_weight >= 0))
  expect_true(all(res$weights$Associated_weight <= 10))
})

test_that("H5 the raw and normalised pair reproduces the pinned parity rule", {
  reference <- utils::read.csv(
    file.path(
      testthat::test_path("fixtures", "python_parity"),
      "stacked_multi_omic_reference.csv"
    ),
    stringsAsFactors = FALSE
  )
  raw <- reference$raw_weight
  expect_equal(reference$weight, raw / sum(raw), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(reference$weight, raw)))

  predictions <- utils::read.csv(
    file.path(
      testthat::test_path("fixtures", "python_parity"),
      "stacking_predictions.csv"
    ),
    stringsAsFactors = FALSE
  )
  outcomes <- utils::read.csv(
    file.path(
      testthat::test_path("fixtures", "python_parity"),
      "stacking_outcomes.csv"
    ),
    stringsAsFactors = FALSE
  )
  ids <- predictions$sample_id
  prediction_matrix <- as.matrix(
    predictions[, setdiff(names(predictions), "sample_id"), drop = FALSE]
  )
  rownames(prediction_matrix) <- ids
  y <- outcomes$outcome[match(ids, outcomes$sample_id)]
  names(y) <- ids

  res <- stacked_multi_omic(
    prediction_matrix, y, task_type = "binary", n_iter = 1L,
    random_state = 1L
  )
  # The reported pair must obey the same rule as the pinned fixture, whatever
  # draw the search happened to pick.
  expect_equal(
    res$weights$Normalised_weight,
    res$weights$Associated_weight / sum(res$weights$Associated_weight)
  )
  expect_identical(
    rownames(res$weights), colnames(prediction_matrix)
  )
})

test_that("H5 the multiclass weights report the same two columns", {
  ids <- paste0("s", seq_len(12L))
  y <- factor(rep(c("A", "B", "C"), each = 4L), levels = c("A", "B", "C"))
  names(y) <- ids
  probs <- matrix(0.15, nrow = 12L, ncol = 3L,
                  dimnames = list(ids, c("A", "B", "C")))
  probs[cbind(seq_len(12L), as.integer(y))] <- 0.7
  probs <- probs / rowSums(probs)
  res <- stacked_multi_omic(
    list(a = probs, b = probs), y, task_type = "multiclass",
    n_iter = 40L, random_state = 65L
  )
  expect_named(res$weights, c("Associated_weight", "Normalised_weight"))
  expect_equal(
    res$weights$Normalised_weight,
    res$weights$Associated_weight / sum(res$weights$Associated_weight)
  )
})

test_that("H6 the multiclass log loss is floored at eps", {
  ids <- "s1"
  truth <- factor("A")
  saturated <- matrix(
    c(1, 0, 0), nrow = 1L, dimnames = list(ids, c("A", "B", "C"))
  )
  expect_equal(stablr:::.multiclass_log_loss(truth, saturated), 0, tolerance = 1e-12)

  # A zero probability on the observed class is floored at eps, so the loss is
  # bounded by -log(eps) and the achievable multiclass score by log(eps).  The
  # bound is part of the documented contract, so pin it.
  floored <- matrix(
    c(0, 0.5, 0.5), nrow = 1L, dimnames = list(ids, c("A", "B", "C"))
  )
  expect_equal(stablr:::.multiclass_log_loss(truth, floored, eps = 1e-15),
               -log(1e-15), tolerance = 1e-9)
  expect_lt(stablr:::.multiclass_log_loss(truth, floored, eps = 1e-3), 10)
})
