test_that("adjusted_similarity does not reach 1 for identical partial sets", {
  # The correction denominator is min(k1,k2) - max(0, k1+k2-d), not the set
  # size, so A == B of size k in a universe of d scores 1 - k/d. The roxygen
  # previously claimed the range was (-1, 1] with "1 means perfect
  # agreement", and its worked example asserted a value of 1. Both were wrong
  # and are pinned here so the documentation cannot drift back.
  for (k in c(1L, 2L, 5L)) {
    set <- as.character(seq_len(k))
    expect_equal(
      adjusted_similarity(set, set, nb_total_elements = 10L),
      1 - k / 10
    )
  }
})

test_that("pearson_similarity does reach 1 for identical sets", {
  # The contrast with adjusted_similarity matters: the two chance-corrected
  # measures in the same file have different normalizations, and only
  # pearson_similarity is normalized so identical sets score exactly 1.
  for (k in c(1L, 2L, 5L)) {
    set <- as.character(seq_len(k))
    expect_equal(pearson_similarity(set, set, d = 10L), 1)
  }
})

test_that("the two chance-corrected measures disagree on identical partial sets", {
  set <- c("f1", "f2", "f3")
  expect_lt(
    adjusted_similarity(set, set, nb_total_elements = 10L),
    pearson_similarity(set, set, d = 10L)
  )
})

test_that("the universe size is validated before the chance correction runs", {
  # A non-positive universe previously reached sqrt(pi*(1-pi)) and produced
  # NaN plus a "NaNs produced" warning from pearson_similarity, and a
  # misleading union-cardinal error from adjusted_similarity.
  for (bad in list(0, -1, 2.5, NA_real_, Inf, c(2L, 3L), "10")) {
    expect_error(
      adjusted_similarity(c("f1"), c("f1"), nb_total_elements = bad),
      "positive finite integer"
    )
    expect_error(
      pearson_similarity(c("f1"), c("f1"), d = bad),
      "positive finite integer"
    )
  }
})

test_that("a non-positive universe is rejected without leaking NaN warnings", {
  # expect_error alone would pass even if the validator ran after the sqrt();
  # this asserts no warning condition is signalled on the way out.
  warnings <- character(0)
  withCallingHandlers(
    tryCatch(pearson_similarity(c("f1"), c("f1"), d = 0), error = function(e) NULL),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 0L)
})

test_that("summarizing a single feature set returns NA without NaN warnings", {
  # With one run there are no pairs, so .pairwise_upper_triangle returns
  # numeric(0). mean(numeric(0)) is NaN, so the mean branch used to emit
  # "NaNs produced" on its way to an NA.
  for (stat in c("median", "mean")) {
    single <- list(c("a", "b", "c"))
    for (summary_fn in list(adjusted_similarity_measure, pearson_similarity_measure)) {
      signalled <- character(0)
      result <- withCallingHandlers(
        summary_fn(single, 8L, stat = stat),
        warning = function(w) {
          signalled <<- c(signalled, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      expect_identical(result$statistic, NA_real_)
      expect_identical(result$err, c(NA_real_, NA_real_))
      expect_length(signalled, 0L)
    }
  }
})

test_that("an empty run list summarizes as NA, like a single run", {
  # Zero runs is the same "no pairs to compare" situation as one run, so it
  # must produce the same NA rather than a different failure mode.
  for (summary_fn in list(adjusted_similarity_measure, pearson_similarity_measure)) {
    result <- summary_fn(list(), 8L)
    expect_identical(result$statistic, NA_real_)
    expect_identical(result$err, c(NA_real_, NA_real_))
  }
})

test_that("the summary statistics are unchanged for two or more sets", {
  # Guards the n == 0 short-circuit against altering the ordinary path.
  sets <- list(c("a", "b", "c"), c("b", "c", "d"), c("a", "d"))
  expect_equal(
    adjusted_similarity_measure(sets, 8L, stat = "mean")$statistic,
    mean(c(
      adjusted_similarity(sets[[1]], sets[[2]], 8L),
      adjusted_similarity(sets[[1]], sets[[3]], 8L),
      adjusted_similarity(sets[[2]], sets[[3]], 8L)
    ))
  )
  expect_equal(
    pearson_similarity_measure(sets, 8L, stat = "median")$statistic,
    stats::median(c(
      pearson_similarity(sets[[1]], sets[[2]], 8L),
      pearson_similarity(sets[[1]], sets[[3]], 8L),
      pearson_similarity(sets[[2]], sets[[3]], 8L)
    ))
  )
})

test_that("set metrics treat a repeated identifier as one feature", {
  # All six set metrics call unique() (directly or via the set-based
  # intersect()/setdiff()), so a duplicated identifier must not inflate the
  # true-positive count. This matches the pinned Python reference, whose
  # fdr/tpr/fscore_similarity all compute set(list).intersection(...).
  repeated <- c("a", "a", "b")
  truth <- c("a")

  expect_equal(jaccard_similarity(repeated, truth), 0.5)
  expect_equal(fdr_similarity(repeated, truth), 0.5)
  expect_equal(tpr_similarity(repeated, truth), 1)
  expect_equal(fscore_similarity(repeated, truth, beta = 1), 2 / 3)
})

test_that("adjusted_similarity is the one measure that still counts duplicates", {
  # Documented divergence from the pinned Python reference, pinned here so it
  # is a decision rather than an accident. Python computes the intersection
  # from set(list1) but the union from len(list1), so a repeated identifier
  # inflates its union and not its overlap. R calls unique() up front, so both
  # sides use deduplicated cardinalities. The two therefore disagree only when
  # an input actually repeats an identifier, which upstream validation of
  # feature names forbids, and which no parity fixture exercises.
  expect_equal(adjusted_similarity(c("a", "a", "b"), c("a"), 10L), 0.8)
  expect_equal(
    adjusted_similarity(c("a", "b"), c("a", "b"), 10L),
    0.8,
    info = "deduplicated identical sets still score 1 - k/d"
  )
})

test_that("multiclass MCC matches a hand-computed confusion matrix", {
  # .classification_metrics computes MCC with the multiclass generalization
  # citing Gorodkin 2004, and it previously had no unit test. The expected
  # value is derived here from the same confusion matrix so the test is
  # self-documenting rather than a hard-coded number.
  truth <- c("a", "a", "a", "b", "b", "b", "c", "c", "c")
  pred <- c("a", "a", "b", "b", "b", "b", "c", "c", "a")

  classes <- c("a", "b", "c")
  cm <- table(factor(truth, classes), factor(pred, classes))
  s <- sum(cm)
  tp <- sum(diag(cm))
  actual_positive <- rowSums(cm)
  predicted_positive <- colSums(cm)

  # Gorodkin's generalization: numerator is the trace of the off-covariance
  # matrix, denominator its Frobenius norm.
  numerator <- s * tp - sum(predicted_positive * actual_positive)
  denominator <- sqrt(
    (s^2 - sum(predicted_positive^2)) * (s^2 - sum(actual_positive^2))
  )
  expected <- numerator / denominator

  expect_equal(
    .classification_metrics(truth, pred)$mcc,
    expected
  )
})

test_that("binary MCC matches its closed form", {
  truth <- c(0, 0, 1, 1, 1, 0, 1, 0)
  pred <- c(0, 1, 1, 1, 0, 0, 1, 0)
  tp <- 3; tn <- 3; fp <- 1; fn <- 1
  expected <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

  expect_equal(
    .classification_metrics(truth, pred)$mcc,
    expected
  )
})

test_that("the FDP+ plot is labelled as a proportion, not a rate", {
  skip_if_not_installed("ggplot2")
  set.seed(4242L)
  x <- matrix(rnorm(60 * 5), nrow = 60, ncol = 5,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:5)))
  y <- rnorm(60)
  y <- stats::setNames(y, rownames(x))
  fit <- stabl_fit(x, y, lambda_grid = "auto", n_bootstraps = 5L,
                   n_lambda = 3L, random_state = 11L)

  plot <- plot_fdr_graph(fit)
  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$plot$labels$y, "FDP+ estimate")
  expect_equal(built$plot$labels$title, "FDP+ Estimate")
})
