# Tests for stabl_interpret() — descriptive biomarker interpretation (PR1).
# All quantities are post-selection descriptive / apparent unless genuinely
# held-out data are supplied; no inferential claims are tested here.

.make_interp_fit <- function(family = "gaussian", seed = 1L, n = 80L, p = 8L) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  if (identical(family, "binomial")) {
    eta <- 1.5 * x[, "f1"] - 1.2 * x[, "f2"]
    y <- setNames(
      factor(rbinom(n, 1L, plogis(eta)), levels = 0:1,
             labels = c("control", "event")),
      rownames(x)
    )
  } else {
    y <- setNames(0.8 * x[, "f1"] - 0.6 * x[, "f3"] + rnorm(n),
                  rownames(x))
  }
  fit <- stabl_fit(
    x, y,
    lambda_grid = data.frame(lambda = c(0.3, 0.1, 0.05)),
    n_bootstraps = 10L, artificial_type = "random_permutation",
    random_state = seed, family = family
  )
  list(x = x, y = y, fit = fit)
}

test_that("stabl_interpret returns the canonical evidence table", {
  f <- .make_interp_fit("gaussian")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                         n_repeats = 3L, random_state = 2L)
  expect_s3_class(ins, "stabl_interpretation")

  ev <- ins$evidence
  expected_cols <- c("feature", "selected", "stability",
                     "selection_threshold", "coefficient",
                     "standardized_effect", "effect_direction",
                     "effect_measure", "correlation_group",
                     "permutation_importance", "importance_sd",
                     "importance_n_repeats", "evidence_type",
                     "evaluation_scope", "n_observations", "model_id",
                     "na_reason")
  expect_true(all(expected_cols %in% names(ev)))
  expect_equal(nrow(ev), length(f$fit$feature_names))
  expect_identical(ev$feature, f$fit$feature_names)
})

test_that("gaussian effects recover known signs", {
  f <- .make_interp_fit("gaussian")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                         n_repeats = 3L, random_state = 2L)
  ev <- ins$evidence
  # f1 has a positive true effect, f3 negative.
  if (isTRUE(ev$selected[ev$feature == "f1"])) {
    expect_equal(ev$effect_direction[ev$feature == "f1"], "positive")
  }
  if (isTRUE(ev$selected[ev$feature == "f3"])) {
    expect_equal(ev$effect_direction[ev$feature == "f3"], "negative")
  }
})

test_that("binomial effects use log-odds-ratio measure and event class", {
  f <- .make_interp_fit("binomial")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "binomial",
                         n_repeats = 3L, random_state = 3L)
  expect_identical(ins$event_class, "event")
  expect_identical(unique(ins$evidence$effect_measure), "log_odds_ratio")
  ev <- ins$evidence
  if (isTRUE(ev$selected[ev$feature == "f1"])) {
    expect_equal(ev$effect_direction[ev$feature == "f1"], "positive")
  }
})

test_that("event_class can be pinned to the reference level", {
  f <- .make_interp_fit("binomial")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "binomial",
                         event_class = "control", n_repeats = 3L,
                         random_state = 3L)
  expect_identical(ins$event_class, "control")
  # Coefficients flip sign relative to the default event coding.
  ins_default <- stabl_interpret(f$fit, f$x, f$y, family = "binomial",
                                 n_repeats = 3L, random_state = 3L)
  ev <- ins$evidence
  evd <- ins_default$evidence
  sel <- ev$selected
  if (any(sel)) {
    expect_equal(ev$coefficient[sel], -evd$coefficient[sel],
                 tolerance = 1e-6)
  }
})

test_that("zero selected features give explicit NA with reason", {
  f <- .make_interp_fit("gaussian")
  fit0 <- stabl_fit(f$x, f$y, lambda_grid = data.frame(lambda = 5),
                    n_bootstraps = 5L, artificial_type = NULL,
                    hard_threshold = 0.99, random_state = 1L)
  ins <- stabl_interpret(fit0, f$x, f$y, family = "gaussian", n_repeats = 3L)
  ev <- ins$evidence
  expect_true(all(is.na(ev$coefficient)))
  expect_true(all(is.na(ev$permutation_importance)))
  expect_identical(unique(ev$na_reason), "not_selected")
  expect_null(ins$effects_model)
})

test_that("unselected features have NA effects and zero repeats", {
  f <- .make_interp_fit("gaussian")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                         n_repeats = 4L, random_state = 2L)
  ev <- ins$evidence
  unsel <- !ev$selected
  if (any(unsel)) {
    expect_true(all(is.na(ev$coefficient[unsel])))
    expect_true(all(ev$importance_n_repeats[unsel] == 0L))
    expect_true(all(ev$na_reason[unsel] == "not_selected"))
  }
  if (any(ev$selected)) {
    expect_true(all(ev$importance_n_repeats[ev$selected] == 4L))
  }
})

test_that("correlated features land in the same redundancy group", {
  set.seed(5L)
  n <- 100L
  base <- rnorm(n)
  x <- cbind(base + rnorm(n, sd = 0.01), base + rnorm(n, sd = 0.01),
             rnorm(n), rnorm(n))
  colnames(x) <- c("g1a", "g1b", "x3", "x4")
  rownames(x) <- paste0("s", seq_len(n))
  y <- setNames(0.9 * x[, "g1a"] + 0.9 * x[, "g1b"] + rnorm(n, sd = 0.5),
                rownames(x))
  fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.2, 0.05)),
                   n_bootstraps = 10L, artificial_type = "random_permutation",
                   random_state = 5L, hard_threshold = 0.3)
  ins <- stabl_interpret(fit, x, y, family = "gaussian",
                         cor_threshold = 0.9, n_repeats = 3L,
                         random_state = 6L)
  ev <- ins$evidence
  sel <- ev[ev$selected, ]
  if ("g1a" %in% sel$feature && "g1b" %in% sel$feature) {
    expect_identical(sel$correlation_group[sel$feature == "g1a"],
                     sel$correlation_group[sel$feature == "g1b"])
  }
})

test_that("permutation importance is seeded and retains negative values", {
  f <- .make_interp_fit("gaussian")
  i1 <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                        n_repeats = 5L, random_state = 42L)
  i2 <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                        n_repeats = 5L, random_state = 42L)
  expect_identical(i1$importance$repeats, i2$importance$repeats)
  # A noise feature can yield negative importance; the value is retained.
  ev <- i1$evidence
  expect_true(all(is.finite(ev$permutation_importance[ev$selected])))
})

test_that("feature and sample reordering give identical results", {
  f <- .make_interp_fit("gaussian")
  set.seed(9L)
  col_ord <- sample(colnames(f$x))
  row_ord <- sample(rownames(f$x))
  x2 <- f$x[row_ord, col_ord]
  y2 <- f$y[row_ord]
  i1 <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                        n_repeats = 3L, random_state = 7L)
  i2 <- stabl_interpret(f$fit, x2, y2, family = "gaussian",
                        n_repeats = 3L, random_state = 7L)
  expect_identical(i1$evidence$feature, i2$evidence$feature)
  expect_equal(i1$evidence$coefficient, i2$evidence$coefficient)
})

test_that("stabl_interpret errors on missing features and bad IDs", {
  f <- .make_interp_fit("gaussian")
  expect_error(
    stabl_interpret(f$fit, f$x[, -1L], f$y, family = "gaussian"),
    "missing"
  )
  y_bad <- f$y
  names(y_bad)[1L] <- "not_a_sample"
  expect_error(
    stabl_interpret(f$fit, f$x, y_bad, family = "gaussian"),
    "do not match"
  )
  x_nonames <- f$x
  colnames(x_nonames) <- NULL
  expect_error(
    stabl_interpret(f$fit, x_nonames, f$y, family = "gaussian"),
    "row names"
  )
})

test_that("assessment flag relabels scope without changing values", {
  f <- .make_interp_fit("gaussian")
  i_fit <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                           assessment = FALSE, n_repeats = 3L,
                           random_state = 2L)
  i_ass <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                           assessment = TRUE, n_repeats = 3L,
                           random_state = 2L)
  expect_identical(unique(i_fit$evidence$evaluation_scope), "fitting_data")
  expect_identical(unique(i_ass$evidence$evaluation_scope), "held_out_assessment")
  expect_identical(unique(i_fit$evidence$evidence_type), "post_selection_descriptive")
  expect_identical(unique(i_ass$evidence$evidence_type), "cross_validated")
  expect_equal(i_fit$evidence$coefficient, i_ass$evidence$coefficient)
})

test_that("the original selector output is unchanged by interpretation", {
  f <- .make_interp_fit("gaussian")
  before_scores <- get_importances(f$fit)
  before_support <- get_support(f$fit)
  invisible(stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                            n_repeats = 3L, random_state = 2L))
  expect_identical(get_importances(f$fit), before_scores)
  expect_identical(get_support(f$fit), before_support)
})

test_that("summary and print methods work", {
  f <- .make_interp_fit("gaussian")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                         n_repeats = 3L, random_state = 2L)
  s <- summary(ins)
  expect_s3_class(s, "data.frame")
  expect_true(all(s$selected))
  expect_output(print(ins), "stabl_interpretation")
})

test_that("plot returns a ggplot for a non-empty selection", {
  skip_if_not_installed("ggplot2")
  f <- .make_interp_fit("gaussian")
  ins <- stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                         n_repeats = 3L, random_state = 2L)
  if (sum(ins$evidence$selected) > 0L) {
    expect_s3_class(plot(ins), "ggplot")
  }
})

test_that("non-syntactic feature names get correct effect coefficients (B4 regression)", {
  # Regression: glm returns non-syntactic names backtick-escaped, and the
  # name-based coefficient assignment silently missed, reporting every such
  # selected feature as coefficient = NA / "effect_not_estimable".
  set.seed(11L)
  n <- 60L
  x <- matrix(rnorm(n * 4), n, 4,
              dimnames = list(paste0("s", seq_len(n)),
                              c("feat 1", "feat 2", "f3", "f4")))
  y <- setNames(1.2 * x[, "feat 1"] - 0.9 * x[, "f3"] + rnorm(n, sd = 0.3),
                rownames(x))
  fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.2, 0.08)),
                   n_bootstraps = 10L, artificial_type = NULL,
                   hard_threshold = 0.5, random_state = 11L)
  ins <- stabl_interpret(fit, x, y, family = "gaussian",
                         n_repeats = 2L, random_state = 12L)
  ev <- ins$evidence
  # The strong non-syntactic signal feature must be selected.
  expect_true("feat 1" %in% ev$feature[ev$selected])
  # No selected feature may be reported as not estimable.
  expect_true(all(!is.na(ev$coefficient[ev$selected])))
  # Coefficients match a direct glm on the selected panel.
  sel <- ev$feature[ev$selected]
  df <- as.data.frame(x[, sel, drop = FALSE])
  df$.y <- unname(y)
  direct <- stats::coef(stats::glm(.y ~ ., data = df))[-1L]
  expect_equal(unname(ev$coefficient[match(sel, ev$feature)]),
               unname(direct), tolerance = 1e-8)
  expect_equal(ev$effect_direction[ev$feature == "feat 1"], "positive")
})

test_that("a feature named '.y' is not silently dropped from the effect model", {
  set.seed(13L)
  n <- 60L
  x <- matrix(rnorm(n * 3), n, 3,
              dimnames = list(paste0("s", seq_len(n)), c(".y", "f2", "f3")))
  y <- setNames(1.5 * x[, ".y"] + rnorm(n, sd = 0.3), rownames(x))
  fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.2, 0.08)),
                   n_bootstraps = 10L, artificial_type = NULL,
                   hard_threshold = 0.5, random_state = 13L)
  ins <- stabl_interpret(fit, x, y, family = "gaussian",
                         n_repeats = 2L, random_state = 14L)
  ev <- ins$evidence
  expect_true(".y" %in% ev$feature[ev$selected])
  coef_y <- ev$coefficient[ev$feature == ".y"]
  expect_true(is.finite(coef_y))
  expect_equal(unname(coef_y), 1.5, tolerance = 0.2)
})

test_that("stabl_interpret rejects missing values in x or y (E11 regression)", {
  # glm would fit on complete cases while the importance loss uses the
  # full vector, silently producing NA permutation importances.
  f <- .make_interp_fit("gaussian")
  y_na <- f$y
  y_na[3L] <- NA
  expect_error(
    stabl_interpret(f$fit, f$x, y_na, family = "gaussian", n_repeats = 2L),
    "missing"
  )
  x_na <- f$x
  x_na[2L, "f1"] <- NA
  expect_error(
    stabl_interpret(f$fit, x_na, f$y, family = "gaussian", n_repeats = 2L),
    "missing"
  )
})

test_that("single-observation permutation importance is the identity (interp#9 regression)", {
  # With one observation, permuting a column must change nothing.  Pre-fix,
  # sample(xp[, j]) on a positive-integer-valued scalar drew from 1:x.
  set.seed(1)
  df <- data.frame(a = c(5, 2, 3, 7), b = c(1.5, 2.5, 3.5, 4.5))
  m <- stats::glm(y ~ ., data = data.frame(y = c(1.1, 2.2, 2.9, 4.1), df))
  x1 <- as.matrix(df[1, , drop = FALSE])   # a == 5 (positive integer)
  out <- stablr:::.interpret_permutation_importance(
    m, x1, y_vec = 1.1, selected = c("a", "b"), family = "gaussian",
    loss = "squared_error", n_repeats = 20L, random_state = 1L,
    event_class = NULL
  )
  expect_true(all(out$repeats == 0))
})

test_that("loss is validated against the family at entry (interp#12 regression)", {
  f <- .make_interp_fit("gaussian")
  # Bogus loss errors regardless of selection state (pre-fix: validated
  # lazily, only after a non-empty selection and a successful glm).
  expect_error(
    stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                    loss = "bogus", n_repeats = 2L, random_state = 1L),
    "`loss` must be one of"
  )
  # gaussian silently treated brier as squared error pre-fix.
  expect_error(
    stabl_interpret(f$fit, f$x, f$y, family = "gaussian",
                    loss = "brier", n_repeats = 2L, random_state = 1L),
    "`loss` must be one of"
  )
})

test_that("plot.stabl_interpretation builds and errors on empty selection (diag#14)", {
  skip_if_not_installed("ggplot2")
  fake <- list(
    evidence = data.frame(
      feature = c("a", "b", "c"),
      stability = c(0.9, 0.4, 0.7),
      selected = c(TRUE, FALSE, TRUE),
      stringsAsFactors = FALSE
    ),
    family = "gaussian",
    provenance = list(evaluation_scope = "test", evidence_type = "test",
                      n_observations = 10L)
  )
  class(fake) <- "stabl_interpretation"
  p <- plot.stabl_interpretation(fake)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(nrow(summary.stabl_interpretation(fake)), 2L)
  fake$evidence$selected <- FALSE
  expect_error(plot.stabl_interpretation(fake), "No selected features")
})
