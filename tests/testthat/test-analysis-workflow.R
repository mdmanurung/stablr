# Tests for stabl_analyze (PR3): explicit exploratory/validated modes,
# fold-safe preprocessing threading, held-out probabilities and importance,
# and the adversarial leakage checks required by Gate V3.

# ── Helpers ──────────────────────────────────────────────────────────────────

.make_binary_data <- function(n = 24L, p = 8L, seed = 42L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  y <- setNames(factor(rep(c("A", "B"), each = n / 2L)), ids)
  x <- matrix(rnorm(n * p), n, p, dimnames = list(ids, paste0("f", 1:p)))
  x[, 1] <- x[, 1] + 1.5 * as.numeric(y == "B")
  x[, 2] <- x[, 2] - 1.2 * as.numeric(y == "B")
  list(x = x, y = y)
}

.validated_call <- function(x, y, ...) {
  args <- list(
    x = x, y = y, family = "binomial", mode = "validated",
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 2L, inner_v = 2L, n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
    random_state = 7L, n_repeats = 3L
  )
  overrides <- list(...)
  args[names(overrides)] <- overrides
  do.call(stabl_analyze, args)
}

# ── Mode and family gating ───────────────────────────────────────────────────

test_that("mode is required and validated explicitly", {
  d <- .make_binary_data()
  expect_error(
    stabl_analyze(d$x, d$y, family = "binomial"),
    "`mode` is required"
  )
  expect_error(
    stabl_analyze(d$x, d$y, family = "binomial", mode = "guess"),
    "should be one of"
  )
})

test_that("validated mode rejects non-binary families and outcomes clearly", {
  d <- .make_binary_data()
  y_num <- setNames(as.numeric(d$y == "B"), names(d$y))
  expect_error(
    stabl_analyze(d$x, y_num, family = "gaussian", mode = "validated"),
    "binary classification only"
  )
  y3 <- setNames(factor(rep(c("A", "B", "C"), each = 8L)), names(d$y))
  expect_error(
    stabl_analyze(d$x, y3, family = "binomial", mode = "validated"),
    "exactly two classes"
  )
  expect_error(
    stabl_analyze(d$x, y3, family = "binomial", mode = "exploratory"),
    "exactly two classes"
  )
})

test_that("a fitted preprocessing object is rejected as leakage", {
  d <- .make_binary_data()
  pp <- stabl_preprocess_fit(d$x, transforms = c("center", "scale"))
  expect_error(
    stabl_analyze(d$x, d$y, family = "binomial", mode = "exploratory",
                  preprocessing = pp),
    "leak"
  )
  expect_error(
    .validated_call(d$x, d$y, preprocessing = pp),
    "leak"
  )
})

# ── Exploratory mode ─────────────────────────────────────────────────────────

test_that("exploratory selection matches an equivalent direct stabl_fit call", {
  d <- .make_binary_data(n = 40L)
  grid <- data.frame(lambda = c(0.2, 0.1))
  direct <- suppressWarnings(stabl_fit(
    d$x, d$y, lambda_grid = grid, family = "binomial",
    n_bootstraps = 5L, artificial_type = "random_permutation",
    random_state = 7L
  ))
  ex <- suppressWarnings(stabl_analyze(
    d$x, d$y, family = "binomial", mode = "exploratory",
    lambda_grid = grid, n_bootstraps = 5L,
    artificial_type = "random_permutation", random_state = 7L, n_repeats = 2L
  ))
  expect_identical(get_support(direct), get_support(ex$fit))
  expect_identical(get_importances(direct), get_importances(ex$fit))
})

test_that("exploratory output is labelled apparent/post-selection", {
  d <- .make_binary_data()
  ex <- suppressWarnings(stabl_analyze(
    d$x, d$y, family = "binomial", mode = "exploratory",
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    n_bootstraps = 3L, artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 7L, n_repeats = 2L
  ))
  expect_s3_class(ex, "stabl_analysis")
  expect_identical(ex$mode, "exploratory")
  expect_identical(ex$evidence_type, "post_selection_descriptive")
  expect_identical(ex$evaluation_scope, "fitting_data")
  expect_s3_class(ex$interpretation, "stabl_interpretation")
  expect_identical(
    unique(ex$interpretation$evidence$evidence_type),
    "post_selection_descriptive"
  )
})

test_that("exploratory mode applies preprocessing before fitting", {
  d <- .make_binary_data()
  d$x[, "f8"] <- 5  # zero-variance feature
  ex <- suppressWarnings(stabl_analyze(
    d$x, d$y, family = "binomial", mode = "exploratory",
    preprocessing = c("zv", "center", "scale"),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    n_bootstraps = 3L, artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 7L, n_repeats = 2L
  ))
  expect_false("f8" %in% ex$preprocessing$retained_features)
  expect_false("f8" %in% ex$fit$feature_names)
})

# ── Validated mode structure ─────────────────────────────────────────────────

test_that("validated mode returns one held-out prediction per sample per repeat", {
  d <- .make_binary_data()
  va <- suppressWarnings(.validated_call(d$x, d$y))
  expect_s3_class(va, "stabl_analysis")
  expect_identical(va$mode, "validated")
  expect_identical(va$evidence_type, "cross_validated")
  expect_identical(va$evaluation_scope, "held_out_assessment")

  preds <- va$predictions
  expect_true(all(c("repeat_id", "fold", "fold_id", "sample_id", "truth",
                    "predicted", "event_probability") %in% names(preds)))
  # Exactly one prediction per participant per repeat.
  expect_equal(
    anyDuplicated(paste(preds$repeat_id, preds$sample_id)),
    0L
  )
  expect_equal(sort(unique(preds$sample_id)), sort(names(d$y)))
  expect_true(all(preds$event_probability >= 0 & preds$event_probability <= 1))
  expect_identical(va$event_class, "B")
  expect_true(all(preds$truth %in% c("A", "B")))
  expect_true(all(preds$predicted %in% c("A", "B")))

  pooled <- va$performance$pooled
  expect_true(all(c("accuracy", "balanced_error_rate", "macro_f1", "mcc",
                    "brier", "log_loss", "auroc") %in% names(pooled)))
  expect_equal(nrow(va$performance$by_repeat), 1L)
})

test_that("event-class mapping is explicit and reversible", {
  d <- .make_binary_data(n = 40L)
  va_b <- suppressWarnings(.validated_call(d$x, d$y, event_class = "B"))
  va_a <- suppressWarnings(.validated_call(d$x, d$y, event_class = "A"))
  expect_equal(va_b$predictions$event_probability,
               1 - va_a$predictions$event_probability)
  # f1 is higher in class B by construction: high f1 rows should have high
  # P(B) on average.
  preds <- va_b$predictions
  f1_vals <- d$x[preds$sample_id, "f1"]
  expect_gt(stats::cor(f1_vals, preds$event_probability), 0)
})

test_that("validated mode is deterministic for a fixed seed", {
  d <- .make_binary_data()
  va1 <- suppressWarnings(.validated_call(d$x, d$y))
  va2 <- suppressWarnings(.validated_call(d$x, d$y))
  expect_identical(va1$predictions, va2$predictions)
  expect_identical(va1$selected_features, va2$selected_features)
})

# ── Adversarial leakage checks (Gate V3) ─────────────────────────────────────

test_that("permuting assessment outcomes leaves fold training untouched", {
  d <- .make_binary_data()
  va <- suppressWarnings(.validated_call(d$x, d$y, stratified = FALSE))

  fold1_valid <- va$nested$outer_folds[[1L]]$valid_ids
  y2 <- d$y
  y2[fold1_valid] <- ifelse(y2[fold1_valid] == "A", "B", "A")
  va2 <- suppressWarnings(.validated_call(d$x, y2, stratified = FALSE))

  # Fold membership is fixed (unstratified folds do not depend on y).
  expect_identical(
    va$nested$outer_folds[[1L]]$valid_ids,
    va2$nested$outer_folds[[1L]]$valid_ids
  )
  # Fold-1 training-side artifacts are unchanged: selection, importances,
  # and the predictions/probabilities for the held-out rows.
  p1 <- va$predictions[va$predictions$fold_id == "Repeat1_Fold1", ]
  p2 <- va2$predictions[va2$predictions$fold_id == "Repeat1_Fold1", ]
  expect_identical(p1$predicted, p2$predicted)
  expect_identical(p1$event_probability, p2$event_probability)
  s1 <- va$selected_features[va$selected_features$fold_id == "Repeat1_Fold1", ]
  s2 <- va2$selected_features[va2$selected_features$fold_id == "Repeat1_Fold1", ]
  expect_identical(s1, s2)
  # The altered truth is only visible at metric evaluation time.
  expect_identical(p1$sample_id, p2$sample_id)
  expect_false(identical(
    factor(p1$truth, levels = c("A", "B")),
    factor(p2$truth, levels = c("A", "B"))
  ))
})

test_that("extreme assessment values do not change training-fitted transforms", {
  d <- .make_binary_data()
  va <- suppressWarnings(.validated_call(
    d$x, d$y, preprocessing = c("center", "scale"), stratified = FALSE
  ))

  fold1_valid <- va$nested$outer_folds[[1L]]$valid_ids
  x2 <- d$x
  x2[fold1_valid, ] <- x2[fold1_valid, ] * 1000
  va2 <- suppressWarnings(.validated_call(
    x2, d$y, preprocessing = c("center", "scale"), stratified = FALSE
  ))

  pp1 <- va$nested$fold_results[[1L]]$fit$preprocessing
  pp2 <- va2$nested$fold_results[[1L]]$fit$preprocessing
  expect_identical(pp1$params$center$means, pp2$params$center$means)
  expect_identical(pp1$params$scale$sds, pp2$params$scale$sds)
  s1 <- va$selected_features[va$selected_features$fold_id == "Repeat1_Fold1", ]
  s2 <- va2$selected_features[va2$selected_features$fold_id == "Repeat1_Fold1", ]
  expect_identical(s1, s2)
})

test_that("zero-variance features are removed inside fold training only", {
  d <- .make_binary_data()
  d$x[, "f8"] <- 3.14  # constant: must never be selected
  va <- suppressWarnings(.validated_call(
    d$x, d$y, preprocessing = c("zv", "center", "scale")
  ))
  expect_false("f8" %in% va$selected_features$feature)
  for (fr in va$nested$fold_results) {
    expect_false("x::f8" %in% fr$fit$preprocessing$retained_features)
  }
})

# ── Held-out importance ──────────────────────────────────────────────────────

test_that("held-out importance uses fold predictors and training groups", {
  d <- .make_binary_data(n = 40L)
  d$x[, 1] <- d$x[, 1] + as.numeric(d$y == "B")  # total shift 2.5
  va <- suppressWarnings(.validated_call(d$x, d$y, n_bootstraps = 5L))
  imp <- va$heldout_importance
  expect_true(all(c("feature", "mean_importance", "sd_importance",
                    "n_folds_selected", "n_repeats", "loss", "event_class",
                    "evidence_type", "evaluation_scope", "na_reason") %in%
                    names(imp$summary)))
  expect_identical(unique(imp$summary$evidence_type), "cross_validated")
  expect_identical(unique(imp$summary$evaluation_scope), "held_out_assessment")
  expect_identical(imp$loss, "brier")
  # Directional check: the two true signals carry more held-out importance
  # than the noise features on average.
  sel <- imp$summary[imp$summary$n_folds_selected > 0L, ]
  expect_true(all(c("f1", "f2") %in% sel$feature))
  is_signal <- sel$feature %in% c("f1", "f2")
  if (any(!is_signal)) {
    expect_gt(mean(sel$mean_importance[is_signal]),
              mean(sel$mean_importance[!is_signal], na.rm = TRUE))
  } else {
    # Only the true signals were ever selected: the purest possible outcome.
    succeed()
  }
  # Fold-level table and training-derived groups exist per fold.
  expect_true(all(c("repeat_id", "fold_id", "feature_set", "set_type",
                    "importance_mean", "importance_sd") %in% names(imp$folds)))
  expect_true(length(imp$groups_by_fold) >= 1L)
  # Importance is computed with the fold's own predictor: the retained
  # predictor objects exist for every fold.
  preds <- lapply(va$nested$fold_results, function(fr) fr$fit$predictor)
  expect_true(all(vapply(preds, function(p) !is.null(p), logical(1L))))
})

test_that("held-out importance supports log_loss and group sets", {
  d <- .make_binary_data()
  # Correlated substitute for f1.
  d$x[, "f6"] <- d$x[, "f1"] + rnorm(nrow(d$x), sd = 0.05)
  va <- suppressWarnings(.validated_call(d$x, d$y, loss = "log_loss",
                                         cor_threshold = 0.5))
  expect_identical(va$heldout_importance$loss, "log_loss")
  folds <- va$heldout_importance$folds
  if (any(folds$set_type == "group")) {
    grp <- folds[folds$set_type == "group", ]
    expect_true(all(startsWith(grp$feature_set, "group:")))
    expect_true(all(is.finite(grp$importance_mean)))
  }
  expect_true(all(is.finite(folds$importance_mean)))
})

# ── Final refit ──────────────────────────────────────────────────────────────

test_that("final_refit is explicit, labelled, and separate from CV evidence", {
  d <- .make_binary_data()
  va <- suppressWarnings(.validated_call(d$x, d$y, final_refit = TRUE))
  expect_false(is.null(va$final_refit))
  expect_identical(va$final_refit$evidence_type, "final_refit")
  expect_identical(
    unique(va$final_refit$interpretation$evidence$evidence_type),
    "final_refit"
  )
  expect_s3_class(va$final_refit$fit, "stabl_fit")
  # CV evidence labels are untouched by the refit.
  expect_identical(va$evidence_type, "cross_validated")
  expect_match(va$final_refit$note, "different model")
})

# ── Study-disjoint validated mode ────────────────────────────────────────────

test_that("validated mode respects study-disjoint isolation", {
  d <- .make_binary_data(n = 24L)
  # Interleave study membership so every study is class-mixed (3 A / 3 B).
  studies <- setNames(rep(paste0("study", 1:4), times = 6L), names(d$y))
  va <- suppressWarnings(.validated_call(
    d$x, d$y, study_labels = studies, outer_v = 4L, inner_v = 2L
  ))
  expect_identical(va$provenance$fold_scheme, "study_disjoint")
  expect_false(is.null(va$performance_by_study))
  expect_equal(nrow(va$performance_by_study), 4L)
  for (fold in va$nested$outer_folds) {
    held_out_studies <- fold$study_ids
    train_studies <- unique(as.character(studies[fold$train_ids]))
    expect_length(intersect(held_out_studies, train_studies), 0L)
  }
})

# ── Binary probability metrics ───────────────────────────────────────────────

test_that("binary probability metrics are exact on known inputs", {
  truth <- factor(c("A", "A", "B", "B"))
  p <- c(0.1, 0.2, 0.8, 0.9)  # P(B)
  m <- .binary_probability_metrics(truth, p, "B")
  expect_equal(m$brier, mean((c(0, 0, 1, 1) - p)^2))
  expect_equal(m$auroc, 1)
  expect_equal(m$n, 4L)
  expect_equal(m$n_events, 2L)

  m_rev <- .binary_probability_metrics(truth, 1 - p, "B")
  expect_equal(m_rev$auroc, 0)

  m_half <- .binary_probability_metrics(truth, rep(0.5, 4L), "B")
  expect_equal(m_half$log_loss, log(2))

  # Single-class truth: AUROC is an explicit diagnostic NA, not an error.
  one_class <- factor(c("A", "A"), levels = c("A", "B"))
  m_na <- .binary_probability_metrics(one_class, c(0.3, 0.4), "B")
  expect_true(is.na(m_na$auroc))
  expect_true(is.finite(m_na$brier))

  # Event-class mapping is validated.
  expect_error(
    .binary_probability_metrics(truth, p, "C"),
    "event_class"
  )
  expect_error(
    .binary_probability_metrics(truth, c(0.2, 1.2, 0.3, 0.4), "B"),
    "\\[0, 1\\]"
  )
})

# ── Selected-feature predictor units ─────────────────────────────────────────

test_that("zero-feature binary fallback is the training event prevalence", {
  y_train <- factor(c("A", "A", "B", "B", "B"))
  out <- .fit_selected_predictor(
    x_train = matrix(numeric(), 5L, 0L),
    y_train = y_train,
    x_valid = matrix(numeric(), 3L, 0L),
    levels = c("A", "B"),
    event_class = "B"
  )
  expect_identical(out$predicted, rep("B", 3L))
  expect_equal(out$probability, rep(0.6, 3L))
  expect_identical(out$fallback, "training_prevalence")
  expect_equal(out$predict_event_prob(matrix(numeric(), 2L, 0L)), rep(0.6, 2L))
})

test_that("one-feature binary predictor matches a direct logit glm", {
  set.seed(1L)
  xtr <- rnorm(30L)
  y_train <- factor(ifelse(xtr + rnorm(30L, sd = 0.5) > 0, "B", "A"))
  xva <- rnorm(7L)
  out <- .fit_selected_predictor(
    x_train = matrix(xtr, ncol = 1L),
    y_train = y_train,
    x_valid = matrix(xva, ncol = 1L),
    levels = c("A", "B"),
    event_class = "B"
  )
  direct <- suppressWarnings(glm(y_train ~ xtr, family = binomial))
  expect_equal(
    out$probability,
    unname(predict(direct, newdata = data.frame(xtr = xva), type = "response"))
  )
  # Event-class flip maps probabilities to 1 - p.
  out_a <- .fit_selected_predictor(
    x_train = matrix(xtr, ncol = 1L),
    y_train = y_train,
    x_valid = matrix(xva, ncol = 1L),
    levels = c("A", "B"),
    event_class = "A"
  )
  expect_equal(out_a$probability, 1 - out$probability)
  # Hard classes are unchanged by the event-class choice.
  expect_identical(out_a$predicted, out$predicted)
})

test_that("multi-feature binary probabilities match cv.glmnet event mapping", {
  skip_if_not_installed("glmnet")
  set.seed(3L)
  xtr <- matrix(rnorm(60 * 3), 60L, 3L)
  y_train <- factor(ifelse(xtr[, 1L] + rnorm(60L, sd = 0.3) > 0, "B", "A"))
  xva <- matrix(rnorm(10 * 3), 10L, 3L)
  set.seed(9L)
  out <- .fit_selected_predictor(
    x_train = xtr, y_train = y_train, x_valid = xva,
    levels = c("A", "B"), event_class = "B"
  )
  set.seed(9L)
  fit <- glmnet::cv.glmnet(xtr, y_train, family = "binomial",
                           type.measure = "class", nfolds = 5L)
  p_direct <- as.numeric(predict(fit, newx = xva, s = "lambda.min",
                                 type = "response")[, 1L])
  expect_equal(out$probability, p_direct)
  # glmnet models P(second level) for a two-level factor: P(A) = 1 - P(B).
  set.seed(9L)
  out_a <- .fit_selected_predictor(
    x_train = xtr, y_train = y_train, x_valid = xva,
    levels = c("A", "B"), event_class = "A"
  )
  expect_equal(out_a$probability, 1 - p_direct)
})

test_that("validated mode tolerates folds with empty selection (IV regression)", {
  # Regression for the IV-benchmark failure (2026-10-15): on null data some
  # outer folds select nothing, and the empty fold-importance branch built a
  # malformed zero-row data.frame ("arguments imply differing number of
  # rows: 1, 0").  These exact settings produced empty-selection folds on
  # seeds 6, 9 and 10 of the IV-0 design.
  set.seed(6L)
  n <- 120L
  p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", 1:p)))
  pr <- rep(0.5, n)
  y <- setNames(
    factor(ifelse(stats::rbinom(n, 1L, pr) == 1L, "B", "A"),
           levels = c("A", "B")),
    rownames(x)
  )
  va <- suppressWarnings(stabl_analyze(
    x, y, family = "binomial", mode = "validated",
    preprocessing = c("center", "scale"),
    lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
    outer_v = 3L, inner_v = 3L, n_bootstraps = 50L,
    artificial_type = "random_permutation", random_state = 6L,
    n_repeats = 10L
  ))
  # The empty-selection branch was exercised: at least one of the 3 outer
  # folds selected nothing.
  expect_lt(length(unique(va$selected_features$fold_id)), 3L)
  # Importance summary stays well-formed: one row per feature, features
  # never selected in any fold are labelled not_selected.
  s <- va$heldout_importance$summary
  expect_equal(nrow(s), p)
  expect_true(all(paste0("f", 1:p) %in% s$feature))
  expect_true(any(s$n_folds_selected == 0L))
  expect_true(all(s$na_reason[s$n_folds_selected == 0L] == "not_selected"))
})
