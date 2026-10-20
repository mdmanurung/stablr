# Tests for stabl_diagnose() (1.0.0 spec, A1 — selection reliability
# diagnostics). Hand-worked score fixtures with exact expected counts and
# FDP+ values; parity preservation against compute_fdp_plus() and the
# byte-pinned Python-parity fixtures; null-saturation signatures; edge
# cases (empty support, ties, no artificials, zero/saturated scores);
# all response families; the diagnostic never changes the fit or support.

# ---- helpers ----------------------------------------------------------------

# Run compute_fdp_plus() on a fixture, then diagnose with its outputs —
# mirroring how stabl_fit() populates the object.
diagnose_fixture <- function(real, art, prop, thresholds,
                             hard_threshold = NULL,
                             learner = "lasso", family = "gaussian",
                             decoy = "knockoff_equi",
                             low_decoy_count = 0L,
                             with_fdp = TRUE) {
  fdp <- if (with_fdp) {
    compute_fdp_plus(real, art, prop, fdr_threshold_range = thresholds)
  } else {
    list(FDRs = NULL, min_fdr = NULL, fdr_min_threshold = NULL)
  }
  .diagnose_stabl_scores(
    stabl_scores            = real,
    stabl_scores_artificial = art,
    fdr_threshold_range     = thresholds,
    FDRs                    = fdp$FDRs,
    min_fdr                 = fdp$min_fdr,
    fdr_min_threshold       = fdp$fdr_min_threshold,
    hard_threshold          = hard_threshold,
    artificial_type         = decoy,
    artificial_proportion   = prop,
    effective_artificial_proportion = prop,
    base_learner            = learner,
    family                  = family,
    n_bootstraps            = 100L,
    feature_names           = rownames(real),
    low_decoy_count         = low_decoy_count
  )
}

# Tiny real fit used by several blocks (fast: 10 bootstraps, 3 lambdas).
make_tiny_fit <- function(family = "gaussian") {
  set.seed(1)
  n <- 60L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- if (family == "binomial") {
    setNames(rbinom(n, 1L, 0.5), rownames(x))
  } else {
    setNames(x[, 1L] * 2 + rnorm(n), rownames(x))
  }
  stabl_fit(
    x, y,
    lambda_grid     = data.frame(lambda = c(0.2, 0.1, 0.05)),
    family          = family,
    n_bootstraps    = 10L,
    artificial_type = "random_permutation",
    random_state    = 1L
  )
}

# ---- 1. Hand-worked fixture: unique interior minimum (exact values) ---------

test_that("hand-worked fixture A: exact threshold path and FDP+ values", {
  # real max scores 0.9, 0.8, 0.7, 0.6, 0.4, 0.2 ; art max scores 0.45, 0.30
  # prop 1/3 ; thresholds 0, 0.25, 0.5, 0.75
  # strict-> counts: real 6,5,4,2 ; art 2,2,0,0
  # FDP+ = (3*n_art + 1)/pmax(1, n_real): 7/6, 7/5, 1/4, 1/2
  real <- cbind(c(0.9, 0.8, 0.7, 0.6, 0.4, 0.2), rep(0.1, 6L))
  art  <- cbind(c(0.45, 0.30), rep(0.1, 2L))
  d <- diagnose_fixture(real, art, 1/3, c(0, 0.25, 0.5, 0.75))

  expect_s3_class(d, "stabl_diagnostics")
  expect_equal(d$threshold_path$n_real, c(6, 5, 4, 2))
  expect_equal(d$threshold_path$n_artificial, c(2, 2, 0, 0))
  expect_equal(d$threshold_path$fdp_plus, c(7/6, 1.4, 0.25, 0.5))
  expect_equal(d$threshold_path$prop_real, c(1, 5/6, 4/6, 2/6))
  expect_equal(d$threshold_path$prop_artificial, c(1, 1, 0, 0))

  expect_equal(d$fdp_curve$min_fdr, 0.25)
  expect_equal(d$fdp_curve$minima_thresholds, 0.5)
  expect_false(d$fdp_curve$is_tied)
  expect_false(d$fdp_curve$is_flat)
  expect_false(d$fdp_curve$is_boundary)

  expect_identical(d$cutoff$selection_rule, "fdp_plus_minimum")
  expect_equal(d$cutoff$selected_threshold, 0.5)
  expect_identical(d$cutoff$comparator, "strictly greater (>)")
  expect_equal(d$cutoff$denominator$artificial_proportion, 1/3)

  # per-lambda view at the selected threshold
  expect_equal(d$per_lambda$n_real, c(4, 0))
  expect_equal(d$per_lambda$n_artificial, c(0, 0))

  # clean fixture: no warnings (evidence rows for lasso x knockoff_equi x
  # gaussian are mixed, so no blanket evidence code)
  expect_identical(d$warning_codes, character(0))
  expect_false(d$support_empty)
  expect_false(d$null_saturation$detected)
})

# ---- 2. Ties, flat minima, boundary solutions --------------------------------

test_that("fixture B: tied minima (non-contiguous) flagged, first-index wins", {
  # FDP+: 1.25, 1.0, 1.5, 1.0 -> min at thresholds 0.25 and 0.75
  real <- cbind(c(0.9, 0.7, 0.4, 0.1), rep(0.05, 4L))
  art  <- cbind(c(0.6, 0.2), rep(0.05, 2L))
  d <- diagnose_fixture(real, art, 0.5, c(0, 0.25, 0.5, 0.75))

  expect_equal(d$fdp_curve$min_fdr, 1.0)
  expect_equal(d$fdp_curve$minima_thresholds, c(0.25, 0.75))
  expect_true(d$fdp_curve$is_tied)
  expect_false(d$fdp_curve$is_flat)
  expect_false(d$fdp_curve$is_boundary)
  expect_equal(d$cutoff$selected_threshold, 0.25)  # first minimum
  expect_true("STABL_DIAG_TIED_MINIMA" %in% d$warning_codes)
  expect_false("STABL_DIAG_FLAT_MINIMUM" %in% d$warning_codes)
  expect_false("STABL_DIAG_BOUNDARY_MINIMUM" %in% d$warning_codes)
})

test_that("fixture C: contiguous plateau flagged as flat", {
  # FDP+: 4/3, 1.25, 1.25, 1.25 -> flat minimum over 0.1..0.5
  real <- cbind(c(0.9, 0.6, 0.05), rep(0.01, 3L))
  art  <- cbind(c(0.55, 0.04), rep(0.01, 2L))
  d <- diagnose_fixture(real, art, 2/3, c(0, 0.1, 0.2, 0.5))

  expect_equal(d$fdp_curve$minima_thresholds, c(0.1, 0.2, 0.5))
  expect_true(d$fdp_curve$is_tied)
  expect_true(d$fdp_curve$is_flat)
  expect_false(d$fdp_curve$is_boundary)
  expect_true("STABL_DIAG_TIED_MINIMA" %in% d$warning_codes)
  expect_true("STABL_DIAG_FLAT_MINIMUM" %in% d$warning_codes)
})

test_that("fixture D: minimum at the range edge flagged as boundary", {
  # FDP+: 4/3, 1.25, 1.0 -> unique minimum at the last threshold
  real <- cbind(c(0.4, 0.3, 0.2), rep(0.01, 3L))
  art  <- cbind(c(0.35, 0.25), rep(0.01, 2L))
  d <- diagnose_fixture(real, art, 2/3, c(0, 0.25, 0.5))

  expect_equal(d$fdp_curve$minima_thresholds, 0.5)
  expect_true(d$fdp_curve$is_boundary)
  expect_true("STABL_DIAG_BOUNDARY_MINIMUM" %in% d$warning_codes)
})

# ---- 3. Null-saturation probe -------------------------------------------------

test_that("fixture E: decoy-leaked scores trigger null saturation", {
  # Importance-leakage signature: decoys score like reals, so the FDP+
  # minimum (0.875 at t = 0.5) sits in a region selecting most of both.
  # FDP+: t0 (2*5+1)/10 = 1.1 ; t0.25 (2*4+1)/8 = 1.125 ;
  #       t0.5 (2*3+1)/8 = 0.875 ; t0.75 (2*1+1)/3 = 1.0
  real <- cbind(c(0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.2, 0.1),
                rep(0.05, 10L))
  art  <- cbind(c(0.8, 0.7, 0.6, 0.5, 0.3), rep(0.05, 5L))
  d <- diagnose_fixture(real, art, 0.5, c(0, 0.25, 0.5, 0.75),
                        learner = "rf")

  expect_equal(d$cutoff$selected_threshold, 0.5)
  expect_true(d$null_saturation$detected)
  expect_equal(d$null_saturation$at_selected$prop_real, 0.8)
  expect_equal(d$null_saturation$at_selected$prop_artificial, 0.6)
  expect_equal(d$null_saturation$near_zero_region$thresholds, 0)
  expect_true("STABL_DIAG_NULL_SATURATION" %in% d$warning_codes)
  expect_output(print(d), "NULL SATURATION")
})

test_that("fixture E2: user hard_threshold in a saturated region is probed off-grid", {
  real <- matrix(1, 3L, 2L)
  art  <- matrix(1, 2L, 2L)
  d <- diagnose_fixture(real, art, 2/3, c(0, 0.5), hard_threshold = 0.01)
  expect_identical(d$cutoff$selection_rule, "user_hard_threshold")
  expect_true(d$null_saturation$detected)
  expect_equal(d$null_saturation$at_selected$prop_real, 1)
  expect_equal(d$null_saturation$at_selected$prop_artificial, 1)
  expect_true("STABL_DIAG_NULL_SATURATION" %in% d$warning_codes)
  expect_true("STABL_DIAG_USER_THRESHOLD" %in% d$warning_codes)
})

test_that("healthy fixture does not trigger null saturation", {
  real <- cbind(c(0.9, 0.8, 0.7, 0.6, 0.4, 0.2), rep(0.1, 6L))
  art  <- cbind(c(0.45, 0.30), rep(0.1, 2L))
  d <- diagnose_fixture(real, art, 1/3, c(0, 0.25, 0.5, 0.75))
  expect_false(d$null_saturation$detected)
  expect_false("STABL_DIAG_NULL_SATURATION" %in% d$warning_codes)
})

# ---- 4. User threshold vs FDP+ rule; empty support ---------------------------

test_that("fixture F: user hard_threshold distinguished; empty support flagged", {
  real <- cbind(c(0.1, 0.2, 0.15), rep(0.05, 3L))
  art  <- cbind(c(0.1, 0.12), rep(0.05, 2L))
  d <- diagnose_fixture(real, art, 2/3, c(0, 0.25, 0.5),
                        hard_threshold = 0.5, decoy = "random_permutation",
                        with_fdp = FALSE)

  expect_identical(d$cutoff$selection_rule, "user_hard_threshold")
  expect_equal(d$cutoff$selected_threshold, 0.5)
  expect_null(d$fdp_curve)
  expect_true(d$support_empty)
  expect_true("STABL_DIAG_USER_THRESHOLD" %in% d$warning_codes)
  expect_true("STABL_DIAG_EMPTY_SUPPORT" %in% d$warning_codes)
  expect_true(any(grepl("hard_threshold", d$assumption_notes)))
})

# ---- 5. No artificials --------------------------------------------------------

test_that("fixture G: no-artificial fit degrades gracefully", {
  real <- cbind(c(0.9, 0.2, 0.1), rep(0.05, 3L))
  d <- .diagnose_stabl_scores(
    stabl_scores = real,
    stabl_scores_artificial = NULL,
    fdr_threshold_range = c(0, 0.25, 0.5),
    hard_threshold = 0.5,
    base_learner = "lasso", family = "gaussian",
    n_bootstraps = 100L, low_decoy_count = 0L
  )
  expect_true(all(is.na(d$threshold_path$n_artificial)))
  expect_true(all(is.na(d$threshold_path$fdp_plus)))
  expect_null(d$fdp_curve)
  expect_null(d$decoy_provenance)
  expect_true("STABL_DIAG_NO_ARTIFICIALS" %in% d$warning_codes)
  expect_true("STABL_DIAG_USER_THRESHOLD" %in% d$warning_codes)
  expect_true(any(grepl("without artificial features", d$calibration_scope)))
  expect_false(d$support_empty)  # 0.9 > 0.5

  d2 <- .diagnose_stabl_scores(
    stabl_scores = real,
    stabl_scores_artificial = NULL,
    fdr_threshold_range = c(0, 0.25, 0.5),
    hard_threshold = NULL,
    base_learner = "lasso", family = "gaussian",
    n_bootstraps = 100L, low_decoy_count = 0L
  )
  expect_identical(d2$cutoff$selection_rule, "none")
  expect_true("STABL_DIAG_NO_THRESHOLD" %in% d2$warning_codes)
  expect_true(d2$support_empty)
})

# ---- 6. Zero and saturated importances ----------------------------------------

test_that("fixture H: all-zero scores run clean with structural flags", {
  real <- matrix(0, 3L, 2L)
  art  <- matrix(0, 2L, 2L)
  d <- diagnose_fixture(real, art, 2/3, c(0, 0.5))
  expect_equal(d$fdp_curve$fdp_plus, c(1, 1))
  expect_true(d$fdp_curve$is_tied)
  expect_true(d$fdp_curve$is_flat)
  expect_true(d$fdp_curve$is_boundary)
  expect_true(d$support_empty)
  expect_true("STABL_DIAG_EMPTY_SUPPORT" %in% d$warning_codes)
})

test_that("fixture I: all-one scores are capped out (refusal), path shows saturation", {
  # FDP+ = 4/3 > 1 everywhere -> fdr_min_threshold capped at 1 -> strict >
  # selects nothing: the machinery refuses a fully saturated path.
  real <- matrix(1, 3L, 2L)
  art  <- matrix(1, 2L, 2L)
  d <- diagnose_fixture(real, art, 2/3, c(0, 0.5))
  expect_equal(d$cutoff$fdr_min_threshold, 1)
  expect_true(d$support_empty)
  expect_true("STABL_DIAG_EMPTY_SUPPORT" %in% d$warning_codes)
  expect_false(d$null_saturation$detected)
  # the saturated near-zero region is still reported as data
  expect_equal(d$null_saturation$near_zero_region$thresholds, 0)
})

# ---- 7. Parity preservation ----------------------------------------------------

test_that("diagnostics reproduce compute_fdp_plus outputs exactly", {
  real <- cbind(c(0.9, 0.8, 0.7, 0.6, 0.4, 0.2), rep(0.1, 6L))
  art  <- cbind(c(0.45, 0.30), rep(0.1, 2L))
  thresholds <- c(0, 0.25, 0.5, 0.75)
  fdp <- compute_fdp_plus(real, art, 1/3, fdr_threshold_range = thresholds)
  d <- diagnose_fixture(real, art, 1/3, thresholds)
  expect_identical(d$fdp_curve$fdp_plus, fdp$FDRs)
  expect_identical(d$fdp_curve$min_fdr, fdp$min_fdr)
  expect_identical(d$cutoff$fdr_min_threshold, fdp$fdr_min_threshold)
})

test_that("legacy Python-parity FDP+ inputs keep their pinned behavior", {
  fixture_dir <- test_path("fixtures", "python_parity")
  skip_if_not(dir.exists(fixture_dir),
              "canonical python-parity fixtures not available")
  scores <- utils::read.csv(file.path(fixture_dir, "stability_scores.csv"),
                            check.names = FALSE, stringsAsFactors = FALSE)
  curve  <- utils::read.csv(file.path(fixture_dir, "fdp_plus_curve.csv"),
                            check.names = FALSE, stringsAsFactors = FALSE)
  score_matrix <- function(kind) {
    sub <- scores[scores$feature_type == kind, , drop = FALSE]
    out <- matrix(NA_real_,
                  length(unique(sub$feature)), length(unique(sub$lambda_id)),
                  dimnames = list(unique(sub$feature), unique(sub$lambda_id)))
    for (i in seq_len(nrow(sub))) {
      out[sub$feature[i], sub$lambda_id[i]] <- sub$stability_score[i]
    }
    out
  }
  real <- score_matrix("real")
  art  <- score_matrix("artificial")
  fdp <- compute_fdp_plus(real, art, 0.5, fdr_threshold_range = curve$threshold)
  d <- diagnose_fixture(real, art, 0.5, curve$threshold)

  # exact replay of the pinned Python curve and first-minimum policy
  expect_equal(d$fdp_curve$fdp_plus, curve$fdp_plus, tolerance = 1e-12)
  expect_equal(d$cutoff$fdr_min_threshold, 0.25, tolerance = 1e-12)
  expect_equal(d$fdp_curve$minima_thresholds[[1L]], 0.25, tolerance = 1e-12)
  expect_identical(d$fdp_curve$fdp_plus, fdp$FDRs)
})

# ---- 8. Evidence profiles ------------------------------------------------------

test_that("evidence profile statuses follow the benchmark-backed table", {
  real <- cbind(c(0.9, 0.8, 0.7, 0.6, 0.4, 0.2), rep(0.1, 6L))
  art  <- cbind(c(0.45, 0.30), rep(0.1, 2L))
  thrs <- c(0, 0.25, 0.5, 0.75)

  d_ok <- diagnose_fixture(real, art, 1/3, thrs,
                           learner = "lasso", decoy = "knockoff_equi",
                           family = "gaussian")
  expect_true("validated_in_stablr" %in% d_ok$evidence_profile$status)
  expect_false("STABL_DIAG_EVIDENCE_CONCERNING" %in% d_ok$warning_codes)
  expect_false("STABL_DIAG_EVIDENCE_NOT_EVALUATED" %in% d_ok$warning_codes)
  expect_true(all(!is.na(d_ok$evidence_profile$provenance)))

  d_rf <- diagnose_fixture(real, art, 1/3, thrs,
                           learner = "rf", decoy = "random_permutation")
  expect_true(all(d_rf$evidence_profile$status == "empirically_concerning"))
  expect_true("STABL_DIAG_EVIDENCE_CONCERNING" %in% d_rf$warning_codes)

  d_ob <- diagnose_fixture(real, art, 1/3, thrs, learner = "oblique_rf")
  expect_true("STABL_DIAG_EVIDENCE_CONCERNING" %in% d_ob$warning_codes)

  # splsda + binomial has MINT-benchmark calibration evidence (validated)
  d_sp <- diagnose_fixture(real, art, 1/3, thrs,
                           learner = "splsda", family = "binomial",
                           decoy = "random_permutation")
  expect_true(all(d_sp$evidence_profile$status == "validated_in_stablr"))
  expect_false("STABL_DIAG_EVIDENCE_NOT_EVALUATED" %in% d_sp$warning_codes)

  # splsda outside the evaluated family still surfaces absence of evidence
  d_ne <- diagnose_fixture(real, art, 1/3, thrs,
                           learner = "splsda", family = "gaussian")
  expect_true(all(d_ne$evidence_profile$status == "not_evaluated"))
  expect_true("STABL_DIAG_EVIDENCE_NOT_EVALUATED" %in% d_ne$warning_codes)

  # lasso + knockoff_equi + binomial matches only the MINT high-dimensional
  # multi-cohort row: concerning, and it is the sole matching row
  d_mint <- diagnose_fixture(real, art, 1/3, thrs,
                             learner = "lasso", family = "binomial",
                             decoy = "knockoff_equi")
  expect_identical(d_mint$evidence_profile$status, "empirically_concerning")
  expect_true("STABL_DIAG_EVIDENCE_CONCERNING" %in% d_mint$warning_codes)

  d_none <- diagnose_fixture(real, art, 1/3, thrs,
                             learner = "nonexistent_learner")
  expect_identical(d_none$evidence_profile$status, "not_evaluated")
  expect_true("STABL_DIAG_EVIDENCE_NOT_EVALUATED" %in% d_none$warning_codes)
})

# ---- 9. Integration: real fits, all families, read-only contract --------------

test_that("stabl_diagnose works on real fits and never changes the fit", {
  fit <- make_tiny_fit("gaussian")
  expect_identical(fit$base_learner, "lasso")
  expect_identical(fit$family, "gaussian")

  support_before <- get_support(fit)
  fit_snapshot   <- fit
  d <- stabl_diagnose(fit)

  expect_s3_class(d, "stabl_diagnostics")
  expect_identical(d$context$family, "gaussian")
  expect_identical(d$context$base_learner, "lasso")
  expect_identical(d$context$artificial_type, "random_permutation")
  expect_false(is.null(d$decoy_provenance))
  expect_identical(d$decoy_provenance$requested_type, "random_permutation")

  # read-only contract
  expect_identical(get_support(fit), support_before)
  expect_identical(fit, fit_snapshot)

  # print / plot smoke
  expect_output(print(d), "STABL selection diagnostics")
  expect_output(print(d), "selection rule")
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot(d)
  expect_type(plots, "list")
  expect_s3_class(plots$fdp_curve, "ggplot")
  expect_s3_class(plots$threshold_path, "ggplot")
  expect_s3_class(plots$score_distributions, "ggplot")
  plots_one <- plot(d, which = "fdp_curve")
  expect_length(plots_one, 1L)
})

test_that("stabl_diagnose works on binomial fits", {
  fit <- make_tiny_fit("binomial")
  d <- stabl_diagnose(fit)
  expect_s3_class(d, "stabl_diagnostics")
  expect_identical(d$context$family, "binomial")
})

test_that("diagnostics cover all response families (workhorse level)", {
  # The diagnostic operates on stability scores only; the response family
  # enters via context/evidence, so family coverage is exact at this level.
  real <- cbind(c(0.9, 0.8, 0.7, 0.6, 0.4, 0.2), rep(0.1, 6L))
  art  <- cbind(c(0.45, 0.30), rep(0.1, 2L))
  for (fam in c("gaussian", "binomial", "multinomial", "ordinal", "cox")) {
    d <- diagnose_fixture(real, art, 1/3, c(0, 0.25, 0.5, 0.75),
                          family = fam)
    expect_s3_class(d, "stabl_diagnostics")
    expect_identical(d$context$family, fam)
  }
})

test_that("stabl_diagnose works on ordinal fits", {
  skip_if_not_installed("ordinalNet")
  set.seed(2)
  n <- 60L; p <- 6L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- factor(sample(c("a", "b", "c"), n, replace = TRUE),
              levels = c("a", "b", "c"))
  names(y) <- rownames(x)
  fit <- suppressWarnings(stabl_fit(
    x, y, lambda_grid = "auto", base_learner = "ordinal",
    family = "ordinal", n_bootstraps = 10L,
    artificial_type = "random_permutation", random_state = 99L
  ))
  d <- stabl_diagnose(fit)
  expect_s3_class(d, "stabl_diagnostics")
  expect_identical(d$context$family, "ordinal")
  expect_identical(d$context$base_learner, "ordinal")
})

test_that("stabl_diagnose works on cox fits", {
  skip_if_not_installed("survival")
  set.seed(3)
  n <- 50L; p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- survival::Surv(time = rexp(n, rate = 0.1),
                      event = rbinom(n, size = 1L, prob = 0.7))
  rownames(y) <- rownames(x)
  fit <- stabl_fit(
    x, y,
    lambda_grid     = data.frame(lambda = seq(0.02, 0.2, length.out = 3L)),
    base_learner    = "lasso", family = "cox",
    n_bootstraps    = 5L, artificial_type = "random_permutation",
    random_state    = 31L
  )
  d <- stabl_diagnose(fit)
  expect_s3_class(d, "stabl_diagnostics")
  expect_identical(d$context$family, "cox")
})

test_that("stabl_diagnose works on stabl() wrapper objects", {
  skip_if_not_installed("knockoff")
  set.seed(4)
  x <- matrix(rnorm(60 * 8), 60, 8,
              dimnames = list(NULL, paste0("f", 1:8)))
  y <- setNames(x[, 1L] * 2 + rnorm(60), rownames(x))
  fit <- stabl(x, y, n_bootstraps = 10L, random_state = 1L)
  d <- stabl_diagnose(fit)
  expect_s3_class(d, "stabl_diagnostics")
  expect_identical(d$context$base_learner, "lasso")
})

# ---- 10. Input validation -------------------------------------------------------

test_that("stabl_diagnose rejects bad inputs", {
  expect_error(stabl_diagnose(list()), "stabl_scores_")
  fit <- make_tiny_fit("gaussian")
  expect_error(stabl_diagnose(fit, bogus = 1), "Unused arguments")
  expect_error(stabl_diagnose(fit, saturation_proportion = 0))
  expect_error(stabl_diagnose(fit, saturation_proportion = 2))
  expect_error(stabl_diagnose(fit, saturation_threshold = -0.1))
  expect_error(stabl_diagnose(fit, low_decoy_count = -1L))
})

test_that("low-decoy-count instability fires at the selected threshold", {
  # fixture A has zero decoys above the selected threshold 0.5
  real <- cbind(c(0.9, 0.8, 0.7, 0.6, 0.4, 0.2), rep(0.1, 6L))
  art  <- cbind(c(0.45, 0.30), rep(0.1, 2L))
  d <- diagnose_fixture(real, art, 1/3, c(0, 0.25, 0.5, 0.75),
                        low_decoy_count = 5L)
  expect_true("STABL_DIAG_LOW_DECOY_COUNT" %in% d$warning_codes)
  expect_true(any(grepl("unstable", d$assumption_notes)))
})

test_that("NULL metadata dimension matches only 'any' rows (B6 regression)", {
  res <- stablr:::.match_diagnostic_evidence("lasso", NULL, "gaussian")
  expect_gt(nrow(res), 0L)
  expect_true(all(res$decoy == "any"))
  expect_true(all(res$learner %in% c("any", "lasso")))
  expect_true(all(res$family %in% c("any", "gaussian")))
})

test_that("assumption_notes are single strings, not fragments (B7 regression)", {
  # selection_rule == "none" note (fixture G setup).
  real <- cbind(c(0.9, 0.2, 0.1), rep(0.05, 3L))
  d <- .diagnose_stabl_scores(
    stabl_scores = real, stabl_scores_artificial = NULL,
    fdr_threshold_range = c(0, 0.25, 0.5), hard_threshold = NULL,
    base_learner = "lasso", family = "gaussian",
    n_bootstraps = 100L, low_decoy_count = 0L
  )
  expect_true(any(grepl("No threshold is available.*empty by contract",
                        d$assumption_notes)))
  # Boundary-minimum note (fixture D setup).
  real2 <- cbind(c(0.4, 0.3, 0.2), rep(0.01, 3L))
  art2  <- cbind(c(0.35, 0.25), rep(0.01, 2L))
  d2 <- diagnose_fixture(real2, art2, 2/3, c(0, 0.25, 0.5))
  expect_true("STABL_DIAG_BOUNDARY_MINIMUM" %in% d2$warning_codes)
  expect_true(any(grepl("edge of fdr_threshold_range.*widening the range",
                        d2$assumption_notes)))
})

test_that("min_fdr not exactly in FDRs yields zero minima, not a crash", {
  # Hand-fixture path only: real fits always have min_fdr == min(FDRs) on
  # the same vector.  Pre-fix this threw "subscript out of bounds".
  d <- .diagnose_stabl_scores(
    stabl_scores            = cbind(c(0.4, 0.3, 0.2), rep(0.01, 3L)),
    stabl_scores_artificial = cbind(c(0.35, 0.25), rep(0.01, 2L)),
    fdr_threshold_range     = c(0, 0.25, 0.5),
    FDRs                    = c(1.2, 0.8, 1.0),
    min_fdr                 = 0.75,  # not exactly any element of FDRs
    fdr_min_threshold       = 0.25,
    artificial_type         = "knockoff_equi",
    artificial_proportion   = 2/3,
    effective_artificial_proportion = 2/3,
    base_learner            = "lasso", family = "gaussian",
    n_bootstraps            = 100L,
    feature_names           = c("f1", "f2"),
    low_decoy_count         = 0L
  )
  expect_identical(d$fdp_curve$n_minima, 0L)
  expect_false(d$fdp_curve$is_boundary)
  expect_false(d$fdp_curve$is_tied)
})
