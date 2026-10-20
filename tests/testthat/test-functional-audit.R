# Functional audit (2026-10-21): close the gaps found by the covr
# line-coverage run (87.18% baseline) and the integration-chain matrix.
#
#   * Chains previously covered only piecewise (no single end-to-end file):
#       - study tools: study_split -> train_validate -> study_diagnostics
#       - artificial features: make_artificial_features -> compute_fdp_plus
#   * Retained internal cooperative validation branch (fenced from the public
#     API in 1.0.0, F-CV-04): non-mse measures (binomial/poisson deviance,
#     auc, class, mae) were untested.
#   * Multinomial late-fusion stacking with a validation set
#     (.apply_multiclass_stack_weights was unreachable by the existing
#     train-only multinomial test).
#   * Registered S3 methods never executed: family.cv.multiview,
#     plot.multiview.
#   * Bootstrap class_weights block and ordinal-learner error branches.
#
# Vendored multiview utilities that no stablr call path reaches (auc,
# auc.mat, cv.elnet, mvglmnet.fit, pen_function, process_exclude, and the
# cooperative-utils index helpers) are intentionally NOT tested here: they
# are kept for upstream file parity with multiview 1.0 and are a documented
# audit exception. The multiview-parity fixtures lock every vendored path
# stablr actually exercises against CRAN multiview 1.0.

# ── Helpers ───────────────────────────────────────────────────────────────────

.make_fa_study_fixture <- function(n_studies = 4L, n_per_study = 6L,
                                   n_features = 6L, seed = 11L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n_studies * n_per_study))
  y <- setNames(factor(rep(c("A", "B"), times = length(ids) / 2L)), ids)
  x <- matrix(
    rnorm(length(ids) * n_features, sd = 0.5), length(ids), n_features,
    dimnames = list(ids, paste0("g", seq_len(n_features)))
  )
  x[, 1] <- x[, 1] + 1.2 * as.numeric(y == "B")
  study_labels <- setNames(
    rep(paste0("study", seq_len(n_studies)), each = n_per_study), ids
  )
  list(x_list = list(rna = x), y = y, study_labels = study_labels)
}

.make_fa_coop_fixture <- function(n = 40L, seed = 31L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  x_a <- matrix(rnorm(n * 5L), n, 5L,
                dimnames = list(ids, paste0("a", seq_len(5L))))
  x_b <- matrix(rnorm(n * 5L), n, 5L,
                dimnames = list(ids, paste0("b", seq_len(5L))))
  signal <- 0.9 * x_a[, 1L] - 0.7 * x_b[, 2L]
  list(
    x = list(omic_a = x_a, omic_b = x_b),
    signal = signal,
    tr = seq_len(30L), vl = 31L:40L
  )
}

.fa_coop_fit <- function(fixture, y, family, type_measure) {
  x_tr <- lapply(fixture$x, function(m) m[fixture$tr, , drop = FALSE])
  x_vl <- lapply(fixture$x, function(m) m[fixture$vl, , drop = FALSE])
  cooperative_args <- .normalize_cooperative_multiomic_args(
    x_list                   = x_tr,
    family                   = family,
    rho                      = c(0, 0.3),
    cooperative_selection    = "validation",
    cooperation_selector     = "lambda.min",
    cooperation_type_measure = type_measure,
    cooperation_nfolds       = 5L,
    x_valid_list             = x_vl,
    y_valid                  = y[fixture$vl]
  )
  suppressWarnings(.cooperative_multiomic_fit(
    x_train_list     = x_tr,
    y_train          = y[fixture$tr],
    x_valid_list     = x_vl,
    y_valid          = y[fixture$vl],
    groups_train     = NULL,
    family           = family,
    random_state     = 16L,
    cooperative_args = cooperative_args
  ))
}

# ── 1. study-tools chain ─────────────────────────────────────────────────────

test_that("study chain composes end to end: split -> train_validate -> diagnostics", {
  fixture <- .make_fa_study_fixture()

  split <- stabl_multiomic_study_split(
    x_list = fixture$x_list, y = fixture$y,
    study_labels = fixture$study_labels, prop_train = 0.5, random_state = 3L
  )
  expect_setequal(c(split$train_ids, split$valid_ids), names(fixture$y))

  fit <- suppressWarnings(stabl_multiomic_train_validate(
    x_train_list = lapply(fixture$x_list, function(m) {
      m[split$train_ids, , drop = FALSE]
    }),
    y_train = fixture$y[split$train_ids],
    x_valid_list = lapply(fixture$x_list, function(m) {
      m[split$valid_ids, , drop = FALSE]
    }),
    y_valid = fixture$y[split$valid_ids],
    lambda_grid = data.frame(lambda = c(0.3, 0.1)),
    family = "binomial", n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 0.5, random_state = 5L
  ))
  expect_true(all(c("fits", "selected_features", "selected_train",
                    "selected_valid") %in% names(fit)))
  expect_true(length(fit$selected_features$rna) >= 1L)

  diag <- stabl_multiomic_study_diagnostics(
    x_list = fixture$x_list, y = fixture$y, study_labels = fixture$study_labels
  )
  expect_equal(nrow(diag$study_summary), 4L)
  expect_equal(diag$study_summary$n_samples, rep(6L, 4L))
  # the split is study-disjoint and partitions exactly the study structure
  # the diagnostics report
  expect_identical(split$fold_scheme, "study_disjoint")
  train_studies <- unique(as.character(fixture$study_labels[split$train_ids]))
  valid_studies <- unique(as.character(fixture$study_labels[split$valid_ids]))
  expect_length(intersect(train_studies, valid_studies), 0L)
  expect_setequal(
    union(train_studies, valid_studies),
    paste0("study", seq_len(4L))
  )
})

# ── 2. artificial-features chain ─────────────────────────────────────────────

test_that("artificial-feature chain composes: make_artificial_features -> scoring -> compute_fdp_plus", {
  set.seed(21L)
  n <- 60L
  p <- 10L
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  # five features share a latent signal so several real features separate
  # cleanly from the artificial floor (FDP+ has a +1 numerator, so its
  # attainable minimum here is 1 / n_real, never 0)
  latent <- rnorm(n)
  for (j in 1:5) {
    x[, j] <- 0.9 * latent + sqrt(1 - 0.9^2) * rnorm(n)
  }
  y <- 2.5 * latent + rnorm(n, sd = 0.3)

  gen <- make_artificial_features(
    x, n_injected = p, type = "random_permutation", random_state = 21L
  )
  expect_equal(dim(gen$x_augmented), c(n, 2L * p))
  expect_length(gen$noise_col_indices, p)

  # minimal stability scorer: bootstrap |correlation| over the augmented design
  B <- 25L
  set.seed(22L)
  score_all <- replicate(B, {
    idx <- sample(seq_len(n), n %/% 2L)
    c(abs(stats::cor(gen$x_augmented[idx, , drop = FALSE], y[idx])))
  })
  rownames(score_all) <- colnames(gen$x_augmented)
  # the artificial block is appended as the LAST n_injected columns of
  # x_augmented (stabl_fit.R splits scores at n_features + seq_len(n_injected));
  # noise_col_indices is the source-column map, not augmented positions
  s_real <- score_all[seq_len(p), , drop = FALSE]
  s_art  <- score_all[p + seq_len(p), , drop = FALSE]

  pi_hat <- length(gen$noise_col_indices) / ncol(gen$x_augmented)
  res <- compute_fdp_plus(
    stabl_scores            = s_real,
    stabl_scores_artificial = s_art,
    artificial_proportion   = pi_hat,
    fdr_threshold_range     = seq(0, 0.99, by = 0.01)
  )
  expect_length(res$FDRs, 100L)
  expect_true(all(is.finite(res$FDRs)))
  expect_equal(res$min_fdr, min(res$FDRs))

  # hand-check the STABL formula on the composed scores at one threshold
  # (aggregation is per-feature row max; counting is strict `>`)
  q_real <- apply(s_real, 1L, max)
  q_art  <- apply(s_art, 1L, max)
  t0 <- 0.5
  expect_equal(
    unname(res$FDRs[abs(seq(0, 0.99, by = 0.01) - t0) < 1e-12]),
    ((1 / pi_hat) * sum(q_art > t0) + 1) / max(1, sum(q_real > t0))
  )

  # the composed pipeline separates the five planted signals from the
  # artificial floor: the min-FDP+ threshold is attained where only real
  # signal passes (attainable floor is 1 / n_real because of the +1)
  expect_lte(res$min_fdr, 0.2)
  expect_true(all(q_real[paste0("f", 1:5)] > res$fdr_min_threshold))

  # the generator is deterministic under random_state
  gen2 <- make_artificial_features(
    x, n_injected = p, type = "random_permutation", random_state = 21L
  )
  expect_identical(gen$x_augmented, gen2$x_augmented)
})

# ── 3. retained cooperative validation branch: non-mse measures ─────────────

test_that("validation-based cooperative selection supports binomial deviance", {
  fixture <- .make_fa_coop_fixture()
  y <- factor(ifelse(fixture$signal > median(fixture$signal), "B", "A"))
  names(y) <- rownames(fixture$x$omic_a)

  cf <- .fa_coop_fit(fixture, y, family = "binomial", type_measure = "deviance")
  expect_equal(cf$selection, "validation")
  expect_true(is.finite(cf$score))
  expect_gte(cf$score, 0) # deviance is non-negative
})

test_that("validation-based cooperative selection supports binomial auc", {
  fixture <- .make_fa_coop_fixture()
  y <- factor(ifelse(fixture$signal > median(fixture$signal), "B", "A"))
  names(y) <- rownames(fixture$x$omic_a)

  cf <- .fa_coop_fit(fixture, y, family = "binomial", type_measure = "auc")
  expect_equal(cf$selection, "validation")
  expect_true(is.finite(cf$score))
  expect_gte(cf$score, 0)
  expect_lte(cf$score, 1)
})

test_that("validation-based cooperative selection supports binomial class error", {
  fixture <- .make_fa_coop_fixture()
  y <- factor(ifelse(fixture$signal > median(fixture$signal), "B", "A"))
  names(y) <- rownames(fixture$x$omic_a)

  cf <- .fa_coop_fit(fixture, y, family = "binomial", type_measure = "class")
  expect_equal(cf$selection, "validation")
  expect_true(is.finite(cf$score))
  expect_gte(cf$score, 0)
  expect_lte(cf$score, 1)
})

test_that("poisson cooperative deviance is fenced publicly and verified directly", {
  fixture <- .make_fa_coop_fixture()
  set.seed(32L)
  y <- setNames(
    rpois(length(fixture$signal), exp(0.3 * scale(fixture$signal)[, 1L])),
    rownames(fixture$x$omic_a)
  )

  # Public fence (F-CV-04): cooperative fusion accepts gaussian/binomial
  # only; the poisson branch is retained internally for upstream parity.
  expect_error(
    .fa_coop_fit(fixture, y, family = "poisson", type_measure = "deviance"),
    "currently supports family"
  )

  # Direct verification of the retained branch against an independent
  # dpois-based deviance (the lgamma terms cancel in the difference).
  eta <- log(c(0.5, 1.5, 2.5, 4))
  yy  <- c(0, 1, 2, 5)
  expected <- mean(2 * (stats::dpois(yy, yy, log = TRUE) -
                          stats::dpois(yy, exp(eta), log = TRUE)))
  expect_equal(.cooperative_poisson_deviance(eta, yy), expected)
  # the internal metric router still dispatches to the poisson branch
  expect_equal(
    .cooperative_validation_metric(yy, eta, family = "poisson",
                                   type_measure = "deviance"),
    expected
  )
})

test_that("cooperative validation metrics match independent computations", {
  # binomial deviance vs dbinom (saturated binary model has zero log-lik)
  p  <- c(0.2, 0.7, 0.9)
  yb <- factor(c("A", "B", "B"))
  expect_equal(
    .cooperative_binomial_deviance(p, yb),
    mean(-2 * stats::dbinom(c(0L, 1L, 1L), size = 1L, prob = p, log = TRUE))
  )

  # binary outcome coercion: factor, 0/1 numeric, and non-{0,1} numeric
  expect_identical(
    .cooperative_coerce_binary_outcome(factor(c("A", "B", "A"))),
    c(0L, 1L, 0L)
  )
  expect_identical(.cooperative_coerce_binary_outcome(c(0, 1, 1)),
                   c(0L, 1L, 1L))
  expect_identical(.cooperative_coerce_binary_outcome(c(1, 2, 2)),
                   c(0L, 1L, 1L))

  # dispatch error branches
  expect_error(
    .cooperative_validation_metric(1:3, 1:3, family = "cox",
                                   type_measure = "deviance"),
    "does not support family 'cox'"
  )
  expect_error(
    .cooperative_validation_metric(1:3, 1:3, family = "gaussian",
                                   type_measure = "bogus"),
    "Unsupported cooperative type measure"
  )
})

test_that("validation-based cooperative selection supports gaussian mae", {
  fixture <- .make_fa_coop_fixture()
  set.seed(33L)
  y <- setNames(fixture$signal + rnorm(length(fixture$signal), sd = 0.3),
                rownames(fixture$x$omic_a))

  cf <- .fa_coop_fit(fixture, y, family = "gaussian", type_measure = "mae")
  expect_equal(cf$selection, "validation")
  expect_true(is.finite(cf$score))
  expect_gte(cf$score, 0)
})

# ── 4. multinomial late fusion with a validation split ───────────────────────

test_that("multinomial late fusion applies stack weights on a validation split", {
  set.seed(41L)
  n <- 36L
  ids <- paste0("s", seq_len(n))
  # interleave classes so both the train and the validation split observe
  # all three levels (the late-fusion validator requires matching levels)
  y <- setNames(factor(rep(c("A", "B", "C"), times = 12L)), ids)
  signal <- model.matrix(~ y - 1)
  x_a <- matrix(rnorm(n * 6L, sd = 0.2), n, 6L,
                dimnames = list(ids, paste0("a", seq_len(6L))))
  x_b <- matrix(rnorm(n * 5L, sd = 0.2), n, 5L,
                dimnames = list(ids, paste0("b", seq_len(5L))))
  x_a[, 1:3] <- x_a[, 1:3] + signal
  x_b[, 1:3] <- x_b[, 1:3] + signal
  tr <- seq_len(24L)
  vl <- 25L:36L

  fit <- suppressWarnings(stabl_multiomic_train_validate(
    x_train_list = list(omic_a = x_a[tr, ], omic_b = x_b[tr, ]),
    y_train = y[tr],
    x_valid_list = list(omic_a = x_a[vl, ], omic_b = x_b[vl, ]),
    y_valid = y[vl],
    lambda_grid = data.frame(lambda = c(0.05, 0.01)),
    artificial_type = NULL, hard_threshold = 1e-9,
    n_bootstraps = 2L, sample_fraction = 1,
    family = "multinomial", random_state = 41L,
    late_fusion = TRUE, n_iter_lf = 100L
  ))

  lf <- fit$late_fusion
  expect_equal(lf$task_type, "multiclass")
  expect_true(!is.null(lf$valid_predictions))
  expect_true(
    all(paste0("prob_", levels(y)) %in% names(lf$valid_predictions))
  )
  prob_cols <- as.matrix(lf$valid_predictions[, paste0("prob_", levels(y))])
  expect_true(all(is.finite(prob_cols)))
  expect_true(all(prob_cols >= 0 & prob_cols <= 1))
})

# ── 5. vendored multiview S3 methods ─────────────────────────────────────────

test_that("family and plot S3 methods work for vendored multiview fits", {
  set.seed(51L)
  n <- 24L
  x_a <- matrix(rnorm(n * 4L), n, 4L,
                dimnames = list(NULL, paste0("a", seq_len(4L))))
  x_b <- matrix(rnorm(n * 3L), n, 3L,
                dimnames = list(NULL, paste0("b", seq_len(3L))))
  y <- 0.8 * x_a[, 1L] + rnorm(n, sd = 0.5)
  x_list <- list(a = x_a, b = x_b)
  lam <- c(0.3, 0.1, 0.03)

  fit <- suppressWarnings(.cooperative_backend_fit(
    x_list, y, family = "gaussian", rho = 0.2, lambda = lam
  ))
  expect_s3_class(fit, "multiview")
  expect_identical(family(fit), "gaussian")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit))

  cvf <- suppressWarnings(.mv_cv_multiview(
    x_list, y, family = gaussian(), rho = 0.2, nfolds = 4L, lambda = lam
  ))
  expect_s3_class(cvf, "cv.multiview")
  expect_identical(family(cvf), "gaussian")
})

# ── 6. bootstrap class_weights ───────────────────────────────────────────────

test_that("sample_bootstrap_indices honours and validates class_weights", {
  y <- factor(rep(c("A", "B"), each = 10L))

  s1 <- sample_bootstrap_indices(
    y, n_subsamples = 50L, class_weights = c(A = 3, B = 1), random_state = 61L
  )
  s2 <- sample_bootstrap_indices(
    y, n_subsamples = 50L, class_weights = c(A = 3, B = 1), random_state = 61L
  )
  expect_identical(s1, s2)

  idx <- as.matrix(s1)
  expect_true(all(idx >= 1L & idx <= length(y)))
  # 3:1 weights must oversample class A beyond the balanced 0.5 baseline
  prop_a <- mean(y[as.vector(idx)] == "A")
  expect_gt(prop_a, 0.6)

  expect_error(
    sample_bootstrap_indices(y, 5L, class_weights = c(2, 1)),
    "named numeric vector"
  )
  expect_error(
    sample_bootstrap_indices(y, 5L, class_weights = c(A = 2)),
    "All observed classes"
  )
  expect_error(
    sample_bootstrap_indices(
      y, 5L, class_weights = c(A = 2, B = 1), stratify = TRUE
    ),
    "cannot be combined"
  )
  expect_error(
    sample_bootstrap_indices(
      y, 5L, class_weights = c(A = 2, B = 1),
      strata = data.frame(site = rep("x", 20L))
    ),
    "cannot be combined"
  )
})

# ── 7. ordinal-learner error branches ────────────────────────────────────────

test_that("ordinal outcome validation rejects non-factor and NA outcomes", {
  expect_error(.validate_ordinal_outcome(c(1, 2, 3)), "factor")
  y <- ordered(rep(c("lo", "mid", "hi"), each = 3L),
               levels = c("lo", "mid", "hi"))
  y[1L] <- NA
  expect_error(.validate_ordinal_outcome(y), "must not contain NA")
})

test_that("ordinal batch adapter rejects malformed lambda grids", {
  adapter <- .make_ordinal_batch_adapter()
  x <- matrix(rnorm(30L * 4L), 30L, 4L)
  y <- ordered(rep(c("lo", "mid", "hi"), each = 10L),
               levels = c("lo", "mid", "hi"))
  expect_error(
    adapter(x, y, data.frame(lambda = c(0.2, NA))),
    "finite numeric"
  )
})

# ── 8. stratified-bootstrap branches ─────────────────────────────────────────
# Second-pass coverage: input-validation and allocation branches in
# bootstrap_helpers.R that the first-pass class_weights tests do not reach.

test_that("stratified bootstrap validates and coerces strata designs", {
  y <- factor(rep(c("A", "B"), each = 10L))

  expect_error(
    sample_bootstrap_indices(y, 5L, stratify = NA),
    "TRUE or FALSE"
  )

  # matrix and list strata are coerced to a design frame
  expect_no_error(sample_bootstrap_indices(
    y, 6L, strata = matrix(rep(c("x", "y"), each = 10L), ncol = 1L),
    random_state = 71L
  ))
  expect_no_error(sample_bootstrap_indices(
    y, 6L, strata = list(site = rep(c("x", "y"), each = 10L)),
    random_state = 71L
  ))

  # malformed strata fail closed
  expect_error(
    sample_bootstrap_indices(y, 6L, strata = c("x", "y")),
    "must have length"
  )
  expect_error(
    sample_bootstrap_indices(y, 6L, strata = data.frame(a = rep("x", 5L))),
    "must have .* rows"
  )
  expect_error(
    sample_bootstrap_indices(y, 6L, strata = data.frame(a = rep("x", 20L))[, 0]),
    "at least one column"
  )
  expect_error(
    sample_bootstrap_indices(y, 6L, strata = c(rep("x", 19L), NA)),
    "cannot contain missing values"
  )
})

test_that("stratified count allocation covers edge branches", {
  # single realised stratum: fast path returns the full count
  expect_identical(
    .stratified_counts(rep("only", 8L), 5L, replace = TRUE),
    c(only = 5L)
  )
  # n_subsamples below the realised-strata minimum
  expect_error(
    .stratified_counts(rep(c("a", "b"), each = 5L), 1L, replace = TRUE),
    "at least the number of realised strata"
  )
  expect_error(
    .stratified_counts(
      rep(c("a", "b"), each = 5L), 3L, replace = TRUE, min_per_stratum = 2L
    ),
    "to allocate"
  )
  # replace = FALSE requires every stratum to meet the minimum
  expect_error(
    .stratified_counts(
      c(rep("a", 5L), "b"), 4L, replace = FALSE, min_per_stratum = 2L
    ),
    "when `replace = FALSE`"
  )
  # largest-remainder bumping can overshoot; the trim loop rebalances
  expect_identical(
    .stratified_counts(
      c(rep("a", 100L), "b", "c", "d"), 4L, replace = TRUE
    ),
    c(a = 1L, b = 1L, c = 1L, d = 1L)
  )
  # replace = FALSE fails closed when a count exceeds its stratum size
  expect_error(
    .stratified_counts(c("a", "a", "b"), 5L, replace = FALSE),
    "exceeds an observed stratum size"
  )
})

test_that("bootstrap strata subsetting matches sample ids fail-closed", {
  ids <- paste0("s", 1:4)
  df <- data.frame(site = c("x", "x", "y", "y"))

  # no row ids -> positional, sample ids assigned
  out <- .subset_bootstrap_strata_by_ids(df, ids)
  expect_equal(rownames(out), ids)
  # matching row ids -> reordered to sample_ids (values follow their id)
  df2 <- df
  rownames(df2) <- rev(ids)
  out2 <- .subset_bootstrap_strata_by_ids(df2, ids)
  expect_equal(rownames(out2), ids)
  expect_equal(out2$site, c("y", "y", "x", "x"))
  # mismatched row ids -> error
  df3 <- df
  rownames(df3) <- c("s1", "s2", "s3", "zzz")
  expect_error(.subset_bootstrap_strata_by_ids(df3, ids), "Sample mismatch")
  # wrong row count without ids -> error
  expect_error(
    .subset_bootstrap_strata_by_ids(df[1:3, , drop = FALSE], ids),
    "one row per sample"
  )
  # matrix and list inputs are coerced to design frames
  expect_s3_class(
    .subset_bootstrap_strata_by_ids(matrix(c("x", "x", "y", "y"), ncol = 1L), ids),
    "data.frame"
  )
  expect_s3_class(
    .subset_bootstrap_strata_by_ids(list(site = c("x", "x", "y", "y")), ids),
    "data.frame"
  )
  # atomic strata: positional without names, matched with names
  expect_identical(
    .subset_bootstrap_strata_by_ids(c("x", "x", "y", "y"), ids),
    stats::setNames(c("x", "x", "y", "y"), ids)
  )
  expect_error(
    .subset_bootstrap_strata_by_ids(c("x", "y"), ids),
    "one value per sample"
  )
  expect_identical(
    .subset_bootstrap_strata_by_ids(c(s4 = "y", s3 = "y", s2 = "x", s1 = "x"), ids),
    c(s1 = "x", s2 = "x", s3 = "y", s4 = "y")
  )
  expect_error(
    .subset_bootstrap_strata_by_ids(c(s1 = "x", s2 = "x", s3 = "y", zz = "y"), ids),
    "Sample mismatch"
  )
  expect_null(.subset_bootstrap_strata_by_ids(NULL, ids))

  # row-id detection: NULL / default / NA ids are not usable ids
  expect_false(.bootstrap_strata_has_row_ids(NULL, ids, 4L))
  expect_false(.bootstrap_strata_has_row_ids(c("1", "2", "3", "4"), ids, 4L))
  expect_false(.bootstrap_strata_has_row_ids(c("s1", NA, "s3", "s4"), ids, 4L))
  expect_true(.bootstrap_strata_has_row_ids(ids, ids, 4L))
})

# ── 9. sparse-group-lasso and coefficient-helper branches ────────────────────

test_that("sparse-group-lasso adapter validates inputs and fits multinomial", {
  skip_if_not_installed("sparsegl")
  set.seed(79L)
  x <- matrix(rnorm(40L * 4L), 40L, 4L,
              dimnames = list(NULL, paste0("g", 1:4)))
  y3 <- factor(rep(c("A", "B", "C"), length.out = 40L))
  grp <- c(1L, 1L, 2L, 2L)

  # factory validation branches
  expect_error(
    make_sgl_adapter(family = "cox", feature_groups = grp),
    "not supported by sparse_group_lasso"
  )
  expect_error(
    make_sgl_adapter(feature_groups = integer(0)),
    "cannot be empty"
  )
  expect_error(
    make_sgl_adapter(feature_groups = c(1, NA, 2, 2)),
    "cannot contain NA"
  )
  expect_error(
    make_sgl_adapter(feature_groups = grp, alpha_fixed = 2),
    "numeric scalar in \\[0, 1\\]"
  )

  ad <- make_sgl_adapter(family = "multinomial", feature_groups = grp)
  # group vector must match the design width at call time
  expect_error(
    ad(x[, 1:3], y3, data.frame(lambda = 0.1)),
    "length must match ncol"
  )
  # multinomial one-vs-rest path, alpha resolved from the grid row and
  # from the package default
  suppressWarnings({
    sel_alpha <- ad(x, y3, data.frame(lambda = 0.1, alpha = 0.5))
    sel_default <- ad(x, y3, data.frame(lambda = 0.1))
  })
  expect_type(sel_alpha, "logical")
  expect_length(sel_alpha, 4L)
  expect_type(sel_default, "logical")
  # invalid alpha inside the grid row fails closed
  expect_error(
    ad(x, y3, data.frame(lambda = 0.1, alpha = 2)),
    "`alpha` in `lambda_val`"
  )
  # alpha_fixed resolves directly
  ad_fixed <- make_sgl_adapter(
    family = "gaussian", feature_groups = grp, alpha_fixed = 0.3
  )
  suppressWarnings(
    sel_fixed <- ad_fixed(x, rnorm(40L), data.frame(lambda = 0.1))
  )
  expect_type(sel_fixed, "logical")
})

test_that("sparse-group-lasso batch adapter covers alpha batching and multinomial", {
  skip_if_not_installed("sparsegl")
  set.seed(80L)
  x <- matrix(rnorm(40L * 4L), 40L, 4L)
  y <- rnorm(40L)
  y3 <- factor(rep(c("A", "B", "C"), length.out = 40L))
  grp <- c(1L, 1L, 2L, 2L)
  grid <- data.frame(lambda = c(0.2, 0.1))

  expect_error(
    .make_sgl_batch_adapter("cox", grp, NULL, 1e-8),
    "not supported by sparse_group_lasso"
  )
  suppressWarnings({
    # alpha_fixed: one batch over the whole grid
    b_fixed <- .make_sgl_batch_adapter("gaussian", grp, 0.5, 1e-8)(x, y, grid)
    # alpha column: one batch per distinct alpha
    b_alpha <- .make_sgl_batch_adapter("gaussian", grp, NULL, 1e-8)(
      x, y, data.frame(lambda = c(0.2, 0.1, 0.15), alpha = c(0.5, 0.5, 0.9))
    )
    # multinomial: per-class one-vs-rest max
    b_multi <- .make_sgl_batch_adapter("multinomial", grp, NULL, 1e-8)(x, y3, grid)
  })
  expect_equal(dim(b_fixed), c(4L, 2L))
  expect_type(b_fixed, "logical")
  expect_equal(dim(b_alpha), c(4L, 3L))
  expect_equal(dim(b_multi), c(4L, 2L))
})

test_that("coefficient helpers cover list-coef and near-lambda branches", {
  skip_if_not_installed("glmnet")
  set.seed(83L)
  x <- matrix(rnorm(30L * 5L), 30L, 5L)
  y3 <- factor(rep(c("A", "B", "C"), each = 10L))
  suppressWarnings(
    fit <- glmnet::glmnet(x, y3, family = "multinomial", lambda = c(0.1, 0.05))
  )
  # multinomial coef() returns a per-class list -> rowMaxs aggregation
  co <- .feature_abs_coefs(fit, s = 0.05, family = "multinomial")
  expect_type(co, "double")
  expect_length(co, 5L)
  expect_true(all(co >= 0))

  # near-lambda matching: exact fast path, tolerance match, and misses
  expect_identical(.near_lambda_match(c(0.1, 0.05), c(0.1, 0.05)), c(1L, 2L))
  expect_identical(.near_lambda_match(0.1 + 1e-15, c(0.1, 0.05)), 1L)
  expect_true(is.na(.near_lambda_match(9, c(0.1, 0.05))))
  expect_true(is.na(.near_lambda_match(NA_real_, c(0.1, 0.05))))
})

test_that("group bootstrap validates its own replace and stratify inputs", {
  y <- factor(rep(c("A", "B"), each = 6L))
  groups <- rep(paste0("G", 1:4), each = 3L)
  expect_error(
    group_bootstrap_indices(y, groups, n_subsamples = 13L, replace = FALSE),
    "cannot exceed sample size"
  )
  expect_error(
    group_bootstrap_indices(y, groups, n_subsamples = 4L, stratify = NA),
    "TRUE or FALSE"
  )
})

test_that("degenerate-outcome guard passes single-class draws through", {
  idx <- .guard_diverse_bootstrap(
    factor(rep("A", 5L)), 1:5, function() 1:5, "probe", "hint"
  )
  expect_identical(idx, 1:5)
  expect_error(
    .stratified_counts(character(0), 3L, replace = TRUE),
    "at least one sample"
  )
})

test_that("sparse-group-lasso multinomial coerces non-factor outcomes", {
  skip_if_not_installed("sparsegl")
  set.seed(84L)
  x <- matrix(rnorm(40L * 4L), 40L, 4L)
  y3 <- factor(rep(c("A", "B", "C"), length.out = 40L))
  grp <- c(1L, 1L, 2L, 2L)

  ad <- make_sgl_adapter(family = "multinomial", feature_groups = grp)
  # character outcomes are coerced to factor inside the multinomial path
  suppressWarnings(sel <- ad(x, as.character(y3), data.frame(lambda = 0.1)))
  expect_type(sel, "logical")
  expect_length(sel, 4L)
  # a single observed class fails closed
  expect_error(
    ad(x, factor(rep("A", 40L)), data.frame(lambda = 0.1)),
    "at least 2 classes"
  )
  # the batch adapter applies the same coercion
  suppressWarnings(
    b <- .make_sgl_batch_adapter("multinomial", grp, NULL, 1e-8)(
      x, as.character(y3), data.frame(lambda = 0.1)
    )
  )
  expect_equal(dim(b), c(4L, 1L))
})
