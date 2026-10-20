# B2/B3 — study-level diagnostics and study-disjoint train/valid split
# (0.3.0 cross-study reliability round; inferred scope, flagged in NEWS).

.make_study_tools_fixture <- function(n_studies = 3L, n_per_study = 4L,
                                      n_features = 4L, seed = 1L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n_studies * n_per_study))
  y <- setNames(factor(rep(c("A", "B"), times = length(ids) / 2L)), ids)
  x <- matrix(
    rnorm(length(ids) * n_features, sd = 0.3),
    length(ids), n_features,
    dimnames = list(ids, paste0("g", seq_len(n_features)))
  )
  study_labels <- setNames(
    rep(paste0("study", seq_len(n_studies)), each = n_per_study), ids
  )
  list(x_list = list(rna = x), y = y, study_labels = study_labels)
}

test_that("study diagnostics report composition and recover a planted shift", {
  fixture <- .make_study_tools_fixture(seed = 2L)
  # plant a location shift on study1's first two features
  fixture$x_list$rna[fixture$study_labels == "study1", 1:2] <-
    fixture$x_list$rna[fixture$study_labels == "study1", 1:2] + 3

  diag <- stabl_multiomic_study_diagnostics(
    x_list = fixture$x_list,
    y = fixture$y,
    study_labels = fixture$study_labels
  )

  expect_named(diag, c("study_summary", "feature_shift"))
  expect_equal(nrow(diag$study_summary), 3L)
  expect_equal(nrow(diag$feature_shift), 3L)
  expect_equal(diag$study_summary$n_samples, rep(4L, 3L))
  expect_equal(diag$study_summary$n_A, rep(2L, 3L))
  expect_equal(diag$study_summary$prop_B, rep(0.5, 3L))

  shifted <- diag$feature_shift[diag$feature_shift$study == "study1", ]
  nullish <- diag$feature_shift[diag$feature_shift$study == "study2", ]
  expect_true(shifted$smd_max_abs > 2)
  expect_true(shifted$smd_max_abs > nullish$smd_max_abs)
  expect_equal(diag$feature_shift$n_zero_sd, rep(0L, 3L))
})

test_that("study diagnostics handle numeric outcomes and fail closed", {
  fixture <- .make_study_tools_fixture(seed = 3L)
  y_num <- setNames(rnorm(length(fixture$y)), names(fixture$y))

  diag <- stabl_multiomic_study_diagnostics(
    x_list = fixture$x_list, y = y_num, study_labels = fixture$study_labels
  )
  expect_true(all(c("outcome_mean", "outcome_sd") %in% names(diag$study_summary)))
  expect_equal(diag$study_summary$n_samples, rep(4L, 3L))

  expect_error(
    stabl_multiomic_study_diagnostics(
      x_list = fixture$x_list, y = fixture$y,
      study_labels = unname(fixture$study_labels)
    ),
    "named vector"
  )
  one_sample <- fixture$study_labels
  one_sample[[1L]] <- "studyX"
  expect_error(
    stabl_multiomic_study_diagnostics(
      x_list = fixture$x_list, y = fixture$y, study_labels = one_sample
    ),
    "at least 2 samples"
  )
})

test_that("study split is study-disjoint, deterministic, and class-balanced", {
  fixture <- .make_study_tools_fixture(n_studies = 5L, n_per_study = 4L, seed = 4L)

  s1 <- stabl_multiomic_study_split(
    x_list = fixture$x_list, y = fixture$y,
    study_labels = fixture$study_labels,
    prop_train = 0.6, random_state = 7L
  )
  s2 <- stabl_multiomic_study_split(
    x_list = fixture$x_list, y = fixture$y,
    study_labels = fixture$study_labels,
    prop_train = 0.6, random_state = 7L
  )
  expect_identical(s1, s2)
  expect_identical(s1$fold_scheme, "study_disjoint")
  expect_length(intersect(s1$train_ids, s1$valid_ids), 0L)
  expect_true(length(s1$train_ids) > 0L)
  expect_true(length(s1$valid_ids) > 0L)
  # sides are study-pure
  expect_length(
    intersect(
      unique(as.character(fixture$study_labels[s1$train_ids])),
      unique(as.character(fixture$study_labels[s1$valid_ids]))
    ),
    0L
  )
  # both sides keep both outcome classes
  expect_true(all(table(fixture$y[s1$train_ids]) > 0L))
  expect_true(all(table(fixture$y[s1$valid_ids]) > 0L))
  # the split realises the target proportion within study granularity
  achieved <- length(s1$train_ids) / (length(s1$train_ids) + length(s1$valid_ids))
  expect_true(abs(achieved - 0.6) <= 0.25)
})

test_that("study split fails closed on degenerate requests", {
  fixture <- .make_study_tools_fixture(n_studies = 2L, n_per_study = 4L, seed = 5L)

  expect_error(
    stabl_multiomic_study_split(
      x_list = fixture$x_list, y = fixture$y,
      study_labels = fixture$study_labels, prop_train = 1
    ),
    "strictly between 0 and 1"
  )
  expect_error(
    stabl_multiomic_study_split(
      x_list = fixture$x_list, y = fixture$y,
      study_labels = unname(fixture$study_labels)
    ),
    "named vector"
  )
  # a single study cannot be split
  solo <- fixture$study_labels
  solo[] <- "study1"
  expect_error(
    stabl_multiomic_study_split(
      x_list = fixture$x_list, y = fixture$y, study_labels = solo
    ),
    "at least 2 studies"
  )
  # extreme prop_train with 2 equal-size studies strands one side
  expect_error(
    stabl_multiomic_study_split(
      x_list = fixture$x_list, y = fixture$y,
      study_labels = fixture$study_labels,
      prop_train = 0.99, random_state = 8L
    ),
    "prop_train"
  )
})

test_that("stabl_multiomic_cv records the fold scheme", {
  fixture <- .make_study_tools_fixture(n_studies = 3L, n_per_study = 4L, seed = 6L)
  y_num <- setNames(as.numeric(fixture$y) - 1L, names(fixture$y))

  cv_iid <- suppressWarnings(stabl_multiomic_cv(
    fixture$x_list, y_num,
    lambda_grid = data.frame(lambda = 0.1), v = 2L,
    n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
    random_state = 9L
  ))
  expect_identical(cv_iid$fold_scheme, "sample_iid")

  cv_grouped <- suppressWarnings(stabl_multiomic_cv(
    fixture$x_list, y_num,
    lambda_grid = data.frame(lambda = 0.1), v = 2L,
    n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
    random_state = 9L, groups = fixture$study_labels
  ))
  expect_identical(cv_grouped$fold_scheme, "group_disjoint")
  # group-disjoint folds: no study spans both sides of any fold
  for (fold in cv_grouped$folds) {
    expect_length(
      intersect(
        unique(as.character(fixture$study_labels[fold$train_ids])),
        unique(as.character(fixture$study_labels[fold$valid_ids]))
      ),
      0L
    )
  }
})

# C2 (inferred scope): pooled vs study-disjoint scheme comparison.
test_that("scheme comparison runs both arms and reports the gap", {
  fixture <- .make_study_tools_fixture(n_studies = 6L, n_per_study = 8L, seed = 8L)
  cmp <- suppressWarnings(stabl_multiomic_scheme_comparison(
    x_list = fixture$x_list,
    y = fixture$y,
    study_labels = fixture$study_labels,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 2L, inner_v = 2L, n_bootstraps = 10L,
    artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 10L
  ))
  expect_identical(cmp$sample_iid$fold_scheme, "sample_iid")
  expect_identical(cmp$study_disjoint$fold_scheme, "study_disjoint")
  expect_equal(nrow(cmp$comparison), 2L)
  expect_identical(cmp$comparison$scheme, c("sample_iid", "study_disjoint"))
  expect_true(all(is.finite(cmp$comparison$balanced_error_rate)))
  expect_equal(
    cmp$ber_gap,
    cmp$comparison$balanced_error_rate[[2L]] -
      cmp$comparison$balanced_error_rate[[1L]]
  )
  expect_error(
    stabl_multiomic_scheme_comparison(
      fixture$x_list, fixture$y, study_labels = NULL
    ),
    "study_labels"
  )
})

# C4 (inferred scope): stability of stability across study-disjoint folds.
test_that("selection stability returns a fold x fold Jaccard matrix", {
  fixture <- .make_study_tools_fixture(n_studies = 6L, n_per_study = 8L, seed = 11L)
  # plant signal so selections are non-degenerate
  fixture$x_list$rna[, 1:2] <-
    fixture$x_list$rna[, 1:2] + 1.5 * (fixture$y == "A")

  nested <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 2L, inner_v = 2L, n_bootstraps = 10L,
    artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 12L,
    study_labels = fixture$study_labels
  ))
  stab <- stabl_multiomic_selection_stability(nested)
  fold_ids <- vapply(nested$outer_folds, `[[`, character(1L), "fold_id")
  expect_identical(dimnames(stab$similarity_matrix), list(fold_ids, fold_ids))
  expect_identical(
    dim(stab$similarity_matrix),
    c(length(fold_ids), length(fold_ids))
  )
  expect_true(all(diag(stab$similarity_matrix) == 1))
  expect_true(all(stab$similarity_matrix >= 0 & stab$similarity_matrix <= 1))
  expect_true(is.finite(stab$mean_pairwise_similarity))
  expect_true(stab$mean_pairwise_similarity >= 0)
  expect_identical(stab$fold_scheme, "study_disjoint")
  expect_true(all(stab$fold_selection_sizes >= 0))

  # fail-closed: sample-iid objects and foreign objects are rejected
  nested_iid <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 2L, inner_v = 2L, n_bootstraps = 10L,
    artificial_type = NULL, hard_threshold = 1e-9,
    sample_fraction = 1, random_state = 12L
  ))
  expect_error(
    stabl_multiomic_selection_stability(nested_iid),
    "study_disjoint"
  )
  expect_error(
    stabl_multiomic_selection_stability(list()),
    "stabl_multiomic_nested_cv"
  )
})

# C5 (inferred scope): grouped-CV single-class guardrail.
test_that("grouped CV fails closed on single-class assessment folds", {
  set.seed(14L)
  ids <- paste0("s", 1:8)
  y <- setNames(factor(c(rep("A", 4), rep("B", 4))), ids)
  x <- matrix(rnorm(8 * 4), 8, 4, dimnames = list(ids, paste0("g", 1:4)))
  groups <- setNames(rep(c("g1", "g2"), each = 4L), ids)

  # v = 2 with two single-class groups: every assignment yields a
  # single-class assessment fold, so the guard fires deterministically.
  expect_error(
    stabl_multiomic_cv(
      x_list = list(rna = x), y = y,
      lambda_grid = data.frame(lambda = 0.1), v = 2L,
      n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
      random_state = 3L, groups = groups, family = "binomial"
    ),
    "single outcome class"
  )

  # non-regression: mixed-composition groups with classification pass the
  # guard and still record the group-disjoint scheme
  y_mixed <- setNames(factor(rep(c("A", "B"), times = 4L)), ids)
  x_mixed <- x
  x_mixed[, 1:2] <- x_mixed[, 1:2] + 1.5 * (y_mixed == "A")
  cv_ok <- suppressWarnings(stabl_multiomic_cv(
    x_list = list(rna = x_mixed), y = y_mixed,
    lambda_grid = data.frame(lambda = 0.1), v = 2L,
    n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
    random_state = 3L, groups = groups, family = "binomial"
  ))
  expect_identical(cv_ok$fold_scheme, "group_disjoint")
})

test_that("zero-variance partitions with different means give Inf SMD (E9 regression)", {
  # Both partitions constant: differing constants are a maximal, perfectly
  # separated shift (true SMD = +-Inf) and must surface in smd_max_abs
  # rather than be masked as NA and counted as n_zero_sd.
  x_study <- matrix(5, 4, 2, dimnames = list(NULL, c("a", "b")))
  x_rest <- matrix(1, 6, 2, dimnames = list(NULL, c("a", "b")))
  x_rest[, "b"] <- 5  # equal constants: no shift, no variance information
  smd <- stablr:::.study_block_smd_abs(x_study, x_rest)
  expect_identical(unname(smd["a"]), Inf)
  expect_true(is.na(smd["b"]))

  x_list <- list(rna = rbind(x_study, x_rest))
  shift <- stablr:::.study_feature_shift(
    x_list, c(rep("s1", 4L), rep("s2", 6L)), c("s1", "s2")
  )
  expect_identical(shift$smd_max_abs[shift$study == "s1"], Inf)
  expect_identical(shift$n_zero_sd[shift$study == "s1"], 1L)  # only "b"
})

test_that(".study_outcome_summary fails closed on matrix-like y (A5 regression)", {
  y <- matrix(rnorm(6), 6, 1, dimnames = list(paste0("s", 1:6), "out"))
  expect_error(
    stablr:::.study_outcome_summary(y, rep(c("A", "B"), each = 3), c("A", "B")),
    "vector or factor"
  )
})

test_that("scheme_comparison rejects ambiguous mixed named/positional calls (A10 regression)", {
  set.seed(1)
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(36), 12, 3, dimnames = list(ids, paste0("f", 1:3)))
  y <- setNames(factor(rep(c("A", "B"), 6)), ids)
  sl <- setNames(rep(c("st1", "st2"), each = 6), ids)
  # (x, y, sl, y = y2): the named `y` wins the formal and the positional `y`
  # would silently re-bind to `study_labels` without the guard (review A10).
  expect_error(
    stabl_multiomic_scheme_comparison(list(m = x), y, sl, y = y),
    "Ambiguous call"
  )
  # A fully named duplicate still errors at match time, before the body runs.
  expect_error(
    stabl_multiomic_scheme_comparison(
      x_list = list(m = x), y = y, study_labels = sl, y = y
    ),
    "matched by multiple"
  )
  # Unambiguous mixed calls are unaffected: unnamed args before the first
  # named reserved formal bind exactly as written.
  expect_error(
    stabl_multiomic_scheme_comparison(list(m = x), y, study_labels = NULL),
    "`study_labels` is required"
  )
})

test_that("scheme_comparison validates x_list before study_labels (A10 regression)", {
  ids <- paste0("s", 1:12)
  y <- setNames(factor(rep(c("A", "B"), 6)), ids)
  sl <- setNames(rep(c("st1", "st2"), each = 6), ids)
  # Pre-fix this errored opaquely via rownames(x_list[[1]]) on an empty list.
  expect_error(
    stabl_multiomic_scheme_comparison(list(), y, sl),
    "`x_list` must be a non-empty named list"
  )
})

test_that(".study_feature_shift counts NA-data SMDs separately from zero-SD (A12 regression)", {
  ids <- paste0("s", 1:6)
  x <- matrix(c(5, 5, 5, 5, 5, 5,        # f1: constant everywhere -> 0/0
                NA, NA, NA, 1, 2, 3,     # f2: all-missing in study A
                1, 2, 3, 2, 3, 4),       # f3: ordinary
              6, 3, dimnames = list(ids, paste0("f", 1:3)))
  shift <- stablr:::.study_feature_shift(
    list(b = x),
    setNames(rep(c("A", "B"), each = 3), ids),
    c("A", "B")
  )
  row_a <- shift[shift$study == "A", ]
  expect_equal(row_a$n_zero_sd, 1)
  expect_equal(row_a$n_missing, 1)
})
