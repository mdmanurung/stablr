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
