# B1 — study-disjoint nested CV (0.3.0 cross-study reliability round).
#
# Acceptance criteria follow derivations/b1-study-disjoint-nested-cv-design.md:
# 1. outer-fold study disjointness and full study coverage per repeat;
# 2. inner folds study-disjoint within each outer training split;
# 3. fail-closed validation errors;
# 4. determinism and study-order permutation stability;
# 5. study_labels = NULL reproduces the sample-i.i.d. folds byte-identically;
# 6. end-to-end smoke on synthetic multi-study data.
#
# Design deviation (recorded): the design sketch proposed a 3-study smoke
# example, but its own validation 3 (each outer-training split needs at least
# inner_v + 1 studies, inner_v >= 2) makes inner CV impossible with fewer
# than 4 studies. The smoke test therefore uses 4 studies.

.make_study_fixture <- function(study_sizes, n_features = 6L, seed = 1L,
                                signal_scale = 2) {
  set.seed(seed)
  ids <- character()
  y_chr <- character()
  study_chr <- character()
  for (i in seq_along(study_sizes)) {
    n_i <- study_sizes[[i]]
    ids_i <- paste0("s", length(ids) + seq_len(n_i))
    ids <- c(ids, ids_i)
    y_chr <- c(y_chr, rep(c("A", "B"), length.out = n_i))
    study_chr <- c(study_chr, rep(paste0("study", i), n_i))
  }
  y <- setNames(factor(y_chr), ids)
  study_labels <- setNames(study_chr, ids)
  x <- matrix(
    rnorm(length(ids) * n_features, sd = 0.3),
    length(ids), n_features,
    dimnames = list(ids, paste0("g", seq_len(n_features)))
  )
  signal <- stats::model.matrix(~ y - 1) * signal_scale
  x[, seq_len(ncol(signal))] <- x[, seq_len(ncol(signal)), drop = FALSE] + signal
  list(x_list = list(rna = x), y = y, study_labels = study_labels)
}

.fold_study_map <- function(folds) {
  out <- list()
  for (i in seq_along(folds)) {
    for (s in folds[[i]]$study_ids) out[[s]] <- i
  }
  out
}

test_that("greedy assignment keeps even fold sizes under prevalence drift", {
  # heterogeneous study compositions: unconstrained imbalance minimisation
  # would concentrate studies into few folds (e.g. 4/1/1 for 6 studies)
  set.seed(17L)
  n_per <- 10L
  ids <- paste0("s", seq_len(6L * n_per))
  prevalences <- c(0.15, 0.3, 0.5, 0.7, 0.85, 0.5)
  y_chr <- unlist(lapply(prevalences, function(p) {
    sample(c("A", "B"), n_per, replace = TRUE, prob = c(1 - p, p))
  }))
  y <- setNames(factor(y_chr, levels = c("A", "B")), ids)
  study_labels <- setNames(
    rep(paste0("study", 1:6), each = n_per), ids
  )
  folds <- .make_repeated_study_folds(
    study_labels, y, v = 3L, repeats = 1L, stratified = TRUE,
    random_state = 61L
  )
  fold_sizes <- vapply(folds, function(f) length(f$study_ids), integer(1L))
  # capacity = ceil(6 / 3) = 2 studies per fold; no empty folds
  expect_true(all(fold_sizes <= 2L))
  expect_true(all(fold_sizes >= 1L))
  expect_equal(sum(fold_sizes), 6L)
  # every outer-training split keeps inner_v + 1 = 3 studies feasible
  for (fold in folds) {
    train_studies <- unique(as.character(study_labels[fold$train_ids]))
    expect_true(length(train_studies) >= 3L)
  }
})

test_that("outer folds are study-disjoint and cover all studies per repeat", {
  fixture <- .make_study_fixture(rep(4L, 5L), seed = 9L)
  res <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 5L,
    outer_repeats = 2L,
    inner_v = 2L,
    n_bootstraps = 5L,
    artificial_type = NULL,
    hard_threshold = 1e-9,
    sample_fraction = 1,
    random_state = 21L,
    study_labels = fixture$study_labels
  ))

  expect_identical(res$fold_scheme, "study_disjoint")
  expect_s3_class(res, "stabl_multiomic_nested_cv")
  expect_length(res$outer_folds, 10L)
  expect_true(all(grepl("^Repeat[0-9]+_Fold[0-9]+$",
                        vapply(res$outer_folds, `[[`, character(1L), "fold_id"))))

  all_studies <- sort(unique(as.character(fixture$study_labels)))
  for (rep_i in 1:2) {
    rep_folds <- res$outer_folds[
      vapply(res$outer_folds, `[[`, integer(1L), "repeat") == rep_i
    ]
    valid_studies <- unlist(lapply(rep_folds, `[[`, "study_ids"))
    expect_setequal(valid_studies, all_studies)
    expect_equal(anyDuplicated(valid_studies), 0L)
    for (fold in rep_folds) {
      # no shared study between the two sides of a split
      train_studies <- unique(as.character(fixture$study_labels[fold$train_ids]))
      valid_studies_fold <- unique(as.character(fixture$study_labels[fold$valid_ids]))
      expect_length(intersect(train_studies, valid_studies_fold), 0L)
      # recorded membership matches the held-out samples
      expect_setequal(fold$study_ids, valid_studies_fold)
      expect_setequal(
        fold$valid_ids,
        names(fixture$y)[as.character(fixture$study_labels) %in% fold$study_ids]
      )
      # per-fold class counts are recorded and match the held-out truth
      expect_equal(
        as.integer(fold$class_counts),
        as.integer(table(fixture$y[fold$valid_ids]))
      )
      # every validation fold holds both outcome classes
      expect_true(sum(fold$class_counts > 0L) >= 2L)
    }
  }
})

test_that("inner folds are study-disjoint within each outer training split", {
  fixture <- .make_study_fixture(rep(4L, 5L), seed = 10L)
  res <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 5L,
    inner_v = 2L,
    n_bootstraps = 5L,
    artificial_type = NULL,
    hard_threshold = 1e-9,
    sample_fraction = 1,
    random_state = 22L,
    study_labels = fixture$study_labels
  ))

  for (fold_result in res$fold_results) {
    inner_folds <- fold_result$inner_folds
    expect_length(inner_folds, 2L)
    inner_valid_studies <- unlist(lapply(inner_folds, `[[`, "study_ids"))
    # inner validation folds partition the outer-training studies
    train_studies <- unique(as.character(
      fixture$study_labels[fold_result$outer_fold$train_ids]
    ))
    expect_setequal(inner_valid_studies, train_studies)
    expect_equal(anyDuplicated(inner_valid_studies), 0L)
    for (fold in inner_folds) {
      train_studies_fold <- unique(as.character(
        fixture$study_labels[fold$train_ids]
      ))
      valid_studies_fold <- unique(as.character(
        fixture$study_labels[fold$valid_ids]
      ))
      expect_length(intersect(train_studies_fold, valid_studies_fold), 0L)
    }
  }
})

test_that("study_labels validation fails closed with actionable errors", {
  fixture <- .make_study_fixture(rep(4L, 4L), seed = 11L)
  y <- fixture$y
  study_labels <- fixture$study_labels

  # unnamed labels
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = fixture$x_list, y = y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 4L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = unname(study_labels)
    )),
    "named vector"
  )

  # names do not cover the sample IDs (one missing)
  incomplete <- study_labels[names(study_labels) != names(study_labels)[1L]]
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = fixture$x_list, y = y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 4L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = incomplete
    )),
    "exactly the sample IDs"
  )

  # names carry an unexpected sample ID
  extra <- c(study_labels, setNames("study1", "not_a_sample"))
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = fixture$x_list, y = y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 4L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = extra
    )),
    "exactly the sample IDs"
  )

  # a 1-sample study is rejected
  one_sample <- study_labels
  one_sample[[1L]] <- "studyX"
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = fixture$x_list, y = y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 4L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = one_sample
    )),
    "at least 2 samples"
  )

  # fewer studies than outer folds
  two_study_fixture <- .make_study_fixture(rep(8L, 2L), seed = 12L)
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = two_study_fixture$x_list, y = two_study_fixture$y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 3L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = two_study_fixture$study_labels
    )),
    "requires at least max\\(outer_v, 2\\)"
  )

  # 2 studies with inner_v = 2: the training split would hold 1 study
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = two_study_fixture$x_list, y = two_study_fixture$y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 2L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = two_study_fixture$study_labels
    )),
    "inner_v \\+ 1"
  )

  # a single-class validation fold is rejected with fold composition
  unbalanced <- .make_study_fixture(rep(4L, 4L), seed = 13L)
  first_study <- unique(as.character(unbalanced$study_labels))[[1L]]
  unbalanced$study_labels[unbalanced$study_labels == first_study] <- first_study
  unbalanced$y[unbalanced$study_labels == first_study] <- "A"
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = unbalanced$x_list, y = unbalanced$y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 4L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = unbalanced$study_labels
    )),
    "single outcome class"
  )

  # sample-level strata cannot be combined with study labels
  strata <- setNames(rep(c("lo", "hi"), length.out = length(y)), names(y))
  expect_error(
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = fixture$x_list, y = y,
      candidates = list(list(name = "rna", blocks = "rna")),
      lambda_grid = data.frame(lambda = c(0.2, 0.1)),
      outer_v = 4L, inner_v = 2L, n_bootstraps = 5L,
      artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
      random_state = 1L, study_labels = study_labels, strata = strata
    )),
    "cannot be combined"
  )
})

test_that("study-disjoint folds are deterministic and permutation-stable", {
  fixture <- .make_study_fixture(rep(4L, 5L), seed = 14L)

  f1 <- .make_repeated_study_folds(
    fixture$study_labels, fixture$y, v = 5L, repeats = 1L,
    stratified = TRUE, random_state = 31L
  )
  f2 <- .make_repeated_study_folds(
    fixture$study_labels, fixture$y, v = 5L, repeats = 1L,
    stratified = TRUE, random_state = 31L
  )
  expect_identical(f1, f2)

  # identical-composition studies: permuting sample rows must not change the
  # study-to-fold assignment
  ord <- sample(names(fixture$y))
  f3 <- .make_repeated_study_folds(
    fixture$study_labels[ord], fixture$y[ord], v = 5L, repeats = 1L,
    stratified = TRUE, random_state = 31L
  )
  expect_identical(.fold_study_map(f1), .fold_study_map(f3))

  # unstratified assignment is deterministic too and stays study-disjoint
  u1 <- .make_repeated_study_folds(
    fixture$study_labels, fixture$y, v = 5L, repeats = 1L,
    stratified = FALSE, random_state = 32L
  )
  u2 <- .make_repeated_study_folds(
    fixture$study_labels, fixture$y, v = 5L, repeats = 1L,
    stratified = FALSE, random_state = 32L
  )
  expect_identical(u1, u2)
  valid_studies <- unlist(lapply(u1, `[[`, "study_ids"))
  expect_setequal(valid_studies, unique(as.character(fixture$study_labels)))
  expect_equal(anyDuplicated(valid_studies), 0L)
})

test_that("study_labels = NULL reproduces the sample-iid folds byte-identically", {
  fixture <- .make_study_fixture(rep(6L, 3L), seed = 15L)
  res <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L,
    inner_v = 3L,
    n_bootstraps = 3L,
    artificial_type = NULL,
    hard_threshold = 1e-9,
    sample_fraction = 1,
    random_state = 41L
  ))

  expected <- .make_repeated_cv_folds(
    y = fixture$y, v = 3L, repeats = 1L, stratified = TRUE, random_state = 41L
  )
  expect_identical(res$outer_folds, expected)
  expect_identical(
    paste(capture.output(str(res$outer_folds)), collapse = "\n"),
    paste(capture.output(str(expected)), collapse = "\n")
  )
  expect_identical(res$fold_scheme, "sample_iid")
  expect_null(res$study_labels)
})

test_that("end-to-end smoke: 4 synthetic studies report the study-disjoint scheme", {
  fixture <- .make_study_fixture(c(8L, 8L, 7L, 7L), seed = 16L)  # n = 30
  res <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 4L,
    inner_v = 2L,
    n_bootstraps = 10L,
    artificial_type = NULL,
    hard_threshold = 1e-9,
    sample_fraction = 1,
    random_state = 51L,
    study_labels = fixture$study_labels
  ))

  expect_identical(res$fold_scheme, "study_disjoint")
  expect_identical(
    sort(unique(as.character(res$study_labels))),
    sort(unique(as.character(fixture$study_labels)))
  )
  expect_true(is.finite(res$performance$balanced_error_rate))
  expect_output(print(res), "study_disjoint")
  expect_output(print(res), "studies:\\s+4")
})

# C1 (inferred scope): per-study performance reporting.
test_that("performance_by_study pools per-study held-out predictions", {
  fixture <- .make_study_fixture(c(8L, 8L, 8L, 8L, 8L, 8L), seed = 16L)  # n = 48
  # study6 is single-class: its pooled held-out BER is undefined (NA), but
  # pairing it with a mixed study in every fold keeps all validation folds
  # two-class so the CV itself runs.
  ids6 <- names(fixture$study_labels)[fixture$study_labels == "study6"]
  fixture$y[ids6] <- factor(rep("A", length(ids6)), levels = levels(fixture$y))

  res <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L,
    inner_v = 2L,
    n_bootstraps = 10L,
    artificial_type = NULL,
    hard_threshold = 1e-9,
    sample_fraction = 1,
    random_state = 51L,
    study_labels = fixture$study_labels
  ))

  pbs <- res$performance_by_study
  expect_false(is.null(pbs))
  expect_equal(nrow(pbs), 6L)
  expect_setequal(pbs$study, paste0("study", 1:6))
  expect_equal(sum(pbs$n_predictions), nrow(res$outer_predictions))
  expect_true(all(pbs$accuracy >= 0 & pbs$accuracy <= 1))
  expect_true(is.na(pbs$balanced_error_rate[pbs$study == "study6"]))
  expect_true(all(is.finite(
    pbs$balanced_error_rate[pbs$study != "study6"]
  )))

  # sample-iid arm: the field stays NULL (additive, zero behavior change)
  res_iid <- suppressWarnings(stabl_multiomic_nested_cv(
    x_list = fixture$x_list,
    y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 3L,
    inner_v = 2L,
    n_bootstraps = 10L,
    artificial_type = NULL,
    hard_threshold = 1e-9,
    sample_fraction = 1,
    random_state = 51L
  ))
  expect_null(res_iid$performance_by_study)
})

# ── PR3: preprocessing under study-disjoint folds ────────────────────────────

test_that("held-out studies never enter fold-fitted preprocessing", {
  fixture <- .make_study_fixture(c(6L, 6L, 6L, 6L), seed = 11L)
  args <- list(
    x_list = fixture$x_list, y = fixture$y,
    candidates = list(list(name = "rna", blocks = "rna")),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    outer_v = 4L, inner_v = 2L, n_bootstraps = 3L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_fraction = 1,
    random_state = 12L,
    study_labels = fixture$study_labels,
    preprocessing = c("center", "scale"),
    keep_fold_predictors = TRUE
  )
  fit <- suppressWarnings(do.call(stabl_multiomic_nested_cv, args))

  # Shift every feature of the first fold's held-out study by a large
  # constant: fold-1's training-fitted transform must be bit-identical.
  fold1_valid <- fit$outer_folds[[1L]]$valid_ids
  x2 <- fixture$x_list
  x2$rna[fold1_valid, ] <- x2$rna[fold1_valid, ] + 1e6
  args2 <- args
  args2$x_list <- x2
  fit2 <- suppressWarnings(do.call(stabl_multiomic_nested_cv, args2))

  pp1 <- fit$fold_results[[1L]]$fit$preprocessing
  pp2 <- fit2$fold_results[[1L]]$fit$preprocessing
  expect_identical(pp1$params$center$means, pp2$params$center$means)
  expect_identical(pp1$params$scale$sds, pp2$params$scale$sds)
  expect_identical(
    fit$fold_results[[1L]]$fit$selected_features,
    fit2$fold_results[[1L]]$fit$selected_features
  )
})
