# Nested cross-validation helpers for multi-omic STABL workflows.

#' Multi-Omic STABL Nested Cross-Validation
#'
#' Runs a repeated outer cross-validation loop for generalisation assessment.
#' Within each outer training split, candidate STABL workflows are compared
#' using an inner cross-validation loop. The selected candidate is then refit
#' on the full outer-training split and evaluated once on the held-out outer
#' fold.
#'
#' @param x_list Named list of omic matrices or data frames with identical
#'   sample row names.
#' @param y Named factor or character vector aligned to `x_list` rows.
#' @param candidates Optional list of candidate definitions. Each candidate is
#'   a list with `name` and `blocks`; `blocks` names one or more omics from
#'   `x_list`. When `NULL`, one candidate is created for each omic plus one
#'   early-fusion candidate using all omics.
#' @param lambda_grid Either `"auto"`, a single lambda-grid data frame, or a
#'   named list of lambda grids. Named entries can match omic names or candidate
#'   names.
#' @param outer_v Number of outer folds per repeat.
#' @param outer_repeats Number of repeated outer CV rounds.
#' @param inner_v Number of inner folds used for candidate selection.
#' @param stratified Logical. When `TRUE` (default), outer and inner folds
#'   preserve class proportions as closely as possible.
#' @param strata Optional named vector used to stratify folds instead of `y`.
#'   Categorical values are used directly. Numeric values are binned into
#'   quantile groups before fold assignment. Cannot be combined with
#'   `study_labels`.
#' @param strata_bins Number of quantile bins for numeric `strata`.
#' @param study_labels Optional named vector assigning each sample to a study
#'   (cohort, site, or batch); names must cover exactly the sample IDs. When
#'   `NULL` (default), folds are assigned sample-wise exactly as in prior
#'   releases (`fold_scheme = "sample_iid"`). When supplied, the CV unit
#'   becomes the study (`fold_scheme = "study_disjoint"`): outer folds assign
#'   whole studies (with `outer_v` equal to the number of studies this is
#'   leave-one-study-out), inner candidate-selection folds are study-disjoint
#'   within each outer training split, and `stratified = TRUE` is reinterpreted
#'   at the study level through a deterministic largest-first greedy pass that
#'   minimises class-proportion imbalance across folds (sample-level
#'   stratification is impossible under study disjointness). Study
#'   disjointness concerns the CV folds only: within-training bootstrap
#'   resampling is untouched, so `bootstrap_strata` (sampling design) and SGL
#'   `feature_groups` (feature grouping) remain orthogonal to `study_labels`
#'   (CV unit).
#' @param preprocessing Optional fold-safe preprocessing specification:
#'   `NULL` (default; identity, legacy behaviour), a character vector of
#'   transforms for [stabl_preprocess_fit()] (`"zv"`, `"impute_median"`,
#'   `"center"`, `"scale"`), or an unprepped [recipes::recipe()] (requires the
#'   `recipes` package).  The transform is fitted anew on every inner- and
#'   outer-training split and applied unchanged to the corresponding
#'   held-out rows; a fitted `stabl_preprocessing` object is rejected because
#'   it would leak information across folds.  Recipe steps see the
#'   block-prefixed candidate columns (`block::feature`).
#' @param event_class For binary outcomes, the level treated as the event
#'   when `return_probabilities = TRUE`; defaults to the second level of `y`.
#' @param return_probabilities Logical; when `TRUE` and `y` is binary, each
#'   outer held-out prediction gains an `event_probability` column giving
#'   P(`event_class`) from the fold-trained selected-feature predictor, and
#'   the object records the `event_class`.  Hard-class columns are unchanged.
#' @param keep_fold_predictors Logical; when `TRUE`, each outer fold's fitted
#'   downstream predictor (with a `predict_event_prob()` closure for binary
#'   outcomes) and its fold-fitted preprocessing transform are retained in
#'   `fold_results`, enabling held-out permutation importance with the same
#'   fold-trained predictor.  Defaults to `FALSE` (legacy behaviour).
#' @param metric Candidate-selection metric, `"ber"` or `"accuracy"`.
#' @param family Passed to [stabl_fit()]. Defaults to `"multinomial"`.
#' @param n_bootstraps Passed to [stabl_fit()].
#' @param artificial_type Passed to [stabl_fit()].
#' @param hard_threshold Passed to [stabl_fit()].
#' @param n_lambda Passed to [stabl_fit()] when `lambda_grid = "auto"`.
#' @inheritParams stabl_fit
#' @param l1_ratio Passed to [stabl_fit()] when `lambda_grid = "auto"`.
#'   Use this with `base_learner = "elastic_net"` to generate alpha-aware
#'   train-fold grids.
#' @param workers Passed to [stabl_fit()] for bootstrap-level parallelism.
#' @param cv_workers Number of outer folds to evaluate in parallel. Uses
#'   `parallel::mclapply()` on Unix-like systems and falls back to sequential
#'   execution on Windows.
#' @param ... Additional arguments passed to [stabl_fit()].
#'
#' @note Cooperative fusion (`cooperative_fusion = TRUE`) is not supported in
#'   nested CV. Use [stabl_multiomic_train_validate()] or
#'   [stabl_multiomic_cv()] for cooperative workflows. The held-out metrics
#'   estimate predictive performance for the declared candidate workflow; they
#'   do not establish dataset-wide FDP control. Late-fusion stacking is not a
#'   candidate type in this helper.
#'
#' @return An object of class `"stabl_multiomic_nested_cv"` containing fold
#'   definitions, inner candidate diagnostics, outer held-out predictions,
#'   selected features, and aggregate performance. The object also records the
#'   sampling unit as `fold_scheme` (`"sample_iid"` or `"study_disjoint"`) and
#'   the validated `study_labels`; under `fold_scheme = "study_disjoint"` each
#'   outer/inner fold additionally carries `study_ids` (per-fold study
#'   membership) and `class_counts` (per-fold held-out class counts). When
#'   `study_labels` is supplied, `performance_by_study` holds per-study
#'   outer held-out accuracy and balanced error rate (pooled across folds and
#'   repeats; `NA` balanced error rate for studies whose pooled held-out
#'   predictions contain a single outcome class).
#'
#' @note Batch harmonization (e.g. ComBat) fitted on the full dataset before
#'   cross-validation leaks study-level information across folds and inflates
#'   performance estimates. Fit harmonization within training folds, or
#'   interpret all reported metrics as conditional on the harmonization step;
#'   study-disjoint folds do not repair leakage from pooled preprocessing.
#'
#' @examples
#' set.seed(3L)
#' ids <- paste0("s", 1:18)
#' y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
#' signal <- stats::model.matrix(~ y - 1)
#' x <- matrix(
#'   rnorm(18 * 6, sd = 0.3), 18, 6,
#'   dimnames = list(ids, paste0("g", 1:6))
#' )
#' x[, 1:3] <- x[, 1:3] + signal
#' nested <- suppressWarnings(stabl_multiomic_nested_cv(
#'   x_list = list(rna = x), y = y,
#'   candidates = list(list(name = "rna", blocks = "rna")),
#'   lambda_grid = data.frame(lambda = c(0.2, 0.1)),
#'   outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
#'   artificial_type = NULL, hard_threshold = 1e-9,
#'   sample_fraction = 1, random_state = 3L
#' ))
#' nested$performance$balanced_error_rate
#'
#' # Study-disjoint folds: each block of four samples is treated as one
#' # study; outer CV is leave-one-study-out.
#' set.seed(4L)
#' ids2 <- paste0("s", 1:20)
#' y2 <- setNames(factor(rep(c("A", "B"), times = 10L)), ids2)
#' x2 <- matrix(
#'   rnorm(20 * 6, sd = 0.3), 20, 6,
#'   dimnames = list(ids2, paste0("g", 1:6))
#' )
#' x2[, 1:2] <- x2[, 1:2] + stats::model.matrix(~ y2 - 1)
#' study_labels <- setNames(rep(paste0("study", 1:5), each = 4L), ids2)
#' nested2 <- suppressWarnings(stabl_multiomic_nested_cv(
#'   x_list = list(rna = x2), y = y2,
#'   candidates = list(list(name = "rna", blocks = "rna")),
#'   lambda_grid = data.frame(lambda = c(0.2, 0.1)),
#'   outer_v = 5L, inner_v = 2L, n_bootstraps = 3L,
#'   artificial_type = NULL, hard_threshold = 1e-9,
#'   sample_fraction = 1, random_state = 4L,
#'   study_labels = study_labels
#' ))
#' nested2$fold_scheme
#' nested2$outer_folds[[1L]]$study_ids
#' @export
stabl_multiomic_nested_cv <- function(
    x_list,
    y,
    candidates      = NULL,
    lambda_grid     = "auto",
    outer_v         = 5L,
    outer_repeats   = 1L,
    inner_v         = 5L,
    stratified      = TRUE,
    strata          = NULL,
    strata_bins     = 5L,
    metric          = c("ber", "accuracy"),
    family          = "multinomial",
    n_bootstraps    = 100L,
    artificial_type = "random_permutation",
    hard_threshold  = NULL,
    random_state    = NULL,
    n_lambda        = 30L,
    l1_ratio        = NULL,
    workers         = 1L,
    cv_workers      = 1L,
    study_labels    = NULL,
    preprocessing   = NULL,
    event_class     = NULL,
    return_probabilities = FALSE,
    keep_fold_predictors = FALSE,
    ...
) {
  metric <- match.arg(metric)
  outer_v <- .validate_scalar_integer_like(outer_v, "outer_v", min = 2L)
  outer_repeats <- .validate_scalar_integer_like(outer_repeats, "outer_repeats", min = 1L)
  inner_v <- .validate_scalar_integer_like(inner_v, "inner_v", min = 2L)
  strata_bins <- .validate_scalar_integer_like(strata_bins, "strata_bins", min = 1L)
  cv_workers <- .validate_scalar_integer_like(cv_workers, "cv_workers", min = 1L)
  workers <- .validate_scalar_integer_like(workers, "workers", min = 1L)
  n_lambda <- .validate_scalar_integer_like(n_lambda, "n_lambda", min = 1L)
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }
  validate_multiomic_inputs(x_list = x_list, y = y)

  y <- .subset_outcome_by_ids(y, rownames(x_list[[1L]]))
  y <- factor(y)
  sample_ids <- names(y)
  if (nlevels(y) < 2L) {
    stop("`y` must contain at least two outcome classes.", call. = FALSE)
  }
  .validate_preprocessing_spec(preprocessing)
  if (!isTRUE(return_probabilities) && !identical(return_probabilities, FALSE)) {
    stop("`return_probabilities` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(keep_fold_predictors) && !identical(keep_fold_predictors, FALSE)) {
    stop("`keep_fold_predictors` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(return_probabilities) && nlevels(y) != 2L) {
    stop(
      paste0(
        "`return_probabilities = TRUE` requires a binary outcome; ",
        "multiclass probability schemas are not part of this helper."
      ),
      call. = FALSE
    )
  }
  if (!is.null(event_class)) {
    if (nlevels(y) != 2L) {
      stop("`event_class` is only defined for binary outcomes.", call. = FALSE)
    }
    event_class <- as.character(event_class)
    if (!event_class %in% levels(y)) {
      stop("`event_class` must be one of the outcome levels.", call. = FALSE)
    }
  }
  strata_labels <- .make_nested_strata(
    strata = strata,
    y = y,
    sample_ids = sample_ids,
    bins = strata_bins
  )

  if (!is.null(study_labels) && !is.null(strata)) {
    stop(
      paste0(
        "`study_labels` and `strata` cannot be combined: study-disjoint ",
        "folds are assigned at the study level, so sample-level ",
        "stratification does not apply."
      ),
      call. = FALSE
    )
  }

  candidates <- .normalize_stabl_nested_candidates(candidates, names(x_list))

  if (is.null(study_labels)) {
    if (isTRUE(stratified) && min(table(strata_labels)) < max(outer_v, inner_v)) {
      stop("Each stratum must have at least `max(outer_v, inner_v)` samples.", call. = FALSE)
    }
    fold_scheme <- "sample_iid"
    outer_folds <- .make_repeated_cv_folds(
      y = strata_labels,
      v = outer_v,
      repeats = outer_repeats,
      stratified = stratified,
      random_state = random_state
    )
    inner_folds_pre <- NULL
  } else {
    study_labels <- .validate_study_labels(study_labels, sample_ids)
    n_studies <- length(unique(as.character(study_labels)))
    if (n_studies < max(outer_v, 2L)) {
      stop(sprintf(
        paste0(
          "`study_labels` defines %d %s; study-disjoint outer CV ",
          "requires at least max(outer_v, 2) = %d %s."
        ),
        n_studies,
        if (n_studies == 1L) "study" else "studies",
        max(outer_v, 2L),
        if (max(outer_v, 2L) == 1L) "study" else "studies"
      ), call. = FALSE)
    }
    fold_scheme <- "study_disjoint"
    outer_folds <- .make_repeated_study_folds(
      study_labels = study_labels,
      y = y,
      v = outer_v,
      repeats = outer_repeats,
      stratified = stratified,
      random_state = random_state
    )
    inner_folds_pre <- .prebuild_study_disjoint_inner_folds(
      outer_folds = outer_folds,
      study_labels = study_labels,
      y = y,
      inner_v = inner_v,
      stratified = stratified,
      random_state = random_state
    )
  }

  fold_runner <- function(outer_i) {
    .run_nested_cv_fold(
      outer_i         = outer_i,
      outer_folds     = outer_folds,
      strata_labels   = strata_labels,
      inner_v         = inner_v,
      stratified      = stratified,
      x_list          = x_list,
      y               = y,
      candidates      = candidates,
      lambda_grid     = lambda_grid,
      metric          = metric,
      family          = family,
      n_bootstraps    = n_bootstraps,
      artificial_type = artificial_type,
      hard_threshold  = hard_threshold,
      random_state    = random_state,
      n_lambda        = n_lambda,
      l1_ratio        = l1_ratio,
      workers         = workers,
      inner_folds     = if (is.null(inner_folds_pre)) NULL else inner_folds_pre[[outer_i]],
      preprocessing   = preprocessing,
      event_class     = event_class,
      return_probabilities = return_probabilities,
      keep_fold_predictors = keep_fold_predictors,
      ...
    )
  }

  if (cv_workers > 1L && .Platform$OS.type != "windows") {
    outer_results <- parallel::mclapply(
      seq_along(outer_folds),
      fold_runner,
      mc.cores = min(cv_workers, length(outer_folds))
    )
  } else {
    outer_results <- lapply(seq_along(outer_folds), fold_runner)
  }

  fold_results <- lapply(outer_results, `[[`, "fold_result")
  names(fold_results) <- vapply(outer_folds, `[[`, character(1L), "fold_id")
  inner_rows <- lapply(outer_results, `[[`, "diagnostics")
  prediction_rows <- lapply(outer_results, `[[`, "predictions")
  feature_rows <- lapply(outer_results, `[[`, "selected_features")

  predictions <- do.call(rbind, prediction_rows)
  diagnostics <- do.call(rbind, inner_rows)
  selected_features <- do.call(rbind, feature_rows)
  performance <- .classification_metrics_from_preds_df(predictions, y)
  performance_by_study <- if (!is.null(study_labels)) {
    .per_study_performance(predictions, study_labels, y)
  } else {
    NULL
  }

  structure(
    list(
      outer_folds = outer_folds,
      fold_results = fold_results,
      diagnostics = diagnostics,
      outer_predictions = predictions,
      selected_features = selected_features,
      performance = performance,
      performance_by_study = performance_by_study,
      levels = levels(y),
      metric = metric,
      candidates = candidates,
      strata = strata_labels,
      stratified = isTRUE(stratified),
      fold_scheme = fold_scheme,
      study_labels = study_labels,
      cv_workers = cv_workers,
      stabl_workers = workers,
      preprocessing = preprocessing,
      event_class = if (isTRUE(return_probabilities)) {
        if (is.null(event_class)) levels(y)[2L] else event_class
      } else {
        NULL
      }
    ),
    class = "stabl_multiomic_nested_cv"
  )
}

.normalize_stabl_nested_candidates <- function(candidates, omic_names) {
  if (is.null(candidates)) {
    out <- lapply(omic_names, function(block) list(name = block, blocks = block))
    names(out) <- omic_names
    if (length(omic_names) > 1L) {
      out$early_fusion <- list(name = "early_fusion", blocks = omic_names)
    }
    return(out)
  }

  if (!is.list(candidates) || length(candidates) == 0L) {
    stop("`candidates` must be a non-empty list.", call. = FALSE)
  }

  out <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    cand <- candidates[[i]]
    if (is.character(cand)) cand <- list(name = cand[[1L]], blocks = cand)
    if (is.null(cand$name) || is.null(cand$blocks)) {
      stop("Each candidate must define `name` and `blocks`.", call. = FALSE)
    }
    cand$name <- as.character(cand$name[[1L]])
    cand$blocks <- as.character(cand$blocks)
    if (!all(cand$blocks %in% omic_names)) {
      stop("Candidate `blocks` must name omics in `x_list`.", call. = FALSE)
    }
    out[[i]] <- cand
  }
  names(out) <- vapply(out, `[[`, character(1L), "name")
  if (anyDuplicated(names(out))) {
    stop("Candidate names must be unique.", call. = FALSE)
  }
  out
}

.make_nested_strata <- function(strata, y, sample_ids, bins = 5L) {
  if (is.null(strata)) {
    out <- factor(y)
    names(out) <- sample_ids
    return(out)
  }

  if (is.null(names(strata))) {
    if (length(strata) != length(sample_ids)) {
      stop("Unnamed `strata` must have the same length as `y`.", call. = FALSE)
    }
    names(strata) <- sample_ids
  }
  if (!all(sample_ids %in% names(strata))) {
    stop("`strata` must contain names for every sample in `y`.", call. = FALSE)
  }
  strata <- strata[sample_ids]

  if (anyNA(strata)) {
    stop("`strata` cannot contain missing values.", call. = FALSE)
  }

  if (is.numeric(strata) || is.integer(strata)) {
    if (bins < 2L) {
      stop("`strata_bins` must be at least 2 for numeric strata.", call. = FALSE)
    }
    probs <- seq(0, 1, length.out = bins + 1L)
    breaks <- unique(stats::quantile(as.numeric(strata), probs = probs,
                                     na.rm = TRUE, names = FALSE,
                                     type = 7))
    if (length(breaks) < 3L) {
      out <- factor(rep("all", length(strata)))
    } else {
      out <- cut(as.numeric(strata), breaks = breaks, include.lowest = TRUE,
                 ordered_result = TRUE)
    }
  } else {
    out <- factor(strata)
  }

  names(out) <- sample_ids
  out
}

.make_repeated_cv_folds <- function(y, v, repeats = 1L, stratified = TRUE,
                                    random_state = NULL) {
  rep_list <- lapply(seq_len(as.integer(repeats)), function(rep_i) {
    rep_seed  <- .derive_nested_seed(random_state, rep_i, 10L)
    rep_folds <- .make_cv_folds(y = y, v = v, stratified = stratified,
                                random_state = rep_seed)
    lapply(seq_along(rep_folds), function(fold_i) {
      fold             <- rep_folds[[fold_i]]
      fold[["repeat"]] <- rep_i
      fold$fold        <- paste0("Fold", fold_i)
      fold$fold_id     <- paste0("Repeat", rep_i, "_Fold", fold_i)
      fold
    })
  })
  unlist(rep_list, recursive = FALSE)
}

.make_cv_folds <- function(y, v, stratified = TRUE, random_state = NULL) {
  if (isTRUE(stratified)) {
    return(.make_stratified_folds(y = y, v = v, random_state = random_state))
  }
  .make_unstratified_folds(y = y, v = v, random_state = random_state)
}

.make_stratified_folds <- function(y, v, random_state = NULL) {
  y <- factor(y)
  v <- as.integer(v)
  if (min(table(y)) < v) {
    stop("Each class must have at least `v` samples for stratified folds.", call. = FALSE)
  }
  .with_local_seed(
    if (!is.null(random_state)) as.integer(random_state),
    {
      ids <- names(y)
      fold_ids <- integer(length(y))
      names(fold_ids) <- ids
      for (lvl in levels(y)) {
        class_ids <- sample(ids[y == lvl])
        fold_ids[class_ids] <- rep(seq_len(v), length.out = length(class_ids))
      }
      .build_fold_list(ids, fold_ids, v)
    }
  )
}

.make_unstratified_folds <- function(y, v, random_state = NULL) {
  ids <- names(y)
  v <- as.integer(v)
  if (length(ids) < v) {
    stop("The number of samples must be at least `v`.", call. = FALSE)
  }
  .with_local_seed(
    if (!is.null(random_state)) as.integer(random_state),
    {
      shuffled <- sample(ids)
      fold_ids <- rep(seq_len(v), length.out = length(shuffled))
      names(fold_ids) <- shuffled
      .build_fold_list(ids, fold_ids, v)
    }
  )
}

.build_fold_list <- function(ids, fold_ids, v) {
  lapply(seq_len(v), function(i) {
    valid_ids <- names(fold_ids)[fold_ids == i]
    list(
      train_ids = setdiff(ids, valid_ids),
      valid_ids = valid_ids,
      fold      = paste0("Fold", i)
    )
  })
}

.validate_study_labels <- function(study_labels, sample_ids) {
  if (is.null(study_labels)) return(NULL)
  if (is.null(names(study_labels))) {
    stop("`study_labels` must be a named vector; its names must be sample IDs.",
         call. = FALSE)
  }
  if (anyNA(study_labels)) {
    stop("`study_labels` cannot contain missing values.", call. = FALSE)
  }
  label_names <- names(study_labels)
  if (length(study_labels) != length(sample_ids) ||
      !setequal(label_names, sample_ids)) {
    missing_ids <- setdiff(sample_ids, label_names)
    unexpected_ids <- setdiff(label_names, sample_ids)
    stop(sprintf(
      paste0(
        "`study_labels` names must cover exactly the sample IDs ",
        "(missing: %s; unexpected: %s)."
      ),
      if (length(missing_ids)) paste(head(missing_ids, 5L), collapse = ", ") else "none",
      if (length(unexpected_ids)) paste(head(unexpected_ids, 5L), collapse = ", ") else "none"
    ), call. = FALSE)
  }
  study_labels <- as.character(study_labels[sample_ids])
  names(study_labels) <- sample_ids
  sizes <- table(study_labels)
  if (any(sizes < 2L)) {
    bad <- names(sizes)[sizes < 2L]
    stop(sprintf(
      paste0(
        "Each study must contain at least 2 samples; %s %s only %s. ",
        "A 1-sample study usually indicates an ID parsing mistake."
      ),
      paste(sprintf("'%s'", bad), collapse = ", "),
      if (length(bad) == 1L) "has" else "have",
      paste(head(as.integer(sizes[sizes < 2L]), 5L), collapse = ", ")
    ), call. = FALSE)
  }
  out <- factor(study_labels)
  names(out) <- sample_ids
  out
}

.make_repeated_study_folds <- function(study_labels, y, v, repeats = 1L,
                                       stratified = TRUE, random_state = NULL) {
  rep_list <- lapply(seq_len(as.integer(repeats)), function(rep_i) {
    rep_seed  <- .derive_nested_seed(random_state, rep_i, 10L)
    rep_folds <- .make_study_disjoint_folds(
      study_labels = study_labels,
      y = y,
      v = v,
      stratified = stratified,
      random_state = rep_seed
    )
    lapply(seq_along(rep_folds), function(fold_i) {
      fold             <- rep_folds[[fold_i]]
      fold[["repeat"]] <- rep_i
      fold$fold        <- paste0("Fold", fold_i)
      fold$fold_id     <- paste0("Repeat", rep_i, "_Fold", fold_i)
      fold
    })
  })
  folds <- unlist(rep_list, recursive = FALSE)
  .assert_study_fold_class_coverage(folds, context = "outer")
  folds
}

.make_study_disjoint_folds <- function(study_labels, y, v, stratified = TRUE,
                                       random_state = NULL) {
  study_labels <- study_labels[names(y)]
  if (anyNA(study_labels)) {
    stop("`study_labels` must be defined for every sample in `y`.", call. = FALSE)
  }
  studies <- sort(unique(as.character(study_labels)))
  n_studies <- length(studies)
  v <- as.integer(v)
  if (n_studies < v) {
    stop("The number of studies must be at least `v`.", call. = FALSE)
  }

  class_tab <- table(
    study = factor(as.character(study_labels), levels = studies),
    class = factor(as.character(y))
  )

  if (isTRUE(stratified)) {
    assignment <- .greedy_balanced_study_assignment(class_tab, v, random_state)
  } else {
    assignment <- .with_local_seed(
      if (!is.null(random_state)) as.integer(random_state),
      {
        shuffled <- sample(studies)
        setNames(rep(seq_len(v), length.out = n_studies), shuffled)
      }
    )
  }

  study_of_sample <- as.character(study_labels)
  lapply(seq_len(v), function(i) {
    fold_studies <- sort(names(assignment)[assignment == i])
    valid_ids <- names(y)[study_of_sample %in% fold_studies]
    list(
      train_ids = setdiff(names(y), valid_ids),
      valid_ids = valid_ids,
      fold      = paste0("Fold", i),
      study_ids = fold_studies,
      class_counts = table(y[valid_ids])
    )
  })
}

.greedy_balanced_study_assignment <- function(class_tab, v, random_state) {
  studies <- rownames(class_tab)
  overall_prop <- colSums(class_tab) / sum(class_tab)

  ordered_studies <- .with_local_seed(
    if (!is.null(random_state)) as.integer(random_state),
    {
      # Largest studies first; the seeded shuffle fixes the order of size
      # ties so the assignment is deterministic for a given `random_state`
      # and invariant to sample-row order.
      shuffled <- sample(studies)
      sizes <- as.integer(rowSums(class_tab))[shuffled]
      shuffled[order(-sizes, method = "radix")]
    }
  )

  fold_class_counts <- matrix(0L, nrow = v, ncol = ncol(class_tab))
  fold_sizes <- integer(v)
  fold_n_studies <- integer(v)
  assignment <- setNames(integer(length(ordered_studies)), ordered_studies)

  assign_study <- function(study, allowed) {
    study_counts <- as.integer(class_tab[study, ])
    study_size <- sum(study_counts)
    imbalance <- vapply(allowed, function(f) {
      new_size <- fold_sizes[f] + study_size
      sum(abs((fold_class_counts[f, ] + study_counts) / new_size - overall_prop))
    }, numeric(1L))
    best <- allowed[imbalance == min(imbalance)]
    if (length(best) > 1L) {
      best <- best[fold_sizes[best] == min(fold_sizes[best])]
      best <- best[fold_n_studies[best] == min(fold_n_studies[best])]
      best <- best[[1L]]
    }
    fold_class_counts[best, ] <<- fold_class_counts[best, ] + study_counts
    fold_sizes[best] <<- fold_sizes[best] + study_size
    fold_n_studies[best] <<- fold_n_studies[best] + 1L
    assignment[study] <<- best
  }

  # Pass 1: seed every fold with one of the largest studies (class-balanced
  # among the empty folds) so no validation fold can end up empty; with
  # `v == n_studies` this is the whole leave-one-study-out assignment.
  for (study in ordered_studies[seq_len(v)]) {
    assign_study(study, allowed = which(fold_n_studies == 0L))
  }
  # Pass 2: greedily assign the remaining studies by class-proportion
  # imbalance across folds, restricted to folds below the even-split
  # capacity (ceil(n_studies / v) studies per fold). Under outcome-
  # prevalence drift, unconstrained imbalance minimisation can concentrate
  # studies into few folds, which distorts fold sizes and can leave an
  # outer-training split with fewer than inner_v + 1 studies.
  if (length(ordered_studies) > v) {
    capacity <- ceiling(nrow(class_tab) / v)
    for (study in ordered_studies[-seq_len(v)]) {
      assign_study(study, allowed = which(fold_n_studies < capacity))
    }
  }
  assignment
}

.assert_study_fold_class_coverage <- function(folds, context = "outer") {
  for (fold in folds) {
    counts <- fold$class_counts
    if (sum(counts > 0L) < 2L) {
      fold_name <- if (is.null(fold$fold_id)) fold$fold else fold$fold_id
      stop(sprintf(
        paste0(
          "Study-disjoint %s validation fold `%s` contains a single outcome ",
          "class (held-out counts: %s; studies: %s); balanced error rate is ",
          "undefined for single-class held-out sets. Adjust `study_labels`, ",
          "reduce `outer_v`/`inner_v`, or keep `stratified = TRUE`."
        ),
        context,
        fold_name,
        paste(sprintf("%s=%d", names(counts), as.integer(counts)),
              collapse = ", "),
        paste(fold$study_ids, collapse = ", ")
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

.prebuild_study_disjoint_inner_folds <- function(outer_folds, study_labels, y,
                                                 inner_v, stratified,
                                                 random_state) {
  lapply(seq_along(outer_folds), function(outer_i) {
    outer <- outer_folds[[outer_i]]
    train_ids <- outer$train_ids
    train_study_labels <- study_labels[train_ids]
    n_train_studies <- length(unique(as.character(train_study_labels)))
    if (n_train_studies < inner_v + 1L) {
      stop(sprintf(
        paste0(
          "Outer fold `%s` trains on %d %s; study-disjoint inner ",
          "candidate selection requires at least inner_v + 1 = %d %s. ",
          "Reduce `inner_v` or supply more studies."
        ),
        outer$fold_id,
        n_train_studies,
        if (n_train_studies == 1L) "study" else "studies",
        inner_v + 1L,
        if (inner_v + 1L == 1L) "study" else "studies"
      ), call. = FALSE)
    }
    inner_folds <- .make_study_disjoint_folds(
      study_labels = train_study_labels,
      y = y[train_ids],
      v = inner_v,
      stratified = stratified,
      random_state = .derive_nested_seed(random_state, outer_i, 1000L)
    )
    inner_folds <- lapply(seq_along(inner_folds), function(fold_i) {
      fold <- inner_folds[[fold_i]]
      fold$fold <- paste0("Fold", fold_i)
      fold
    })
    .assert_study_fold_class_coverage(
      inner_folds,
      context = sprintf("inner (within outer fold `%s`)", outer$fold_id)
    )
    inner_folds
  })
}

.derive_nested_seed <- function(random_state, index, offset) {
  if (is.null(random_state)) return(NULL)
  # Compute in double before the modulo to avoid integer overflow when
  # `index` exceeds ~271500 (as.integer(index) * 7919L would return NA).
  as.integer((as.double(random_state) + as.double(offset) + as.double(index) * 7919) %%
               .Machine$integer.max)
}

.evaluate_stabl_candidates_inner <- function(x_list, y, train_ids, candidates,
                                             inner_folds, lambda_grid, metric,
                                             family, n_bootstraps,
                                             artificial_type, hard_threshold,
                                             random_state, n_lambda, l1_ratio,
                                             workers, preprocessing = NULL,
                                             ...) {
  candidate_rows <- list()
  k <- 1L

  for (cand_name in names(candidates)) {
    candidate <- candidates[[cand_name]]
    pred_rows <- list()

    for (fold_i in seq_along(inner_folds)) {
      fold <- inner_folds[[fold_i]]
      inner_train <- fold$train_ids
      inner_valid <- fold$valid_ids
      fit <- .fit_stabl_nested_candidate(
        x_list = x_list,
        y = y,
        train_ids = inner_train,
        valid_ids = inner_valid,
        candidate = candidate,
        lambda_grid = lambda_grid,
        family = family,
        n_bootstraps = n_bootstraps,
        artificial_type = artificial_type,
        hard_threshold = hard_threshold,
        random_state = .derive_nested_seed(random_state, fold_i, k * 100L),
        n_lambda = n_lambda,
        l1_ratio = l1_ratio,
        workers = workers,
        preprocessing = preprocessing,
        ...
      )
      pred_rows[[fold_i]] <- data.frame(
        truth = as.character(y[inner_valid]),
        predicted = fit$predicted,
        stringsAsFactors = FALSE
      )
    }

    preds <- do.call(rbind, pred_rows)
    metrics <- .classification_metrics_from_preds_df(preds, y)
    candidate_rows[[cand_name]] <- data.frame(
      candidate = cand_name,
      accuracy = metrics$accuracy,
      balanced_error_rate = metrics$balanced_error_rate,
      macro_f1 = metrics$macro_f1,
      n_inner_predictions = nrow(preds),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }

  list(summary = do.call(rbind, candidate_rows))
}

.select_nested_candidate <- function(summary, metric) {
  if (metric == "accuracy") {
    return(summary$candidate[which.max(summary$accuracy)])  # ties: first-index wins
  }
  summary$candidate[which.min(summary$balanced_error_rate)]  # ties: first-index wins
}

.fit_stabl_nested_candidate <- function(x_list, y, train_ids, valid_ids,
                                        candidate, lambda_grid, family,
                                        n_bootstraps, artificial_type,
                                        hard_threshold, random_state,
                                        n_lambda, l1_ratio, workers,
                                        preprocessing = NULL,
                                        event_class = NULL,
                                        keep_predictor = FALSE,
                                        ...) {
  x_train <- .candidate_matrix(x_list, candidate, train_ids)
  x_valid <- .candidate_matrix(x_list, candidate, valid_ids)
  y_train <- y[train_ids]

  # Fold-safe preprocessing: the transform is fitted on THIS fold's training
  # matrix only, then applied unchanged to the held-out rows.  Selection,
  # auto lambda grids and decoy construction all see transformed training
  # data; assessment values never influence the fitted parameters.
  pp <- .fit_fold_preprocessing(preprocessing, x_train)
  if (!is.null(pp)) {
    x_train <- stabl_preprocess_apply(pp, x_train)
    x_valid <- stabl_preprocess_apply(pp, x_valid)
  }

  fit <- stabl_fit(
    x = x_train,
    y = y_train,
    lambda_grid = .resolve_nested_lambda_grid(
      lambda_grid,
      candidate,
      x_train,
      y_train,
      family,
      n_lambda,
      l1_ratio
    ),
    family = family,
    n_bootstraps = n_bootstraps,
    artificial_type = artificial_type,
    hard_threshold = hard_threshold,
    random_state = random_state,
    n_lambda = n_lambda,
    l1_ratio = l1_ratio,
    workers = workers,
    ...
  )

  support <- get_support(fit)
  selected <- names(support)[support]
  importances <- get_importances(fit)
  predictor <- .fit_selected_predictor(
    x_train = x_train[, selected, drop = FALSE],
    y_train = y_train,
    x_valid = x_valid[, selected, drop = FALSE],
    levels = levels(y),
    event_class = event_class
  )

  list(
    stabl_fit = fit,
    selected_features = .split_prefixed_features(selected),
    importances = importances,
    predicted = predictor$predicted,
    probability = predictor$probability,
    predictor = if (isTRUE(keep_predictor)) predictor else NULL,
    preprocessing = pp,
    candidate = candidate$name
  )
}

.candidate_matrix <- function(x_list, candidate, ids) {
  mats <- lapply(candidate$blocks, function(block) {
    x <- as.matrix(x_list[[block]][ids, , drop = FALSE])
    colnames(x) <- paste(block, colnames(x), sep = "::")
    x
  })
  out <- do.call(cbind, mats)
  rownames(out) <- ids
  out
}

.resolve_nested_lambda_grid <- function(lambda_grid, candidate, x_train, y_train,
                                        family, n_lambda, l1_ratio = NULL) {
  if (identical(lambda_grid, "auto")) {
    return("auto")
  }
  if (is.data.frame(lambda_grid)) {
    return(lambda_grid)
  }
  if (is.list(lambda_grid)) {
    if (!is.null(lambda_grid[[candidate$name]])) {
      return(lambda_grid[[candidate$name]])
    }
    if (length(candidate$blocks) == 1L && !is.null(lambda_grid[[candidate$blocks[[1L]]]])) {
      return(lambda_grid[[candidate$blocks[[1L]]]])
    }
  }
  auto_lambda_grid(
    x_train,
    y_train,
    family = family,
    n_lambda = n_lambda,
    l1_ratio = l1_ratio
  )
}

.majority_class_prediction <- function(y_train, n) {
  rep(names(which.max(table(y_train))), n)
}

.classification_metrics_from_preds_df <- function(predictions, y) {
  .classification_metrics(
    truth     = factor(predictions$truth,     levels = levels(y)),
    predicted = factor(predictions$predicted, levels = levels(y))
  )
}

# Per-study outer held-out performance (C1): pools a study's predictions
# across all outer folds and repeats.  Balanced error rate is computed over
# the outcome classes actually present in the study's pooled predictions and
# is NA when only one class is present.
.per_study_performance <- function(predictions, study_labels, y) {
  if (is.null(predictions) || nrow(predictions) == 0L) {
    return(NULL)
  }
  study_of_prediction <- unname(study_labels[as.character(predictions$sample_id)])
  studies <- sort(unique(study_of_prediction))
  rows <- lapply(studies, function(study_id) {
    idx <- which(study_of_prediction == study_id)
    truth <- factor(predictions$truth[idx], levels = levels(y))
    predicted <- factor(predictions$predicted[idx], levels = levels(y))
    tab <- table(truth = truth, predicted = predicted)
    n <- sum(tab)
    accuracy <- unname(sum(diag(tab)) / n)
    recalls <- diag(tab) / pmax(rowSums(tab), 1L)
    present <- rowSums(tab) > 0L
    ber <- if (sum(present) >= 2L) {
      unname(1 - mean(recalls[present]))
    } else {
      NA_real_
    }
    data.frame(
      study = as.character(study_id),
      n_predictions = as.integer(n),
      accuracy = accuracy,
      balanced_error_rate = ber,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Fit the downstream selected-feature predictor: hard classes plus, for
# binary outcomes, event probabilities and the retained predictor itself (a
# model object and a `predict_event_prob(newx)` closure) so validated
# workflows can compute held-out permutation importance with the same
# fold-trained predictor.  `.predict_selected_multinomial()` remains the
# legacy hard-class wrapper with unchanged behaviour.
.fit_selected_predictor <- function(x_train, y_train, x_valid, levels,
                                    event_class = NULL) {
  y_train <- factor(y_train, levels = levels)
  n_sel <- ncol(x_train)
  binary <- length(levels) == 2L
  if (!is.null(event_class)) {
    if (!binary) {
      stop("`event_class` is only defined for binary outcomes.", call. = FALSE)
    }
    event_class <- as.character(event_class)
    if (!event_class %in% levels) {
      stop("`event_class` must be one of the outcome levels.", call. = FALSE)
    }
  }
  # Both downstream engines (stats::glm and glmnet) model P(second level)
  # for a two-level factor response; map to the requested event class.
  event <- if (binary) (event_class %||% levels[2L]) else NULL
  to_event <- if (binary && identical(event, levels[1L])) {
    function(p) 1 - p
  } else {
    function(p) p
  }

  # Training-derived fallback: majority hard class; binary event prevalence
  # (never a fixed 0.5) for the probability.
  prevalence_fallback <- function() {
    prevalence <- if (binary) mean(y_train == event) else NULL
    list(
      predicted = .majority_class_prediction(y_train, nrow(x_valid)),
      probability = if (binary) rep(prevalence, nrow(x_valid)) else NULL,
      model = NULL,
      predict_event_prob = if (binary) {
        function(newx) rep(prevalence, nrow(newx))
      } else {
        NULL
      },
      fallback = "training_prevalence"
    )
  }

  # Zero selected features, or a degenerate one-class training outcome.
  if (n_sel == 0L || length(unique(y_train)) < 2L) {
    return(prevalence_fallback())
  }

  # Exactly one selected feature: glmnet requires at least two columns, so
  # the downstream predictor is a genuine one-feature model.  Binary uses a
  # logit glm; multinomial uses a nearest-class-centroid discriminant rule
  # (per-class training means of the feature, midpoint decision boundaries),
  # which is the one-predictor LDA direction and needs no iterative fit.
  # Separation or other numerical failures fall back to the training-derived
  # majority class / event prevalence.
  if (n_sel == 1L) {
    return(tryCatch(
      {
        if (binary) {
          # suppressWarnings: complete separation warns ("fitted probabilities
          # numerically 0 or 1"); the resulting fit is still a valid
          # one-feature predictor and is preferred over the majority class.
          fit <- suppressWarnings(.with_linear_algebra_numerical_infeasibility(
            stats::glm(
              y_train ~ .x,
              data   = data.frame(.x = as.numeric(x_train[, 1L])),
              family = stats::binomial(link = "logit"),
              control = stats::glm.control(maxit = 50L)
            ),
            source = "one-feature binomial selected-feature fit",
            subclass = "stablr_downstream_numerical_infeasibility"
          ))
          prob <- unname(stats::predict(
            fit,
            newdata = data.frame(.x = as.numeric(x_valid[, 1L])),
            type = "response"
          ))
          list(
            predicted = ifelse(prob >= 0.5, levels[2L], levels[1L]),
            probability = to_event(prob),
            model = fit,
            predict_event_prob = function(newx) {
              to_event(unname(stats::predict(
                fit,
                newdata = data.frame(.x = as.numeric(newx[, 1L])),
                type = "response"
              )))
            },
            fallback = NA_character_
          )
        } else {
          xv <- as.numeric(x_train[, 1L])
          centroids <- vapply(
            levels,
            function(lv) mean(xv[y_train == lv]),
            numeric(1L)
          )
          if (any(!is.finite(centroids))) {
            .abort_numerical_infeasibility(
              "stablr_downstream_numerical_infeasibility",
              "One-feature multinomial centroid is not finite.",
              centroids = centroids
            )
          }
          x_new <- as.numeric(x_valid[, 1L])
          dist <- abs(outer(x_new, centroids, "-"))
          list(
            predicted = levels[max.col(-dist, ties.method = "first")],
            probability = NULL,
            model = list(centroids = centroids),
            predict_event_prob = NULL,
            fallback = NA_character_
          )
        }
      },
      stablr_numerical_infeasibility = function(e) {
        prevalence_fallback()
      }
    ))
  }

  tryCatch(
    {
      class_counts <- table(y_train)
      if (any(class_counts < 3L)) {
        .abort_numerical_infeasibility(
          "stablr_downstream_numerical_infeasibility",
          paste0(
            "Selected-feature classification requires at least three ",
            "training observations per class."
          ),
          training_class_counts = class_counts
        )
      }
      fit <- .with_glmnet_numerical_infeasibility(
        glmnet::cv.glmnet(
          x = x_train,
          y = y_train,
          family = if (binary) "binomial" else "multinomial",
          type.measure = "class",
          nfolds = min(5L, min(class_counts))
        ),
        source = "glmnet selected-feature classification fit"
      )
      predicted <- as.character(stats::predict(
        fit,
        newx = x_valid,
        s = "lambda.min",
        type = "class"
      )[, 1L])
      probability <- NULL
      predict_event_prob <- NULL
      if (binary) {
        probability <- to_event(as.numeric(stats::predict(
          fit,
          newx = x_valid,
          s = "lambda.min",
          type = "response"
        )[, 1L]))
        predict_event_prob <- function(newx) {
          to_event(as.numeric(stats::predict(
            fit,
            newx = newx,
            s = "lambda.min",
            type = "response"
          )[, 1L]))
        }
      }
      list(
        predicted = predicted,
        probability = probability,
        model = fit,
        predict_event_prob = predict_event_prob,
        fallback = NA_character_
      )
    },
    stablr_numerical_infeasibility = function(e) {
      prevalence_fallback()
    }
  )
}

.predict_selected_multinomial <- function(x_train, y_train, x_valid, levels) {
  .fit_selected_predictor(x_train, y_train, x_valid, levels)$predicted
}

.split_prefixed_features <- function(features) {
  out <- list()
  if (length(features) == 0L) return(out)
  parts <- strsplit(features, "::", fixed = TRUE)
  blocks <- vapply(parts, `[`, character(1L), 1L)
  names_only <- vapply(parts, function(x) paste(x[-1L], collapse = "::"), character(1L))
  split(names_only, blocks)
}

.stabl_nested_feature_table <- function(selected_features, importances, repeat_id,
                                        fold, fold_id, method, candidate) {
  block_names <- names(selected_features)
  if (length(block_names) == 0L) {
    return(data.frame(
      method = character(), candidate = character(), repeat_id = integer(),
      fold = character(), fold_id = character(), block = character(),
      feature = character(), score = numeric(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(block_names, function(block) {
    feats    <- selected_features[[block]]
    prefixed <- paste(block, feats, sep = "::")
    data.frame(
      method = method, candidate = candidate, repeat_id = repeat_id,
      fold = fold, fold_id = fold_id, block = block,
      feature = feats, score = unname(importances[prefixed]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.run_nested_cv_fold <- function(outer_i, outer_folds, strata_labels, inner_v,
                               stratified, x_list, y, candidates, lambda_grid,
                               metric, family, n_bootstraps, artificial_type,
                               hard_threshold, random_state, n_lambda, l1_ratio,
                               workers, inner_folds = NULL, preprocessing = NULL,
                               event_class = NULL, return_probabilities = FALSE,
                               keep_fold_predictors = FALSE, ...) {
  outer <- outer_folds[[outer_i]]
  train_ids <- outer$train_ids
  valid_ids <- outer$valid_ids

  if (is.null(inner_folds)) {
    inner_seed <- .derive_nested_seed(random_state, outer_i, 1000L)
    inner_folds <- .make_cv_folds(
      y = strata_labels[train_ids],
      v = inner_v,
      stratified = stratified,
      random_state = inner_seed
    )
  }

  inner_eval <- .evaluate_stabl_candidates_inner(
    x_list = x_list,
    y = y,
    train_ids = train_ids,
    candidates = candidates,
    inner_folds = inner_folds,
    lambda_grid = lambda_grid,
    metric = metric,
    family = family,
    n_bootstraps = n_bootstraps,
    artificial_type = artificial_type,
    hard_threshold = hard_threshold,
    random_state = .derive_nested_seed(random_state, outer_i, 2000L),
    n_lambda = n_lambda,
    l1_ratio = l1_ratio,
    workers = workers,
    preprocessing = preprocessing,
    ...
  )

  selected_name <- .select_nested_candidate(inner_eval$summary, metric)
  selected_candidate <- candidates[[selected_name]]

  outer_fit <- .fit_stabl_nested_candidate(
    x_list = x_list,
    y = y,
    train_ids = train_ids,
    valid_ids = valid_ids,
    candidate = selected_candidate,
    lambda_grid = lambda_grid,
    family = family,
    n_bootstraps = n_bootstraps,
    artificial_type = artificial_type,
    hard_threshold = hard_threshold,
    random_state = .derive_nested_seed(random_state, outer_i, 3000L),
    n_lambda = n_lambda,
    l1_ratio = l1_ratio,
    workers = workers,
    preprocessing = preprocessing,
    event_class = event_class,
    keep_predictor = keep_fold_predictors,
    ...
  )

  pred_df <- data.frame(
    repeat_id = outer[["repeat"]],
    fold = outer$fold,
    fold_id = outer$fold_id,
    sample_id = valid_ids,
    truth = as.character(y[valid_ids]),
    predicted = outer_fit$predicted,
    selected_candidate = selected_name,
    stringsAsFactors = FALSE
  )
  # Binary event probability for the named event class; rows stay aligned
  # with `sample_id`.  Only added on explicit request for binary outcomes.
  if (isTRUE(return_probabilities) && !is.null(outer_fit$probability)) {
    pred_df$event_probability <- outer_fit$probability
  }
  feature_df <- .stabl_nested_feature_table(
    selected_features = outer_fit$selected_features,
    importances = outer_fit$importances,
    repeat_id = outer[["repeat"]],
    fold = outer$fold,
    fold_id = outer$fold_id,
    method = "stablr",
    candidate = selected_name
  )

  inner_diag <- inner_eval$summary
  inner_diag$repeat_id <- outer[["repeat"]]
  inner_diag$fold <- outer$fold
  inner_diag$fold_id <- outer$fold_id

  list(
    fold_result = list(
      outer_fold = outer,
      inner_folds = inner_folds,
      inner_results = inner_eval,
      selected_candidate = selected_name,
      fit = outer_fit
    ),
    diagnostics = inner_diag,
    predictions = pred_df,
    selected_features = feature_df
  )
}

#' @describeIn stabl_multiomic_nested_cv Print a concise summary of a
#'   `stabl_multiomic_nested_cv` object; invisibly returns `x`.
#' @export
print.stabl_multiomic_nested_cv <- function(x, ...) {
  cat("<stabl_multiomic_nested_cv>\n")
  cat("  outer folds: ", length(x$outer_folds), "\n", sep = "")
  cat("  fold scheme: ",
      if (is.null(x$fold_scheme)) "sample_iid" else x$fold_scheme, "\n",
      sep = "")
  if (!is.null(x$study_labels)) {
    cat("  studies:     ", length(unique(as.character(x$study_labels))), "\n", sep = "")
  }
  cat("  candidates:  ", paste(names(x$candidates), collapse = ", "), "\n", sep = "")
  cat("  metric:      ", x$metric, "\n", sep = "")
  cat("  accuracy:    ", sprintf("%.3f", x$performance$accuracy), "\n", sep = "")
  cat("  BER:         ", sprintf("%.3f", x$performance$balanced_error_rate), "\n", sep = "")
  invisible(x)
}
