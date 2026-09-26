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
#'   quantile groups before fold assignment.
#' @param strata_bins Number of quantile bins for numeric `strata`.
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
#'   execution on Windows.  Must be `1L` unless `random_state` is supplied:
#'   forked workers inherit an identical `.Random.seed` snapshot, so with no
#'   seed every fold would draw from one common stream and the result would be
#'   unreproducible.
#' @param ... Additional arguments passed to [stabl_fit()].
#'
#' @section Candidate selection is optimistic:
#' The held-out numbers are the score of the *declared candidate workflow*, not
#' of a candidate chosen without looking at the data.  Within each outer fold
#' the winner is picked by a single inner-CV number (`which.min` balanced error
#' rate, or `which.max` accuracy -- ties go to the first candidate), so the
#' outer estimate inherits the winner's-curse bias of that selection step; no
#' one-standard-error rule is applied.  This is a scientific decision for the
#' analysis author, not a defect to be smoothed over here.
#'
#' `outer_repeats` defaults to `1L`, which leaves the outer estimate with the
#' variance of a single repeated-CV split.  Use `outer_repeats > 1L` whenever a
#' generalisation number from this helper is going to be reported or compared.
#'
#' The inner-CV comparison is *paired*: every candidate is fitted on the same
#' inner folds with the same per-fold seeds, so a difference in inner score
#' reflects a difference in the candidate rather than a difference in the
#' bootstrap draws.
#'
#' @section RNG:
#' When `random_state` is an integer, `stabl_multiomic_nested_cv()` does not
#' mutate the caller's `.Random.seed`.  Every seeded stage -- fold construction,
#' each inner fit, each outer refit, and the downstream `cv.glmnet` -- runs
#' inside a `.with_local_seed()` scope.  With `random_state = NULL` the fold
#' draws come from ambient RNG state and therefore do advance the caller's
#' stream, as any unseeded stochastic routine does.
#'
#' Reproducibility is conditional on `RNGkind()` and on `sample.kind`; neither is
#' pinned by this function, so changing them between runs changes the folds.
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
#'   selected features, and aggregate performance.
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
  strata_labels <- .make_nested_strata(
    strata = strata,
    y = y,
    sample_ids = sample_ids,
    bins = strata_bins
  )

  if (isTRUE(stratified) && min(table(strata_labels)) < max(outer_v, inner_v)) {
    stop("Each stratum must have at least `max(outer_v, inner_v)` samples.", call. = FALSE)
  }

  candidates <- .normalize_stabl_nested_candidates(candidates, names(x_list))
  outer_folds <- .make_repeated_cv_folds(
    y = strata_labels,
    v = outer_v,
    repeats = outer_repeats,
    stratified = stratified,
    random_state = random_state
  )

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
      ...
    )
  }

  if (cv_workers > 1L && .Platform$OS.type != "windows") {
    if (is.null(random_state)) {
      .abort_stablr(
        "stablr_nested_cv_unseeded_parallel",
        paste0(
          "`cv_workers > 1` requires `random_state`. Forked workers inherit an ",
          "identical `.Random.seed` snapshot, so without a seed every outer ",
          "fold would draw from one common RNG stream: the folds would not be ",
          "independent and the result would not be reproducible. Supply ",
          "`random_state`, or use `cv_workers = 1L`."
        ),
        cv_workers = cv_workers
      )
    }
    outer_results <- parallel::mclapply(
      seq_along(outer_folds),
      fold_runner,
      mc.cores = min(cv_workers, length(outer_folds))
    )
  } else {
    outer_results <- lapply(seq_along(outer_folds), fold_runner)
  }

  # mclapply returns a `try-error` object for a failed worker instead of
  # propagating the error.  Without this check the failure surfaces as
  # `subscript out of bounds` from the `[[` below, which names neither the fold
  # nor the cause.
  .check_nested_cv_fold_results(outer_results, outer_folds)

  fold_results <- lapply(outer_results, `[[`, "fold_result")
  names(fold_results) <- vapply(outer_folds, `[[`, character(1L), "fold_id")
  inner_rows <- lapply(outer_results, `[[`, "diagnostics")
  prediction_rows <- lapply(outer_results, `[[`, "predictions")
  feature_rows <- lapply(outer_results, `[[`, "selected_features")

  predictions <- do.call(rbind, prediction_rows)
  diagnostics <- do.call(rbind, inner_rows)
  selected_features <- do.call(rbind, feature_rows)
  performance <- .classification_metrics_from_preds_df(predictions, y)

  structure(
    list(
      outer_folds = outer_folds,
      fold_results = fold_results,
      diagnostics = diagnostics,
      outer_predictions = predictions,
      selected_features = selected_features,
      performance = performance,
      levels = levels(y),
      metric = metric,
      candidates = candidates,
      strata = strata_labels,
      stratified = isTRUE(stratified),
      cv_workers = cv_workers,
      stabl_workers = workers
    ),
    class = "stabl_multiomic_nested_cv"
  )
}

# Raise a typed condition naming the fold that failed, instead of letting
# `lapply(results, "[[", "fold_result")` fail with an unrelated
# "subscript out of bounds" when a worker returned a `try-error`.
.check_nested_cv_fold_results <- function(results, folds) {
  failed <- vapply(results, function(result) inherits(result, "try-error"),
                   logical(1L))
  if (!any(failed)) {
    return(invisible(results))
  }
  index <- which(failed)[[1L]]
  condition <- results[[index]]
  fold_id <- tryCatch(
    vapply(folds, `[[`, character(1L), "fold_id")[[index]],
    error = function(e) NA_character_
  )
  # mclapply stores the worker's `try-error` as a character string, not as a
  # condition, so conditionMessage() only applies to the second case.
  detail <- if (inherits(condition, "condition")) {
    conditionMessage(condition)
  } else {
    paste(as.character(condition), collapse = " ")
  }
  .abort_stablr(
    "stablr_nested_cv_fold_failed",
    paste0(
      "Nested-CV outer fold ", index, " (", fold_id, ") failed: ", detail
    ),
    parent = if (inherits(condition, "condition")) condition else NULL,
    fold_index = index,
    fold_id = fold_id
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

.derive_nested_seed <- function(random_state, index, offset) {
  if (is.null(random_state)) return(NULL)
  # Compute in double before the modulo to avoid integer overflow when
  # `index` exceeds ~271500 (as.integer(index) * 7919L would return NA).
  as.integer((as.double(random_state) + as.double(offset) + as.double(index) * 7919) %%
               .Machine$integer.max)
}

# Seed offsets for the four derived streams in one nested-CV run.  They must stay
# distinct: 1000L drives inner-fold construction, 2000L drives every candidate's
# inner fit (independently of which candidate it is -- common random numbers),
# 3000L drives the outer refit, and 4000L drives the downstream `cv.glmnet`.
.NESTED_INNER_FOLD_SEED_OFFSET <- 1000L
.NESTED_INNER_FIT_SEED_OFFSET  <- 2000L
.NESTED_OUTER_FIT_SEED_OFFSET  <- 3000L
.NESTED_CV_GLMNET_SEED_OFFSET  <- 4000L

# `glmnet::cv.glmnet()` generates its own internal fold assignment, so it
# consumes ambient RNG state and its result changes with whatever ran before it.
# Wrapping the call in `.with_local_seed()` gives it a deterministic fold
# assignment derived from the caller's `random_state` and leaves the caller's
# `.Random.seed` untouched.  With `random_state = NULL` the call is left
# unseeded, which is the honest behaviour for an unseeded run.
#
# This helper exists so the three unseeded `cv.glmnet()` sites in the package
# converge on one implementation.  Sites outside `R/nested_cv.R` should call it
# with the seed already threaded through their own arguments.
.seeded_cv_glmnet <- function(x, y, family, type.measure, nfolds,
                              random_state = NULL, ...) {
  .with_local_seed(
    random_state,
    glmnet::cv.glmnet(
      x = x,
      y = y,
      family = family,
      type.measure = type.measure,
      nfolds = nfolds,
      ...
    )
  )
}

.evaluate_stabl_candidates_inner <- function(x_list, y, train_ids, candidates,
                                             inner_folds, lambda_grid, metric,
                                             family, n_bootstraps,
                                             artificial_type, hard_threshold,
                                             random_state, n_lambda, l1_ratio,
                                             workers, ...) {
  candidate_rows <- list()
  grid_sources <- character(0)

  for (cand_name in names(candidates)) {
    candidate <- candidates[[cand_name]]
    pred_rows <- list()

    for (fold_i in seq_along(inner_folds)) {
      fold <- inner_folds[[fold_i]]
      inner_train <- fold$train_ids
      inner_valid <- fold$valid_ids
      # The seed depends on the inner fold only, never on the candidate index:
      # all candidates are fitted on the same folds with the same bootstrap
      # draws (common random numbers).  Deriving it from the candidate index
      # injected sampling noise into a single-number selection rule, and
      # `k * 100L` also collided with the 1000L and 3000L offsets used for the
      # inner-fold and outer-refit seeds once a run had ten candidates.
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
        random_state = .derive_nested_seed(
          random_state, fold_i, .NESTED_INNER_FIT_SEED_OFFSET
        ),
        n_lambda = n_lambda,
        l1_ratio = l1_ratio,
        workers = workers,
        ...
      )
      grid_sources <- c(grid_sources, fit$lambda_grid_source)
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
  }

  summary <- do.call(rbind, candidate_rows)
  # A named-list `lambda_grid` with no entry for some candidate falls back to an
  # auto grid on that candidate's training data.  The per-candidate outcome is
  # recorded here so the substitution is visible in the returned diagnostics
  # rather than only in a warning.  The resolution depends on the grid
  # specification and the candidate, not on the fold, so the first fold's label
  # is the candidate's label.
  summary$lambda_grid_source <- names(candidate_rows)
  summary$lambda_grid_resolved <- unname(grid_sources[
    seq.int(1L, length(grid_sources), by = length(inner_folds))
  ])
  list(summary = summary)
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
                                        n_lambda, l1_ratio, workers, ...) {
  x_train <- .candidate_matrix(x_list, candidate, train_ids)
  x_valid <- .candidate_matrix(x_list, candidate, valid_ids)
  y_train <- y[train_ids]

  # Shared resolver: the named-list lookup and its `auto_lambda_grid` fallback
  # are defined once in learner_adapters.R, and the fallback now warns.
  resolved_grid <- .resolve_lambda_grid(
    lambda_grid,
    keys = .nested_lambda_grid_keys(candidate),
    x = x_train, y = y_train,
    family = family,
    n_lambda = n_lambda,
    l1_ratio = l1_ratio,
    on_missing_key = "warn"
  )

  fit <- stabl_fit(
    x = x_train,
    y = y_train,
    lambda_grid = resolved_grid$grid,
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
  predicted <- .predict_selected_multinomial(
    x_train = x_train[, selected, drop = FALSE],
    y_train = y_train,
    x_valid = x_valid[, selected, drop = FALSE],
    levels = levels(y),
    random_state = .derive_nested_seed(
      random_state, 1L, .NESTED_CV_GLMNET_SEED_OFFSET
    )
  )

  list(
    stabl_fit = fit,
    selected_features = .split_prefixed_features(selected),
    importances = importances,
    predicted = predicted,
    candidate = candidate$name,
    lambda_grid_source = resolved_grid$source
  )
}

# Ordered lookup keys for a named-list `lambda_grid`: the candidate name first,
# then its single block name (the two cases the resolver can satisfy).  A
# multi-block candidate only matches on its own name.
.nested_lambda_grid_keys <- function(candidate) {
  keys <- candidate$name
  if (length(candidate$blocks) == 1L) {
    keys <- c(keys, candidate$blocks[[1L]])
  }
  keys
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

# `.resolve_nested_lambda_grid()` used to live here.  The named-list lookup, the
# "auto" expansion and the `auto_lambda_grid()` fallback now have one
# implementation, `.resolve_lambda_grid()` in `R/learner_adapters.R`, which also
# validates the resolved grid and warns on the fallback path.  The fallback
# outcome is recorded per candidate in the returned `diagnostics` table.

.majority_class_prediction <- function(y_train, n) {
  rep(names(which.max(table(y_train))), n)
}

.classification_metrics_from_preds_df <- function(predictions, y) {
  .classification_metrics(
    truth     = factor(predictions$truth,     levels = levels(y)),
    predicted = factor(predictions$predicted, levels = levels(y))
  )
}

.predict_selected_multinomial <- function(x_train, y_train, x_valid, levels,
                                         random_state = NULL) {
  y_train <- factor(y_train, levels = levels)
  # cv.glmnet needs at least two predictor columns, so a STABL support of size
  # 0 or 1 has no downstream model to fit; fall back to the majority class.
  if (ncol(x_train) < 2L || length(unique(y_train)) < 2L) {
    return(.majority_class_prediction(y_train, nrow(x_valid)))
  }

  pred <- tryCatch(
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
        .seeded_cv_glmnet(
          x = x_train,
          y = y_train,
          family = if (length(levels) == 2L) "binomial" else "multinomial",
          type.measure = "class",
          nfolds = min(5L, min(class_counts)),
          random_state = random_state
        ),
        source = "glmnet selected-feature classification fit"
      )
      as.character(stats::predict(
        fit,
        newx = x_valid,
        s = "lambda.min",
        type = "class"
      )[, 1L])
    },
    stablr_numerical_infeasibility = function(e) {
      .majority_class_prediction(y_train, nrow(x_valid))
    }
  )

  pred
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
                               workers, ...) {
  outer <- outer_folds[[outer_i]]
  train_ids <- outer$train_ids
  valid_ids <- outer$valid_ids

  inner_seed <- .derive_nested_seed(
    random_state, outer_i, .NESTED_INNER_FOLD_SEED_OFFSET
  )
  inner_folds <- .make_cv_folds(
    y = strata_labels[train_ids],
    v = inner_v,
    stratified = stratified,
    random_state = inner_seed
  )

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
    random_state = .derive_nested_seed(
      random_state, outer_i, .NESTED_INNER_FIT_SEED_OFFSET
    ),
    n_lambda = n_lambda,
    l1_ratio = l1_ratio,
    workers = workers,
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
    random_state = .derive_nested_seed(
      random_state, outer_i, .NESTED_OUTER_FIT_SEED_OFFSET
    ),
    n_lambda = n_lambda,
    l1_ratio = l1_ratio,
    workers = workers,
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
  cat("  candidates:  ", paste(names(x$candidates), collapse = ", "), "\n", sep = "")
  cat("  metric:      ", x$metric, "\n", sep = "")
  cat("  accuracy:    ", sprintf("%.3f", x$performance$accuracy), "\n", sep = "")
  cat("  BER:         ", sprintf("%.3f", x$performance$balanced_error_rate), "\n", sep = "")
  invisible(x)
}
