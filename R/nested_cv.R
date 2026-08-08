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
#'   selected features, and aggregate performance.
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

.evaluate_stabl_candidates_inner <- function(x_list, y, train_ids, candidates,
                                             inner_folds, lambda_grid, metric,
                                             family, n_bootstraps,
                                             artificial_type, hard_threshold,
                                             random_state, n_lambda, l1_ratio,
                                             workers, ...) {
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
                                        n_lambda, l1_ratio, workers, ...) {
  x_train <- .candidate_matrix(x_list, candidate, train_ids)
  x_valid <- .candidate_matrix(x_list, candidate, valid_ids)
  y_train <- y[train_ids]

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
  predicted <- .predict_selected_multinomial(
    x_train = x_train[, selected, drop = FALSE],
    y_train = y_train,
    x_valid = x_valid[, selected, drop = FALSE],
    levels = levels(y)
  )

  list(
    stabl_fit = fit,
    selected_features = .split_prefixed_features(selected),
    importances = importances,
    predicted = predicted,
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

.predict_selected_multinomial <- function(x_train, y_train, x_valid, levels) {
  y_train <- factor(y_train, levels = levels)
  if (ncol(x_train) == 0L || length(unique(y_train)) < 2L) {
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
        glmnet::cv.glmnet(
          x = x_train,
          y = y_train,
          family = if (length(levels) == 2L) "binomial" else "multinomial",
          type.measure = "class",
          nfolds = min(5L, min(class_counts))
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

  inner_seed <- .derive_nested_seed(random_state, outer_i, 1000L)
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
    random_state = .derive_nested_seed(random_state, outer_i, 2000L),
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
    random_state = .derive_nested_seed(random_state, outer_i, 3000L),
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
