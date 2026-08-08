# Leakage-safe late fusion.  The historical implementation remains in
# multiomic_workflows.R because it is also the public Python-parity path.

.late_fusion_oof_fit <- function(x_train_list, y_train, x_valid_list, y_valid,
                                 groups_train, bootstrap_strata_train,
                                 lambda_by_omic, fit_params,
                                 omic_names, family, n_iter_lf, nfolds,
                                 random_state, ...) {
  task_type <- .family_to_task_type(family)
  binary_mapping <- NULL
  y_downstream <- y_train
  if (identical(task_type, "binary")) {
    mapped <- .late_fusion_binary_outcome(y_train)
    y_downstream <- mapped$y
    binary_mapping <- mapped$mapping
  }
  sample_ids <- rownames(x_train_list[[1L]])
  if (length(sample_ids) < nfolds) {
    stop("`late_fusion_nfolds` cannot exceed the number of training samples.",
         call. = FALSE)
  }
  folds <- .make_multiomic_cv_folds(
    sample_ids = sample_ids, groups = groups_train, v = nfolds,
    random_state = .derive_nested_seed(random_state, 1L, 29001L),
    strata = if (task_type %in% c("binary", "multiclass")) y_downstream else NULL
  )
  fold_id <- setNames(integer(length(sample_ids)), sample_ids)
  fold_details <- vector("list", length(folds))
  warnings_seen <- character()

  if (identical(task_type, "multiclass")) {
    classes <- levels(factor(y_train))
    oof <- setNames(lapply(omic_names, function(.) {
      matrix(NA_real_, length(sample_ids), length(classes),
             dimnames = list(sample_ids, classes))
    }), omic_names)
  } else {
    oof <- matrix(NA_real_, length(sample_ids), length(omic_names),
                  dimnames = list(sample_ids, omic_names))
  }

  for (i in seq_along(folds)) {
    fold <- folds[[i]]
    train_ids <- fold$train_ids
    assess_ids <- fold$valid_ids
    fold_id[assess_ids] <- i
    fold_seed <- .derive_nested_seed(random_state, i, 31001L)
    fold_params <- fit_params
    fold_params$groups <- if (is.null(groups_train)) NULL else groups_train[train_ids]
    fold_params$bootstrap_strata <- .subset_bootstrap_strata_by_ids(
      bootstrap_strata_train, train_ids, "bootstrap_strata_train"
    )
    fold_params$random_state <- fold_seed

    captured <- withCallingHandlers(
      .late_fusion_select_per_omic_safe(
        x_train_list = .subset_multiomic_rows(x_train_list, train_ids),
        y_train = .subset_outcome_by_ids(y_train, train_ids),
        lambda_by_omic = lambda_by_omic,
        x_valid_list = .subset_multiomic_rows(x_train_list, assess_ids),
        omic_names = omic_names, fit_params = fold_params, ...
      ),
      warning = function(w) {
        warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )

    fallback <- setNames(character(length(omic_names)), omic_names)
    for (omic in omic_names) {
      pred <- .late_fusion_fit_omic_safe(
        x_train_sel = captured$selected_train[[omic]],
        y_train = .subset_outcome_by_ids(y_downstream, train_ids),
        x_valid_sel = captured$selected_valid[[omic]],
        task_type = task_type,
        levels = if (identical(task_type, "multiclass")) classes else NULL,
        selection_error = captured$selection_errors[[omic]]
      )
      warnings_seen <- unique(c(warnings_seen, pred$warnings))
      fallback[[omic]] <- pred$fallback_reason %||% NA_character_
      if (identical(task_type, "multiclass")) {
        oof[[omic]][assess_ids, ] <- pred$valid_preds
      } else {
        oof[assess_ids, omic] <- pred$valid_preds
      }
    }
    fold_artificial_provenance <- lapply(
      captured$fits,
      function(x) x$artificial_provenance %||% NULL
    )
    fold_details[[i]] <- list(
      fold = fold$fold, train_ids = train_ids, assessment_ids = assess_ids,
      seed = fold_seed,
      selected_features = captured$selected_features,
      artificial_provenance = fold_artificial_provenance,
      artificial_feature_provenance = fold_artificial_provenance,
      fallback_reasons = fallback,
      selector_errors = captured$selection_errors
    )
  }

  if (any(fold_id == 0L)) stop("Internal error: incomplete OOF coverage.", call. = FALSE)
  stacked <- stacked_multi_omic(
    predictions = oof, y = y_downstream, task_type = task_type,
    n_iter = n_iter_lf,
    random_state = .derive_nested_seed(random_state, 1L, 47001L)
  )

  full_captured <- withCallingHandlers(
    .fit_multiomic_per_omic(
      x_train_list = x_train_list, y_train = y_train,
      lambda_by_omic = lambda_by_omic, x_valid_list = x_valid_list,
      omic_names = omic_names, fit_params = fit_params, ...
    ),
    warning = function(w) {
      warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )

  full <- .late_fusion_full_refit_predictions(
    full_per_omic = full_captured, y_train = y_downstream,
    x_valid_list = x_valid_list, task_type = task_type,
    omic_names = omic_names,
    levels = if (identical(task_type, "multiclass")) levels(factor(y_train)) else NULL,
    weights = stacked$weights$Associated_weight
  )
  result <- list(
    weights = stacked$weights,
    train_predictions = stacked$predictions,
    valid_predictions = full$valid_predictions,
    score = stacked$score,
    full_per_omic = full_captured,
    provenance = list(
      training_mode = "oof", fold_id = fold_id,
      seeds = list(folds = vapply(fold_details, `[[`, integer(1L), "seed"),
                   stacking = .derive_nested_seed(random_state, 1L, 47001L),
                   full_refit = random_state),
      folds = fold_details, warnings = warnings_seen,
      binary_event_mapping = binary_mapping,
      full_refit = full$diagnostics
    )
  )
  if (identical(task_type, "multiclass")) {
    result$task_type <- "multiclass"
    result$levels <- stacked$levels
    result$log_loss <- stacked$log_loss
    result$train_metrics <- .classification_metrics(
      factor(y_train, levels = stacked$levels),
      factor(stacked$predictions$predicted_class, levels = stacked$levels)
    )
    if (!is.null(y_valid) && !is.null(full$valid_predictions)) {
      result$valid_metrics <- .classification_metrics(
        factor(y_valid, levels = stacked$levels),
        factor(full$valid_predictions$predicted_class, levels = stacked$levels)
      )
    }
  }
  result
}

.late_fusion_fit_omic_safe <- function(x_train_sel, y_train, x_valid_sel,
                                       task_type, levels = NULL,
                                       selection_error = NULL) {
  fallback_reason <- if (!is.null(selection_error)) {
    paste0("selector_fit_error: ", selection_error)
  } else if (ncol(x_train_sel) == 0L) {
    "no_selected_features"
  } else {
    NULL
  }
  y_mean <- if (identical(task_type, "regression")) mean(unname(y_train)) else NA_real_
  captured_warnings <- character()
  out <- tryCatch(
    withCallingHandlers(
      .late_fusion_fit_omic(
        x_train_sel, y_train, x_valid_sel, y_mean, task_type, levels
      ),
      warning = function(w) {
        captured_warnings <<- c(captured_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      fallback_reason <<- paste0("downstream_fit_error: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(fallback_reason) && identical(task_type, "binary")) {
    prior <- mean(as.numeric(y_train))
    out$train_preds <- rep(prior, nrow(x_train_sel))
    if (!is.null(x_valid_sel)) out$valid_preds <- rep(prior, nrow(x_valid_sel))
  }
  if (is.null(out)) {
    empty <- x_train_sel[, FALSE, drop = FALSE]
    empty_valid <- if (is.null(x_valid_sel)) NULL else
      x_valid_sel[, FALSE, drop = FALSE]
    out <- .late_fusion_fit_omic(
      empty, y_train, empty_valid, y_mean, task_type, levels
    )
  } else if (is.null(out$model) && is.null(fallback_reason)) {
    fallback_reason <- "downstream_fit_fallback"
  }
  valid_values <- c(
    unlist(out$train_preds, use.names = FALSE),
    unlist(out$valid_preds, use.names = FALSE)
  )
  if (length(valid_values) && any(!is.finite(valid_values))) {
    stop("Late-fusion downstream predictions must be finite.", call. = FALSE)
  }
  out$fallback_reason <- fallback_reason
  out$warnings <- captured_warnings
  out
}

.late_fusion_select_per_omic_safe <- function(x_train_list, y_train,
                                               lambda_by_omic, x_valid_list,
                                               omic_names, fit_params, ...) {
  out <- list(
    fits = setNames(vector("list", length(omic_names)), omic_names),
    selected_features = setNames(vector("list", length(omic_names)), omic_names),
    selected_train = setNames(vector("list", length(omic_names)), omic_names),
    selected_valid = if (is.null(x_valid_list)) NULL else
      setNames(vector("list", length(omic_names)), omic_names),
    selection_errors = setNames(vector("list", length(omic_names)), omic_names)
  )
  for (omic in omic_names) {
    one <- tryCatch(
      .fit_multiomic_per_omic(
        x_train_list = x_train_list[omic], y_train = y_train,
        lambda_by_omic = lambda_by_omic[omic],
        x_valid_list = if (is.null(x_valid_list)) NULL else x_valid_list[omic],
        omic_names = omic, fit_params = fit_params, ...
      ),
      error = function(e) e
    )
    if (inherits(one, "error")) {
      out$selection_errors[[omic]] <- conditionMessage(one)
      out$selected_features[[omic]] <- character()
      out$selected_train[[omic]] <- as.matrix(x_train_list[[omic]])[, FALSE, drop = FALSE]
      if (!is.null(x_valid_list)) {
        out$selected_valid[[omic]] <- as.matrix(x_valid_list[[omic]])[, FALSE, drop = FALSE]
      }
    } else {
      out$fits[[omic]] <- one$fits[[omic]]
      out$selected_features[[omic]] <- one$selected_features[[omic]]
      out$selected_train[[omic]] <- one$selected_train[[omic]]
      if (!is.null(x_valid_list)) out$selected_valid[[omic]] <- one$selected_valid[[omic]]
    }
  }
  out
}

.late_fusion_binary_outcome <- function(y) {
  ids <- names(y)
  if (is.logical(y)) {
    values <- as.integer(y)
    labels <- c("FALSE", "TRUE")
  } else if (is.numeric(y) && setequal(unique(unname(y)), c(0, 1))) {
    values <- as.integer(y)
    labels <- c("0", "1")
  } else {
    f <- factor(y)
    if (nlevels(f) != 2L || anyNA(f)) {
      stop("Binary outcomes must contain exactly two complete event levels.",
           call. = FALSE)
    }
    values <- as.integer(f) - 1L
    labels <- levels(f)
  }
  names(values) <- ids
  list(
    y = values,
    mapping = setNames(c(0L, 1L), labels)
  )
}

.late_fusion_full_refit_predictions <- function(full_per_omic, y_train,
                                                x_valid_list, task_type,
                                                omic_names, levels, weights) {
  fallbacks <- setNames(character(length(omic_names)), omic_names)
  downstream_warnings <- setNames(vector("list", length(omic_names)), omic_names)
  if (identical(task_type, "multiclass")) {
    valid <- if (is.null(x_valid_list)) NULL else .named_omic_list(omic_names)
  } else {
    valid <- if (is.null(x_valid_list)) NULL else matrix(
      NA_real_, nrow(x_valid_list[[1L]]), length(omic_names),
      dimnames = list(rownames(x_valid_list[[1L]]), omic_names)
    )
  }
  for (omic in omic_names) {
    pred <- .late_fusion_fit_omic_safe(
      full_per_omic$selected_train[[omic]], y_train,
      if (is.null(x_valid_list)) NULL else full_per_omic$selected_valid[[omic]],
      task_type, levels
    )
    fallbacks[[omic]] <- pred$fallback_reason %||% NA_character_
    downstream_warnings[[omic]] <- pred$warnings
    if (!is.null(valid)) {
      if (identical(task_type, "multiclass")) valid[[omic]] <- pred$valid_preds
      else valid[, omic] <- pred$valid_preds
    }
  }
  combined <- if (is.null(valid)) NULL else if (identical(task_type, "multiclass")) {
    .apply_multiclass_stack_weights(valid, weights)
  } else {
    .weighted_masked_mean(valid, weights)
  }
  artificial_provenance <- lapply(
    full_per_omic$fits,
    function(x) x$artificial_provenance %||% NULL
  )
  list(
    valid_predictions = combined,
    diagnostics = list(
      selected_features = full_per_omic$selected_features,
      fallback_reasons = fallbacks,
      warnings = downstream_warnings,
      artificial_provenance = artificial_provenance,
      artificial_feature_provenance = artificial_provenance
    )
  )
}
