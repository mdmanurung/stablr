# analysis_workflow.R — explicit exploratory/validated analysis modes (PR3).
#
# `stabl_analyze()` runs STABL feature selection under an explicitly declared
# evidence mode.  The mode is a required argument: stablr never implies
# validation.
#
#   * "exploratory": full-data selection plus descriptive interpretation
#     (PR1).  All metrics are apparent / post-selection and are labelled as
#     such; nothing here is held-out evidence.
#   * "validated": binary classification only.  The complete procedure
#     (fold-fitted preprocessing -> STABL selection -> tuning -> downstream
#     predictor) is isolated inside each training partition of the existing
#     nested-CV machinery, yielding genuinely held-out event probabilities,
#     hard labels, selected panels and held-out permutation importance.
#     Gaussian, ordinal, Cox, multiclass-probability and fusion/stacking
#     variants fail clearly.
#
# An optional `final_refit = TRUE` adds an explicitly labelled
# full-development refit for future inference; it is a different model from
# the cross-validation models and plays no role in performance estimation.

#' STABL analysis with an explicit evidence mode
#'
#' Runs STABL feature selection in one of two explicitly declared modes:
#' `"exploratory"` (full-data selection with descriptive, apparent
#' interpretation) or `"validated"` (binary classification with the whole
#' pipeline — preprocessing, selection, tuning and prediction — isolated
#' inside nested cross-validation training partitions).
#'
#' @param x Numeric feature matrix (samples × features) with row names
#'   (sample IDs) and column names (feature names).
#' @param y Outcome aligned to the rows of `x`.  Named vector or factor whose
#'   names cover exactly the row names of `x`.  For `family = "binomial"`,
#'   exactly two classes; for `family = "gaussian"`, numeric.
#' @param family `"binomial"` or `"gaussian"`.  Validated mode currently
#'   supports `"binomial"` only.
#' @param mode **Required.** `"exploratory"` or `"validated"`.  There is no
#'   default: the evidence mode must be chosen explicitly.
#' @param preprocessing Optional fold-safe preprocessing specification:
#'   `NULL` (identity), a character vector of transforms for
#'   [stabl_preprocess_fit()] (`"zv"`, `"impute_median"`, `"center"`,
#'   `"scale"`), or an unprepped [recipes::recipe()].  In validated mode the
#'   transform is fitted anew on every inner- and outer-training split; in
#'   exploratory mode it is fitted once on the full data (an apparent,
#'   fitting-data transform).  Fitted `stabl_preprocessing` objects are
#'   rejected.
#' @param resampling Resampling scheme for validated mode; currently only
#'   `"nested_cv"` (via [stabl_multiomic_nested_cv()]).  Ignored in
#'   exploratory mode.
#' @param event_class For binary outcomes, the level treated as the event for
#'   probabilities, odds ratios and probability metrics; defaults to the
#'   second level of the outcome factor.
#' @param outer_v,outer_repeats,inner_v,stratified,study_labels,lambda_grid,n_bootstraps,artificial_type,hard_threshold,metric
#'   Validated-mode resampling and selection settings, passed to
#'   [stabl_multiomic_nested_cv()].
#' @param cor_method,cor_threshold Correlation method and |correlation|
#'   threshold for redundancy groups.  In validated mode, groups are derived
#'   on each fold's transformed training data only.
#' @param n_repeats Number of permutation repeats for importance (per fold in
#'   validated mode).
#' @param loss Loss for permutation importance: `"squared_error"`
#'   (gaussian), `"brier"` or `"log_loss"` (binomial).  Positive importance
#'   means `permuted_loss - baseline_loss > 0`; negative values are retained.
#' @param random_state Seed.  Fold assignment, selection, and per-fold
#'   importance seeds are all derived from it deterministically.
#' @param final_refit Logical; validated mode only.  When `TRUE`, an
#'   explicitly labelled full-development refit (selection + interpretation
#'   on all data) is added for future inference.  It is a different model
#'   from the cross-validation models and is not part of performance
#'   estimation.
#' @param ... Additional arguments passed to [stabl_fit()] (and, in validated
#'   mode, to [stabl_multiomic_nested_cv()]).
#'
#' @return An object of class `"stabl_analysis"` with a `mode` field.
#'   Exploratory mode returns the full-data `fit` (a `"stabl_fit"`), its
#'   `interpretation` (a `"stabl_interpretation"`, evidence type
#'   `"post_selection_descriptive"`), the fitted `preprocessing` transform,
#'   and `provenance`.  Validated mode returns the underlying
#'   `"stabl_multiomic_nested_cv"` object in `nested`; fold-level held-out
#'   `predictions` (one row per sample per repeat, with `event_probability`
#'   for the named `event_class`); `performance` (pooled and per-repeat
#'   accuracy, BER, macro-F1, MCC, Brier score, clipped log loss and AUROC);
#'   fold-level `selected_features`; `heldout_importance` (per-feature and
#'   per-group permutation importance computed with each fold's own
#'   predictor on that fold's held-out rows, with correlation groups derived
#'   from training data only); and optionally `final_refit`.  Validated
#'   evidence is labelled `"cross_validated"` with evaluation scope
#'   `"held_out_assessment"`.
#'
#' @examples
#' set.seed(2L)
#' ids <- paste0("s", 1:80)
#' x <- matrix(rnorm(80 * 8), 80, 8,
#'             dimnames = list(ids, paste0("f", 1:8)))
#' y <- setNames(factor(rep(c("A", "B"), each = 40L)), ids)
#' x[, 1] <- x[, 1] + 1.5 * as.numeric(y == "B")
#' x[, 2] <- x[, 2] - 1.2 * as.numeric(y == "B")
#'
#' exploratory <- suppressWarnings(stabl_analyze(
#'   x, y, family = "binomial", mode = "exploratory",
#'   lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
#'   n_bootstraps = 10L, random_state = 2L, n_repeats = 3L
#' ))
#' summary(exploratory)
#'
#' validated <- suppressWarnings(stabl_analyze(
#'   x, y, family = "binomial", mode = "validated",
#'   preprocessing = c("center", "scale"),
#'   lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
#'   outer_v = 3L, inner_v = 3L, n_bootstraps = 10L,
#'   random_state = 2L, n_repeats = 3L
#' ))
#' validated$performance$pooled
#' @export
stabl_analyze <- function(x,
                          y,
                          family = c("binomial", "gaussian"),
                          mode,
                          preprocessing = NULL,
                          resampling = c("nested_cv"),
                          event_class = NULL,
                          outer_v = 5L,
                          outer_repeats = 1L,
                          inner_v = 5L,
                          stratified = TRUE,
                          study_labels = NULL,
                          lambda_grid = "auto",
                          n_bootstraps = 100L,
                          artificial_type = "random_permutation",
                          hard_threshold = NULL,
                          metric = c("ber", "accuracy"),
                          cor_method = c("pearson", "spearman"),
                          cor_threshold = 0.7,
                          n_repeats = 20L,
                          loss = NULL,
                          random_state = NULL,
                          final_refit = FALSE,
                          ...) {
  if (missing(mode) || is.null(mode)) {
    stop(
      paste0(
        "`mode` is required and must be given explicitly as ",
        "\"exploratory\" or \"validated\"; stablr never implies validation."
      ),
      call. = FALSE
    )
  }
  mode <- match.arg(mode, c("exploratory", "validated"))
  family <- match.arg(family)
  if (identical(mode, "validated") && !identical(family, "binomial")) {
    stop(
      paste0(
        "Validated mode currently supports binary classification only ",
        "(family = \"binomial\").  Validated gaussian, ordinal, Cox, ",
        "multiclass-probability and fusion/stacking workflows are not ",
        "implemented; use mode = \"exploratory\" for descriptive full-data ",
        "analysis."
      ),
      call. = FALSE
    )
  }
  metric <- match.arg(metric)
  cor_method <- match.arg(cor_method)
  cor_threshold <- .validate_proportion(cor_threshold, "cor_threshold",
                                        allow_one = TRUE)
  n_repeats <- .validate_scalar_integer_like(n_repeats, "n_repeats", min = 1L)
  if (!isTRUE(final_refit) && !identical(final_refit, FALSE)) {
    stop("`final_refit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }
  .validate_preprocessing_spec(preprocessing)

  aligned <- .analyze_align_inputs(x, y, family)
  x <- aligned$x
  y <- aligned$y
  event_class <- .analyze_event_class(y, family, event_class)

  call <- match.call()
  provenance <- list(
    call = call,
    mode = mode,
    family = family,
    event_class = event_class,
    random_state = random_state,
    preprocessing = if (is.null(preprocessing)) {
      "identity"
    } else if (is.character(preprocessing)) {
      paste(preprocessing, collapse = ", ")
    } else {
      "recipes"
    },
    n_observations = nrow(x),
    n_features = ncol(x),
    outcome_levels = if (is.factor(y)) levels(y) else NULL,
    stablr_version = tryCatch(
      as.character(utils::packageVersion("stablr")),
      error = function(e) NA_character_
    )
  )

  if (identical(mode, "exploratory")) {
    if (!missing(resampling)) {
      message("`resampling` is ignored in exploratory mode.")
    }
    if (isTRUE(final_refit)) {
      message(
        "`final_refit` is ignored in exploratory mode: the exploratory fit ",
        "is already a full-data fit."
      )
    }
    return(.analyze_exploratory(
      x = x, y = y, family = family, preprocessing = preprocessing,
      lambda_grid = lambda_grid, n_bootstraps = n_bootstraps,
      artificial_type = artificial_type, hard_threshold = hard_threshold,
      random_state = random_state, cor_method = cor_method,
      cor_threshold = cor_threshold, n_repeats = n_repeats, loss = loss,
      event_class = event_class, provenance = provenance, ...
    ))
  }

  .analyze_validated(
    x = x, y = y, preprocessing = preprocessing, resampling = resampling,
    event_class = event_class, outer_v = outer_v, outer_repeats = outer_repeats,
    inner_v = inner_v, stratified = stratified, study_labels = study_labels,
    lambda_grid = lambda_grid, n_bootstraps = n_bootstraps,
    artificial_type = artificial_type, hard_threshold = hard_threshold,
    metric = metric, cor_method = cor_method, cor_threshold = cor_threshold,
    n_repeats = n_repeats, loss = loss, random_state = random_state,
    final_refit = final_refit, provenance = provenance, ...
  )
}

# ── Input validation ─────────────────────────────────────────────────────────

.analyze_align_inputs <- function(x, y, family) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric after as.matrix().", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("`x` must have row names (sample IDs) and column names (features).",
         call. = FALSE)
  }
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop("`x` row names and column names must be unique.", call. = FALSE)
  }
  if (is.null(names(y))) {
    if (length(y) != nrow(x)) {
      stop("Unnamed `y` must have the same length as nrow(x).", call. = FALSE)
    }
    names(y) <- rownames(x)
  }
  missing_ids <- setdiff(rownames(x), names(y))
  unexpected_ids <- setdiff(names(y), rownames(x))
  if (length(missing_ids) || length(unexpected_ids)) {
    stop(sprintf(
      "`y` names must cover exactly the row names of `x` (missing: %s; unexpected: %s).",
      if (length(missing_ids)) paste(head(missing_ids, 5L), collapse = ", ") else "none",
      if (length(unexpected_ids)) paste(head(unexpected_ids, 5L), collapse = ", ") else "none"
    ), call. = FALSE)
  }
  y <- y[rownames(x)]

  if (identical(family, "gaussian")) {
    if (is.factor(y) || is.character(y)) {
      stop("`y` must be numeric for family = \"gaussian\".", call. = FALSE)
    }
    y <- setNames(as.numeric(y), rownames(x))
    if (anyNA(y)) {
      stop("`y` cannot contain missing values.", call. = FALSE)
    }
  } else {
    y <- factor(y)
    names(y) <- rownames(x)
    if (anyNA(y)) {
      stop("`y` cannot contain missing values.", call. = FALSE)
    }
    if (nlevels(y) != 2L) {
      stop(sprintf(
        paste0(
          "`y` must have exactly two classes for family = \"binomial\" ",
          "(found %d).  Multiclass analysis is available through ",
          "stabl_multiomic_nested_cv() with family = \"multinomial\"."
        ),
        nlevels(y)
      ), call. = FALSE)
    }
  }
  list(x = x, y = y)
}

.analyze_event_class <- function(y, family, event_class) {
  if (!identical(family, "binomial")) {
    if (!is.null(event_class)) {
      stop("`event_class` is only defined for family = \"binomial\".",
           call. = FALSE)
    }
    return(NULL)
  }
  if (is.null(event_class)) return(levels(y)[2L])
  event_class <- as.character(event_class)
  if (!event_class %in% levels(y)) {
    stop(sprintf(
      "`event_class` must be one of the outcome levels (%s).",
      paste(levels(y), collapse = ", ")
    ), call. = FALSE)
  }
  event_class
}

# ── Exploratory mode ─────────────────────────────────────────────────────────

.analyze_exploratory <- function(x, y, family, preprocessing, lambda_grid,
                                 n_bootstraps, artificial_type, hard_threshold,
                                 random_state, cor_method, cor_threshold,
                                 n_repeats, loss, event_class, provenance,
                                 ...) {
  # Full-data transform: apparent by construction and labelled as such.
  pp <- .fit_fold_preprocessing(preprocessing, x)
  x_pp <- if (is.null(pp)) x else stabl_preprocess_apply(pp, x)

  fit <- stabl_fit(
    x = x_pp,
    y = y,
    lambda_grid = lambda_grid,
    family = family,
    n_bootstraps = n_bootstraps,
    artificial_type = artificial_type,
    hard_threshold = hard_threshold,
    random_state = random_state,
    ...
  )
  insight <- stabl_interpret(
    fit, x_pp, y,
    family = family,
    assessment = FALSE,
    cor_method = cor_method,
    cor_threshold = cor_threshold,
    n_repeats = n_repeats,
    loss = loss,
    random_state = random_state,
    event_class = event_class
  )

  structure(
    list(
      mode = "exploratory",
      family = family,
      event_class = event_class,
      fit = fit,
      interpretation = insight,
      preprocessing = pp,
      evidence_type = "post_selection_descriptive",
      evaluation_scope = "fitting_data",
      provenance = c(provenance, list(resampling = NA_character_)),
      final_refit = NULL
    ),
    class = "stabl_analysis"
  )
}

# ── Validated mode ───────────────────────────────────────────────────────────

.analyze_validated <- function(x, y, preprocessing, resampling, event_class,
                               outer_v, outer_repeats, inner_v, stratified,
                               study_labels, lambda_grid, n_bootstraps,
                               artificial_type, hard_threshold, metric,
                               cor_method, cor_threshold, n_repeats, loss,
                               random_state, final_refit, provenance, ...) {
  resampling <- match.arg(resampling, c("nested_cv"))
  # y is already a validated two-class factor (see .analyze_align_inputs).

  nested <- stabl_multiomic_nested_cv(
    x_list = list(x = x),
    y = y,
    candidates = list(list(name = "x", blocks = "x")),
    lambda_grid = lambda_grid,
    outer_v = outer_v,
    outer_repeats = outer_repeats,
    inner_v = inner_v,
    stratified = stratified,
    metric = metric,
    family = "binomial",
    n_bootstraps = n_bootstraps,
    artificial_type = artificial_type,
    hard_threshold = hard_threshold,
    random_state = random_state,
    study_labels = study_labels,
    preprocessing = preprocessing,
    event_class = event_class,
    return_probabilities = TRUE,
    keep_fold_predictors = TRUE,
    ...
  )

  preds <- nested$outer_predictions
  performance <- list(
    pooled = .analyze_performance_row(preds, y, event_class),
    by_repeat = do.call(rbind, lapply(
      split(preds, preds$repeat_id),
      .analyze_performance_row,
      y = y,
      event_class = event_class
    ))
  )
  rownames(performance$by_repeat) <- NULL

  loss <- loss %||% "brier"
  heldout <- .analyze_fold_importance(
    nested = nested,
    x = x,
    y = y,
    event_class = event_class,
    loss = loss,
    n_repeats = n_repeats,
    cor_method = cor_method,
    cor_threshold = cor_threshold,
    random_state = random_state
  )

  fr <- NULL
  if (isTRUE(final_refit)) {
    # Explicitly labelled full-development refit for future inference.  It is
    # a different model from the CV models and plays no role in the
    # performance estimates above.
    pp_full <- .fit_fold_preprocessing(preprocessing, x)
    x_full <- if (is.null(pp_full)) x else stabl_preprocess_apply(pp_full, x)
    fit_full <- stabl_fit(
      x = x_full,
      y = y,
      lambda_grid = lambda_grid,
      family = "binomial",
      n_bootstraps = n_bootstraps,
      artificial_type = artificial_type,
      hard_threshold = hard_threshold,
      random_state = random_state,
      ...
    )
    interp_full <- stabl_interpret(
      fit_full, x_full, y,
      family = "binomial",
      assessment = FALSE,
      cor_method = cor_method,
      cor_threshold = cor_threshold,
      n_repeats = n_repeats,
      loss = loss,
      random_state = random_state,
      event_class = event_class
    )
    interp_full$evidence$evidence_type <- "final_refit"
    interp_full$provenance$evidence_type <- "final_refit"
    fr <- list(
      fit = fit_full,
      interpretation = interp_full,
      preprocessing = pp_full,
      evidence_type = "final_refit",
      note = paste0(
        "Full-development refit for future inference; a different model ",
        "from the cross-validation models and not part of performance ",
        "estimation."
      )
    )
  }

  structure(
    list(
      mode = "validated",
      family = "binomial",
      event_class = event_class,
      resampling = resampling,
      nested = nested,
      predictions = preds,
      performance = performance,
      performance_by_study = nested$performance_by_study,
      selected_features = nested$selected_features,
      heldout_importance = heldout,
      evidence_type = "cross_validated",
      evaluation_scope = "held_out_assessment",
      provenance = c(provenance, list(
        resampling = resampling,
        fold_scheme = nested$fold_scheme,
        loss = loss
      )),
      final_refit = fr
    ),
    class = "stabl_analysis"
  )
}

# One row of hard-class + probability performance for a predictions frame.
.analyze_performance_row <- function(preds, y, event_class) {
  truth <- factor(preds$truth, levels = levels(y))
  predicted <- factor(preds$predicted, levels = levels(y))
  hard <- .classification_metrics(truth, predicted)
  prob <- .binary_probability_metrics(truth, preds$event_probability,
                                      event_class)
  data.frame(
    n = nrow(preds),
    n_events = prob$n_events,
    accuracy = hard$accuracy,
    balanced_error_rate = hard$balanced_error_rate,
    macro_f1 = hard$macro_f1,
    mcc = hard$mcc,
    brier = prob$brier,
    log_loss = prob$log_loss,
    auroc = prob$auroc,
    stringsAsFactors = FALSE
  )
}

# Rebuild the single-block candidate matrix used inside the nested-CV
# machinery (block name "x", columns prefixed "x::").
.analyze_candidate_matrix <- function(x, ids) {
  xm <- as.matrix(x[ids, , drop = FALSE])
  colnames(xm) <- paste("x", colnames(xm), sep = "::")
  rownames(xm) <- ids
  xm
}

# Held-out permutation importance: for every outer fold, use THAT fold's
# trained predictor and fold-fitted preprocessing on the fold's held-out
# rows.  Correlation groups are derived on the fold's transformed training
# data only.  Positive importance means permuted loss exceeds baseline loss;
# negative values are retained.
.analyze_fold_importance <- function(nested, x, y, event_class, loss,
                                     n_repeats, cor_method, cor_threshold,
                                     random_state) {
  loss_fun <- .interpret_loss(loss, "binomial")
  fold_results <- nested$fold_results
  fold_rows <- list()
  groups_by_fold <- list()

  for (i in seq_along(fold_results)) {
    fr <- fold_results[[i]]
    fit_rec <- fr$fit
    predictor <- fit_rec$predictor
    fold_id <- fr$outer_fold$fold_id
    if (is.null(predictor) || is.null(predictor$predict_event_prob)) next

    support <- get_support(fit_rec$stabl_fit)
    selected <- names(support)[support]
    train_ids <- fr$outer_fold$train_ids
    valid_ids <- fr$outer_fold$valid_ids

    if (length(selected) > 0L) {
      x_tr <- .analyze_candidate_matrix(x, train_ids)
      x_va <- .analyze_candidate_matrix(x, valid_ids)
      pp <- fit_rec$preprocessing
      if (!is.null(pp)) {
        x_tr <- stabl_preprocess_apply(pp, x_tr)
        x_va <- stabl_preprocess_apply(pp, x_va)
      }
      xs_tr <- x_tr[, selected, drop = FALSE]
      xs_va <- x_va[, selected, drop = FALSE]
      y_vec <- as.integer(y[valid_ids] == event_class)

      groups <- .interpret_correlation_groups(xs_tr, selected, cor_method,
                                              cor_threshold)
      groups_by_fold[[fold_id]] <- lapply(groups$groups, function(g) {
        sub("^x::", "", g)
      })

      # Individual features plus multi-member groups (one shared row
      # permutation across all columns of a group).
      sets <- as.list(selected)
      names(sets) <- selected
      multi <- groups$groups[lengths(groups$groups) > 1L]
      if (length(multi)) {
        names(multi) <- paste0("group:", names(multi))
        sets <- c(sets, multi)
      }

      imp <- .permutation_importance_sets(
        predict_fn = predictor$predict_event_prob,
        x_eval = xs_va,
        y_vec = y_vec,
        sets = sets,
        loss_fun = loss_fun,
        n_repeats = n_repeats,
        random_state = .derive_nested_seed(random_state, i, 7000L)
      )
      fold_rows[[fold_id]] <- data.frame(
        repeat_id = fr$outer_fold[["repeat"]],
        fold_id = fold_id,
        feature_set = names(sets),
        set_type = ifelse(startsWith(names(sets), "group:"), "group", "feature"),
        importance_mean = imp$mean,
        importance_sd = imp$sd,
        n_repeats = n_repeats,
        n_selected = length(selected),
        stringsAsFactors = FALSE
      )
    } else {
      groups_by_fold[[fold_id]] <- list()
      # Zero-row placeholder: every column must be length 0 so the frame
      # binds cleanly with non-empty fold rows.
      fold_rows[[fold_id]] <- data.frame(
        repeat_id = integer(),
        fold_id = character(),
        feature_set = character(),
        set_type = character(),
        importance_mean = numeric(),
        importance_sd = numeric(),
        n_repeats = integer(),
        n_selected = integer(),
        stringsAsFactors = FALSE
      )
    }
  }

  folds_df <- do.call(rbind, fold_rows)
  rownames(folds_df) <- NULL

  # Per-feature summary across the folds in which the feature was selected.
  feat_rows <- folds_df[folds_df$set_type == "feature", , drop = FALSE]
  feats <- colnames(x)
  summary_df <- data.frame(
    feature = feats,
    mean_importance = NA_real_,
    sd_importance = NA_real_,
    n_folds_selected = 0L,
    n_repeats = n_repeats,
    loss = loss,
    event_class = event_class,
    evidence_type = "cross_validated",
    evaluation_scope = "held_out_assessment",
    stringsAsFactors = FALSE
  )
  if (nrow(feat_rows)) {
    feat_rows$feature <- sub("^x::", "", feat_rows$feature_set)
    for (f in unique(feat_rows$feature)) {
      vals <- feat_rows$importance_mean[feat_rows$feature == f]
      idx <- which(summary_df$feature == f)
      summary_df$mean_importance[idx] <- mean(vals)
      summary_df$sd_importance[idx] <- if (length(vals) > 1L) {
        stats::sd(vals)
      } else {
        NA_real_
      }
      summary_df$n_folds_selected[idx] <- length(vals)
    }
  }
  summary_df$na_reason <- ifelse(
    summary_df$n_folds_selected == 0L, "not_selected", NA_character_
  )

  list(
    summary = summary_df,
    folds = folds_df,
    groups_by_fold = groups_by_fold,
    loss = loss,
    event_class = event_class
  )
}

# Permutation importance for named column sets with a prediction callback.
# Each set is permuted jointly with one shared row permutation per repeat;
# importance is permuted_loss - baseline_loss (negative values retained).
.permutation_importance_sets <- function(predict_fn, x_eval, y_vec, sets,
                                         loss_fun, n_repeats, random_state) {
  out <- list(
    mean = setNames(rep(NA_real_, length(sets)), names(sets)),
    sd = setNames(rep(NA_real_, length(sets)), names(sets)),
    repeats = NULL
  )
  if (length(sets) == 0L || nrow(x_eval) < 2L) return(out)

  base_loss <- loss_fun(y_vec, predict_fn(x_eval))
  if (!is.null(random_state)) set.seed(random_state)
  rep_mat <- matrix(NA_real_, nrow = n_repeats, ncol = length(sets),
                    dimnames = list(NULL, names(sets)))
  for (r in seq_len(n_repeats)) {
    for (s in names(sets)) {
      cols <- sets[[s]]
      xp <- x_eval
      idx <- sample(nrow(xp))
      xp[, cols] <- xp[idx, cols, drop = FALSE]
      rep_mat[r, s] <- loss_fun(y_vec, predict_fn(xp)) - base_loss
    }
  }
  out$mean <- colMeans(rep_mat)
  out$sd <- apply(rep_mat, 2L, stats::sd)
  out$repeats <- rep_mat
  out
}

# ── S3 methods ───────────────────────────────────────────────────────────────

#' @describeIn stabl_analyze Print a concise summary of a `"stabl_analysis"`
#'   object; invisibly returns `x`.
#' @export
print.stabl_analysis <- function(x, ...) {
  cat("<stabl_analysis:", x$mode, ">\n")
  cat("  family:      ", x$family, "\n", sep = "")
  if (!is.null(x$event_class)) {
    cat("  event class: ", x$event_class, "\n", sep = "")
  }
  cat("  preprocessing: ", x$provenance$preprocessing, "\n", sep = "")
  if (identical(x$mode, "exploratory")) {
    cat("  selected:    ", x$provenance$n_features -
          sum(!x$interpretation$evidence$selected), "of",
        x$provenance$n_features, "features\n")
    cat("  evidence:    ", x$evidence_type, " (", x$evaluation_scope,
        "; apparent)\n", sep = "")
  } else {
    p <- x$performance$pooled
    cat("  resampling:  ", x$resampling, " (", x$provenance$fold_scheme, ")\n",
        sep = "")
    cat("  held-out:    ", p$n, " predictions; acc ",
        sprintf("%.3f", p$accuracy), ", BER ",
        sprintf("%.3f", p$balanced_error_rate), ", brier ",
        sprintf("%.3f", p$brier), ", log-loss ",
        sprintf("%.3f", p$log_loss), ", AUROC ",
        if (is.na(p$auroc)) "NA (single-class)" else sprintf("%.3f", p$auroc),
        "\n", sep = "")
    cat("  evidence:    ", x$evidence_type, " (", x$evaluation_scope, ")\n",
        sep = "")
    if (!is.null(x$final_refit)) {
      cat("  final refit: yes (separate full-development model)\n")
    }
  }
  invisible(x)
}

#' @describeIn stabl_analyze Return the key evidence table invisibly: the
#'   interpretation evidence for exploratory mode, or the held-out
#'   performance and importance summary for validated mode.
#' @param object A `"stabl_analysis"` object (for `summary`).
#' @export
summary.stabl_analysis <- function(object, ...) {
  print(object)
  if (identical(object$mode, "exploratory")) {
    ev <- object$interpretation$evidence
    sel <- ev[ev$selected, , drop = FALSE]
    cat("\nSelected features (post-selection descriptive):\n")
    if (nrow(sel)) {
      print(sel[, c("feature", "stability", "coefficient", "effect_direction",
                    "correlation_group", "permutation_importance"),
                drop = FALSE])
    } else {
      cat("  (none)\n")
    }
    return(invisible(ev))
  }
  cat("\nHeld-out performance (pooled):\n")
  print(object$performance$pooled)
  imp <- object$heldout_importance$summary
  imp <- imp[!is.na(imp$mean_importance), , drop = FALSE]
  cat("\nHeld-out permutation importance (selected features):\n")
  if (nrow(imp)) {
    print(imp[, c("feature", "mean_importance", "sd_importance",
                  "n_folds_selected"), drop = FALSE])
  } else {
    cat("  (none)\n")
  }
  invisible(list(
    performance = object$performance,
    heldout_importance = object$heldout_importance$summary
  ))
}
