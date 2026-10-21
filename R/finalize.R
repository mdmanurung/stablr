# finalize.R — portable final predictive models (1.0.0 spec A6).
#
# `stabl_finalize()` turns a completed development result (a `"stabl_fit"`,
# or a `"stabl_analysis"` carrying one) into a `"stabl_predictor"`: the
# selected panel, the training-fitted preprocessing transform, the event /
# level mapping and a downstream predictive engine refitted on the full
# development data, packed into a single serializable object.
#
# Design rules (spec A6):
#   * No training observations are stored.  One-feature glm engines are
#     reduced to their coefficient vector (predict.glm would need the qr /
#     residuals slots, which encode the training data); cv.glmnet objects
#     store no raw data; recipes preprocessing has its retained training
#     `template` stripped.  Fallback engines store only training-derived
#     scalars (event prevalence / majority class / outcome mean).
#   * Predict-time preprocessing applies the STORED transform; imputation
#     and normalization parameters are never recomputed from new data.
#   * `validate_newdata_schema()` enforces unique, named, numeric columns:
#     missing or duplicated required columns are always an error; extra
#     columns are an error unless explicitly dropped.

#' Finalize a STABL selection into a portable predictive model
#'
#' Takes a completed development result — a [stabl_fit()] object, or a
#' [stabl_analyze()] result — and refits a downstream predictive engine on
#' the selected panel using the full development data.  The result is a
#' `"stabl_predictor"` that can be saved with [saveRDS()] and applied to
#' new data with [predict()][predict.stabl_predictor()] without carrying
#' any training observations.
#'
#' @param development_result A `"stabl_fit"` from [stabl_fit()], or a
#'   `"stabl_analysis"` from [stabl_analyze()].  For a validated analysis,
#'   the object must have been produced with `final_refit = TRUE` (the
#'   full-development refit supplies the selected panel); the
#'   cross-validated performance summary is carried into the predictor's
#'   `development` field, explicitly labelled as distinct from the
#'   full-data refit stored here.  For an exploratory analysis the
#'   full-data fit is used and the development summary is labelled
#'   apparent / post-selection.
#' @param x Development feature matrix (samples x features) with row and
#'   column names; the data on which the downstream engine and the
#'   preprocessing transform are fitted.  Must contain every selected
#'   feature.
#' @param y Development outcome aligned to the rows of `x` (named, or
#'   length `nrow(x)`).  Factor with exactly two levels for
#'   `family = "binomial"`, numeric for `family = "gaussian"`.  The family
#'   is taken from the fitted object, not re-declared.
#' @param preprocessing Optional preprocessing specification: `NULL`
#'   (identity), a character vector of transforms for
#'   [stabl_preprocess_fit()] (`"zv"`, `"impute_median"`, `"center"`,
#'   `"scale"`), or an unprepped [recipes::recipe()].  Fitted here on the
#'   full development data — the correct choice for a final model, since
#'   no held-out claim is derived from this fit.  The stored transform is
#'   applied unchanged at predict time.
#' @param event_class For binary outcomes, the level treated as the event
#'   for probability predictions; defaults to the second level of the
#'   outcome factor.
#' @param random_state Seed for the downstream engine's internal
#'   cross-validation (glmnet path).  `NULL` leaves the ambient stream;
#'   set it for a bitwise-reproducible final model.
#'
#' @details
#' The downstream engine mirrors the nested-CV machinery: an empty panel
#' or degenerate outcome falls back to a deterministic intercept-only
#' model (training event prevalence / majority class for classification,
#' training mean for regression); a one-feature panel uses a glm
#' (stored as bare coefficients); two or more features use a seeded
#' [glmnet::cv.glmnet()] fit.  Only binomial and gaussian fits are
#' supported; other families fail clearly.
#'
#' When `preprocessing` is `NULL`, new data must supply exactly the
#' selected panel.  When a transform is fitted, new data must supply all
#' development columns, because the stored transform may depend on any of
#' them; the panel is selected after transformation.
#'
#' @return A `"stabl_predictor"` object with the packed `model`, `engine`
#'   tag, `family`, outcome `levels` / `event_class`, `selected_features`,
#'   `required_columns`, the fitted `preprocessing` transform, the
#'   `learner` specification, `feature_provenance` (selection scores),
#'   a `development` summary, the stablr version and the call.  No
#'   training observations are stored.
#'
#' @examples
#' set.seed(2L)
#' ids <- paste0("s", 1:60)
#' x <- matrix(rnorm(60 * 6), 60, 6,
#'             dimnames = list(ids, paste0("f", 1:6)))
#' y <- setNames(2 * x[, 1] - x[, 2] + rnorm(60), ids)
#'
#' fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
#'                  family = "gaussian", n_bootstraps = 20L, random_state = 2L)
#' final <- stabl_finalize(fit, x, y, preprocessing = c("center", "scale"),
#'                         random_state = 2L)
#' final
#' predict(final, x[1:3, ], type = "response")
#' @export
stabl_finalize <- function(development_result, x, y, preprocessing = NULL,
                           event_class = NULL, random_state = NULL) {
  call <- match.call()
  extracted <- .finalize_extract_fit(development_result)
  fit <- extracted$fit
  development <- extracted$development
  .check_fitted_stabl(fit)

  family <- fit$family
  if (!family %in% c("binomial", "gaussian")) {
    stop(sprintf(
      paste0(
        "stabl_finalize() supports binomial and gaussian fits; got ",
        "family = \"%s\".  Other families are not yet portable."
      ),
      family
    ), call. = FALSE)
  }
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }
  .validate_preprocessing_spec(preprocessing)

  aligned <- .analyze_align_inputs(x, y, family)
  x <- aligned$x
  y <- aligned$y

  levels <- .modeled_class_levels(y, family)
  event_class <- .analyze_event_class(y, family, event_class)

  support <- get_support(fit)
  selected <- names(support)[support]

  # Fit the transform on the full development data.  This is the final
  # model: no held-out claim is derived from this fit, so a full-data
  # transform is correct (and required for deployment).
  pp <- .fit_fold_preprocessing(preprocessing, x)
  x_tr <- if (is.null(pp)) x else stabl_preprocess_apply(pp, x)

  missing_panel <- setdiff(selected, colnames(x_tr))
  if (length(missing_panel)) {
    stop(sprintf(
      paste0(
        "Selected feature(s) not present after preprocessing: %s.  ",
        "Pass the development matrix and preprocessing specification ",
        "used at selection time."
      ),
      paste(head(missing_panel, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  # A prepped recipe retains the baked training data in `template`; the
  # portable object must not carry training observations.  bake() on new
  # data does not need it.
  if (!is.null(pp) && inherits(pp, "stabl_preprocessing_recipes") &&
      !is.null(pp$prepped$template)) {
    pp$prepped$template <- NULL
  }

  # Downstream engine on the selected panel, reusing the nested-CV
  # engines.  A one-row dummy stands in for `x_valid`: the engines also
  # return immediate predictions, which are discarded here.
  x_sel <- x_tr[, selected, drop = FALSE]
  x_dummy <- x_sel[seq_len(min(1L, nrow(x_sel))), , drop = FALSE]
  engine_result <- if (identical(family, "gaussian")) {
    .fit_selected_regressor(x_sel, y, x_dummy, random_state = random_state)
  } else {
    .fit_selected_predictor(x_sel, y, x_dummy, levels = levels,
                            event_class = event_class,
                            random_state = random_state)
  }
  packed <- .finalize_pack_model(engine_result, family)

  required_columns <- if (is.null(pp)) selected else colnames(x)

  structure(
    list(
      model = packed$model,
      engine = packed$engine,
      family = family,
      levels = levels,
      event_class = event_class,
      selected_features = selected,
      required_columns = required_columns,
      preprocessing = pp,
      learner = packed$learner,
      feature_provenance = data.frame(
        feature = fit$feature_names,
        stabl_score = unname(get_importances(fit)[fit$feature_names]),
        selected = as.logical(support[fit$feature_names]),
        row.names = NULL,
        stringsAsFactors = FALSE
      ),
      development = development,
      stablr_version = tryCatch(
        as.character(utils::packageVersion("stablr")),
        error = function(e) NA_character_
      ),
      call = call
    ),
    class = "stabl_predictor"
  )
}

# Extract the full-data stabl_fit and build the development summary from a
# supported development result.
.finalize_extract_fit <- function(development_result) {
  if (inherits(development_result, "stabl_fit")) {
    return(list(
      fit = development_result,
      development = list(
        source = "stabl_fit",
        evidence_type = "full_data_refit",
        performance = NULL,
        note = paste0(
          "Full-data refit for deployment.  No validation summary was ",
          "supplied; apparent performance on the development data would ",
          "be optimistic and is not reported."
        )
      )
    ))
  }
  if (inherits(development_result, "stabl_analysis")) {
    mode <- development_result$mode
    if (identical(mode, "validated")) {
      fr <- development_result$final_refit
      if (is.null(fr) || is.null(fr$fit)) {
        stop(
          paste0(
            "A validated \"stabl_analysis\" must be produced with ",
            "final_refit = TRUE to be finalized (the full-development ",
            "refit supplies the selected panel).  Rerun stabl_analyze() ",
            "with final_refit = TRUE, or pass a \"stabl_fit\" directly."
          ),
          call. = FALSE
        )
      }
      return(list(
        fit = fr$fit,
        development = list(
          source = "stabl_analysis_validated",
          evidence_type = "full_data_refit",
          validation_evidence_type = development_result$evidence_type,
          performance = development_result$performance,
          performance_by_study = development_result$performance_by_study,
          note = paste0(
            "Full-data refit for deployment.  The `performance` field ",
            "above is the nested cross-validation estimate from model ",
            "development; it describes the development PROCEDURE, not ",
            "this refitted model, whose apparent performance on the ",
            "development data would be optimistic."
          )
        )
      ))
    }
    if (identical(mode, "exploratory")) {
      return(list(
        fit = development_result$fit,
        development = list(
          source = "stabl_analysis_exploratory",
          evidence_type = "full_data_refit",
          validation_evidence_type = development_result$evidence_type,
          performance = NULL,
          note = paste0(
            "Full-data refit for deployment.  Development was ",
            "exploratory (post-selection descriptive); no held-out ",
            "performance estimate exists."
          )
        )
      ))
    }
  }
  stop(
    paste0(
      "`development_result` must be a \"stabl_fit\" or a \"stabl_analysis\" ",
      "(see stabl_fit() and stabl_analyze())."
    ),
    call. = FALSE
  )
}

# Pack an engine result into a minimal, training-data-free model store.
# One-feature glm fits are reduced to their coefficient vector: predict.glm
# requires the qr and residuals slots, which encode the training data, so
# the glm object itself is never stored.  cv.glmnet objects store no raw
# data and are kept whole.
.finalize_pack_model <- function(engine_result, family) {
  fallback <- engine_result$fallback
  model <- engine_result$model

  if (!is.na(fallback) && identical(fallback, "training_prevalence")) {
    return(list(
      engine = "training_prevalence",
      model = list(
        prevalence = unname(engine_result$probability[1L]),
        majority_class = as.character(engine_result$predicted[1L])
      ),
      learner = list(engine = "training_prevalence", family = family,
                     fallback = fallback)
    ))
  }
  if (!is.na(fallback) && identical(fallback, "training_mean")) {
    return(list(
      engine = "training_mean",
      model = list(intercept = unname(model$intercept)),
      learner = list(engine = "training_mean", family = family,
                     fallback = fallback)
    ))
  }
  if (inherits(model, "glm")) {
    tag <- if (identical(family, "gaussian")) "glm_gaussian" else "glm_logit"
    return(list(
      engine = tag,
      model = list(coefficients = unname(stats::coef(model))),
      learner = list(engine = tag, family = family,
                     link = if (identical(family, "gaussian")) {
                       "identity"
                     } else {
                       "logit"
                     },
                     fallback = NA_character_)
    ))
  }
  if (inherits(model, "cv.glmnet")) {
    return(list(
      engine = "cvglmnet",
      model = model,
      learner = list(engine = "cvglmnet", family = family, s = "lambda.min",
                     type_measure = if (identical(family, "gaussian")) {
                       "mse"
                     } else {
                       "class"
                     },
                     fallback = NA_character_)
    ))
  }
  if (is.list(model) && !is.null(model$centroids)) {
    return(list(
      engine = "centroid",
      model = list(centroids = model$centroids),
      learner = list(engine = "centroid", family = family,
                     fallback = NA_character_)
    ))
  }
  stop("Unrecognized downstream engine result; cannot pack model.",
       call. = FALSE)
}

#' Validate new data against a `"stabl_predictor"` schema
#'
#' Checks `newdata` against the column schema recorded at
#' [stabl_finalize()] time and returns the required columns as a numeric
#' matrix in training order.  Column names must be present, non-empty and
#' unique; missing or duplicated required columns are always an error;
#' non-numeric columns are an error; extra columns are an error unless
#' `extra_columns = "drop"` is given explicitly.
#'
#' @param object A `"stabl_predictor"` from [stabl_finalize()].
#' @param newdata A data frame or numeric matrix with column names.
#' @param extra_columns Policy for columns not required by the model:
#'   `"error"` (default) or `"drop"` (explicit opt-in; the extra columns
#'   are discarded).
#'
#' @return A numeric matrix with exactly the required columns in training
#'   order.  For an intercept-only model (empty selected panel, no
#'   preprocessing) the matrix has zero columns and `nrow(newdata)` rows:
#'   the model ignores all feature columns.
#'
#' @examples
#' set.seed(2L)
#' ids <- paste0("s", 1:60)
#' x <- matrix(rnorm(60 * 6), 60, 6,
#'             dimnames = list(ids, paste0("f", 1:6)))
#' y <- setNames(2 * x[, 1] - x[, 2] + rnorm(60), ids)
#' fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
#'                  family = "gaussian", n_bootstraps = 20L, random_state = 2L)
#' final <- stabl_finalize(fit, x, y, random_state = 2L)
#' req <- final$required_columns
#' stablr:::validate_newdata_schema(final, x[1:3, rev(req)])  # order does not matter
#' @keywords internal
validate_newdata_schema <- function(object, newdata,
                                    extra_columns = c("error", "drop")) {
  if (!inherits(object, "stabl_predictor")) {
    stop("`object` must be a \"stabl_predictor\" (see stabl_finalize()).",
         call. = FALSE)
  }
  extra_columns <- match.arg(extra_columns)
  if (is.data.frame(newdata)) {
    non_numeric <- names(newdata)[!vapply(newdata, is.numeric, logical(1L))]
    if (length(non_numeric)) {
      stop(sprintf(
        "`newdata` columns must be numeric; non-numeric: %s.",
        paste(head(non_numeric, 5L), collapse = ", ")
      ), call. = FALSE)
    }
    newdata <- as.matrix(newdata)
  }
  if (!is.matrix(newdata)) {
    stop("`newdata` must be a data frame or matrix with column names.",
         call. = FALSE)
  }
  if (!is.numeric(newdata)) {
    stop("`newdata` must be numeric.", call. = FALSE)
  }
  cn <- colnames(newdata)
  if (is.null(cn) || any(!nzchar(cn))) {
    stop("`newdata` must have non-empty column (feature) names.",
         call. = FALSE)
  }
  dup <- unique(cn[duplicated(cn)])
  if (length(dup)) {
    stop(sprintf(
      "`newdata` column names must be unique; duplicated: %s.",
      paste(head(dup, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  required <- object$required_columns
  if (length(required) == 0L) {
    # Intercept-only model: every column is ignored by definition.
    return(matrix(NA_real_, nrow = nrow(newdata), ncol = 0L,
                  dimnames = list(rownames(newdata), character(0L))))
  }
  missing_feats <- setdiff(required, cn)
  if (length(missing_feats)) {
    stop(sprintf(
      "`newdata` is missing %d required feature(s): %s.",
      length(missing_feats),
      paste(head(missing_feats, 5L), collapse = ", ")
    ), call. = FALSE)
  }
  extra <- setdiff(cn, required)
  if (length(extra) && identical(extra_columns, "error")) {
    stop(sprintf(
      paste0(
        "`newdata` has %d unexpected column(s): %s.  ",
        "Use extra_columns = \"drop\" to discard them explicitly."
      ),
      length(extra), paste(head(extra, 5L), collapse = ", ")
    ), call. = FALSE)
  }
  newdata[, required, drop = FALSE]
}

#' @describeIn stabl_finalize Predict from a portable `"stabl_predictor"`.
#'   The stored preprocessing transform is applied unchanged (imputation
#'   and normalization parameters are never recomputed from new data) and
#'   columns are matched by name, so input column order does not matter.
#' @param object A `"stabl_predictor"` from [stabl_finalize()].
#' @param newdata Development-schema new data (see
#'   [validate_newdata_schema()]).
#' @param type Prediction type: `"response"` (numeric; the only type for
#'   gaussian), `"class"` (hard class; binomial default) or
#'   `"probability"` (event probability for the stored `event_class`;
#'   binomial only).
#' @param extra_columns Policy for unexpected columns; see
#'   [validate_newdata_schema()].
#' @param ... Ignored.
#' @export
predict.stabl_predictor <- function(object, newdata, type = NULL,
                                    extra_columns = c("error", "drop"),
                                    ...) {
  if (!inherits(object, "stabl_predictor")) {
    stop("`object` must be a \"stabl_predictor\" (see stabl_finalize()).",
         call. = FALSE)
  }
  family <- object$family
  gaussian <- identical(family, "gaussian")
  type <- if (is.null(type)) {
    if (gaussian) "response" else "class"
  } else {
    match.arg(type, if (gaussian) "response" else c("class", "probability"))
  }
  if (identical(type, "probability") && gaussian) {
    stop("`type = \"probability\"` is only defined for binomial models.",
         call. = FALSE)
  }

  x_new <- validate_newdata_schema(object, newdata,
                                   extra_columns = extra_columns)
  if (!is.null(object$preprocessing)) {
    x_new <- stabl_preprocess_apply(object$preprocessing, x_new)
  }
  n <- nrow(x_new)
  x_sel <- x_new[, object$selected_features, drop = FALSE]

  engine <- object$engine
  model <- object$model
  levels <- object$levels

  # p2 is always P(second level); the event probability maps through the
  # stored event class.
  p2 <- NULL
  pred_class <- NULL
  pred_response <- NULL

  if (identical(engine, "training_mean")) {
    pred_response <- rep(model$intercept, n)
  } else if (identical(engine, "training_prevalence")) {
    p_event_const <- model$prevalence
    pred_class <- rep(model$majority_class, n)
    # Reconstruct P(levels[2]) from the stored EVENT probability.
    p2 <- if (identical(object$event_class, levels[1L])) {
      1 - p_event_const
    } else {
      p_event_const
    }
    p2 <- rep(p2, n)
  } else if (identical(engine, "glm_gaussian")) {
    pred_response <- model$coefficients[1L] +
      model$coefficients[2L] * as.numeric(x_sel[, 1L])
  } else if (identical(engine, "glm_logit")) {
    p2 <- stats::plogis(
      model$coefficients[1L] + model$coefficients[2L] * as.numeric(x_sel[, 1L])
    )
  } else if (identical(engine, "centroid")) {
    d <- abs(outer(as.numeric(x_sel[, 1L]), model$centroids, "-"))
    pred_class <- names(model$centroids)[max.col(-d, ties.method = "first")]
  } else if (identical(engine, "cvglmnet")) {
    s <- object$learner$s
    if (gaussian) {
      pred_response <- as.numeric(stats::predict(
        model, newx = x_sel, s = s
      )[, 1L])
    } else {
      p2 <- as.numeric(stats::predict(
        model, newx = x_sel, s = s, type = "response"
      )[, 1L])
    }
  } else {
    stop("Unknown engine in \"stabl_predictor\" object.", call. = FALSE)
  }

  if (gaussian) {
    return(unname(pred_response))
  }
  if (is.null(pred_class)) {
    pred_class <- ifelse(p2 >= 0.5, levels[2L], levels[1L])
  }
  if (identical(type, "probability")) {
    p_event <- if (identical(object$event_class, levels[1L])) 1 - p2 else p2
    return(unname(p_event))
  }
  factor(pred_class, levels = levels)
}

#' @describeIn stabl_finalize Print a concise summary of a
#'   `"stabl_predictor"`; invisibly returns `x`.
#' @export
print.stabl_predictor <- function(x, ...) {
  cat("<stabl_predictor:", x$family, ">\n")
  cat("  engine:      ", x$engine, "\n", sep = "")
  cat("  panel:       ", length(x$selected_features), " of ",
      nrow(x$feature_provenance), " features\n", sep = "")
  if (length(x$selected_features)) {
    cat("  selected:    ",
        paste(utils::head(x$selected_features, 5L), collapse = ", "),
        if (length(x$selected_features) > 5L) ", ..." else "", "\n",
        sep = "")
  }
  cat("  requires:    ", length(x$required_columns),
      " input column(s)\n", sep = "")
  if (!is.null(x$event_class)) {
    cat("  event class: ", x$event_class, "\n", sep = "")
  }
  cat("  preprocessing: ", if (is.null(x$preprocessing)) {
    "identity"
  } else {
    paste(x$preprocessing$transforms, collapse = ", ")
  }, "\n", sep = "")
  cat("  development: ", x$development$source,
      " (full-data refit; not cross-validated performance)\n", sep = "")
  cat("  stablr:      ", x$stablr_version, "\n", sep = "")
  invisible(x)
}
