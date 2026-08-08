#' Get the Feature Selection Mask from a Fitted STABL Object
#'
#' Returns a named logical vector that is `TRUE` for every feature whose
#' maximum stability score exceeds the effective threshold.  This is the
#' primary accessor for downstream use of a fitted STABL model: index your
#' data matrix with the returned mask, or pass the object to
#' [get_feature_names_out()] to obtain the names directly.
#'
#' @details
#' **Threshold resolution order:**
#' 1. `new_hard_threshold` (if supplied to this call).
#' 2. `object$hard_threshold` (if a hard threshold was given to [stabl_fit()]).
#' 3. `object$fdr_min_threshold_` (the FDP+-optimal threshold computed from
#'    artificial features during fitting).
#'
#' An error is raised when none of these is available — which can only happen
#' if `stabl_fit()` was called with both `artificial_type = NULL` and no
#' `hard_threshold`.
#'
#' **Explore fallback:** When `explore = TRUE` was set during fitting and no
#' feature's score exceeds the threshold, the function returns the top
#' `n_explore` features instead of an all-`FALSE` vector.  This is useful in
#' exploratory analyses where you want at least some candidates even when the
#' signal is weak.
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param new_hard_threshold Numeric in `(0, 1]` or `NULL`.  When supplied,
#'   overrides the threshold stored in `object` for this call only.
#'
#' @return Named logical vector of length `object$n_features_in_`.  Names are
#'   the original feature names from the training matrix.
#'
#' @seealso [get_feature_names_out()] to get names directly,
#'   [get_importances()] to inspect raw stability scores.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(30 * 5), 30, 5,
#'             dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
#' y <- setNames(rnorm(30), rownames(x))
#' fit <- stabl_fit(x, y,
#'                  lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
#'                  n_bootstraps   = 4L,
#'                  hard_threshold = 0.3,
#'                  random_state   = 1L)
#' get_support(fit)                        # named logical vector
#' get_support(fit, new_hard_threshold = 0.5)  # stricter threshold
#' @export
get_support <- function(object, new_hard_threshold = NULL) {
  UseMethod("get_support")
}

#' @export
get_support.stabl_fit <- function(object, new_hard_threshold = NULL) {
  .check_fitted_stabl(object)

  threshold  <- .resolve_stabl_threshold(object, new_hard_threshold)
  max_scores <- get_importances(object)
  mask       <- max_scores > threshold

  # explore fallback: top n_explore features if nothing passes threshold.
  # Use direct index selection to avoid over-selecting on tied scores
  # (the old `sort()[n_exp] - 0.01` approach selected ALL features tied at
  # the n_exp-th score, e.g. all-zero scores produced mask = all TRUE).
  if (!any(mask) && isTRUE(object$explore)) {
    n_exp       <- min(object$n_explore, length(max_scores))
    top_idx     <- order(max_scores, decreasing = TRUE)[seq_len(n_exp)]
    mask[top_idx] <- TRUE
  }

  mask
}

#' Get the Full Stability Score Matrix from a Fitted STABL Object
#'
#' Returns the raw stability-score matrix accumulated over all bootstrap
#' iterations.  Inspecting this matrix is useful for diagnostics: you can
#' check how stable features are across the regularisation path, identify
#' features that are consistently selected at many lambda values (robustness),
#' and spot features that peak only at one extreme of the path (fragility).
#'
#' The stability score for feature \eqn{j} at regularisation strength
#' \eqn{\lambda_k} is defined as the fraction of bootstrap subsamples in
#' which feature \eqn{j} received a non-zero coefficient when the model was
#' fitted at \eqn{\lambda_k}.  Values lie in \eqn{[0, 1]}.
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#'
#' @return Numeric matrix with one row per original feature and one column per
#'   lambda in the fitted grid.  Row names are the feature names from the
#'   training matrix; column names are not set (use
#'   `object$fitted_lambda_grid` to map column indices to lambda values).
#'
#' @seealso [get_importances()] for the per-feature maximum score (scalar
#'   summary), [get_support()] for the binary selection mask.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(30 * 5), 30, 5,
#'             dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
#' y <- setNames(rnorm(30), rownames(x))
#' fit <- stabl_fit(x, y,
#'                  lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
#'                  n_bootstraps   = 4L,
#'                  hard_threshold = 0.3,
#'                  random_state   = 1L)
#' scores <- get_stabl_scores(fit)  # features x lambdas matrix
#' dim(scores)
#' @export
get_stabl_scores <- function(object) {
  UseMethod("get_stabl_scores")
}

#' @export
get_stabl_scores.stabl_fit <- function(object) {
  .check_fitted_stabl(object)
  object$stabl_scores_
}

#' Get the Names of Selected Features from a Fitted STABL Object
#'
#' Convenience wrapper around [get_support()] that returns only the names of
#' the features that pass the stability threshold, ready for use as column
#' selectors in downstream modelling (e.g.
#' `x[, get_feature_names_out(fit)]`).
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param new_hard_threshold Numeric in `(0, 1]` or `NULL`.  Forwarded to
#'   [get_support()]; see that function for the full threshold resolution
#'   order.
#'
#' @return Character vector of selected feature names.  An empty character
#'   vector is returned when no feature passes the threshold (and
#'   `explore = FALSE` was used during fitting).
#'
#' @seealso [get_support()] for the binary mask, [get_importances()] for
#'   ranked scores.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(30 * 5), 30, 5,
#'             dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
#' y <- setNames(rnorm(30), rownames(x))
#' fit <- stabl_fit(x, y,
#'                  lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
#'                  n_bootstraps   = 4L,
#'                  hard_threshold = 0.3,
#'                  random_state   = 1L)
#' sel <- get_feature_names_out(fit)   # character vector
#' x_sel <- x[, sel, drop = FALSE]    # subset to selected columns
#' @export
get_feature_names_out <- function(object, new_hard_threshold = NULL) {
  UseMethod("get_feature_names_out")
}

#' @export
get_feature_names_out.stabl_fit <- function(object, new_hard_threshold = NULL) {
  mask <- get_support(object, new_hard_threshold = new_hard_threshold)
  object$feature_names[mask]
}

#' Transform New Data to the Selected STABL Feature Set
#'
#' Subsets a matrix or data frame to the features selected by a fitted
#' [stabl_fit()] object. Columns are returned in fitted feature order, rows are
#' preserved, and the output remains two-dimensional even when no features are
#' selected.
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param x A matrix or `data.frame` with feature names in columns.
#' @param new_hard_threshold Numeric in `(0, 1]` or `NULL`. When supplied,
#'   overrides the fitted threshold for this transformation only.
#'
#' @return A matrix when `x` is a matrix, or a `data.frame` when `x` is a
#'   `data.frame`, containing selected columns in fitted feature order.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(
#'   rnorm(30 * 5), 30, 5,
#'   dimnames = list(paste0("s", 1:30), paste0("f", 1:5))
#' )
#' y <- setNames(rnorm(30), rownames(x))
#' fit <- stabl_fit(
#'   x, y, lambda_grid = data.frame(lambda = c(0.2, 0.1)),
#'   n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
#'   random_state = 1L
#' )
#' transformed <- transform_stabl(fit, x)
#' dim(transformed)
#' try(transform_stabl(fit, unname(x)))
#'
#' @seealso [get_feature_names_out()], [get_support()]
#' @export
transform_stabl <- function(object, x, new_hard_threshold = NULL) {
  .check_fitted_stabl(object)
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("`x` must be a matrix or data.frame.", call. = FALSE)
  }
  feature_names <- colnames(x)
  if (is.null(feature_names) || anyNA(feature_names) || any(feature_names == "")) {
    stop("`x` must have non-empty column names.", call. = FALSE)
  }
  if (anyDuplicated(feature_names)) {
    stop("`x` column names must not contain duplicates.", call. = FALSE)
  }

  selected <- get_feature_names_out(object, new_hard_threshold = new_hard_threshold)
  missing <- setdiff(selected, feature_names)
  if (length(missing) > 0L) {
    stop(
      "`x` is missing selected feature columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x[, selected, drop = FALSE]
}

#' Get Cooperative-Fusion Selected Features
#'
#' Returns the feature names selected by the cooperative-fusion branch of a
#' fitted multi-omic workflow.  This accessor gives downstream code a stable
#' public surface instead of reaching into `$cooperative_fusion` directly.
#'
#' @param object A fitted `"stabl_multiomic_fit"` object returned by
#'   [stabl_multiomic_train_validate()] with `cooperative_fusion = TRUE`, or a
#'   `"stabl_multiomic_cv"` object returned by [stabl_multiomic_cv()] with
#'   `cooperative_fusion = TRUE`.
#' @param view Optional character scalar naming one omic view.  When supplied
#'   for a `"stabl_multiomic_fit"`, only that view's cooperative feature names
#'   are returned.
#'
#' @return For `"stabl_multiomic_fit"`, a named list of character vectors, or
#'   a character vector when `view` is supplied.  For `"stabl_multiomic_cv"`, a
#'   named list keyed by fold, where each element has the same structure as the
#'   `"stabl_multiomic_fit"` method.
#'
#' @seealso [get_cooperative_diagnostics()],
#'   [stabl_multiomic_train_validate()], [stabl_multiomic_cv()]
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' n <- 40L
#' x_list <- list(
#'   omic1 = matrix(rnorm(n * 5), n, 5,
#'                  dimnames = list(paste0("s", 1:n), paste0("g", 1:5))),
#'   omic2 = matrix(rnorm(n * 4), n, 4,
#'                  dimnames = list(paste0("s", 1:n), paste0("p", 1:4)))
#' )
#' y <- setNames(rnorm(n), paste0("s", 1:n))
#' train <- 1:30
#' valid <- 31:40
#' x_train <- lapply(x_list, function(x) x[train, , drop = FALSE])
#' x_valid <- lapply(x_list, function(x) x[valid, , drop = FALSE])
#' lam <- data.frame(lambda = c(0.3, 0.1))
#' fit <- stabl_multiomic_train_validate(
#'   x_train, y[train], lambda_grid = lam,
#'   x_valid_list = x_valid, y_valid = y[valid],
#'   family = "gaussian",
#'   n_bootstraps = 4L,
#'   hard_threshold = 0.2, cooperative_fusion = TRUE,
#'   random_state = 1L
#' )
#' get_cooperative_features(fit)           # named list of feature vectors
#' get_cooperative_features(fit, "omic1") # one view only
#' }
#' @export
get_cooperative_features <- function(object, view = NULL) {
  UseMethod("get_cooperative_features")
}

#' @export
get_cooperative_features.stabl_multiomic_fit <- function(object, view = NULL) {
  cf <- .check_cooperative_branch(object)
  features <- cf$selected_features

  if (is.null(view)) {
    return(features)
  }

  if (length(view) != 1L || is.na(view) || !(view %in% names(features))) {
    stop(
      "`view` must name one of the cooperative-fusion omic views.",
      call. = FALSE
    )
  }

  features[[view]]
}

#' @export
get_cooperative_features.stabl_multiomic_cv <- function(object, view = NULL) {
  out <- lapply(object$fold_results, get_cooperative_features, view = view)
  names(out) <- names(object$fold_results)
  out
}

#' Get Cooperative-Fusion Tuning Diagnostics
#'
#' Returns the cooperative-fusion tuning diagnostics from a fitted multi-omic
#' workflow.  For train/validation fits this is the per-candidate tuning table;
#' for outer cross-validation fits this is the fold diagnostics table restricted
#' to cooperative diagnostic columns.
#'
#' @param object A `"stabl_multiomic_fit"` or `"stabl_multiomic_cv"` object
#'   with cooperative fusion enabled.
#'
#' @return A `data.frame` of cooperative tuning diagnostics.
#'
#' @seealso [get_cooperative_features()],
#'   [stabl_multiomic_train_validate()], [stabl_multiomic_cv()]
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' n <- 40L
#' x_list <- list(
#'   omic1 = matrix(rnorm(n * 5), n, 5,
#'                  dimnames = list(paste0("s", 1:n), paste0("g", 1:5))),
#'   omic2 = matrix(rnorm(n * 4), n, 4,
#'                  dimnames = list(paste0("s", 1:n), paste0("p", 1:4)))
#' )
#' y <- setNames(rnorm(n), paste0("s", 1:n))
#' train <- 1:30
#' valid <- 31:40
#' x_train <- lapply(x_list, function(x) x[train, , drop = FALSE])
#' x_valid <- lapply(x_list, function(x) x[valid, , drop = FALSE])
#' lam <- data.frame(lambda = c(0.3, 0.1))
#' fit <- stabl_multiomic_train_validate(
#'   x_train, y[train], lambda_grid = lam,
#'   x_valid_list = x_valid, y_valid = y[valid],
#'   family = "gaussian",
#'   n_bootstraps = 4L,
#'   hard_threshold = 0.2, cooperative_fusion = TRUE,
#'   random_state = 1L
#' )
#' get_cooperative_diagnostics(fit)   # data.frame of tuning results
#' }
#' @export
get_cooperative_diagnostics <- function(object) {
  UseMethod("get_cooperative_diagnostics")
}

#' @export
get_cooperative_diagnostics.stabl_multiomic_fit <- function(object) {
  cf <- .check_cooperative_branch(object)
  cf$diagnostics
}

#' @export
get_cooperative_diagnostics.stabl_multiomic_cv <- function(object) {
  diagnostics <- object$diagnostics
  cooperative_cols <- grep("^cooperative_", names(diagnostics), value = TRUE)

  if (length(cooperative_cols) == 0L) {
    stop(
      "No cooperative-fusion diagnostics available; fit with `cooperative_fusion = TRUE`.",
      call. = FALSE
    )
  }

  diagnostics[, c("fold", "omic", cooperative_cols), drop = FALSE]
}

#' Get Per-Feature Importance Scores (Maximum Stability Score)
#'
#' Returns a scalar summary of how stably each feature is selected across the
#' entire regularisation path.  The importance of feature \eqn{j} is defined
#' as \eqn{\max_{k} q_{jk}}, i.e. the highest selection frequency it achieved
#' at any lambda.  This is the score compared against the stability threshold
#' in [get_support()].
#'
#' Because the maximum is taken across lambdas, the importance measure is
#' lenient: a feature qualifies even if it is stable only at one particular
#' penalty strength.  For a more conservative view, use [get_stabl_scores()]
#' and inspect the full path.
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#'
#' @return Named numeric vector of length `object$n_features_in_`, with
#'   values in \eqn{[0, 1]}.  Higher values indicate more stable features.
#'   Names are the original feature names from the training matrix.
#'
#' @seealso [get_support()] to convert importances to a binary selection mask,
#'   [get_stabl_scores()] for the full path matrix.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(30 * 5), 30, 5,
#'             dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
#' y <- setNames(rnorm(30), rownames(x))
#' fit <- stabl_fit(x, y,
#'                  lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
#'                  n_bootstraps   = 4L,
#'                  hard_threshold = 0.3,
#'                  random_state   = 1L)
#' imp <- get_importances(fit)  # named numeric, max stability per feature
#' sort(imp, decreasing = TRUE)
#' @export
get_importances <- function(object) {
  UseMethod("get_importances")
}

#' @export
get_importances.stabl_fit <- function(object) {
  .check_fitted_stabl(object)
  scores        <- rowMaxs(object$stabl_scores_)
  names(scores) <- object$feature_names
  scores
}

#' @describeIn stabl_fit Print a concise summary of a fitted `stabl_fit` object;
#'   invisibly returns `x`.
#' @param ... Ignored; present for S3 `print` method compatibility.
#' @export
print.stabl_fit <- function(x, ...) {
  cat("<stabl_fit>\n")
  cat("  Features in:     ", x$n_features_in_, "\n")
  tryCatch(
    {
      mask <- get_support(x)
      cat("  Features selected:", sum(mask), "\n")
    },
    error = function(e) {
      cat("  Features selected: [threshold not yet resolved]\n")
    }
  )
  if (!is.null(x$min_fdr_)) {
    cat("  Min FDP+:        ", round(x$min_fdr_, 4L), "\n")
    cat("  FDP+ threshold:  ", round(x$fdr_min_threshold_, 4L), "\n")
  }
  if (!is.null(x$hard_threshold)) {
    cat("  Hard threshold:  ", x$hard_threshold, "\n")
  }
  cat("  Artificial:      ", if (is.null(x$artificial_type)) "none"
                             else x$artificial_type, "\n")
  invisible(x)
}

#' @describeIn stabl_multiomic_train_validate Print a concise summary of a
#'   `stabl_multiomic_fit` object; invisibly returns `x`.
#' @export
print.stabl_multiomic_fit <- function(x, ...) {
  omic_names <- names(x$fits)
  cat("<stabl_multiomic_fit>\n")
  cat("  Omics:           ", length(omic_names), " (", paste(omic_names, collapse = ", "), ")\n", sep = "")
  cat("  Per-omic selected features:\n")
  for (omic in omic_names) {
    n_sel <- length(x$selected_features[[omic]])
    cat("    ", omic, ": ", n_sel, "\n", sep = "")
  }
  has_valid <- !is.null(x$selected_valid)
  cat("  Validation data: ", if (has_valid) "yes" else "no", "\n")
  cat("  Early fusion:    ", if (!is.null(x$early_fusion)) {
    paste0("yes (", length(x$early_fusion$selected_features), " features selected)")
  } else "no", "\n")
  cat("  Late fusion:     ", if (!is.null(x$late_fusion)) {
    paste0("yes (score = ", round(x$late_fusion$score, 4L), ")")
  } else "no", "\n")
  if (!is.null(x$cooperative_fusion)) {
    cf <- x$cooperative_fusion
    n_cf_sel <- sum(vapply(cf$selected_features, length, integer(1L)))
    cat("  Cooperative fusion:\n")
    cat("    selection:     ", cf$selection, "\n", sep = "")
    cat("    rho (chosen):  ", cf$rho, "\n", sep = "")
    cat("    selector:      ", cf$selector, "\n", sep = "")
    cat("    type.measure:  ", cf$type_measure, "\n", sep = "")
    cat("    score:         ", round(cf$score, 4L), "\n", sep = "")
    cat("    selected feats: ", n_cf_sel, " (across views)\n", sep = "")
  }
  invisible(x)
}

#' @describeIn stabl_multiomic_cv Print a concise summary of a
#'   `stabl_multiomic_cv` object; invisibly returns `x`.
#' @export
print.stabl_multiomic_cv <- function(x, ...) {
  cat("<stabl_multiomic_cv>\n")
  cat("  Folds:           ", length(x$folds), "\n")
  if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
    omic_names <- unique(x$diagnostics$omic)
    cat("  Omics:           ", length(omic_names), " (",
        paste(omic_names, collapse = ", "), ")\n", sep = "")
    cat("  Mean selected features per omic (across folds):\n")
    for (omic in omic_names) {
      sub <- x$diagnostics[x$diagnostics$omic == omic, , drop = FALSE]
      cat("    ", omic, ": ", round(mean(sub$n_selected), 1L), "\n", sep = "")
    }
  }
  invisible(x)
}

# ---- Internal helpers --------------------------------------------------------

# Resolve the effective threshold for a fitted stabl_fit object.
# Precedence: new_hard_threshold > object$hard_threshold > object$fdr_min_threshold_.
# on_missing controls what happens when all three are NULL:
#   "error"  — stop() with the canonical accessor error message (default)
#   "null"   — return NULL silently (used when threshold is optional, e.g. plot annotations)
#   "na"     — return NA_real_ (used in multiomic threshold helpers)
.resolve_stabl_threshold <- function(object, new_hard_threshold = NULL,
                                     on_missing = c("error", "null", "na")) {
  on_missing <- match.arg(on_missing)
  threshold  <- if (!is.null(new_hard_threshold)) {
    new_hard_threshold
  } else if (!is.null(object$hard_threshold)) {
    object$hard_threshold
  } else {
    object$fdr_min_threshold_
  }

  if (!is.null(threshold)) {
    return(threshold)
  }

  switch(
    on_missing,
    error = stop(
      "No threshold available: fit with `artificial_type` or supply ",
      "`hard_threshold`.",
      call. = FALSE
    ),
    null = NULL,
    na   = NA_real_
  )
}

.check_fitted_stabl <- function(object) {
  if (is.null(object$stabl_scores_)) {
    stop(
      "Object has no `stabl_scores_`: ensure `stabl_fit()` completed.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.check_cooperative_branch <- function(object) {
  if (is.null(object$cooperative_fusion)) {
    stop(
      "No cooperative-fusion branch available; fit with `cooperative_fusion = TRUE`.",
      call. = FALSE
    )
  }

  if (is.null(object$cooperative_fusion$selected_features) ||
      is.null(object$cooperative_fusion$diagnostics)) {
    stop(
      "Malformed cooperative-fusion branch: expected `selected_features` and `diagnostics`.",
      call. = FALSE
    )
  }

  object$cooperative_fusion
}

# Null-coalescing helper (package-internal)
`%||%` <- function(a, b) if (is.null(a)) b else a
