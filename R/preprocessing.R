# preprocessing.R — minimal fold-safe preprocessing (PR2).
#
# Contract: every transform is a fit/apply pair.  `stabl_preprocess_fit()`
# learns parameters from TRAINING data only; `stabl_preprocess_apply()`
# applies them to training or assessment data without refitting.  A transform
# fitted on training data must be completely insensitive to assessment
# values: an outlier or shift in the assessment set must not change the
# fitted parameters.  `NULL` preprocessing leaves the input unchanged.
#
# MVP transforms: numeric centering/scaling, median imputation, and
# training-only zero-variance filtering.  An optional `recipes::prep()` /
# `recipes::bake()` path is provided when the `recipes` package is available
# (Suggests, not Imports).

# ── Fit ──────────────────────────────────────────────────────────────────────

#' Fit a fold-safe preprocessing transform on training data
#'
#' Learns preprocessing parameters from a training matrix only.  The result
#' is a `"stabl_preprocessing"` object that [stabl_preprocess_apply()] applies
#' to training or assessment data without refitting.
#'
#' @param x_train Numeric training matrix (samples × features) with column
#'   names.
#' @param transforms Character vector naming the transforms to fit, in order.
#'   Supported: `"zv"` (zero-variance filter), `"impute_median"` (median
#'   imputation), `"center"`, `"scale"`.  `NULL` or `character(0)` fits an
#'   identity transform that leaves input unchanged.
#' @param zv_tol Non-negative numeric; a feature whose training variance is
#'   `<= zv_tol` is removed by `"zv"`.
#' @param ... Reserved for future use.
#'
#' @return A `"stabl_preprocessing"` object recording the fitted parameters,
#'   the retained feature names, and the feature mapping.
#'
#' @examples
#' x <- matrix(rnorm(40 * 5), 40, 5,
#'             dimnames = list(paste0("s", 1:40), paste0("f", 1:5)))
#' pp <- stabl_preprocess_fit(x, transforms = c("center", "scale"))
#' xt <- stabl_preprocess_apply(pp, x)
#' @export
stabl_preprocess_fit <- function(x_train,
                                 transforms = c("center", "scale"),
                                 zv_tol = 0,
                                 ...) {
  x_train <- as.matrix(x_train)
  if (!is.numeric(x_train)) {
    stop("`x_train` must be numeric.", call. = FALSE)
  }
  if (is.null(colnames(x_train))) {
    stop("`x_train` must have column (feature) names.", call. = FALSE)
  }
  transforms <- transforms %||% character()
  supported <- c("zv", "impute_median", "center", "scale")
  bad <- setdiff(transforms, supported)
  if (length(bad)) {
    stop(sprintf("Unsupported transform(s): %s. Supported: %s",
                 paste(bad, collapse = ", "),
                 paste(supported, collapse = ", ")),
         call. = FALSE)
  }
  zv_tol <- .validate_scalar_numeric(zv_tol, "zv_tol", min = 0)

  feats <- colnames(x_train)
  params <- list()

  # Zero-variance filter (training-only decision).
  if ("zv" %in% transforms) {
    vars <- apply(x_train, 2L, stats::var, na.rm = TRUE)
    keep <- which(is.finite(vars) & vars > zv_tol)
    params$zv <- list(retained = feats[keep])
  } else {
    params$zv <- list(retained = feats)
  }
  retained <- params$zv$retained
  xr <- x_train[, retained, drop = FALSE]

  # Median imputation (training medians).
  if ("impute_median" %in% transforms) {
    meds <- apply(xr, 2L, stats::median, na.rm = TRUE)
    params$impute_median <- list(medians = meds)
  }

  # Centering (training means, after imputation).
  if ("center" %in% transforms) {
    xc <- .preprocess_impute_training(xr, params)
    params$center <- list(means = colMeans(xc, na.rm = TRUE))
  }

  # Scaling (training SDs, after imputation and centering).
  if ("scale" %in% transforms) {
    xs <- .preprocess_impute_training(xr, params)
    if ("center" %in% transforms) {
      xs <- sweep(xs, 2L, params$center$means, "-")
    }
    sds <- apply(xs, 2L, stats::sd, na.rm = TRUE)
    # Guard against zero-SD features (e.g. constant after imputation): use 1
    # so scaling is a no-op rather than a division by zero.
    sds[!is.finite(sds) | sds <= .Machine$double.eps] <- 1
    params$scale <- list(sds = sds)
  }

  out <- list(
    transforms = transforms,
    params = params,
    retained_features = retained,
    feature_mapping = setNames(seq_along(retained), retained),
    n_features_in = length(feats),
    zv_tol = zv_tol
  )
  class(out) <- "stabl_preprocessing"
  out
}

# Fill training NAs with fitted medians (used to derive center/scale params).
.preprocess_impute_training <- function(xr, params) {
  if (is.null(params$impute_median)) return(xr)
  meds <- params$impute_median$medians
  for (j in seq_len(ncol(xr))) {
    idx <- is.na(xr[, j])
    if (any(idx) && is.finite(meds[j])) xr[idx, j] <- meds[j]
  }
  xr
}

# ── Apply (S3 generic) ───────────────────────────────────────────────────────

#' Apply a fitted preprocessing transform
#'
#' Applies a transform fitted by [stabl_preprocess_fit()] to new data.  The
#' fitted parameters are used unchanged; assessment values never influence
#' them.
#'
#' @param object A `"stabl_preprocessing"` object.
#' @param x Numeric matrix (samples × features) with column names.  Must
#'   contain every retained feature; extra columns are dropped.
#' @param ... Reserved for future use.
#'
#' @return A numeric matrix with the retained features, transformed.
#'
#' @examples
#' x <- matrix(rnorm(40 * 5), 40, 5,
#'             dimnames = list(paste0("s", 1:40), paste0("f", 1:5)))
#' pp <- stabl_preprocess_fit(x, transforms = c("center", "scale"))
#' xt <- stabl_preprocess_apply(pp, x)
#' @export
stabl_preprocess_apply <- function(object, x, ...) {
  UseMethod("stabl_preprocess_apply")
}

#' @export
stabl_preprocess_apply.stabl_preprocessing <- function(object, x, ...) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    stop("`x` must have column (feature) names.", call. = FALSE)
  }
  retained <- object$retained_features
  missing_feats <- setdiff(retained, colnames(x))
  if (length(missing_feats)) {
    stop(sprintf("`x` is missing %d required feature(s): %s",
                 length(missing_feats),
                 paste(head(missing_feats, 5L), collapse = ", ")),
         call. = FALSE)
  }
  x <- x[, retained, drop = FALSE]
  p <- object$params

  if ("impute_median" %in% object$transforms) {
    meds <- p$impute_median$medians
    for (j in seq_along(retained)) {
      # Fail closed whenever training had no information for a feature:
      # imputation cannot invent a training median, and sweeping the
      # downstream center/scale parameters (themselves NA) would silently
      # emit an all-NA column.
      if (!is.finite(meds[j])) {
        stop(sprintf(
          "Feature %s is entirely missing in training and cannot be imputed.",
          retained[j]),
          call. = FALSE)
      }
      idx <- is.na(x[, j])
      if (any(idx)) x[idx, j] <- meds[j]
    }
  }
  if ("center" %in% object$transforms) {
    x <- sweep(x, 2L, p$center$means, "-")
  }
  if ("scale" %in% object$transforms) {
    x <- sweep(x, 2L, p$scale$sds, "/")
  }
  x
}

# ── Optional recipes integration (Suggests) ──────────────────────────────────

#' Fit a preprocessing transform from a recipes recipe
#'
#' Wraps [recipes::prep()] on training data only.  Requires the `recipes`
#' package (in `Suggests`, not `Imports`).
#'
#' @param recipe A `recipes::recipe()` object.
#' @param x_train Training data frame or matrix.
#' @param ... Passed to [recipes::prep()].
#'
#' @return A `"stabl_preprocessing_recipes"` object wrapping the prepped
#'   recipe.
#'
#' @examples
#' if (requireNamespace("recipes", quietly = TRUE)) {
#'   x <- matrix(rnorm(40 * 4), 40, 4,
#'               dimnames = list(paste0("s", 1:40), paste0("f", 1:4)))
#'   rec <- recipes::recipe(~ ., data = as.data.frame(x))
#'   rec <- recipes::step_center(rec, recipes::all_predictors())
#'   pp <- stabl_preprocess_fit_recipes(rec, x)
#'   xt <- stabl_preprocess_apply(pp, x)
#' }
#' @export
stabl_preprocess_fit_recipes <- function(recipe, x_train, ...) {
  if (!requireNamespace("recipes", quietly = TRUE)) {
    stop("The `recipes` package is required for recipe-based preprocessing. ",
         "Install it or use stabl_preprocess_fit().", call. = FALSE)
  }
  prepped <- recipes::prep(recipe, training = as.data.frame(x_train), ...)
  out <- list(
    prepped = prepped,
    retained_features = names(recipes::bake(prepped, new_data = NULL)),
    transforms = "recipes"
  )
  class(out) <- c("stabl_preprocessing_recipes", "stabl_preprocessing")
  out
}

#' @export
stabl_preprocess_apply.stabl_preprocessing_recipes <- function(object, x, ...) {
  if (!requireNamespace("recipes", quietly = TRUE)) {
    stop("The `recipes` package is required.", call. = FALSE)
  }
  baked <- recipes::bake(object$prepped, new_data = as.data.frame(x))
  as.matrix(baked)
}

# ── Fold-level specification handling (used by nested CV / workflows) ───────

# Validate a user-supplied preprocessing SPECIFICATION.  A specification is
# `NULL` (identity), a character vector of MVP transforms, or an unprepped
# `recipes::recipe`.  A fitted `stabl_preprocessing` object is rejected: a
# transform fitted outside the fold loop would leak information across folds.
.validate_preprocessing_spec <- function(preprocessing, arg = "`preprocessing`") {
  if (is.null(preprocessing)) return(invisible(NULL))
  if (inherits(preprocessing, "stabl_preprocessing")) {
    stop(sprintf(
      paste0(
        "%s must be an unfitted specification (a character vector of ",
        "transforms or an unprepped recipes::recipe), not a fitted ",
        "`stabl_preprocessing` object: a transform fitted outside the fold ",
        "loop would leak information across folds."
      ),
      arg
    ), call. = FALSE)
  }
  if (is.character(preprocessing)) {
    supported <- c("zv", "impute_median", "center", "scale")
    bad <- setdiff(preprocessing, supported)
    if (length(bad)) {
      stop(sprintf(
        "%s names unsupported transform(s): %s. Supported: %s",
        arg, paste(bad, collapse = ", "), paste(supported, collapse = ", ")
      ), call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (inherits(preprocessing, "recipe")) {
    if (isTRUE(preprocessing$trained)) {
      stop(sprintf(
        paste0(
          "%s must be an unprepped recipe; a prepped recipe carries ",
          "full-data parameters and would leak information across folds."
        ),
        arg
      ), call. = FALSE)
    }
    return(invisible(NULL))
  }
  stop(sprintf(
    paste0(
      "%s must be NULL, a character vector of transforms ",
      "(\"zv\", \"impute_median\", \"center\", \"scale\"), or an unprepped ",
      "recipes::recipe."
    ),
    arg
  ), call. = FALSE)
}

# Fit a preprocessing specification on one fold's training matrix.  Returns
# NULL for the identity (NULL) specification.
.fit_fold_preprocessing <- function(preprocessing, x_train) {
  .validate_preprocessing_spec(preprocessing)
  if (is.null(preprocessing)) return(NULL)
  if (is.character(preprocessing)) {
    return(stabl_preprocess_fit(x_train, transforms = preprocessing))
  }
  stabl_preprocess_fit_recipes(preprocessing, x_train)
}
