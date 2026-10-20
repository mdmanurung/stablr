#' Validate Sample Alignment Across Inputs
#'
#' Checks that the predictor matrix, outcome vector, and optional group vector
#' all refer to exactly the same set of samples (by name), and that they can
#' be safely aligned before modelling.
#'
#' STABL enforces strict name-based alignment (rather than positional
#' alignment) to mirror the pandas index-alignment semantics of the Python
#' reference implementation and to prevent silent sample-order bugs that could
#' corrupt stability scores or introduce outcome leakage.
#'
#' This function is called automatically by [stabl_fit()] and
#' [stabl_multiomic_train_validate()]; you only need to call it directly when
#' building a custom pre-processing step that receives `x` and `y` separately.
#'
#' @param x A `data.frame` or numeric matrix with non-empty, non-`NA` row
#'   names used as sample IDs.
#' @param y A named vector (or matrix-like object such as `survival::Surv`
#'   with row names) whose names identify the samples.  The set of names in
#'   `y` must be identical to the set of row names in `x`; order does not
#'   need to match.
#' @param groups Optional named vector where names are sample IDs and values
#'   are group memberships.  When supplied, its name set must also match
#'   `rownames(x)`.  Pass `NULL` to skip group validation.
#'
#' @return Invisibly returns `TRUE` when all checks pass.  Raises an
#'   informative error as soon as the first violation is found.
#'
#' @examples
#' x <- matrix(
#'   1:12, 4, 3,
#'   dimnames = list(paste0("s", 1:4), paste0("f", 1:3))
#' )
#' y <- setNames(c(0, 1, 0, 1), c("s3", "s1", "s4", "s2"))
#' validate_sample_alignment(x, y)
#' try(validate_sample_alignment(x, setNames(y, paste0("v", 1:4))))
#'
#' @seealso [validate_multiomic_inputs()] for multi-omic list inputs.
#' @export
validate_sample_alignment <- function(x, y, groups = NULL) {
  if (!(is.data.frame(x) || is.matrix(x))) {
    .abort_stablr("stablr_input_error", "`x` must be a data.frame or matrix.")
  }

  sample_ids <- rownames(x)
  if (is.null(sample_ids) || anyNA(sample_ids) || any(sample_ids == "")) {
    .abort_stablr("stablr_input_error", "`x` must have non-empty row names used as sample ids.")
  }
  .validate_unique_names(sample_ids, "`x` row names")
  .validate_feature_names(x, "`x`")

  y_ids <- .outcome_sample_ids(y)
  if (is.null(y_ids)) {
    .abort_stablr("stablr_input_error", "`y` must be a named vector or matrix-like outcome with row names as sample ids.")
  }
  .validate_unique_names(y_ids, "`y` sample ids")

  if (!setequal(sample_ids, y_ids)) {
    .abort_stablr("stablr_input_error", "Sample mismatch between `x` row names and `y` names.")
  }

  y_ordered <- .subset_outcome_by_ids(y, sample_ids)
  if (anyNA(y_ordered)) {
    .abort_stablr("stablr_input_error", "`y` cannot contain missing sample mappings after alignment.")
  }

  if (!is.null(groups)) {
    if (is.null(names(groups))) {
      .abort_stablr("stablr_input_error", "`groups` must be a named vector where names are sample ids.")
    }
    .validate_unique_names(names(groups), "`groups` sample ids")
    if (!setequal(sample_ids, names(groups))) {
      .abort_stablr("stablr_input_error", "Sample mismatch between `x` row names and `groups` names.")
    }
    groups_ordered <- groups[sample_ids]
    if (anyNA(groups_ordered)) {
      .abort_stablr("stablr_input_error", "`groups` cannot contain missing sample mappings after alignment.")
    }
  }

  invisible(TRUE)
}

# Resolve sample identifiers from supported outcome types.
.outcome_sample_ids <- function(y) {
  # Matrix-like outcomes first: names() on a data.frame returns COLUMN
  # names (always non-NULL), which made these branches unreachable and
  # turned sample IDs into column names (review E10).
  if (is.matrix(y) || is.data.frame(y) || inherits(y, "Surv")) {
    return(rownames(y))
  }

  if (!is.null(names(y))) {
    return(names(y))
  }

  NULL
}

# Subset supported outcome types by sample id while preserving their shape.
.subset_outcome_by_ids <- function(y, sample_ids) {
  y_ids <- .outcome_sample_ids(y)
  idx <- match(sample_ids, y_ids)

  if (anyNA(idx)) {
    .abort_stablr("stablr_input_error", "`y` cannot contain missing sample mappings after alignment.")
  }

  if (is.matrix(y) || is.data.frame(y) || inherits(y, "Surv")) {
    return(y[idx, , drop = FALSE])
  }

  if (!is.null(names(y))) {
    return(y[sample_ids])
  }

  y[idx]
}

#' Validate Multi-Omic Input Contract
#'
#' Enforces the canonical `stablr` input contract for multi-omic analyses: a
#' named list of omic tables with identical sample IDs and row order across
#' all views, plus strict alignment with the outcome and optional group
#' vectors.
#'
#' Multi-omic analyses require that all omic matrices have exactly the same
#' rows in exactly the same order so that row-wise operations (bootstrap
#' sampling, cooperative learning) are coherent.  This function catches
#' mismatches early and provides informative error messages naming the
#' offending omic layer, which is much easier to diagnose than silent
#' misalignment detected later in the modelling pipeline.
#'
#' @param x_list Named list of `data.frame` or numeric matrix omic tables.
#'   Each table must have non-empty row names identifying samples.  All tables
#'   must have the same set of row names **in the same order**.
#' @param y Named outcome vector aligned to the samples in `x_list`.  The
#'   name set must match the row names of every omic table.
#' @param groups Optional named group vector.  When supplied its names must
#'   match the sample IDs in `x_list`.  Pass `NULL` to skip group validation.
#'
#' @return Invisibly returns `TRUE` when all checks pass.  Raises an
#'   informative error as soon as the first violation is found.
#'
#' @examples
#' ids <- paste0("s", 1:4)
#' x_list <- list(
#'   rna = matrix(1:12, 4, 3, dimnames = list(ids, paste0("g", 1:3))),
#'   protein = matrix(1:8, 4, 2, dimnames = list(ids, paste0("p", 1:2)))
#' )
#' y <- setNames(c(0, 1, 0, 1), ids)
#' validate_multiomic_inputs(x_list, y)
#' try(validate_multiomic_inputs(x_list[c("rna", "rna")], y))
#'
#' @seealso [validate_sample_alignment()] for single-omic inputs,
#'   [stabl_multiomic_train_validate()] which calls this automatically.
#' @export
validate_multiomic_inputs <- function(x_list, y, groups = NULL) {
  if (!is.list(x_list) || length(x_list) == 0L) {
    .abort_stablr("stablr_input_error", "`x_list` must be a non-empty named list.")
  }

  if (is.null(names(x_list)) || anyNA(names(x_list)) || any(names(x_list) == "")) {
    .abort_stablr("stablr_input_error", "`x_list` must have non-empty names for each omic table.")
  }
  .validate_unique_names(names(x_list), "omic names")

  for (name in names(x_list)) {
    x <- x_list[[name]]
    if (!(is.data.frame(x) || is.matrix(x))) {
      .abort_stablr("stablr_input_error", sprintf("Omic '%s' must be a data.frame or matrix.", name))
    }
    .validate_feature_names(x, sprintf("omic '%s'", name))
    validate_sample_alignment(x = x, y = y, groups = groups)
  }

  base_ids <- rownames(x_list[[1]])
  for (name in names(x_list)[-1]) {
    if (!identical(base_ids, rownames(x_list[[name]]))) {
      .abort_stablr("stablr_input_error", sprintf("All omic tables must have identical sample order; mismatch at omic '%s'.", name))
    }
  }

  invisible(TRUE)
}

.supported_cooperation_type_measures <- function(family) {
  switch(
    family,
    gaussian = c("default", "mse", "deviance", "mae"),
    binomial = c("default", "mse", "deviance", "class", "auc", "mae"),
    poisson = c("default", "mse", "deviance", "mae"),
    cox = c("default", "deviance", "C"),
    character(0)
  )
}

.resolve_cooperation_type_measure <- function(family, type_measure) {
  allowed <- .supported_cooperation_type_measures(family)

  if (!(type_measure %in% allowed)) {
    .abort_stablr("stablr_input_error", sprintf(
        "`cooperation_type_measure` must be one of %s for family '%s'.",
        paste(sprintf("'%s'", allowed), collapse = ", "),
        family
      ))
  }

  if (!identical(type_measure, "default")) {
    return(type_measure)
  }

  switch(
    family,
    gaussian = "mse",
    binomial = "deviance",
    poisson = "deviance",
    cox = "deviance",
    .abort_stablr("stablr_input_error", sprintf("No default type_measure for unsupported family '%s'.", family))
  )
}

.normalize_cooperative_multiomic_args <- function(x_list,
                                                  family,
                                                  rho = NULL,
                                                  cooperative_selection = c("cv", "validation"),
                                                  cooperation_selector = c("lambda.min", "lambda.1se"),
                                                  cooperation_type_measure = "default",
                                                  cooperation_nfolds = 5L,
                                                  x_valid_list = NULL,
                                                  y_valid = NULL) {
  cooperative_selection <- match.arg(cooperative_selection)
  cooperation_selector <- match.arg(cooperation_selector)

  if (!.has_cooperative_backend()) {
    .abort_stablr("stablr_input_error", "Cooperative fusion backend is unavailable.")
  }

  if (length(x_list) < 2L) {
    .abort_stablr("stablr_input_error", "`cooperative_fusion = TRUE` requires at least two omic views.")
  }

  if (!(family %in% c("gaussian", "binomial"))) {
    .abort_stablr("stablr_input_error", "`cooperative_fusion = TRUE` currently supports family = 'gaussian' or 'binomial'.")
  }

  if (is.null(rho)) {
    rho <- 0
  }

  if (!is.numeric(rho) || length(rho) == 0L || anyNA(rho) || any(!is.finite(rho))) {
    .abort_stablr("stablr_input_error", "`rho` must be a non-empty numeric scalar or vector.")
  }

  if (any(rho < 0)) {
    .abort_stablr("stablr_input_error", "`rho` must contain only non-negative values.")
  }

  cooperation_nfolds <- .validate_scalar_integer_like(
    cooperation_nfolds,
    "cooperation_nfolds",
    min = 3L
  )
  if (cooperation_nfolds < 3L) {
    .abort_stablr("stablr_input_error", "`cooperation_nfolds` must be greater than or equal to 3.")
  }

  if (identical(cooperative_selection, "validation")) {
    if (is.null(x_valid_list) || is.null(y_valid)) {
      .abort_stablr("stablr_input_error", "`cooperation_selection = 'validation'` requires both `x_valid_list` and `y_valid`.")
    }

    if (identical(cooperation_selector, "lambda.1se")) {
      .abort_stablr("stablr_input_error", "`cooperation_selector = 'lambda.1se'` is only available when `cooperation_selection = 'cv'`.")
    }
  }

  type_measure <- .resolve_cooperation_type_measure(
    family = family,
    type_measure = cooperation_type_measure
  )

  list(
    rho = unique(as.numeric(rho)),
    cooperative_selection = cooperative_selection,
    cooperation_selector = cooperation_selector,
    cooperation_type_measure = type_measure,
    cooperation_nfolds = cooperation_nfolds
  )
}

.validate_unique_names <- function(x, what) {
  if (anyDuplicated(x)) {
    .abort_stablr("stablr_input_error", sprintf("%s must not contain duplicate values.", what))
  }
  invisible(TRUE)
}

.validate_feature_names <- function(x, what = "`x`") {
  feature_names <- colnames(x)
  if (is.null(feature_names)) {
    return(invisible(TRUE))
  }
  if (anyNA(feature_names) || any(feature_names == "")) {
    .abort_stablr("stablr_input_error", sprintf("%s feature names must be non-empty.", what))
  }
  if (anyDuplicated(feature_names)) {
    .abort_stablr("stablr_input_error", sprintf("%s feature names must not contain duplicates.", what))
  }
  invisible(TRUE)
}

.validate_scalar_numeric <- function(x, arg, min = -Inf, max = Inf,
                                     min_open = FALSE, max_open = FALSE,
                                     allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    .abort_stablr("stablr_input_error", sprintf("`%s` must be a finite numeric scalar.", arg))
  }
  if ((min_open && x <= min) || (!min_open && x < min)) {
    .abort_stablr("stablr_input_error", sprintf("`%s` must be %s %s.", arg, if (min_open) ">" else ">=", min))
  }
  if ((max_open && x >= max) || (!max_open && x > max)) {
    .abort_stablr("stablr_input_error", sprintf("`%s` must be %s %s.", arg, if (max_open) "<" else "<=", max))
  }
  x
}

.validate_scalar_integer_like <- function(x, arg, min = -Inf, max = Inf,
                                          allow_null = FALSE) {
  x <- .validate_scalar_numeric(
    x,
    arg,
    min = min,
    max = max,
    allow_null = allow_null
  )
  if (is.null(x)) {
    return(NULL)
  }
  if (x != floor(x)) {
    .abort_stablr("stablr_input_error", sprintf("`%s` must be an integer-like scalar.", arg))
  }
  as.integer(x)
}

.validate_proportion <- function(x, arg, allow_one = TRUE) {
  .validate_scalar_numeric(
    x,
    arg,
    min = 0,
    max = 1,
    min_open = TRUE,
    max_open = !allow_one
  )
}
