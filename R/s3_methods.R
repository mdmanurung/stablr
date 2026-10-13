# S3 convenience methods for fitted STABL objects (plan items C4).
#
# `plot.stabl_fit()` delegates to plot_stabl_path(); `summary.stabl_fit()`
# collects the headline results of a fit into a compact list, printed by
# `print.summary.stabl_fit()` in the same terse style as print.stabl_fit().

#' Plot the Stability Path of a Fitted STABL Object
#'
#' S3 `plot` method for `"stabl_fit"` objects: delegates to
#' [plot_stabl_path()], drawing the stability path of every real feature
#' across the fitted lambda grid, with selected features highlighted and the
#' effective stability threshold marked.
#'
#' @param x A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param ... Forwarded to [plot_stabl_path()]: `new_hard_threshold` (numeric
#'   in `(0, 1]` or `NULL`) and `title` (character scalar).
#'
#' @return The `ggplot` object, invisibly.
#'
#' @seealso [plot_stabl_path()], [plot_fdr_graph()]
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
#' p <- plot(fit)   # same as plot_stabl_path(fit)
#' @export
plot.stabl_fit <- function(x, ...) {
  invisible(plot_stabl_path(x, ...))
}

#' Summarize a Fitted STABL Object
#'
#' S3 `summary` method for `"stabl_fit"` objects: collects the headline
#' results of a fit into a compact list without touching the (potentially
#' large) stability-score matrices.
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param ... Ignored; present for S3 `summary` method compatibility.
#'
#' @return A list of class `"summary.stabl_fit"` with elements:
#'   \describe{
#'     \item{`selected_features`}{Character vector of feature names whose
#'       maximum stability score exceeds the effective threshold (empty when
#'       no feature passes; `NULL` when no threshold is available).}
#'     \item{`n_selected`}{Number of selected features.}
#'     \item{`stabl_threshold`}{Effective stability threshold (hard threshold
#'       or FDP+-optimal threshold), or `NULL` when none is available.}
#'     \item{`min_fdr`}{Minimum FDP+ diagnostic estimate achieved, or `NULL`
#'       when the fit used no artificial features.  Diagnostic estimate, not
#'       an FDR guarantee.}
#'     \item{`lambda_grid_used`}{Integer vector of length 2 giving the
#'       dimensions (rows x columns) of the `data.frame` of lambda
#'       combinations actually fitted; see [get_lambda_grid()].}
#'     \item{`n_features_in`}{Number of original features.}
#'     \item{`n_bootstraps`}{Number of bootstrap iterations used.}
#'   }
#'
#' @seealso [print.summary.stabl_fit()], [get_support()], [get_fdrs()],
#'   [get_lambda_grid()]
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
#' s <- summary(fit)
#' s$n_selected
#' s$selected_features
#' @export
summary.stabl_fit <- function(object, ...) {
  .check_fitted_stabl(object)
  threshold <- .resolve_stabl_threshold(object, on_missing = "null")
  selected  <- tryCatch(
    get_feature_names(object),
    error = function(e) NULL
  )
  structure(
    list(
      selected_features = selected,
      n_selected        = length(selected),
      stabl_threshold   = threshold,
      min_fdr           = object$min_fdr_,
      lambda_grid_used  = dim(object$fitted_lambda_grid),
      n_features_in     = object$n_features_in_,
      n_bootstraps      = object$n_bootstraps %||% NA_integer_
    ),
    class = "summary.stabl_fit"
  )
}

#' Print a Summary of a Fitted STABL Object
#'
#' Prints the compact fit summary produced by [summary.stabl_fit()] in the
#' same terse style as [print.stabl_fit()].
#'
#' @rdname summary.stabl_fit
#' @param x A `"summary.stabl_fit"` object returned by [summary.stabl_fit()].
#' @param ... Ignored; present for S3 `print` method compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [summary.stabl_fit()]
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
#' print(summary(fit))
#' @export
print.summary.stabl_fit <- function(x, ...) {
  cat("<stabl_fit summary>\n")
  cat("  Features in:      ", x$n_features_in, "\n")
  cat("  Features selected:", x$n_selected, "\n")
  if (!is.null(x$stabl_threshold)) {
    cat("  Threshold:        ", round(x$stabl_threshold, 4L), "\n")
  } else {
    cat("  Threshold:        [none available]\n")
  }
  if (!is.null(x$min_fdr)) {
    cat("  Min FDP+:         ", round(x$min_fdr, 4L), "\n")
  }
  cat("  Lambda grid:      ", x$lambda_grid_used[1L], " x ",
      x$lambda_grid_used[2L], "\n", sep = "")
  cat("  Bootstraps:       ", x$n_bootstraps, "\n")
  if (x$n_selected > 0L) {
    cat("  Selected:         ", paste(x$selected_features, collapse = ", "),
        "\n")
  }
  invisible(x)
}
