# 0.3.0 API overhaul: the beginner one-stop-shop entry point.
#
# Two-tier contract:
#   * `stabl()` — for first-time and application users.  Accepts a
#     data.frame or matrix, turns on the statistically safe defaults
#     (model-X knockoff artificial features, FDP+ calibration, strict
#     fallback), and returns an object that prints a plain-language
#     answer.  Roughly eight arguments cover the common cases.
#   * `stabl_fit()` — the expert interface with full control over all
#     27 arguments (lambda grids, bootstrap design, thresholds, adapters).
#     `stabl()` is a thin, opinionated wrapper around it; every `stabl()`
#     result is also a valid `stabl_fit` object, so all accessors and the
#     multiomic tooling keep working.

#' STABL Stability Selection (Beginner Interface)
#'
#' One-stop-shop stability selection for biomarker discovery.  Give it a
#' table of features and an outcome; it runs STABL with statistically safe
#' defaults and returns an object that tells you, in plain language, which
#' features were selected and how reliable the selection is.
#'
#' This is the recommended entry point for new users.  It wraps the expert
#' interface [stabl_fit()] with opinionated defaults:
#' \itemize{
#'   \item \code{artificial_type = "knockoff_equi"} — model-X knockoff
#'     artificial features, the pairing validated for FDP+ calibration
#'     (the package default `"random_permutation"` is retained on the
#'     expert interface for Python-parity workflows but is
#'     anti-conservative for multivariate learners).
#'   \item FDP+ calibration on (no `hard_threshold`), strict knockoff
#'     fallback (`knockoff_fallback = "stop"`).
#'   \item Automatic lambda grid when none is supplied.
#' }
#' Everything else — other learners, bootstrap design, thresholds, groups —
#' is one `...` away: unnamed arguments in `...` are forwarded to
#' [stabl_fit()], which documents the full expert interface.
#'
#' @param x A `data.frame` or numeric matrix of predictors (samples \eqn{\times}
#'   features).  A `data.frame` may contain the outcome as a column (see
#'   `y`).  Row names become sample IDs when present; otherwise IDs are
#'   generated as `s1, s2, ...`.
#' @param y The outcome.  One of: a vector (numeric, character, or factor)
#'   aligned to `x`'s rows (by names when both are named, otherwise by
#'   position); or, when `x` is a `data.frame`, a length-1 character naming
#'   the outcome column, which is then removed from the predictors.
#' @param family Character; `"gaussian"` (default), `"binomial"`,
#'   `"multinomial"`, or `"cox"`.  Use `"binomial"` for two-class outcomes,
#'   `"multinomial"` for more.
#' @param base_learner Character; `"lasso"` (default), `"elastic_net"`,
#'   `"adaptive_lasso"`, `"sparse_group_lasso"`, `"rf"`, or `"oblique_rf"`.
#' @param n_bootstraps Integer; bootstrap iterations per lambda.  Default
#'   `1000L` (the published STABL default); lower it only for quick
#'   exploration.
#' @param random_state Integer or `NULL`; seed for reproducibility.
#' @param ... Additional arguments forwarded to [stabl_fit()] (the expert
#'   interface), e.g. `artificial_proportion`, `sample_fraction`,
#'   `sample_groups`, `knockoff_fallback`.
#'
#' @return An object of class `c("stabl", "stabl_fit")`.  Print it for a
#'   plain-language answer; use [summary()], [plot()], [get_support()],
#'   [get_feature_names()], and [get_fdrs()] for details.  All
#'   `stabl_fit` accessors and workflows apply unchanged.
#'
#' @seealso [stabl_fit()] for the expert interface; [get_feature_names()]
#'   for the selected features; [compute_fdp_plus()] for the calibration
#'   machinery.
#'
#' @examples
#' set.seed(1L)
#' df <- data.frame(
#'   age  = rnorm(40),
#'   m1   = rnorm(40), m2 = rnorm(40), m3 = rnorm(40),
#'   m4   = rnorm(40), m5 = rnorm(40)
#' )
#' df$y <- 0.8 * df$m1 + 0.8 * df$m2 + rnorm(40, sd = 0.5)
#' \donttest{
#' fit <- stabl(df, y = "y", n_bootstraps = 50L, random_state = 1L)
#' fit                      # plain-language print
#' feature_names(fit)       # selected features
#' }
#' @export
stabl <- function(x, y, family = "gaussian", base_learner = "lasso",
                  n_bootstraps = 1000L, random_state = NULL, ...) {
  # ---- x: data.frame or matrix ------------------------------------------------
  if (is.data.frame(x)) {
    # Beginner convenience: y may name a column of x (a plain, unnamed
    # single string that matches a column).
    if (is.character(y) && length(y) == 1L && is.null(names(y)) &&
        y %in% colnames(x)) {
      outcome   <- x[[y]]
      x[[y]]    <- NULL
      y         <- outcome
    }
    if (!all(vapply(x, function(col) is.numeric(col), logical(1L)))) {
      stop(
        "`x` contains non-numeric columns. Remove identifier, factor, and ",
        "character columns (including the outcome column if you passed ",
        "`y` as a column name) before calling stabl().",
        call. = FALSE
      )
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a data.frame or numeric matrix.", call. = FALSE)
  }

  # ---- sample IDs -------------------------------------------------------------
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("s", seq_len(nrow(x)))
  }
  if (is.null(names(y)) && length(y) == nrow(x)) {
    names(y) <- rownames(x)
  }

  fit <- stabl_fit(
    x             = x,
    y             = y,
    lambda_grid   = "auto",
    base_learner  = base_learner,
    family        = family,
    n_bootstraps  = n_bootstraps,
    artificial_type = "knockoff_equi",
    random_state  = random_state,
    ...
  )
  class(fit) <- c("stabl", class(fit))
  fit
}

#' Print a Fitted STABL Object (Beginner Style)
#'
#' Prints a plain-language answer: what was fitted, which features were
#' selected, and how reliable the selection is estimated to be.  The
#' underlying terse `stabl_fit` print remains available by calling
#' `print.stabl_fit()` explicitly or via [stabl_fit()] objects.
#'
#' @param x A `"stabl"` object returned by [stabl()].
#' @param ... Ignored; present for S3 `print` method compatibility.
#'
#' @return The object, invisibly.
#' @seealso [stabl()], [summary.stabl()], [plot.stabl()]
#' @export
print.stabl <- function(x, ...) {
  .check_fitted_stabl(x)
  selected <- tryCatch(get_feature_names(x), error = function(e) NULL)
  min_fdr  <- x$min_fdr_
  cat("<stabl stability selection>\n")
  cat(sprintf(
    "  outcome family: %s | base learner: %s | bootstraps: %s\n",
    x$family %||% "gaussian",
    x$base_learner %||% "lasso",
    x$n_bootstraps %||% "?"
  ))
  if (is.null(selected)) {
    cat("  no features passed the stability threshold\n")
  } else {
    cat(sprintf(
      "  %d feature(s) selected at FDP+ estimate %s:\n",
      length(selected),
      if (is.null(min_fdr)) "n/a" else sprintf("%.3f", min_fdr)
    ))
    cat(" ", paste(utils::head(selected, 10L), collapse = ", "))
    if (length(selected) > 10L) cat(", ...")
    cat("\n")
  }
  cat("  use summary() for details, plot() for stability paths,\n")
  cat("  feature_names() / support() for the selection, get_fdrs() for the\n")
  cat("  FDP+ curve.\n")
  invisible(x)
}

#' Summarize a Fitted STABL Object (Beginner Style)
#'
#' Returns and prints a structured summary with the selection, the FDP+
#' estimate, and the chosen threshold, in plain language.
#'
#' @param object A `"stabl"` object returned by [stabl()].
#' @param ... Ignored; present for S3 `summary` method compatibility.
#'
#' @return A `"summary.stabl_fit"` object (see [summary.stabl_fit()]).
#' @seealso [stabl()], [print.stabl()]
#' @export
summary.stabl <- function(object, ...) {
  s <- summary.stabl_fit(object, ...)
  cat("<stabl stability selection summary>\n")
  cat(sprintf("  features considered : %s\n", s$n_features_in))
  cat(sprintf("  bootstraps          : %s\n", s$n_bootstraps))
  cat(sprintf("  stability threshold : %s\n",
              if (is.null(s$stabl_threshold)) "FDP+-chosen" else
                sprintf("%.3f", s$stabl_threshold)))
  cat(sprintf("  FDP+ estimate       : %s\n",
              if (is.null(s$min_fdr)) "n/a" else sprintf("%.3f", s$min_fdr)))
  cat(sprintf("  selected (%d)       : %s\n", s$n_selected,
              if (s$n_selected) paste(utils::head(s$selected_features, 20L),
                                      collapse = ", ") else "(none)"))
  invisible(s)
}

#' Plot Stability Paths for a Fitted STABL Object (Beginner Style)
#'
#' Delegates to [plot.stabl_fit()] (stability path plus FDR+ curve).
#'
#' @param x A `"stabl"` object returned by [stabl()].
#' @param ... Passed to [plot.stabl_fit()].
#'
#' @return The plotted object, invisibly.
#' @seealso [plot.stabl_fit()], [stabl()]
#' @export
plot.stabl <- function(x, ...) {
  plot.stabl_fit(x, ...)
}

#' Support and Feature Names of a STABL Selection
#'
#' Beginner-friendly generics for the selection: `support()` returns the
#' selected features' names or mask (delegating to [get_support()]) and
#' `feature_names()` returns the selected feature names at the chosen
#' threshold (delegating to [get_feature_names()]).  Methods are provided
#' for `stabl`/`stabl_fit` objects.
#'
#' @param object A fitted `stabl`/`stabl_fit` object.
#' @param ... Passed through (e.g. `new_hard_threshold`).
#'
#' @return See [get_support()] / [get_feature_names()].
#' @seealso [get_support()], [get_feature_names()], [stabl()]
#' @name support
#' @examples
#' \donttest{
#' set.seed(1L)
#' df <- data.frame(m1 = rnorm(40), m2 = rnorm(40), m3 = rnorm(40))
#' df$y <- 0.8 * df$m1 + rnorm(40, sd = 0.5)
#' fit <- stabl(df, y = "y", n_bootstraps = 50L, random_state = 1L)
#' feature_names(fit)
#' }
NULL

#' @rdname support
#' @export
support <- function(object, ...) {
  UseMethod("support")
}

#' @rdname support
#' @export
feature_names <- function(object, ...) {
  UseMethod("feature_names")
}

#' @rdname support
#' @export
support.stabl_fit <- function(object, ...) {
  get_support(object, ...)
}

#' @rdname support
#' @export
support.stabl <- function(object, ...) {
  get_support(object, ...)
}

#' @rdname support
#' @export
feature_names.stabl_fit <- function(object, ...) {
  get_feature_names(object, ...)
}

#' @rdname support
#' @export
feature_names.stabl <- function(object, ...) {
  get_feature_names(object, ...)
}
