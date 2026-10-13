# environment_stability.R — Environment-stability workflow (improvement plan D3)
#
# The user-facing appeal of invariant causal prediction (ICP; Peters, Bühlmann
# & Meinshausen) is stability ACROSS environments: a variable whose estimated
# relationship with the outcome is the same in every environment is a candidate
# invariant (causal) feature. This file delivers the *stability* half of that
# idea without any causal assumptions: fit STABL separately per environment and
# aggregate per-feature selection frequencies by the MINIMUM across
# environments. A feature is only as reliable as its weakest environment.
#
# No new dependencies: environments are plain labels (batch / cohort / site /
# time block), and all machinery comes from stabl_fit() and its accessors.

#' Environment-Stability Analysis Across Pre-Defined Environments
#'
#' Fits [stabl_fit()] separately within each environment (batch, cohort, site,
#' ...) and aggregates the per-feature stability scores across environments.
#' The headline quantity is the **minimum** over environments of each feature's
#' max-over-lambda stability score — the ICP-intersection analog: a feature is
#' only as reliable as its weakest environment.
#'
#' @details
#' **Motivation.** Invariant causal prediction (ICP) selects the variables
#' whose relationship with the outcome is *invariant* across environments
#' (Peters, Bühlmann & Meinshausen, 2016). The inferential target of ICP is
#' causal parenthood, and its classical implementations require an
#' environment-valued *variable* and strong hidden-confounding assumptions.
#' This function deliberately makes **no causal claims**: it reuses the
#' observation that stability across environments is the property users care
#' about, and measures it with the STABL machinery stablr already has
#' (per-feature selection frequencies). Environments are nothing more than
#' labels — batch, cohort, site, sequencing lane, time block — so the workflow
#' applies to ordinary multi-batch omic data where an environment *variable*
#' would not exist. No new dependencies are introduced.
#'
#' **Aggregation rule.** For each environment `e`, the per-feature
#' max-over-lambda stability score is read from the fit with
#' [get_importances()] (the same score [get_support()] thresholds). For each
#' feature, `min_scores` is the **minimum** of these values across
#' environments and `mean_scores` their mean. The minimum is the
#' ICP-intersection analog: a feature that is perfectly stable in one batch but
#' noise-level in another gets a low `min_score`, flagging it as
#' non-transportable. Note that max-over-lambda scores are *lenient* — a
#' lambda grid that extends into the saturated small-lambda regime (where the
#' lasso selects almost every feature) inflates all scores toward 1, so prefer
#' an explicit, top-heavy `lambda_grid` (or a per-environment `"auto"` grid
#' interpreted with this leniency in mind).
#'
#' **Reproducibility.** When `random_state` is supplied, the fit for the
#' i-th environment level (in the order of `levels(environments)`; sorted
#' order for character input) receives the deterministically derived seed
#' `random_state + i` (i = 1, 2, ...). With `random_state = NULL` the
#' per-environment fits are unseeded, matching [stabl_fit()] behavior.
#'
#' **Consistent feature handling.** All environments are fit on the same
#' feature columns of `x` (rows are subset, columns are not), so per-feature
#' scores are directly comparable across environments; this is verified
#' explicitly after fitting.
#'
#' @param x A numeric matrix or `data.frame` (samples × features) with row
#'   names used as sample IDs, as in [stabl_fit()].
#' @param y A named numeric/factor vector whose names are sample IDs, or a
#'   matrix-like outcome (for example `survival::Surv`) with row names as
#'   sample IDs. Only the rows belonging to an environment are used for that
#'   environment's fit.
#' @param environments Factor or character vector of length `nrow(x)` with one
#'   label per sample. `NA` labels and labels that are empty strings are
#'   rejected (fail closed), as are factor levels with zero observations and
#'   designs with fewer than two non-empty environments.
#' @param lambda_grid Passed to [stabl_fit()]. Defaults to `"auto"`, which
#'   derives a data-driven grid *separately within each environment* (see
#'   [auto_lambda_grid()]); supply an explicit grid for full control.
#' @param random_state Integer or `NULL`. Top-level seed; the i-th environment
#'   is fit with the derived seed `random_state + i` (see Details).
#' @param ... Additional arguments forwarded to [stabl_fit()] (for example
#'   `n_bootstraps`, `base_learner`, `family`, `artificial_type`,
#'   `hard_threshold`, `groups`, `bootstrap_strata`). Sample-keyed arguments
#'   such as `groups` are aligned by sample ID inside [stabl_fit()], so full
#'   data-length vectors may be supplied.
#'
#' @return An S3 object of class `"stabl_env_stability"` containing:
#'   \describe{
#'     \item{`env_fits`}{Named list of per-environment `"stabl_fit"` objects.}
#'     \item{`env_scores`}{Numeric matrix (features × environments) of
#'       max-over-lambda stability scores per environment.}
#'     \item{`min_scores`}{Named numeric vector: minimum of `env_scores`
#'       across environments (the transportability summary).}
#'     \item{`mean_scores`}{Named numeric vector: mean of `env_scores` across
#'       environments.}
#'     \item{`environments`}{Character vector of the environment levels used,
#'       in fit order.}
#'     \item{`pairwise_adjusted_similarity`}{Environment × environment matrix
#'       of chance-corrected agreement ([adjusted_similarity()]) between the
#'       features selected in each pair of environments (diagonal = 1 by
#'       definition), or `NULL` when no stability threshold can be resolved
#'       for the fits.}
#'     \item{`call`}{The call that produced this object.}
#'   }
#'
#' @seealso [stabl_fit()], [get_importances()], [get_min_scores()],
#'   [adjusted_similarity()], [auto_lambda_grid()]
#'
#' @references
#' Peters, J., Bühlmann, P., & Meinshausen, N. (2016). Causal inference by
#' using invariant prediction: identification and confidence intervals.
#' *Journal of the Royal Statistical Society: Series B*, 78(5), 947–1012.
#' (Conceptual motivation only; no causal inference is performed here and no
#' ICP software is required.)
#'
#' @examples
#' set.seed(1L)
#' n <- 60L
#' env <- factor(rep(c("A", "B"), each = n / 2L))
#' x <- matrix(rnorm(n * 6L), n, 6L,
#'             dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(6L))))
#' y <- setNames(0.8 * x[, 1] + 0.6 * x[, 2] + rnorm(n, sd = 0.3), rownames(x))
#' res <- stabl_environment_stability(
#'   x, y, env,
#'   lambda_grid   = data.frame(lambda = c(0.5, 0.3, 0.15)),
#'   n_bootstraps  = 10L,
#'   hard_threshold = 0.5,
#'   random_state  = 1L
#' )
#' res$min_scores                      # transportability summary
#' get_min_scores(res)                 # same, via accessor
#' @export
stabl_environment_stability <- function(x,
                                        y,
                                        environments,
                                        lambda_grid = "auto",
                                        random_state = NULL,
                                        ...) {
  call <- match.call()

  # ---- Validate x shape and environment labels (fail closed) ----------------
  if (!(is.data.frame(x) || is.matrix(x))) {
    .abort_stablr(
      "stablr_invalid_environments",
      "`x` must be a data.frame or matrix."
    )
  }
  n_samples <- nrow(x)
  env_factor <- .validate_environment_labels(environments, n_samples)
  env_levels <- levels(env_factor)

  # ---- Derived per-environment seeds (documented: random_state + level i) ---
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(
      random_state,
      "random_state",
      max = .Machine$integer.max - length(env_levels)
    )
  }
  env_seeds <- if (!is.null(random_state)) {
    random_state + seq_along(env_levels)
  } else {
    NULL
  }

  # ---- Fit STABL separately per environment ---------------------------------
  env_fits <- vector("list", length(env_levels))
  names(env_fits) <- env_levels
  for (i in seq_along(env_levels)) {
    lv <- env_levels[[i]]
    idx <- which(env_factor == lv)
    x_env <- x[idx, , drop = FALSE] # rows subset; feature columns shared
    y_env <- .subset_outcome_by_ids(y, rownames(x_env))
    env_fits[[lv]] <- stabl_fit(
      x = x_env,
      y = y_env,
      lambda_grid = lambda_grid,
      ...,
      random_state = env_seeds[[i]]
    )
  }

  # ---- Consistent feature handling across environments ----------------------
  feature_names_per_env <- lapply(env_fits, `[[`, "feature_names")
  if (length(unique(feature_names_per_env)) != 1L) {
    .abort_stablr(
      "stablr_invalid_environments",
      paste0(
        "Per-environment fits disagree on feature names; environment ",
        "stability requires all environments to be fit on the same ",
        "feature columns."
      )
    )
  }

  # ---- Aggregate: min / mean of max-over-lambda scores across environments --
  env_scores <- do.call(cbind, lapply(env_fits, get_importances))
  colnames(env_scores) <- env_levels
  min_scores  <- apply(env_scores, 1L, min)
  mean_scores <- rowMeans(env_scores)

  structure(
    list(
      env_fits                     = env_fits,
      env_scores                   = env_scores,
      min_scores                   = min_scores,
      mean_scores                  = mean_scores,
      environments                 = env_levels,
      pairwise_adjusted_similarity = .pairwise_env_similarity(env_fits),
      call                         = call
    ),
    class = "stabl_env_stability"
  )
}

# ---- Internal helpers --------------------------------------------------------

# Validate the `environments` argument: factor or character, one label per
# sample, no NAs, no empty-string labels, no factor levels with zero
# observations, and at least two non-empty environments. Returns the labels as
# a factor in fit order.
.validate_environment_labels <- function(environments, n_samples) {
  if (is.null(environments) ||
      length(environments) == 0L ||
      (!is.factor(environments) && !is.character(environments))) {
    .abort_stablr(
      "stablr_invalid_environments",
      paste0(
        "`environments` must be a factor or character vector with one ",
        "label per sample."
      )
    )
  }
  if (length(environments) != n_samples) {
    .abort_stablr(
      "stablr_invalid_environments",
      paste0(
        "`environments` must have one label per sample: expected ",
        n_samples, ", got ", length(environments), "."
      )
    )
  }
  if (anyNA(environments)) {
    .abort_stablr(
      "stablr_invalid_environments",
      paste0(
        "`environments` contains NA labels; environment stability fails ",
        "closed on missing environment assignments."
      )
    )
  }
  env_chr <- as.character(environments)
  if (any(!nzchar(env_chr))) {
    .abort_stablr(
      "stablr_invalid_environments",
      "`environments` contains empty-string labels."
    )
  }
  # Keep the caller's factor levels as-is (including empty ones, which are
  # rejected below); for character input, factor() derives sorted levels.
  env_factor <- if (is.factor(environments)) environments else factor(environments)
  counts <- table(env_factor)
  empty_levels <- names(counts)[counts == 0L]
  if (length(empty_levels) > 0L) {
    .abort_stablr(
      "stablr_invalid_environments",
      paste0(
        "`environments` has levels with no observations: ",
        paste(empty_levels, collapse = ", "), "."
      )
    )
  }
  if (length(env_levels <- levels(env_factor)) < 2L) {
    .abort_stablr(
      "stablr_invalid_environments",
      paste0(
        "At least two non-empty environments are required; got ",
        length(env_levels), "."
      )
    )
  }
  env_factor
}

# Pairwise chance-corrected agreement (adjusted_similarity) between the
# features selected in each pair of environments. Returns NULL when no
# threshold can be resolved for the fits (no artificial features and no hard
# threshold anywhere), rather than failing the whole workflow.
.pairwise_env_similarity <- function(env_fits) {
  supports <- tryCatch(
    lapply(env_fits, get_feature_names),
    error = function(e) NULL
  )
  if (is.null(supports)) {
    return(NULL)
  }
  k <- length(supports)
  sim <- matrix(
    NA_real_,
    nrow = k,
    ncol = k,
    dimnames = list(names(env_fits), names(env_fits))
  )
  n_features <- length(env_fits[[1L]]$feature_names)
  for (i in seq_len(k)) {
    sim[i, i] <- 1 # trivial self-agreement
    for (j in seq_len(k)[-i]) {
      sim[i, j] <- adjusted_similarity(
        supports[[i]],
        supports[[j]],
        n_universe = n_features
      )
    }
  }
  sim
}

#' Get the Cross-Environment Minimum Stability Scores
#'
#' Accessor for the `min_scores` element of a `"stabl_env_stability"` object:
#' for each feature, the minimum across environments of its max-over-lambda
#' stability score. This is the transportability summary computed by
#' [stabl_environment_stability()] — a feature is only as reliable as its
#' weakest environment.
#'
#' @param object A `"stabl_env_stability"` object returned by
#'   [stabl_environment_stability()].
#'
#' @return Named numeric vector with one entry per feature, values in
#'   \eqn{[0, 1]}. Higher values indicate features that are stable in *every*
#'   environment.
#'
#' @seealso [stabl_environment_stability()]
#'
#' @examples
#' set.seed(1L)
#' n <- 60L
#' env <- factor(rep(c("A", "B"), each = n / 2L))
#' x <- matrix(rnorm(n * 6L), n, 6L,
#'             dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(6L))))
#' y <- setNames(0.8 * x[, 1] + 0.6 * x[, 2] + rnorm(n, sd = 0.3), rownames(x))
#' res <- stabl_environment_stability(
#'   x, y, env,
#'   lambda_grid   = data.frame(lambda = c(0.5, 0.3, 0.15)),
#'   n_bootstraps  = 10L,
#'   hard_threshold = 0.5,
#'   random_state  = 1L
#' )
#' get_min_scores(res)
#' @export
get_min_scores <- function(object) {
  UseMethod("get_min_scores")
}

#' @export
get_min_scores.stabl_env_stability <- function(object) {
  object$min_scores
}

#' Print a Concise Summary of an Environment-Stability Analysis
#'
#' @describeIn stabl_environment_stability Print a concise summary of a
#'   `"stabl_env_stability"` object; invisibly returns `x`.
#' @param x A `"stabl_env_stability"` object returned by
#'   [stabl_environment_stability()].
#' @param ... Ignored; present for S3 `print` method compatibility.
#' @export
print.stabl_env_stability <- function(x, ...) {
  cat("<stabl_env_stability>\n")
  cat("  Environments:    ", length(x$environments), " (",
      paste(x$environments, collapse = ", "), ")\n", sep = "")
  cat("  Features in:     ", nrow(x$env_scores), "\n")
  cat("  Min score >= 0.5:", sum(x$min_scores >= 0.5), "feature(s)\n")
  top_n <- min(5L, length(x$min_scores))
  top <- names(sort(x$min_scores, decreasing = TRUE))[seq_len(top_n)]
  cat("  Top by min score:", paste(top, collapse = ", "), "\n")
  invisible(x)
}
