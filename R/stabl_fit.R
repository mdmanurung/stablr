#' Fit STABL (Stability-Penalized Feature Selection)
#'
#' Performs the STABL bootstrap stability-selection procedure over a grid of
#' regularisation parameters, optionally injecting artificial features for
#' data-adaptive FDP+ threshold diagnostics. The fitting rule chooses the FDP+
#' minimizer; it is not a universal false-discovery guarantee. This is the R counterpart of
#' `Stabl.fit()` in the Python STABL library.
#'
#' @param x A numeric matrix or `data.frame` (samples \eqn{\times} features)
#'   with row names used as sample IDs.
#' @param y A named numeric/factor vector whose names are sample IDs, or a
#'   matrix-like outcome (for example `survival::Surv`) with row names as
#'   sample IDs.
#' @param lambda_grid A pre-expanded `data.frame` (each row = one parameter
#'   combination) with at least a `lambda` column.  Pass `"auto"` to derive a
#'   data-driven grid via [auto_lambda_grid()] (uses `family` and `n_lambda`).
#' @param base_learner Character; `"lasso"` (alpha = 1), `"elastic_net"`
#'   (reads `alpha` from the `lambda_grid` `alpha` column), or
#'   `"adaptive_lasso"`, or `"sparse_group_lasso"`. Default: `"lasso"`.
#' @param family Character; `glmnet` response family (`"gaussian"`,
#'   `"binomial"`, `"multinomial"`, or `"cox"`). Default: `"gaussian"`.
#' @param n_bootstraps Positive integer; bootstrap iterations per lambda.
#'   Default: `1000L`.
#' @param artificial_type Character or `NULL`; `"random_permutation"`,
#'   `"knockoff"`, `"knockoff_equi"`, `"knockoff_mvr"`, or `NULL` (no
#'   artificial features — requires `hard_threshold`).  Default:
#'   `"random_permutation"`.
#' @param artificial_proportion Numeric in `(0, 1]`; fraction of original
#'   features to inject as artificial noise.  Default: `1.0`.
#' @param sample_fraction Positive numeric; fraction of samples drawn per
#'   bootstrap.  Default: `0.5`.
#' @param replace Logical; sample with replacement?  Default: `FALSE`.
#' @param hard_threshold Numeric in `(0, 1]` or `NULL`.  When supplied, FDP+
#'   control is bypassed and this value is the stability-score cut-off.
#'   Default: `NULL`.
#' @param fdr_threshold_range Numeric vector swept when computing FDP+.
#'   Default: `seq(0, 0.99, by = 0.01)`, matching Python STABL's
#'   `np.arange(0., 1., .01)`.
#' @param explore Logical; if `TRUE` and no features pass the threshold, fall
#'   back to the top `n_explore` features.  Default: `FALSE`.
#' @param n_explore Positive integer; fallback feature count.  Default: `5L`.
#' @param groups Named vector of group IDs (same names as `rownames(x)`) or
#'   `NULL`.  When supplied, [group_bootstrap_indices()] is used instead of
#'   [classic_bootstrap_indices()].  Default: `NULL`.
#' @param stratify_bootstrap Logical; if `TRUE`, bootstrap subsamples are
#'   stratified by the outcome class.  This is intended for classification
#'   tasks with small or imbalanced classes.  Default `FALSE` preserves the
#'   original STABL sampling behavior.
#' @param bootstrap_strata Optional categorical stratification design for
#'   bootstrap sampling.  Provide a named vector, matrix, `data.frame`, or
#'   list with one row/value per sample.  Multiple columns are combined as a
#'   joint interaction stratum, for example outcome class by study group.  When
#'   supplied, this overrides `stratify_bootstrap`.
#' @param n_lambda Integer; number of lambda values when `lambda_grid = "auto"`.
#'   Ignored otherwise.  Default: `30L`.
#' @param l1_ratio Numeric scalar, numeric vector, or `NULL`. Passed to
#'   [auto_lambda_grid()] when `lambda_grid = "auto"`. Use this for
#'   `base_learner = "elastic_net"` so the generated grid contains the
#'   elastic-net `alpha` column. Default `NULL` preserves the lasso-like auto
#'   grid used by existing workflows.
#' @param verbose Logical; emit progress messages.  Default: `FALSE`.
#' @param workers Positive integer; parallel workers for
#'   `furrr::future_map()`.  Actual parallelism requires the caller to invoke
#'   `future::plan(multisession, workers = N)` first.  Default: `1L`.
#' @param random_state Integer or `NULL`; top-level seed.  Default: `NULL`.
#' @param adaptive_gamma Positive numeric scalar. Used only when
#'   `base_learner = "adaptive_lasso"`.
#' @param adaptive_epsilon Positive numeric scalar. Used only when
#'   `base_learner = "adaptive_lasso"`.
#' @param feature_groups Optional feature-group definition for
#'   `base_learner = "sparse_group_lasso"`. Provide a vector/factor of length
#'   `ncol(x)` assigning each original feature to a group.
#' @param corr_group_threshold Optional numeric percentile in `(0, 100]` used
#'   to derive sparse-group feature groups from pairwise correlations of
#'   original features. Used only when
#'   `base_learner = "sparse_group_lasso"`.
#'
#' @return An S3 object of class `"stabl_fit"` containing:
#'   \describe{
#'     \item{`stabl_scores_`}{Numeric matrix (features \eqn{\times} lambdas)
#'       of per-bootstrap selection frequencies (stability scores).}
#'     \item{`stabl_scores_artificial_`}{Analogous matrix for artificial
#'       features, or `NULL` when `artificial_type = NULL`.}
#'     \item{`fitted_lambda_grid`}{The `data.frame` of lambda combinations
#'       actually used.}
#'     \item{`fdr_min_threshold_`}{FDP+-optimal stability threshold, or `NULL`.}
#'     \item{`FDRs_`}{Numeric vector of FDP+ estimates per threshold, or `NULL`.}
#'     \item{`min_fdr_`}{Minimum FDP+ achieved, or `NULL`.}
#'     \item{`fdrs_table_`}{Per-lambda FDP+ matrix, or `NULL`.}
#'     \item{`hard_threshold`}{As supplied.}
#'     \item{`artificial_type`}{As supplied.}
#'     \item{`artificial_provenance`}{List returned by
#'       [make_artificial_features()] recording actual artificial-feature
#'       generation modes and fallback counts/reasons, or `NULL` when
#'       `artificial_type = NULL`.}
#'     \item{`artificial_proportion`}{As supplied.}
#'     \item{`explore`}{As supplied.}
#'     \item{`n_explore`}{As supplied.}
#'     \item{`feature_names`}{Character vector of original feature names.}
#'     \item{`n_features_in_`}{Integer number of original features.}
#'   }
#'
#' @examples
#' set.seed(42)
#' n <- 50L
#' p <- 10L
#' x <- matrix(
#'   rnorm(n * p),
#'   nrow = n,
#'   dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p)))
#' )
#' y <- setNames(0.8 * x[, 1] - 0.5 * x[, 2] + rnorm(n, sd = 0.3), rownames(x))
#' lambda_grid <- data.frame(lambda = seq(0.02, 0.30, length.out = 4L))
#'
#' # Default lasso path
#' fit_lasso <- stabl_fit(
#'   x = x,
#'   y = y,
#'   lambda_grid = lambda_grid,
#'   n_bootstraps = 6L,
#'   artificial_type = "random_permutation",
#'   random_state = 1L
#' )
#' get_feature_names_out(fit_lasso)
#'
#' # Adaptive lasso path
#' fit_adaptive <- stabl_fit(
#'   x = x,
#'   y = y,
#'   lambda_grid = lambda_grid,
#'   base_learner = "adaptive_lasso",
#'   adaptive_gamma = 1.0,
#'   adaptive_epsilon = 1e-6,
#'   n_bootstraps = 6L,
#'   artificial_type = "random_permutation",
#'   random_state = 2L
#' )
#' get_importances(fit_adaptive)
#'
#' # Sparse-group lasso path (optional dependency)
#' if (requireNamespace("sparsegl", quietly = TRUE)) {
#'   feature_groups <- rep(seq_len(5L), each = 2L)
#'   fit_sgl <- stabl_fit(
#'     x = x,
#'     y = y,
#'     lambda_grid = lambda_grid,
#'     base_learner = "sparse_group_lasso",
#'     feature_groups = feature_groups,
#'     n_bootstraps = 5L,
#'     artificial_type = "random_permutation",
#'     random_state = 3L
#'   )
#'   get_support(fit_sgl)
#' }
#'
#' # Multinomial path with three classes
#' y_multi <- setNames(
#'   factor(sample(c("A", "B", "C"), n, replace = TRUE)),
#'   rownames(x)
#' )
#' fit_multi <- suppressWarnings(stabl_fit(
#'   x = x,
#'   y = y_multi,
#'   lambda_grid = lambda_grid,
#'   family = "multinomial",
#'   n_bootstraps = 5L,
#'   artificial_type = "random_permutation",
#'   random_state = 4L
#' ))
#' get_feature_names_out(fit_multi)
#' @export
stabl_fit <- function(
    x,
    y,
    lambda_grid,
    base_learner          = "lasso",
    family                = "gaussian",
    n_bootstraps          = 1000L,
    artificial_type       = "random_permutation",
    artificial_proportion = 1.0,
    sample_fraction       = 0.5,
    replace               = FALSE,
    hard_threshold        = NULL,
    fdr_threshold_range   = seq(0, 0.99, by = 0.01),
    explore               = FALSE,
    n_explore             = 5L,
    groups                = NULL,
    stratify_bootstrap    = FALSE,
    bootstrap_strata      = NULL,
    n_lambda              = 30L,
    l1_ratio              = NULL,
    verbose               = FALSE,
    workers               = 1L,
    random_state          = NULL,
    adaptive_gamma        = 1.0,
    adaptive_epsilon      = 1e-6,
    feature_groups        = NULL,
    corr_group_threshold  = NULL
) {
  # ---- Input coercion -------------------------------------------------------
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix or data.frame.", call. = FALSE)
  }

  validate_sample_alignment(x, y, groups)

  # Align y (and groups) to row order of x
  y <- .subset_outcome_by_ids(y, rownames(x))
  if (!is.null(groups)) groups <- groups[rownames(x)]
  if (!is.null(bootstrap_strata)) {
    bootstrap_strata <- .subset_bootstrap_strata_by_ids(
      bootstrap_strata,
      sample_ids = rownames(x)
    )
  } else if (isTRUE(stratify_bootstrap)) {
    bootstrap_strata <- y
  }
  bootstrap_strata_ids <- .bootstrap_strata_ids(
    strata = bootstrap_strata,
    n = nrow(x),
    arg = "bootstrap_strata"
  )

  # ---- Validate scalar params -----------------------------------------------
  .validate_stabl_params(n_bootstraps, sample_fraction, replace,
                         hard_threshold, artificial_type,
                         artificial_proportion, stratify_bootstrap,
                         n_explore, n_lambda, workers, random_state)
  n_bootstraps <- .validate_scalar_integer_like(n_bootstraps, "n_bootstraps", min = 1L)
  sample_fraction <- .validate_scalar_numeric(sample_fraction, "sample_fraction", min = 0, min_open = TRUE)
  if (!is.null(hard_threshold)) {
    hard_threshold <- .validate_proportion(hard_threshold, "hard_threshold")
  }
  if (!is.null(artificial_type)) {
    artificial_proportion <- .validate_proportion(artificial_proportion, "artificial_proportion")
  }
  n_explore <- .validate_scalar_integer_like(n_explore, "n_explore", min = 1L)
  n_lambda <- .validate_scalar_integer_like(n_lambda, "n_lambda", min = 1L)
  workers <- .validate_scalar_integer_like(workers, "workers", min = 1L)
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }

  n_samples   <- nrow(x)
  n_features  <- ncol(x)
  feat_names  <- colnames(x)
  if (is.null(feat_names)) feat_names <- paste0("x.", seq_len(n_features))

  # ---- Lambda grid ----------------------------------------------------------
  if (identical(lambda_grid, "auto")) {
    if (verbose) message("Computing auto lambda grid...")
    lambda_grid <- auto_lambda_grid(
      x,
      y,
      family = family,
      n_lambda = n_lambda,
      l1_ratio = l1_ratio
    )
  }
  if (!is.data.frame(lambda_grid)) {
    stop("`lambda_grid` must be a data.frame or \"auto\".", call. = FALSE)
  }
  if (!"lambda" %in% names(lambda_grid)) {
    stop("`lambda_grid` must contain a `lambda` column.", call. = FALSE)
  }
  n_lambdas <- nrow(lambda_grid)

  # ---- Derived RNG streams (audit M-5, V-7) ---------------------------------
  # When `random_state` is set we draw three independent integer seeds from a
  # single seeded RNG: one for artificial-feature generation, one for the
  # bootstrap-index draws, and one used as the base of per-iteration seeds
  # passed into the bootstrap loop.  This keeps the three RNG streams
  # independent (artificial features cannot consume RNG that boot indices
  # need) and makes the parallel and sequential loops bit-identical because
  # each iteration seeds its own RNG before calling the learner adapter.
  art_seed  <- NULL
  boot_seed <- NULL
  iter_seed_base <- NULL
  if (!is.null(random_state)) {
    set.seed(random_state)
    rng_seeds      <- sample.int(.Machine$integer.max, size = 3L)
    art_seed       <- rng_seeds[[1L]]
    boot_seed      <- rng_seeds[[2L]]
    iter_seed_base <- rng_seeds[[3L]]
  }

  # ---- Artificial features --------------------------------------------------
  n_injected        <- if (!is.null(artificial_type)) {
    min(n_features, max(1L, as.integer(round(n_features * artificial_proportion))))
  } else {
    0L
  }
  effective_artificial_proportion <- if (!is.null(artificial_type)) {
    n_injected / n_features
  } else {
    NULL
  }
  x_fit             <- x
  noise_col_indices <- NULL
  artificial_provenance <- NULL

  if (!is.null(artificial_type)) {
    if (verbose) message("Generating artificial features (type=",
                         artificial_type, ")...")
    art_result        <- make_artificial_features(x, n_injected,
                                                  artificial_type, art_seed)
    x_fit             <- art_result$x_augmented
    noise_col_indices <- art_result$noise_col_indices
    artificial_provenance <- art_result$artificial_provenance
  }

  # ---- Sparse-group feature groups -----------------------------------------
  sgl_feature_groups <- .resolve_sgl_feature_groups(
    base_learner         = base_learner,
    x_original           = x,
    x_augmented          = x_fit,
    feature_groups       = feature_groups,
    corr_group_threshold = corr_group_threshold,
    noise_col_indices    = noise_col_indices
  )

  # ---- Learner adapter (batch: one model path call per bootstrap) -----------
  # Each batch adapter accepts the full lambda_grid and returns a logical
  # matrix (n_total_features × n_lambdas), reducing n_bootstraps × n_lambdas
  # glmnet calls to just n_bootstraps calls.
  batch_adapter <- .select_batch_adapter(
    base_learner       = base_learner,
    family             = family,
    adaptive_gamma     = adaptive_gamma,
    adaptive_epsilon   = adaptive_epsilon,
    sgl_feature_groups = sgl_feature_groups
  )

  n_total_features <- ncol(x_fit)
  n_subsamples     <- as.integer(floor(sample_fraction * n_samples))
  modeled_classes <- .modeled_class_levels(y, family)

  # Fix 7: reject impossible configuration before entering bootstrap
  if (!replace && n_subsamples > n_samples) {
    stop(
      "`sample_fraction` (= ", sample_fraction, ") implies n_subsamples = ",
      n_subsamples, " > n_samples = ", n_samples,
      " but `replace = FALSE`. Reduce `sample_fraction` or set `replace = TRUE`.",
      call. = FALSE
    )
  }

  # ---- Bootstrap index lists ------------------------------------------------
  if (!is.null(boot_seed)) set.seed(boot_seed)

  boot_sampler <- if (!is.null(modeled_classes)) {
    .classification_bootstrap_sampler(
      y = y,
      groups = groups,
      n_subsamples = n_subsamples,
      replace = replace,
      strata_ids = bootstrap_strata_ids,
      family = family,
      class_levels = modeled_classes
    )
  } else if (is.null(groups)) {
    function(.) classic_bootstrap_indices(
      y = y,
      n_subsamples = n_subsamples,
      replace = replace,
      strata = bootstrap_strata
    )
  } else {
    function(.) group_bootstrap_indices(
      y = y,
      groups = groups,
      n_subsamples = n_subsamples,
      replace = replace,
      strata = bootstrap_strata
    )
  }
  boot_indices <- lapply(seq_len(n_bootstraps), boot_sampler)

  # ---- Per-iteration seeds for adapter calls (audit V-7) -------------------
  # Pre-generate one integer seed per bootstrap iteration so that the learner
  # adapter is reproducible and bit-identical across sequential and parallel
  # execution paths.
  iter_seeds <- if (!is.null(iter_seed_base)) {
    set.seed(iter_seed_base)
    sample.int(.Machine$integer.max, size = n_bootstraps)
  } else NULL

  art_rows <- if (!is.null(artificial_type)) {
    n_features + seq_len(n_injected)
  } else NULL

  if (verbose) message("Running STABL bootstrap loop (",
                       n_bootstraps, " bootstraps, ",
                       n_lambdas, " lambdas, 1 path call per bootstrap)...")

  # ---- Stability accumulation -----------------------------------------------
  boot_result  <- .run_bootstrap_accumulation(
    batch_adapter  = batch_adapter,
    boot_indices   = boot_indices,
    iter_seeds     = iter_seeds,
    x_fit          = x_fit,
    y              = y,
    lambda_grid    = lambda_grid,
    n_features     = n_features,
    feat_names     = feat_names,
    art_rows       = art_rows,
    artificial_type = artificial_type,
    workers        = workers
  )
  stabl_scores_    <- boot_result$stabl_scores_
  stabl_scores_art <- boot_result$stabl_scores_art

  # ---- FDP+ diagnostic thresholding ----------------------------------------
  fdp <- NULL
  if (!is.null(artificial_type)) {
    fdp <- compute_fdp_plus(
      stabl_scores            = stabl_scores_,
      stabl_scores_artificial = stabl_scores_art,
      artificial_proportion   = effective_artificial_proportion,
      fdr_threshold_range     = fdr_threshold_range
    )
  }

  # ---- Assemble S3 result ---------------------------------------------------
  structure(
    list(
      stabl_scores_         = stabl_scores_,
      stabl_scores_artificial_ = stabl_scores_art,
      fitted_lambda_grid    = lambda_grid,
      fdr_min_threshold_    = if (!is.null(fdp)) fdp$fdr_min_threshold else NULL,
      FDRs_                 = if (!is.null(fdp)) fdp$FDRs              else NULL,
      min_fdr_              = if (!is.null(fdp)) fdp$min_fdr            else NULL,
      fdrs_table_           = if (!is.null(fdp)) fdp$fdrs_table         else NULL,
      fdr_threshold_range   = fdr_threshold_range,
      hard_threshold        = hard_threshold,
      artificial_type       = artificial_type,
      artificial_provenance = artificial_provenance,
      artificial_proportion = artificial_proportion,
      effective_artificial_proportion = effective_artificial_proportion,
      stratify_bootstrap    = !is.null(bootstrap_strata_ids),
      bootstrap_strata_levels = if (!is.null(bootstrap_strata_ids)) {
        sort(unique(as.character(bootstrap_strata_ids)))
      } else NULL,
      explore               = explore,
      n_explore             = as.integer(n_explore),
      feature_names         = feat_names,
      n_features_in_        = n_features
    ),
    class = "stabl_fit"
  )
}

.modeled_class_levels <- function(y, family) {
  if (!family %in% c("binomial", "multinomial")) {
    return(NULL)
  }
  values <- as.character(y)
  if (is.factor(y)) {
    return(levels(droplevels(y)))
  }
  sort(unique(values))
}

.classification_bootstrap_is_feasible <- function(y, indices, class_levels,
                                                  minimum_per_class = 2L) {
  observed <- table(factor(as.character(y[indices]), levels = class_levels))
  all(observed >= minimum_per_class)
}

.abort_bootstrap_infeasible <- function(family, n_subsamples, class_counts,
                                        reason, minimum_per_class = 2L) {
  counts_text <- paste0(names(class_counts), "=", as.integer(class_counts),
                        collapse = ", ")
  .abort_stablr(
    "stablr_bootstrap_infeasible",
    paste0(
      "Learner-feasible ", family,
      " bootstrap sampling is impossible: ", reason,
      " Modeled class counts are ", counts_text,
      "; each bootstrap requires at least ", minimum_per_class,
      " observations from every class."
    ),
    family = family,
    n_subsamples = n_subsamples,
    class_counts = class_counts,
    minimum_per_class = minimum_per_class
  )
}

.classification_bootstrap_sampler <- function(y, groups, n_subsamples,
                                              replace, strata_ids, family,
                                              class_levels,
                                              minimum_per_class = 2L,
                                              max_redraws = 1000L) {
  class_counts <- table(factor(as.character(y), levels = class_levels))
  if (any(class_counts < minimum_per_class)) {
    .abort_bootstrap_infeasible(
      family,
      n_subsamples,
      class_counts,
      "the training population has too few observations in at least one class.",
      minimum_per_class
    )
  }
  if (is.null(groups) &&
      n_subsamples < length(class_levels) * minimum_per_class) {
    .abort_bootstrap_infeasible(
      family,
      n_subsamples,
      class_counts,
      paste0("`n_subsamples` is only ", n_subsamples, "."),
      minimum_per_class
    )
  }

  group_levels <- if (is.null(groups)) NULL else unique(groups)
  draw_once <- if (is.null(groups) && is.null(strata_ids)) {
    function() sample.int(
      n = length(y),
      size = n_subsamples,
      replace = replace
    )
  } else if (is.null(groups)) {
    function() .stratified_bootstrap_indices(
      strata_ids,
      n_subsamples,
      replace
    )
  } else if (is.null(strata_ids)) {
    function() .unstratified_group_bootstrap_indices(
      groups,
      group_levels,
      n_subsamples,
      replace
    )
  } else {
    function() .stratified_group_bootstrap_indices(
      strata_ids,
      groups,
      group_levels,
      n_subsamples,
      replace
    )
  }

  function(.) {
    indices <- draw_once()
    if (.classification_bootstrap_is_feasible(
      y,
      indices,
      class_levels,
      minimum_per_class
    )) {
      return(indices)
    }

    if (is.null(groups)) {
      return(.stratified_bootstrap_indices(
        factor(as.character(y), levels = class_levels),
        n_subsamples,
        replace,
        min_per_stratum = minimum_per_class
      ))
    }

    group_is_class_pure <- all(vapply(group_levels, function(group) {
      length(unique(as.character(y[groups == group]))) == 1L
    }, logical(1L)))
    if (group_is_class_pure) {
      repair_target <- max(
        n_subsamples,
        length(class_levels) * minimum_per_class
      )
      return(.stratified_group_bootstrap_indices(
        factor(as.character(y), levels = class_levels),
        groups,
        group_levels,
        repair_target,
        replace,
        min_per_stratum = minimum_per_class
      ))
    }

    for (attempt in seq_len(max_redraws)) {
      indices <- draw_once()
      if (.classification_bootstrap_is_feasible(
        y,
        indices,
        class_levels,
        minimum_per_class
      )) {
        return(indices)
      }
    }
    .abort_bootstrap_infeasible(
      family,
      n_subsamples,
      class_counts,
      paste0("no feasible grouped draw was found after ", max_redraws,
             " deterministic redraws."),
      minimum_per_class
    )
  }
}

# ---- Package-level constants -------------------------------------------------

# Minimum absolute coefficient to count a feature as "selected" in one
# bootstrap replicate. Matches the Python STABL default (1e-5) and is applied
# consistently across all base-learner adapters.
.BOOTSTRAP_COEF_THRESHOLD <- 1e-5

# Default sparse group lasso mixing parameter (asparse in sparsegl).
# Controls L1 vs group-lasso balance; 0.05 = mostly group-lasso.
.SGL_DEFAULT_ALPHA <- 0.05

# Breslow/Efron tie-breaking method for glmnet Cox models.
.COX_TIES <- "efron"

# ---- Internal param validator ------------------------------------------------
.validate_stabl_params <- function(n_bootstraps, sample_fraction, replace,
                                   hard_threshold, artificial_type,
                                   artificial_proportion,
                                   stratify_bootstrap,
                                   n_explore, n_lambda, workers,
                                   random_state) {
  .validate_scalar_integer_like(n_bootstraps, "n_bootstraps", min = 1L)
  .validate_scalar_numeric(sample_fraction, "sample_fraction", min = 0, min_open = TRUE)
  .validate_scalar_integer_like(n_explore, "n_explore", min = 1L)
  .validate_scalar_integer_like(n_lambda, "n_lambda", min = 1L)
  .validate_scalar_integer_like(workers, "workers", min = 1L)
  if (!is.null(random_state)) {
    .validate_scalar_integer_like(random_state, "random_state")
  }
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("`replace` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(stratify_bootstrap) || length(stratify_bootstrap) != 1L ||
      is.na(stratify_bootstrap)) {
    stop("`stratify_bootstrap` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!replace && sample_fraction > 1) {
    stop(
      "`sample_fraction` cannot exceed 1 when `replace = FALSE`.",
      call. = FALSE
    )
  }
  if (!is.null(hard_threshold)) {
    .validate_proportion(hard_threshold, "hard_threshold")
  }
  if (is.null(hard_threshold) && is.null(artificial_type)) {
    stop(
      "Either `hard_threshold` or `artificial_type` must be specified.",
      call. = FALSE
    )
  }
  if (!is.null(artificial_type)) {
    if (!is.character(artificial_type) || length(artificial_type) != 1L ||
        is.na(artificial_type)) {
      stop("`artificial_type` must be a single character string or NULL.", call. = FALSE)
    }
    .validate_proportion(artificial_proportion, "artificial_proportion")
  }
  invisible(NULL)
}

.resolve_sgl_feature_groups <- function(base_learner,
                                        x_original,
                                        x_augmented,
                                        feature_groups,
                                        corr_group_threshold,
                                        noise_col_indices) {
  if (base_learner != "sparse_group_lasso") {
    return(NULL)
  }

  has_explicit <- !is.null(feature_groups)
  has_corr <- !is.null(corr_group_threshold)

  if (has_explicit && has_corr) {
    stop(
      "For `base_learner = \"sparse_group_lasso\"`, provide either `feature_groups` or `corr_group_threshold`, not both.",
      call. = FALSE
    )
  }
  if (!has_explicit && !has_corr) {
    stop(
      "For `base_learner = \"sparse_group_lasso\"`, provide `feature_groups` or `corr_group_threshold`.",
      call. = FALSE
    )
  }

  if (has_explicit) {
    grp <- as.integer(as.factor(feature_groups))
    if (length(grp) != ncol(x_original)) {
      stop("`feature_groups` must have length ncol(x).", call. = FALSE)
    }
    if (anyNA(grp)) {
      stop("`feature_groups` cannot contain NA values.", call. = FALSE)
    }
  } else {
    grp <- .build_corr_groups(x_original, corr_group_threshold)
  }

  .append_noise_groups(grp, noise_col_indices, ncol(x_augmented))
}

.build_corr_groups <- function(x, percentile) {
  if (!is.numeric(percentile) || length(percentile) != 1L ||
      percentile <= 0 || percentile > 100) {
    stop("`corr_group_threshold` must be a numeric scalar in (0, 100].", call. = FALSE)
  }

  p <- ncol(x)
  if (p <= 1L) {
    return(rep.int(1L, p))
  }

  corr <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs"))
  corr[is.na(corr)] <- 0
  corr_vals <- corr[upper.tri(corr, diag = FALSE)]
  # Subtract 0.1 to match Python reference: threshold = np.percentile(corr_val, perc) - 0.1
  # (stabl/stabl.py line 1142).  Without the offset the R port uses a stricter
  # threshold than Python, producing fewer correlation groups than expected.
  cutoff <- as.numeric(stats::quantile(corr_vals, probs = percentile / 100,
                                       names = FALSE, na.rm = TRUE)) - 0.1

  parent <- seq_len(p)
  find_root <- function(i) {
    while (parent[[i]] != i) {
      parent[[i]] <<- parent[[parent[[i]]]]
      i <- parent[[i]]
    }
    i
  }

  # Extract only supra-threshold pairs from the upper triangle — O(p²) memory
  # for the full correlation matrix is unavoidable, but iterating only over the
  # k << p² edges that exceed the cutoff is much faster for sparse graphs.
  edges <- which(corr > cutoff & upper.tri(corr), arr.ind = TRUE)
  for (k in seq_len(nrow(edges))) {
    ri <- find_root(edges[[k, 1L]])
    rj <- find_root(edges[[k, 2L]])
    if (ri != rj) parent[[rj]] <- ri
  }

  roots <- vapply(seq_len(p), find_root, integer(1L))
  as.integer(as.factor(roots))
}

.append_noise_groups <- function(groups, noise_col_indices, total_p) {
  if (is.null(noise_col_indices)) {
    return(groups)
  }

  # Preallocate the full output vector to avoid quadratic c()-in-loop growth.
  n_orig   <- length(groups)
  n_noise  <- length(noise_col_indices)
  out      <- integer(n_orig + n_noise)
  out[seq_len(n_orig)] <- as.integer(groups)
  next_gid <- max(as.integer(groups))

  for (i in seq_len(n_noise)) {
    src <- as.integer(noise_col_indices[[i]])
    if (!is.na(src) && src >= 1L && src <= n_orig) {
      out[n_orig + i] <- groups[[src]]
    } else {
      next_gid        <- next_gid + 1L
      out[n_orig + i] <- next_gid
    }
  }

  if (length(out) != total_p) {
    stop("Sparse-group feature-group construction failed due to length mismatch.", call. = FALSE)
  }

  out
}

.select_batch_adapter <- function(base_learner, family,
                                  adaptive_gamma, adaptive_epsilon,
                                  sgl_feature_groups) {
  switch(
    base_learner,
    lasso = .make_glmnet_batch_adapter(
      family              = family,
      alpha_fixed         = 1.0,
      bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
    ),
    elastic_net = .make_glmnet_batch_adapter(
      family              = family,
      alpha_fixed         = NULL,
      bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
    ),
    adaptive_lasso = .make_adaptive_lasso_batch_adapter(
      family              = family,
      gamma               = adaptive_gamma,
      epsilon             = adaptive_epsilon,
      bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
    ),
    sparse_group_lasso = .make_sgl_batch_adapter(
      family              = family,
      feature_groups      = sgl_feature_groups,
      alpha_fixed         = NULL,
      bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
    ),
    stop(
      "`base_learner` must be one of: \"lasso\", \"elastic_net\", \"adaptive_lasso\", \"sparse_group_lasso\".",
      call. = FALSE
    )
  )
}

.run_bootstrap_accumulation <- function(
    batch_adapter, boot_indices, iter_seeds,
    x_fit, y, lambda_grid,
    n_features, feat_names, art_rows, artificial_type, workers
) {
  n_bootstraps <- length(boot_indices)
  n_lambdas    <- nrow(lambda_grid)

  scores_ <- matrix(0.0, nrow = n_features, ncol = n_lambdas,
                    dimnames = list(feat_names, NULL))
  art_    <- if (!is.null(artificial_type)) {
    matrix(0.0, nrow = length(art_rows), ncol = n_lambdas)
  } else NULL

  process_one <- function(i) {
    idx <- boot_indices[[i]]
    if (!is.null(iter_seeds)) {
      .with_local_seed(iter_seeds[[i]],
                       batch_adapter(x_fit[idx, , drop = FALSE], y[idx], lambda_grid))
    } else {
      batch_adapter(x_fit[idx, , drop = FALSE], y[idx], lambda_grid)
    }
  }

  use_furrr <- (workers > 1L) &&
    requireNamespace("furrr",  quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)

  real_rows <- seq_len(n_features)

  if (use_furrr) {
    # Parallel path: collect all results then stream-accumulate.
    # `seed = TRUE` guards against RNG leaks outside `.with_local_seed`.
    result_list <- furrr::future_map(
      seq_along(boot_indices), process_one,
      .options = furrr::furrr_options(seed = TRUE)
    )
    for (r in result_list) {
      scores_ <- scores_ + r[real_rows, , drop = FALSE]
      if (!is.null(art_)) art_ <- art_ + r[art_rows, , drop = FALSE]
    }
  } else {
    # Sequential path: stream-accumulate one bootstrap at a time (Fix 5).
    for (i in seq_along(boot_indices)) {
      r <- process_one(i)
      scores_ <- scores_ + r[real_rows, , drop = FALSE]
      if (!is.null(art_)) art_ <- art_ + r[art_rows, , drop = FALSE]
    }
  }

  scores_ <- scores_ / n_bootstraps
  if (!is.null(art_)) art_ <- art_ / n_bootstraps

  list(stabl_scores_ = scores_, stabl_scores_art = art_)
}

# Internal: scoped seed helper.  Saves & restores `.Random.seed` so that the
# wrapped expression has its own deterministic RNG state without polluting the
# caller's RNG.  Used throughout the package to make per-bootstrap and
# per-fold calls reproducible.  When `seed` is NULL, `expr` is evaluated
# without any seeding (the caller's RNG state is left untouched).
.with_local_seed <- function(seed, expr) {
  if (!is.null(seed)) {
    has_old <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (has_old) {
      old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        },
        add = TRUE
      )
    }
    set.seed(seed)
  }
  expr
}
