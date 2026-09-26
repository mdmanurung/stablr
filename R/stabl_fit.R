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
#'   threshold control is **suppressed** and this value is the stability-score
#'   cut-off.  Supplying it together with `artificial_type` is allowed -- the
#'   artificial features are still generated and their stability scores are
#'   still returned -- but the FDP+ curve is *not* computed, so `FDRs_`,
#'   `min_fdr_`, `fdr_min_threshold_` and `fdrs_table_` are all `NULL` and
#'   [plot_fdr_graph()] raises an error.  The suppression is recorded in
#'   `fdp_threshold_suppressed`; see Details.  Default: `NULL`.
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
#'   `ncol(x)` assigning each original feature to a group.  Artificial decoys
#'   are *not* placed in their source feature's group for permutation decoys
#'   (see the "Sparse-group artificial features" section).
#' @param corr_group_threshold Optional numeric percentile in `(0, 100]` used
#'   to derive sparse-group feature groups from pairwise correlations of
#'   original features. Used only when
#'   `base_learner = "sparse_group_lasso"`.
#'
#' @section Sparse-group artificial features:
#' FDP+ compares the stability of the real block against the stability of the
#' artificial block, which is only a like-for-like comparison if both blocks
#' face the same group-lasso penalty.  When `artificial_type` generates
#' *permutation* decoys, each decoy therefore receives its own fresh group
#' rather than joining the group of the original feature it was permuted from;
#' otherwise the decoy's selection is coupled to the real block's group
#' structure and its stability score partly measures the real block instead of
#' the decoy's own noise.  For `knockoff*` decoys the decoy is exchangeable
#' with its source by construction, so co-grouping is the honest
#' representation and inheritance stays on.  The resolved choice is recorded in
#' `sgl_group_provenance$noise_group_inheritance`, and the group vector actually
#' used is in `sgl_feature_groups_`.
#'
#' @section RNG and reproducibility:
#' When `random_state` is an integer, `stabl_fit()` does not mutate the
#' caller's `.Random.seed`: it derives three independent sub-seeds (one for
#' artificial-feature generation, one for the bootstrap index draws, one as the
#' base of the per-iteration seeds) and seeds each bootstrap iteration
#' individually, every one of them inside a `.with_local_seed()` scope.  Results
#' are then identical for `workers = 1L` and for a `furrr` parallel run.  That
#' guarantee is conditional on `random_state`; `workers > 1L` with
#' `random_state = NULL` warns, because there is no per-iteration seed to make
#' the parallel run reproducible.
#'
#' With `random_state = NULL` the bootstrap draws come from ambient RNG state
#' and therefore do advance the caller's stream, exactly as any unseeded
#' stochastic routine does.  Seed the fit if you need the caller protected.
#'
#' Reproducibility is also conditional on `RNGkind()` and on
#' `sample.kind`.  `stabl_fit()` does not set either, so changing them between
#' runs changes the sampled indices.  `RNGkind()` in force at fit time is
#' recorded in the returned `rng_kind_` field.
#'
#' @section Threshold resolution and `hard_threshold`:
#' Exactly one of `hard_threshold` and `artificial_type` drives selection.
#' When both are supplied, `hard_threshold` wins (as the accessor
#' documentation has always specified) and the FDP+ threshold search is skipped
#' rather than computed and discarded.  Nothing is lost except the FDP+
#' curve: `stabl_scores_artificial_` and `artificial_provenance` are still
#' returned.  The outcome is recorded in `fdp_threshold_suppressed`, so a
#' consumer can tell "no artificial features" from "artificial features
#' present, FDP+ not used".
#'
#' @section Bootstrap subsample size:
#' `sample_fraction` sets a *lower bound* on the rows in a bootstrap replicate,
#' not an exact count, and the realised size is recorded in the result so the
#' difference is visible.  Two mechanisms legitimately raise it.  Grouped
#' sampling never trims a partially drawn group, so the realised count is the
#' smallest group-aligned value at or above the target.  For `family =
#' "binomial"`/`"multinomial"` with class-pure groups, guaranteeing the minimum
#' number of observations per modeled class requires raising the row target to
#' `n_classes * minimum_per_class` and redrawing stratified; the raised target
#' is reported in `bootstrap_stratified_repair_`.  A replicate that ends up
#' *shorter* than the target is an error, because a fit on fewer rows than
#' requested must not be reported against the requested `sample_fraction`.
#'
#' @details
#' **Cost of the bootstrap loop.** Each bootstrap fits the full lambda path in
#' a single solver call for `base_learner = "lasso"`, `"adaptive_lasso"`, and
#' for `"elastic_net"`/`"sparse_group_lasso"` when a single `alpha` (resp.
#' `asparse`) value is in play.  A grid with \eqn{k} distinct `alpha` values
#' costs \eqn{k} path calls per bootstrap, and `"sparse_group_lasso"` pays
#' `n_lambda` `coef()` extractions per path call because the sparse
#' coefficient helper has no on-grid fast path.
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
#'     \item{`fdp_threshold_suppressed`}{Logical; `TRUE` when
#'       `hard_threshold` was supplied and therefore the FDP+ threshold search
#'       was skipped, `FALSE` when FDP+ drove selection, and `NA` when no
#'       artificial features were generated.  See the "Threshold resolution"
#'       section above.}
#'     \item{`hard_threshold`}{As supplied.}
#'     \item{`artificial_type`}{As supplied.}
#'     \item{`artificial_provenance`}{List returned by
#'       [make_artificial_features()] recording `actual_type`, all selected
#'       generator modes, and the complete ordered fallback history. This is
#'       the scientific source of truth for the generator used; or `NULL` when
#'       `artificial_type = NULL`.}
#'     \item{`artificial_proportion`}{As supplied.}
#'     \item{`explore`}{As supplied.}
#'     \item{`n_explore`}{As supplied.}
#'     \item{`explore_features_`}{Character vector of the `n_explore` highest
#'       importance features resolved at fit time, or `NULL` when
#'       `explore = FALSE`.  This is the fixed set the `explore` fallback
#'       selects for *any* threshold, so every accessor reports the same
#'       fallback rather than recomputing it per call.}
#'     \item{`explore_applied_`}{Logical; `TRUE` when the resolved threshold
#'       admitted no feature and the `explore` fallback therefore determined
#'       the returned support.}
#'     \item{`stratify_bootstrap`}{Logical; `TRUE` when stratified bootstrap
#'       sampling was requested, i.e. when `stratify_bootstrap = TRUE` **or**
#'       `bootstrap_strata` was supplied.  See
#'       `stratify_bootstrap_requested` and `bootstrap_strata_source` for the
#'       literal argument and the origin of the design.}
#'     \item{`stratify_bootstrap_requested`}{The literal `stratify_bootstrap`
#'       argument as supplied.}
#'     \item{`bootstrap_strata_source`}{`"bootstrap_strata"` when an explicit
#'       design was supplied, `"stratify_bootstrap"` when the design was
#'       derived from `y`, or `NULL` when no stratification was requested.}
#'     \item{`bootstrap_strata_applied`}{Logical; `TRUE` when a resolved
#'       stratification design actually reached the bootstrap sampler.}
#'     \item{`bootstrap_subsample_target_`}{Integer; the requested
#'       `floor(sample_fraction * nrow(x))` subsample size.  This is a *lower
#'       bound* on the realised draw, not an exact count: whole groups are never
#'       trimmed by grouped sampling, and the classification feasibility repair
#'       below can raise the target further.}
#'     \item{`bootstrap_subsample_sizes_`}{Integer vector of the realised
#'       subsample size of every bootstrap replicate, so any overshoot against
#'       `bootstrap_subsample_target_` is readable from the result.  A replicate
#'       shorter than the target is an error, not a silently smaller fit.}
#'     \item{`bootstrap_stratified_repair_`}{`NULL` unless the grouped
#'       classification sampler had to raise its row target to guarantee
#'       `minimum_per_class` observations of every modeled class; otherwise a
#'       list with `target_raised`, `requested_target`, `repair_target`,
#'       `n_classes`, `minimum_per_class` and `n_repairs`.}
#'     \item{`sgl_feature_groups_`}{Integer sparse-group assignment actually
#'       used, of length `ncol(x)` plus `n_injected`, or `NULL` when
#'       `base_learner != "sparse_group_lasso"`.}
#'     \item{`sgl_group_provenance`}{For `base_learner = "sparse_group_lasso"`,
#'       a list recording `group_source` (`"feature_groups"` or
#'       `"corr_group_threshold"`), `noise_group_inheritance` (`"on"`, `"off"`,
#'       or `NA` when there are no artificial features), `n_groups`, and
#'       `n_noise_features`; otherwise `NULL`.}
#'     \item{`rng_kind_`}{The `RNGkind()` in force when the fit ran.  Results
#'       are reproducible only for the same `RNGkind()` and `sample.kind`.}
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
  # The validation-side omics are held to a finiteness contract by
  # `.validate_multiomic_validation_inputs()`; the training side needs the same
  # rule here.  It used to be enforced only as a side effect of
  # `make_artificial_features()`, so with `artificial_type = NULL` a `NaN`
  # surfaced as an opaque `glmnet` "x has missing values" message and an `Inf`
  # reached the solver unchallenged.
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite, non-missing values.", call. = FALSE)
  }
  # Guard the degenerate shapes here: with ncol(x) == 0 the artificial block
  # would be empty, `n_injected` would be 0, and
  # `effective_artificial_proportion` would be 0/0 = NaN, surfacing much later
  # as a confusing `.validate_proportion` failure inside compute_fdp_plus().
  if (nrow(x) == 0L) {
    stop("`x` must contain at least one sample row.", call. = FALSE)
  }
  if (ncol(x) == 0L) {
    stop("`x` must contain at least one feature column.", call. = FALSE)
  }

  validate_sample_alignment(x, y, groups)

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
  family <- .validate_stabl_family(family)
  l1_ratio <- .validate_l1_ratio(l1_ratio)
  fdr_threshold_range <- .validate_fdr_threshold_range(fdr_threshold_range)

  stratify_requested <- isTRUE(stratify_bootstrap)
  stratify_design <- if (!is.null(bootstrap_strata)) {
    "bootstrap_strata"
  } else if (stratify_requested) {
    "stratify_bootstrap"
  } else {
    NULL
  }

  # Align y (and groups) to row order of x
  y <- .subset_outcome_by_ids(y, rownames(x))
  if (!is.null(groups)) {
    groups <- groups[rownames(x)]
    .validate_group_ids(groups)
  }
  if (!is.null(bootstrap_strata)) {
    bootstrap_strata <- .subset_bootstrap_strata_by_ids(
      bootstrap_strata,
      sample_ids = rownames(x)
    )
  } else if (stratify_requested) {
    if (identical(family, "cox") || (is.matrix(y) && ncol(y) > 1L)) {
      .abort_stablr(
        "stablr_bootstrap_stratification_infeasible",
        paste0(
          "`stratify_bootstrap = TRUE` derives the stratification design from ",
          "`y`, but for `family = \"cox\"` `y` is a two-column `Surv` object. ",
          "Its `interaction(time, status)` strata are effectively singletons, so ",
          "the stratified sampler cannot allocate a subsample. Supply an ",
          "explicit `bootstrap_strata` design (one categorical value per ",
          "sample) or set `stratify_bootstrap = FALSE`."
        ),
        family = family
      )
    }
    bootstrap_strata <- y
  }
  bootstrap_strata_ids <- .bootstrap_strata_ids(
    strata = bootstrap_strata,
    n = nrow(x),
    arg = "bootstrap_strata"
  )

  n_samples   <- nrow(x)
  n_features  <- ncol(x)
  feat_names  <- colnames(x)
  if (is.null(feat_names)) feat_names <- paste0("x.", seq_len(n_features))

  # ---- Lambda grid ----------------------------------------------------------
  if (identical(lambda_grid, "auto") && verbose) {
    message("Computing auto lambda grid...")
  }
  resolved_grid <- .resolve_lambda_grid(
    lambda_grid,
    keys = character(0),
    x = x, y = y,
    family = family,
    n_lambda = n_lambda,
    l1_ratio = l1_ratio
  )
  lambda_grid <- resolved_grid$grid

  # ---- Derived RNG streams (audit M-5, V-7) ---------------------------------
  # When `random_state` is set we draw three independent integer seeds from a
  # single seeded RNG: one for artificial-feature generation, one for the
  # bootstrap-index draws, and one used as the base of per-iteration seeds
  # passed into the bootstrap loop.  This keeps the three RNG streams
  # independent (artificial features cannot consume RNG that boot indices
  # need) and makes the parallel and sequential loops bit-identical because
  # each iteration seeds its own RNG before calling the learner adapter.
  # That parallel equivalence is conditional on `random_state`: without it
  # there are no per-iteration seeds, so a parallel run cannot reproduce a
  # sequential one.  Every seeding here goes through `.with_local_seed`, so
  # the caller's `.Random.seed` survives the call.
  rng_kind <- RNGkind()
  if (is.null(random_state) && workers > 1L) {
    warning(
      "`workers > 1` with `random_state = NULL` cannot be reproduced: there ",
      "is no per-iteration seed, so the parallel bootstrap draws depend on ",
      "ambient RNG state. Set `random_state` to make the run reproducible.",
      call. = FALSE
    )
  }
  parallel_backend_available <- (workers > 1L) &&
    requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)
  if (workers > 1L && !parallel_backend_available) {
    warning(
      "`workers > 1` was requested but the `furrr` and `future` packages are ",
      "not both available; the bootstrap loop runs sequentially instead. ",
      "Install them (or pre-load `library(future)`) for real parallelism.",
      call. = FALSE
    )
  }

  art_seed  <- NULL
  boot_seed <- NULL
  iter_seed_base <- NULL
  if (!is.null(random_state)) {
    rng_seeds <- .with_local_seed(
      random_state,
      sample.int(.Machine$integer.max, size = 3L)
    )
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
  sgl_groups <- .resolve_sgl_feature_groups(
    base_learner         = base_learner,
    x_original           = x,
    x_augmented          = x_fit,
    feature_groups       = feature_groups,
    corr_group_threshold = corr_group_threshold,
    noise_col_indices    = noise_col_indices,
    artificial_type      = artificial_type,
    actual_type          = artificial_provenance$actual_type
  )

  # ---- Learner adapter ------------------------------------------------------
  # Each batch adapter accepts the full lambda_grid and returns a logical
  # matrix (n_features × n_lambdas).  For a single alpha this is one solver
  # path call per bootstrap; a grid with k distinct alpha values costs k
  # calls, and sparse_group_lasso additionally pays n_lambda coef() calls per
  # path call (see the "Cost of the bootstrap loop" section in the docs).
  batch_adapter <- .select_batch_adapter(
    base_learner       = base_learner,
    family             = family,
    adaptive_gamma     = adaptive_gamma,
    adaptive_epsilon   = adaptive_epsilon,
    sgl_feature_groups = sgl_groups$groups
  )

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
  # The three samplers below are deliberately not the exported ones for the
  # classification branch; see `.classification_bootstrap_sampler()` for why the
  # feasibility contract differs.  All of them ultimately call the same draw
  # primitives (`.stratified_bootstrap_indices()`,
  # `.unstratified_group_bootstrap_indices()`, `.sample_groups_until()`), so the
  # RNG consumption of a given configuration is identical either way.
  stratified_repair <- new.env(parent = emptyenv())
  boot_sampler <- if (!is.null(modeled_classes)) {
    .classification_bootstrap_sampler(
      y = y,
      groups = groups,
      n_subsamples = n_subsamples,
      replace = replace,
      strata_ids = bootstrap_strata_ids,
      family = family,
      class_levels = modeled_classes,
      repair_log = stratified_repair
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

  # The bootstrap draws and the per-iteration seeds are generated inside one
  # `.with_local_seed(boot_seed, ...)` scope: the inner `set.seed(iter_seed_base)`
  # only affects the iteration seeds, and the caller's `.Random.seed` is
  # restored on exit either way.
  draws <- .with_local_seed(boot_seed, {
    list(
      boot_indices = lapply(seq_len(n_bootstraps), boot_sampler),
      iter_seeds = if (!is.null(iter_seed_base)) {
        set.seed(iter_seed_base)
        sample.int(.Machine$integer.max, size = n_bootstraps)
      } else {
        NULL
      }
    )
  })
  boot_indices <- draws$boot_indices
  iter_seeds   <- draws$iter_seeds

  # `n_subsamples` is a lower bound on the realised draw: whole groups are never
  # trimmed, and the class-feasibility repair below can raise the target further.
  # Recording the realised length per bootstrap keeps that visible, and a short
  # draw is an error rather than a smaller-than-requested fit reported against
  # the requested `sample_fraction`.
  subsample_sizes <- vapply(boot_indices, length, integer(1L))
  if (any(subsample_sizes < n_subsamples)) {
    .abort_stablr(
      "stablr_bootstrap_infeasible",
      paste0(
        "A bootstrap replicate was drawn with ", min(subsample_sizes),
        " rows, fewer than the requested n_subsamples = ", n_subsamples,
        "; ", sum(subsample_sizes < n_subsamples), " of ", n_bootstraps,
        " replicates are affected."
      ),
      n_subsamples = n_subsamples,
      realised = min(subsample_sizes)
    )
  }

  art_rows <- if (!is.null(artificial_type)) {
    n_features + seq_len(n_injected)
  } else NULL

  if (verbose) {
    message(
      "Running STABL bootstrap loop (", n_bootstraps, " bootstraps, ",
      nrow(lambda_grid), " lambdas, ",
      .stabl_path_calls_per_bootstrap(base_learner, lambda_grid),
      " solver path call(s) per bootstrap)..."
    )
  }

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
    workers        = workers,
    use_furrr      = parallel_backend_available
  )
  stabl_scores_    <- boot_result$stabl_scores_
  stabl_scores_art <- boot_result$stabl_scores_art

  # ---- FDP+ diagnostic thresholding ----------------------------------------
  # `hard_threshold` and the FDP+ minimiser are two different selection rules
  # and only one of them can win.  When the caller supplied a hard threshold the
  # FDP+ search is skipped rather than computed and discarded, and the fact is
  # recorded in `fdp_threshold_suppressed`.  The artificial-feature scores and
  # provenance are still produced.
  fdp <- NULL
  fdp_threshold_suppressed <- NA
  if (!is.null(artificial_type)) {
    if (is.null(hard_threshold)) {
      fdp <- compute_fdp_plus(
        stabl_scores            = stabl_scores_,
        stabl_scores_artificial = stabl_scores_art,
        artificial_proportion   = effective_artificial_proportion,
        fdr_threshold_range     = fdr_threshold_range
      )
      fdp_threshold_suppressed <- FALSE
    } else {
      fdp_threshold_suppressed <- TRUE
    }
  }

  # ---- Fit-time threshold resolution and explore fallback -------------------
  # The explore fallback is a function of the stability scores and n_explore
  # only, so it is resolved once here and every accessor reads the same answer
  # instead of recomputing (and potentially disagreeing about) it.
  fit_time <- .resolve_fit_time_selection(
    max_scores   = matrixStats::rowMaxs(stabl_scores_),
    feat_names   = feat_names,
    threshold    = if (!is.null(hard_threshold)) hard_threshold else
      if (!is.null(fdp)) fdp$fdr_min_threshold else NULL,
    explore      = explore,
    n_explore    = n_explore
  )

  # ---- Assemble S3 result ---------------------------------------------------
  structure(
    list(
      stabl_scores_         = stabl_scores_,
      stabl_scores_artificial_ = stabl_scores_art,
      fitted_lambda_grid    = lambda_grid,
      lambda_grid_source    = resolved_grid$source,
      fdr_min_threshold_    = if (!is.null(fdp)) fdp$fdr_min_threshold else NULL,
      FDRs_                 = if (!is.null(fdp)) fdp$FDRs              else NULL,
      min_fdr_              = if (!is.null(fdp)) fdp$min_fdr            else NULL,
      fdrs_table_           = if (!is.null(fdp)) fdp$fdrs_table         else NULL,
      fdr_threshold_range   = fdr_threshold_range,
      fdp_threshold_suppressed = fdp_threshold_suppressed,
      hard_threshold        = hard_threshold,
      artificial_type       = artificial_type,
      artificial_provenance = artificial_provenance,
      artificial_proportion = artificial_proportion,
      effective_artificial_proportion = effective_artificial_proportion,
      stratify_bootstrap    = !is.null(stratify_design),
      stratify_bootstrap_requested = stratify_requested,
      bootstrap_strata_source = stratify_design,
      bootstrap_strata_applied = !is.null(bootstrap_strata_ids),
      bootstrap_strata_levels = if (!is.null(bootstrap_strata_ids)) {
        sort(unique(as.character(bootstrap_strata_ids)))
      } else NULL,
      bootstrap_subsample_target_ = n_subsamples,
      bootstrap_subsample_sizes_ = subsample_sizes,
      bootstrap_stratified_repair_ = if (is.null(stratified_repair$target_raised)) {
        NULL
      } else {
        list(
          target_raised = TRUE,
          requested_target = stratified_repair$requested_target,
          repair_target = stratified_repair$repair_target,
          n_classes = stratified_repair$n_classes,
          minimum_per_class = stratified_repair$minimum_per_class,
          n_repairs = stratified_repair$n_repairs
        )
      },
      explore               = explore,
      n_explore             = as.integer(n_explore),
      explore_features_     = fit_time$explore_features,
      explore_applied_      = fit_time$explore_applied,
      support_              = fit_time$support,
      sgl_feature_groups_   = sgl_groups$groups,
      sgl_group_provenance  = sgl_groups$provenance,
      rng_kind_             = rng_kind,
      feature_names         = feat_names,
      n_features_in_        = n_features
    ),
    class = "stabl_fit"
  )
}

# Solver path calls per bootstrap, for the progress message.  A single-alpha
# grid costs one call; a grid with k distinct `alpha` values costs k, because
# the batch adapters fit one path per distinct mixing parameter.
.stabl_path_calls_per_bootstrap <- function(base_learner, lambda_grid) {
  has_alpha_col <- "alpha" %in% names(lambda_grid)
  if (!has_alpha_col) {
    return(1L)
  }
  if (base_learner %in% c("lasso", "adaptive_lasso")) {
    return(1L)
  }
  as.integer(length(unique(lambda_grid[["alpha"]])))
}

# Resolve the fit-time support mask and the explore fallback candidate set.
#
# `explore_features_` is threshold-independent: it is simply the `n_explore`
# highest-importance features, so every accessor reports the same fallback
# whatever threshold it is asked about.  `support_` is the mask for the
# threshold that actually drove this fit, and is `NULL` when no threshold was
# resolvable (an object whose thresholds were removed after fitting).
.resolve_fit_time_selection <- function(max_scores, feat_names, threshold,
                                         explore, n_explore) {
  names(max_scores) <- feat_names
  n_exp <- min(as.integer(n_explore), length(max_scores))
  explore_features <- if (isTRUE(explore) && n_exp > 0L) {
    names(sort(max_scores, decreasing = TRUE))[seq_len(n_exp)]
  } else {
    NULL
  }

  if (is.null(threshold)) {
    return(list(
      support = NULL,
      explore_features = explore_features,
      explore_applied = NA
    ))
  }

  mask <- max_scores > threshold
  explore_applied <- FALSE
  if (!any(mask) && !is.null(explore_features)) {
    mask[explore_features] <- TRUE
    explore_applied <- TRUE
  }
  list(
    support = mask,
    explore_features = explore_features,
    explore_applied = explore_applied
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

# Learner-feasible bootstrap sampler for `family = "binomial"`/`"multinomial"`.
#
# This deliberately does not call `classic_bootstrap_indices()` /
# `group_bootstrap_indices()` even though both are exported.  The exported
# samplers guarantee only "at least two distinct classes" via
# `.guard_diverse_bootstrap()`; a classification fit additionally needs a
# *minimum number of observations* in every modeled class, and when the groups
# are class-pure the only way to guarantee that is to raise the row target above
# `n_subsamples` and redraw stratified.  Routing this through the exported
# samplers would silently weaken the contract (a 1-vs-19 draw is "diverse") and
# would lose the typed `stablr_bootstrap_infeasible` aborts below.  The draw
# primitives are the same ones the exported samplers use, so RNG consumption is
# unchanged; only the acceptance policy differs.  The raised target and the
# realised overshoot are reported through `repair_log` and land in the fit as
# `bootstrap_stratified_repair_` / `bootstrap_subsample_sizes_`.
.classification_bootstrap_sampler <- function(y, groups, n_subsamples,
                                              replace, strata_ids, family,
                                              class_levels,
                                              minimum_per_class = 2L,
                                              max_redraws = 1000L,
                                              repair_log = NULL) {
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
      if (!is.null(repair_log)) {
        repair_log$target_raised <- TRUE
        repair_log$requested_target <- n_subsamples
        repair_log$repair_target <- repair_target
        repair_log$n_classes <- length(class_levels)
        repair_log$minimum_per_class <- minimum_per_class
        repair_log$n_repairs <- (repair_log$n_repairs %||% 0L) + 1L
      }
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

# The four `glmnet` response families this package documents and supports for
# `stabl_fit()`.  Validating here turns a family typo ("gaussain", "binomal")
# into a typed condition instead of an opaque `glmnet` internal error.  The set
# is deliberately not "whatever glmnet accepts": `.modeled_class_levels()`,
# `.feature_abs_coefs()` and `.stop_if_cox_sgl()` all branch on these four.
.stabl_family_levels <- c("gaussian", "binomial", "multinomial", "cox")

.validate_stabl_family <- function(family) {
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("`family` must be a single character string.", call. = FALSE)
  }
  if (!family %in% .stabl_family_levels) {
    stop(
      "`family` must be one of: ",
      paste(sprintf("\"%s\"", .stabl_family_levels), collapse = ", "),
      "; got \"", family, "\".",
      call. = FALSE
    )
  }
  family
}

# `fdr_threshold_range` is swept by compute_fdp_plus(); validate it here so an
# out-of-range or empty range is reported against the fitting call.
.validate_fdr_threshold_range <- function(fdr_threshold_range) {
  if (!is.numeric(fdr_threshold_range) || !is.null(dim(fdr_threshold_range)) ||
      length(fdr_threshold_range) == 0L ||
      any(!is.finite(fdr_threshold_range))) {
    stop("`fdr_threshold_range` must be a non-empty finite numeric vector.",
         call. = FALSE)
  }
  if (any(fdr_threshold_range < 0 | fdr_threshold_range > 1)) {
    stop("`fdr_threshold_range` must lie in [0, 1].", call. = FALSE)
  }
  as.numeric(fdr_threshold_range)
}

.resolve_sgl_feature_groups <- function(base_learner,
                                        x_original,
                                        x_augmented,
                                        feature_groups,
                                        corr_group_threshold,
                                        noise_col_indices,
                                        artificial_type = NULL,
                                        actual_type = NULL) {
  if (base_learner != "sparse_group_lasso") {
    return(list(groups = NULL, provenance = NULL))
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

  generator <- if (!is.null(actual_type)) actual_type else artificial_type
  inherit_groups <- .sgl_noise_group_inheritance(generator)
  groups <- .append_noise_groups(
    grp, noise_col_indices, ncol(x_augmented),
    inherit_groups = identical(inherit_groups, "on")
  )

  list(
    groups = groups,
    provenance = list(
      group_source = if (has_explicit) "feature_groups" else "corr_group_threshold",
      noise_group_inheritance = if (is.null(noise_col_indices)) {
        NA_character_
      } else {
        inherit_groups
      },
      n_groups = length(unique(groups)),
      n_noise_features = length(noise_col_indices)
    )
  )
}

# Should an artificial decoy inherit the sparse-group assignment of the original
# feature it was derived from?
#
# Default OFF for permutation decoys.  FDP+ compares the stability of the real
# block against the artificial block, and it is only a like-for-like comparison
# when both blocks face the same kind of group-lasso penalty.  Putting a
# permutation decoy into its source feature's group couples the decoy's
# selection to the real block's group structure, so the decoy's stability score
# partly measures the real block rather than the decoy's own noise.
#
# Knockoff decoys are the exception: a knockoff is *exchangeable* with its
# source by construction, so co-grouping them is the honest representation and
# inheritance is on.  The resolved value is recorded in
# `sgl_group_provenance$noise_group_inheritance`; pass `"on"`/`"off"` to
# `.append_noise_groups()` directly to override.
.sgl_noise_group_inheritance <- function(generator) {
  if (is.null(generator)) {
    return("off")
  }
  if (is.character(generator) && length(generator) == 1L &&
      startsWith(generator, "knockoff")) {
    return("on")
  }
  "off"
}

# `.build_corr_groups()` collapses features into groups using a
# correlation-threshold graph.  Two behaviours a reader would otherwise have to
# reverse engineer:
#
# 1. `corr[is.na(corr)] <- 0` replaces every undefined correlation with a
#    numerical zero.  Undefined arises for zero-variance columns and for pairs
#    with no complete overlap under `use = "pairwise.complete.obs"`.  A
#    constant column has no estimable correlation, so there is no evidence for
#    grouping it with anything; zero is the "no association" value, so such a
#    column is then treated exactly like any genuinely uncorrelated column --
#    neither force-joined nor treated as maximally correlated.  (Whether a zero
#    ends up above the cutoff is the same question every uncorrelated pair faces,
#    and is answered by the percentile below.)  Propagating `NA` into the edge
#    test instead would make the resulting partition depend on `NA` handling
#    rather than on the data.
# 2. The `- 0.1` subtracted from the percentile is a port of the Python
#    reference (`stabl/stabl.py`:
#    `threshold = np.percentile(corr_val, perc) - 0.1`).  Correlation
#    thresholds are otherwise not comparable across implementations; the offset
#    is what makes the R port reproduce the Python partition, and it is a fixed
#    constant rather than a tuned parameter.
.build_corr_groups <- function(x, percentile) {
  if (!is.numeric(percentile) || length(percentile) != 1L ||
      is.na(percentile) || !is.finite(percentile) ||
      percentile <= 0 || percentile > 100) {
    stop("`corr_group_threshold` must be a numeric scalar in (0, 100].", call. = FALSE)
  }

  p <- ncol(x)
  if (p < 1L) {
    stop("`x` must contain at least one feature column.", call. = FALSE)
  }
  if (p == 1L) {
    return(1L)
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

# Extend an original-feature group assignment over the artificial block.
#
# `inherit_groups = TRUE` (the default, preserved for direct callers of this
# helper) gives each decoy the group of the original feature recorded in
# `noise_col_indices`.  `inherit_groups = FALSE` gives every decoy a fresh group
# id, so the artificial block is penalised independently of the real block and
# the two blocks stay comparable under the group-lasso penalty.  `stabl_fit()`
# passes `inherit_groups = FALSE` for permutation decoys; see
# `.sgl_noise_group_inheritance()`.
.append_noise_groups <- function(groups, noise_col_indices, total_p,
                                inherit_groups = TRUE) {
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
    if (isTRUE(inherit_groups) && !is.na(src) && src >= 1L && src <= n_orig) {
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
    n_features, feat_names, art_rows, artificial_type, workers,
    use_furrr = NULL
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

  if (is.null(use_furrr)) {
    use_furrr <- (workers > 1L) &&
      requireNamespace("furrr",  quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)
  }

  if (use_furrr) {
    batch_adapter <- .make_parallel_safe_batch_adapter(batch_adapter)
  }

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

# Internal: make a factory-created learner closure self-contained before it is
# serialized by future/furrr. Automatic global discovery sees `batch_adapter`
# but does not recursively export the package-private helpers referenced inside
# that nested closure. Binding the exact transitive helpers avoids lookup in an
# installed (and potentially different) stablr namespace on worker processes.
.make_parallel_safe_batch_adapter <- function(batch_adapter) {
  factory_env <- environment(batch_adapter)
  captured_names <- ls(factory_env, all.names = TRUE)
  captured <- if (length(captured_names)) {
    mget(captured_names, envir = factory_env, inherits = FALSE)
  } else {
    list()
  }
  worker_env <- list2env(captured, parent = baseenv())
  helpers <- list(
    .with_glmnet_numerical_infeasibility = .with_glmnet_numerical_infeasibility,
    .with_typed_numerical_infeasibility = .with_typed_numerical_infeasibility,
    .is_glmnet_numerical_error = .is_glmnet_numerical_error,
    .abort_numerical_infeasibility = .abort_numerical_infeasibility,
    .abort_stablr = .abort_stablr,
    .new_stablr_error = .new_stablr_error,
    .feature_abs_coefs_batch = .feature_abs_coefs_batch,
    .feature_abs_coefs = .feature_abs_coefs,
    .feature_abs_coefs_sparsegl_batch = .feature_abs_coefs_sparsegl_batch,
    .feature_abs_coefs_sparsegl = .feature_abs_coefs_sparsegl,
    .near_lambda_match = .near_lambda_match,
    rowMaxs = matrixStats::rowMaxs
  )
  for (name in names(helpers)) {
    helper <- helpers[[name]]
    if (is.function(helper) && !identical(environment(helper), baseenv())) {
      environment(helper) <- worker_env
    }
    assign(name, helper, envir = worker_env)
  }
  environment(batch_adapter) <- worker_env
  batch_adapter
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
