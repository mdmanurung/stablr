# Oblique random forest (aorsf) learner adapter — experimental extension
# (plan item D2).
#
# Statistical design (see the improvement plan, items D2/F5):
# * Negation importance, NOT permutation importance: in an oblique forest
#   every split uses a linear combination of predictors, and negation
#   importance (re-scoring out-of-bag cases with one predictor's coefficients
#   negated) discriminates relevant from irrelevant predictors better than
#   permutation importance does in this setting (Jaeger et al. 2023,
#   DOI:10.1080/10618600.2023.2231048; aorsf documentation).
# * The regularisation path is an importance-threshold path, exactly like the
#   axis-aligned RF adapter (plan D1): thresholds are quantiles of the pooled
#   negation-importance distribution, swept from the sparse end (high
#   threshold) to the dense end (low threshold) — the same sparse-to-dense
#   nesting the glmnet lambda path provides, so max-over-grid stability
#   frequencies and the FDP+ machinery work unchanged.
# * Null pairing: identical to the axis-aligned RF (plan D1).  Oblique random
#   forests are multivariate learners, so real null features correlated with
#   true signals receive leaked importance while permuted nulls
#   (artificial_type = "random_permutation") do not; FDP+ is therefore
#   anti-conservative with permutation nulls and model-X knockoff nulls
#   (knockoff_equi / knockoff_mvr) must be preferred.  Combining
#   "oblique_rf" with "random_permutation" warns.
# * D2 caveat: importance is spread across the linear combinations used in
#   each split, so per-feature selection frequencies are noisier than for
#   axis-aligned random forests; plan for more bootstraps (or larger forests)
#   before trusting small stability-score differences between features.

#' Build an Oblique Random-Forest Learner Adapter (Fenced)
#'
#' Returns a closure that fits ONE `aorsf` oblique random forest on a single
#' bootstrap subsample and returns a logical selection mask over the features,
#' mirroring the role of [make_glmnet_adapter()] and [make_ranger_adapter()]
#' for tree-based learners.
#'
#' @section Selection rule:
#' Instead of a penalty path, the oblique-RF adapter sweeps an
#' **importance-threshold path**.  The grid (see
#' [auto_oblique_threshold_grid()]) contains quantiles of the pooled
#' negation-importance distribution; given a 1-row grid `lambda_val`, a
#' feature is selected when its negation importance is `>=` the row's
#' threshold.  Lower thresholds therefore give supersets of the selections at
#' higher thresholds — the same sparse-to-dense nesting the glmnet lambda path
#' provides, so the stability-accumulation and FDP+ logic in [stabl_fit()]
#' work unchanged.
#'
#' @section Statistical guardrails (read before using):
#' * **Negation importance, not permutation importance.**  In oblique forests
#'   negation importance discriminates relevant from irrelevant predictors
#'   better than permutation importance (Jaeger et al. 2023), so it is the
#'   default here.
#' * **Pair with model-X knockoff artificial features.**  An oblique random
#'   forest is a multivariate learner: real null features correlated with true
#'   signals receive leaked importance, and permuted null features
#'   (`artificial_type = "random_permutation"`) do **not** reproduce that
#'   leakage, so the FDP+ estimate is anti-conservative.  Model-X knockoff
#'   nulls (`"knockoff_equi"` / `"knockoff_mvr"`) preserve the correlation
#'   structure and are the only statistically sound pairing.  Combining
#'   `base_learner = "oblique_rf"` with
#'   `artificial_type = "random_permutation"` in [stabl_fit()] emits a package
#'   warning to this effect.
#' * **Noisier per-feature selection frequencies (D2 caveat).**  Importance is
#'   spread across the linear combinations used in each split, so
#'   per-feature selection frequencies are noisier than for axis-aligned
#'   random forests ([make_ranger_adapter()]).  Small stability-score
#'   differences between features are less reliable; prefer more bootstraps.
#' * **Fenced (WS-3).**  No predeclared remedy reached the calibration
#'   criterion (mean FDP <= 0.15 on both designs; best remedy 0.18 / 0.39).
#'   Knockoff decoy importance is inflated for oblique forests, so FDP+ is
#'   anti-conservative.  Use [make_ranger_adapter()] with knockoff_equi for
#'   knockoff-calibrated selection instead.
#'
#' @param num_trees Positive integer; number of trees per forest (`aorsf`'s
#'   `n_tree`).  Default `500L`.
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`
#'   for use inside [stabl_fit()], where `lambda_val` is a 1-row `data.frame`
#'   from the threshold grid (see [auto_oblique_threshold_grid()]).
#'
#' @details
#' This adapter factory is primarily an internal backend used by [stabl_fit()]
#' (`base_learner = "oblique_rf"`).  The response family is inferred by
#' `aorsf` from `y` (numeric -> regression, factor -> classification,
#' `Surv` -> survival); [stabl_fit()] coerces `y` for the binomial/multinomial
#' families.  The per-tree seeds are drawn from the caller's RNG stream so
#' that results are reproducible under [stabl_fit()]'s `random_state` and
#' identical across sequential and parallel execution.
#'
#' @seealso [auto_oblique_threshold_grid()], [make_ranger_adapter()],
#'   [make_glmnet_adapter()], [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("aorsf", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(200L), 20L, 10L,
#'               dimnames = list(NULL, paste0("f", 1:10)))
#'   y <- x[, 1] * x[, 2] + rnorm(20L)
#'   grid  <- auto_oblique_threshold_grid(x, y, n_grid = 5L, num_trees = 100L,
#'                                        seed = 1L)
#'   adapter <- make_aorsf_adapter(num_trees = 100L)
#'   mask <- adapter(x, y, grid[1L, , drop = FALSE])
#'   cat("selected at sparsest threshold:", sum(mask), "of", ncol(x), "\n")
#' }
#' @export
make_aorsf_adapter <- function(num_trees = 500L) {
  .require_pkg("aorsf", "for base_learner = \"oblique_rf\"")
  num_trees <- .validate_scalar_integer_like(num_trees, "num_trees", min = 1L)
  batch_adapter <- .make_aorsf_batch_adapter(num_trees = num_trees)

  function(x, y, lambda_val) {
    batch_adapter(x, y, lambda_val)[, 1L]
  }
}

#' Build a Data-Driven Importance-Threshold Grid for the Oblique-RF Learner
#' (Fenced)
#'
#' Fits one `aorsf` oblique random forest with negation importance on the full
#' training data and returns the importance-threshold grid as a pre-expanded
#' `data.frame` suitable for the `lambda_grid` argument of [stabl_fit()] with
#' `base_learner = "oblique_rf"`.  This is the oblique-forest analog of
#' [auto_rf_threshold_grid()].
#'
#' Thresholds are quantiles of the pooled negation-importance distribution at
#' probabilities `seq(0.5, 0.999, length.out = n_grid)`, sorted in decreasing
#' order so the path runs from the sparse end (high threshold) to the dense
#' end (low threshold), mirroring the glmnet lambda-path convention.
#'
#' The grid is computed once on the full data (exactly like
#' [auto_lambda_grid()]); each bootstrap subsample is subsequently scored
#' against these fixed thresholds with its own forest.
#'
#' @param x Numeric matrix of predictors (samples by features).  Column names
#'   are recommended; unnamed matrices are named `x1`, `x2`, ... internally.
#' @param y Outcome vector (numeric for regression, factor for
#'   classification, `survival::Surv` for survival).
#' @param n_grid Positive integer; number of threshold values.  Default
#'   `30L` (matches [auto_lambda_grid()]'s default `n_lambda`).
#' @param num_trees Positive integer; number of trees in the grid-defining
#'   forest.  Default `500L`.
#' @param seed Integer scalar or `NULL`; seed for the grid-defining forest.
#'   When `NULL` (default) the forest is seeded from the caller's RNG stream.
#'
#' @return A `data.frame` with columns `lambda` and `threshold` (identical
#'   values: the importance thresholds).  The `lambda` column follows the
#'   grid convention required by [stabl_fit()]; `threshold` documents the
#'   oblique-RF-specific meaning.
#'
#' @seealso [make_aorsf_adapter()], [auto_rf_threshold_grid()],
#'   [auto_lambda_grid()], [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("aorsf", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(200L), 20L, 10L,
#'               dimnames = list(NULL, paste0("f", 1:10)))
#'   y <- x[, 1] * x[, 2] + rnorm(20L)
#'   grid <- auto_oblique_threshold_grid(x, y, n_grid = 5L, num_trees = 100L,
#'                                       seed = 1L)
#'   head(grid)
#' }
#' @export
auto_oblique_threshold_grid <- function(x, y, n_grid = 30L, num_trees = 500L,
                                        seed = NULL) {
  .require_pkg("aorsf", "for base_learner = \"oblique_rf\"")
  n_grid    <- .validate_scalar_integer_like(n_grid, "n_grid", min = 1L)
  num_trees <- .validate_scalar_integer_like(num_trees, "num_trees", min = 1L)
  seed      <- .validate_scalar_integer_like(seed, "seed", allow_null = TRUE)

  # Draw the per-tree seeds from the caller's RNG stream (deterministic under
  # `seed`, or under whatever seed the caller has set when `seed = NULL`).
  tree_seeds <- .with_local_seed(
    seed,
    as.integer(sample.int(.Machine$integer.max, num_trees))
  )
  importance <- .aorsf_fit_negation_importance(x, y, num_trees, tree_seeds)
  .rf_threshold_grid_from_importance(importance, n_grid)
}

# ---- Internal helpers ---------------------------------------------------------

# Fit one aorsf oblique forest and return its negation importance as a named
# numeric vector ordered like the columns of `x`.  aorsf fits through a
# data.frame + formula interface, so the (possibly non-syntactic) predictor
# names are made syntactic for the fit and mapped back afterwards.
.aorsf_fit_negation_importance <- function(x, y, num_trees, tree_seeds) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  raw_names  <- colnames(x)
  safe_names <- make.names(raw_names, unique = TRUE)

  outcome_name <- ".stabl_outcome"
  while (outcome_name %in% safe_names) outcome_name <- paste0(outcome_name, "x")

  data <- as.data.frame(x)
  names(data) <- safe_names
  data[[outcome_name]] <- y
  data <- data[, c(outcome_name, safe_names), drop = FALSE]

  fit <- aorsf::orsf(
    data        = data,
    formula     = stats::as.formula(paste(outcome_name, "~ .")),
    n_tree      = num_trees,
    importance  = "negate",
    tree_seeds  = tree_seeds,
    attach_data = FALSE
  )
  # With importance = "negate" the fit stores the negation importances in
  # `fit$importance` (a named numeric vector; aorsf orders it by its own
  # internal variable order, not by the formula's predictor order).
  vi <- fit$importance

  # Reorder to the caller's column order and restore the original names.
  vi <- vi[safe_names]
  if (length(vi) != ncol(x) || anyNA(vi) || any(!is.finite(vi))) {
    stop(
      "aorsf did not return usable negation importances for this ",
      "bootstrap subsample.",
      call. = FALSE
    )
  }
  stats::setNames(as.numeric(vi), raw_names)
}

# ---- Internal batch adapter ---------------------------------------------------
# Signature matches the other batch adapters:
#   function(x, y, lambda_grid) -> logical matrix (n_features x n_lambdas)
# ONE aorsf forest is fitted per bootstrap subsample; masks for the whole
# grid are extracted from that single negation-importance vector, so the
# bootstrap loop costs n_bootstraps forests total (not n_bootstraps x
# n_lambdas).
.make_aorsf_batch_adapter <- function(family = "gaussian", num_trees = 500L) {
  num_trees <- .validate_scalar_integer_like(num_trees, "num_trees", min = 1L)
  function(x, y, lambda_grid) {
    # aorsf infers the task from `y`: factor -> classification.  Coerce
    # 0/1-numeric or character class labels so binomial/multinomial fits are
    # classification forests, mirroring the glmnet family semantics.
    if (identical(family, "binomial") || identical(family, "multinomial")) {
      if (!is.factor(y)) y <- factor(y)
    }
    n_lambdas  <- nrow(lambda_grid)
    n_features <- ncol(x)
    thresholds <- .rf_thresholds_from_grid(lambda_grid)

    # aorsf grows each tree with its own seed (tree_seeds) and does not
    # consume R's stream otherwise.  Draw the seeds from the caller's (already
    # locally seeded) RNG stream — stabl_fit() wraps every adapter call in
    # .with_local_seed(iter_seeds[[i]], ...) — so each bootstrap iteration is
    # reproducible under `random_state` and bit-identical across sequential
    # and parallel execution.
    tree_seeds <- as.integer(sample.int(.Machine$integer.max, num_trees))

    importance <- .aorsf_fit_negation_importance(x, y, num_trees, tree_seeds)

    masks <- vapply(
      thresholds,
      function(threshold) importance >= threshold,
      logical(n_features)
    )
    dimnames(masks) <- NULL
    masks
  }
}
