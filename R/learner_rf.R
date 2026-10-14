# Random-forest (ranger) learner adapter — experimental extension (plan item D1).
#
# Statistical design (see the improvement plan, items D1/F5):
# * Permutation importance (mean decrease in permutation-based prediction
#   error), NOT Gini importance: Gini importance is biased towards
#   high-cardinality features and its correlation leak is harder to reason
#   about; permutation importance is the better-studied choice (Boruta is
#   likewise built on mean-decrease-accuracy).
# * The regularisation path is an importance-threshold path: thresholds are
#   quantiles of the pooled permutation-importance distribution, swept from
#   the sparse end (high threshold) to the dense end (low threshold) — the
#   same sparse-to-dense nesting the glmnet lambda path provides, so
#   max-over-grid stability frequencies and the FDP+ machinery work unchanged.
# * Null pairing: random forests are multivariate learners, so real null
#   features correlated with true signals receive leaked importance while
#   permuted nulls (artificial_type = "random_permutation") do not.  FDP+ is
#   therefore anti-conservative with permutation nulls; model-X knockoff
#   nulls (knockoff_equi / knockoff_mvr) preserve the correlation structure
#   and must be preferred.  Combining "rf" with "random_permutation" warns.

#' Build a Random-Forest Learner Adapter (Calibrated, Scoped)
#'
#' Returns a closure that fits ONE `ranger` random forest on a single
#' bootstrap subsample and returns a logical selection mask over the features,
#' mirroring the role of [make_glmnet_adapter()] for tree-based learners.
#'
#' @section Selection rule:
#' Instead of a penalty path, the RF adapter sweeps an **importance-threshold
#' path**.  The grid (see [auto_rf_threshold_grid()]) contains quantiles of
#' the pooled permutation-importance distribution; given a 1-row grid
#' `lambda_val`, a feature is selected when its permutation importance is
#' `>= ` the row's threshold.  Lower thresholds therefore give supersets of
#' the selections at higher thresholds — the same sparse-to-dense nesting the
#' glmnet lambda path provides, so the stability-accumulation and FDP+ logic
#' in [stabl_fit()] work unchanged.
#'
#' @section Statistical guardrails (read before using):
#' * **Permutation importance, not Gini.**  Gini importance is biased towards
#'   high-cardinality features; permutation importance is the better-studied
#'   choice for null calibration (Boruta likewise builds on
#'   mean-decrease-accuracy).
#' * **Pair with model-X knockoff artificial features.**  A random forest is a
#'   multivariate learner: real null features correlated with true signals
#'   receive leaked importance, and permuted null features
#'   (`artificial_type = "random_permutation"`) do **not** reproduce that
#'   leakage, so the FDP+ estimate is anti-conservative.  Model-X knockoff
#'   nulls (`"knockoff_equi"` / `"knockoff_mvr"`) preserve the correlation
#'   structure and are the only statistically sound pairing.  Combining
#'   `base_learner = "rf"` with `artificial_type = "random_permutation"` in
#'   [stabl_fit()] emits a package warning to this effect.
#' * **Calibrated within a measured scope (WS-2).**  With knockoff_equi
#'   decoys the mean empirical FDP is 0.025 (linear design) and 0.065
#'   (low-dimensional nonlinear design), meeting the predeclared criterion;
#'   it is NOT met at p = 50 with strong correlation (0.44) or for binary
#'   outcomes (0.26).  Simulation evidence, not a guarantee.
#'
#' @param num_trees Positive integer; number of trees per forest
#'   (`ranger`'s `num.trees`).  Default `500L`.
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`
#'   for use inside [stabl_fit()], where `lambda_val` is a 1-row `data.frame`
#'   from the threshold grid (see [auto_rf_threshold_grid()]).
#'
#' @details
#' This adapter factory is primarily an internal backend used by [stabl_fit()]
#' (`base_learner = "rf"`).  The response family is inferred by `ranger` from
#' `y` (numeric -> regression, factor -> classification, `Surv` -> survival);
#' [stabl_fit()] coerces `y` for the binomial/multinomial families.  The
#' forest seed is drawn from the caller's RNG stream so that results are
#' reproducible under [stabl_fit()]'s `random_state` and identical across
#' sequential and parallel execution.
#'
#' @seealso [auto_rf_threshold_grid()], [make_glmnet_adapter()], [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(200L), 20L, 10L,
#'               dimnames = list(NULL, paste0("f", 1:10)))
#'   y <- x[, 1] * x[, 2] + rnorm(20L)
#'   grid  <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L,
#'                                   seed = 1L)
#'   adapter <- make_ranger_adapter(num_trees = 100L)
#'   mask <- adapter(x, y, grid[1L, , drop = FALSE])
#'   cat("selected at sparsest threshold:", sum(mask), "of", ncol(x), "\n")
#' }
#' @export
make_ranger_adapter <- function(num_trees = 500L) {
  .require_pkg("ranger", "for base_learner = \"rf\"")
  num_trees <- .validate_scalar_integer_like(num_trees, "num_trees", min = 1L)
  batch_adapter <- .make_ranger_batch_adapter(num_trees = num_trees)

  function(x, y, lambda_val) {
    batch_adapter(x, y, lambda_val)[, 1L]
  }
}

#' Build a Data-Driven Importance-Threshold Grid for the RF Learner
#' (Experimental)
#'
#' Fits one `ranger` forest with permutation importance on the full training
#' data and returns the importance-threshold grid as a pre-expanded
#' `data.frame` suitable for the `lambda_grid` argument of [stabl_fit()] with
#' `base_learner = "rf"`.  This is the random-forest analog of
#' [auto_lambda_grid()].
#'
#' Thresholds are quantiles of the pooled permutation-importance distribution
#' at probabilities `seq(0.5, 0.999, length.out = n_grid)`, sorted in
#' decreasing order so the path runs from the sparse end (high threshold) to
#' the dense end (low threshold), mirroring the glmnet lambda-path convention.
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
#'   RF-specific meaning.
#'
#' @seealso [make_ranger_adapter()], [auto_lambda_grid()], [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(200L), 20L, 10L,
#'               dimnames = list(NULL, paste0("f", 1:10)))
#'   y <- x[, 1] * x[, 2] + rnorm(20L)
#'   grid <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L,
#'                                  seed = 1L)
#'   head(grid)
#' }
#' @export
auto_rf_threshold_grid <- function(x, y, n_grid = 30L, num_trees = 500L,
                                   seed = NULL) {
  .require_pkg("ranger", "for base_learner = \"rf\"")
  n_grid    <- .validate_scalar_integer_like(n_grid, "n_grid", min = 1L)
  num_trees <- .validate_scalar_integer_like(num_trees, "num_trees", min = 1L)
  seed      <- .validate_scalar_integer_like(seed, "seed", allow_null = TRUE)

  importance <- .rf_fit_permutation_importance(x, y, num_trees, seed)
  .rf_threshold_grid_from_importance(importance, n_grid)
}

# ---- Internal helpers ---------------------------------------------------------

# Fit one ranger forest and return its permutation importance as a named
# numeric vector of length ncol(x).  ranger requires column names on `x`
# ("No covariates found" otherwise), so unnamed matrices are named here.
.rf_fit_permutation_importance <- function(x, y, num_trees, seed) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  # On very small subsamples (a few dozen rows), ranger's permutation
  # importance can contain NaN for a few features — a numerical artifact of
  # degenerate out-of-bag splits, not a data problem.  Retry with offset
  # seeds (large-prime steps avoid colliding with other iterations' seeds)
  # before failing.
  n_attempts <- 3L
  for (attempt in seq_len(n_attempts)) {
    fit <- ranger::ranger(
      x                         = x,
      y                         = y,
      num.trees                 = num_trees,
      importance                = "permutation",
      respect.unordered.factors = TRUE,
      seed                      = seed + (attempt - 1L) * 104729L
    )
    importance <- fit$variable.importance
    if (!is.null(importance) && length(importance) == ncol(x) &&
        !anyNA(importance)) {
      return(importance)
    }
  }
  NULL
}

# Importance quantiles -> threshold grid, sparse end first (decreasing), the
# same direction as a glmnet lambda path.
.rf_threshold_grid_from_importance <- function(importance, n_grid) {
  probs      <- seq(0.5, 0.999, length.out = n_grid)
  thresholds <- sort(
    stats::quantile(importance, probs = probs, names = FALSE),
    decreasing = TRUE
  )
  data.frame(lambda = thresholds, threshold = thresholds)
}

# Resolve the importance thresholds from a grid row-set.  The `threshold`
# column is authoritative when present; the `lambda` column (required by
# stabl_fit's grid validation) is the fallback for user-supplied grids.
.rf_thresholds_from_grid <- function(lambda_grid) {
  col <- if ("threshold" %in% names(lambda_grid)) "threshold" else "lambda"
  thresholds <- as.numeric(lambda_grid[[col]])
  if (length(thresholds) == 0L || anyNA(thresholds) ||
      any(!is.finite(thresholds))) {
    stop(
      "The RF threshold grid must contain finite numeric `threshold` ",
      "(or `lambda`) values.",
      call. = FALSE
    )
  }
  thresholds
}

# Package warning for the statistically unsound rf + permutation-null pairing.
.warn_rf_random_permutation_nulls <- function() {
  warning(
    "base_learner = \"rf\" with artificial_type = \"random_permutation\": ",
    "a random forest is a multivariate learner, so real null features ",
    "correlated with true signals receive leaked importance while permuted ",
    "null features do not; the FDP+ estimate is therefore ",
    "anti-conservative. Pair the RF learner with model-X knockoff ",
    "artificial features (artificial_type = \"knockoff_equi\" or ",
    "\"knockoff_mvr\"), which preserve the correlation structure.",
    call. = FALSE
  )
}

# Audit A2 (0.2.2): confirmatory guard for the fixed-X knockoff pairing.
# Fixed-X knockoffs are generated once on the full design; their
# exchangeability does not survive the bootstrap row subsampling, so the
# FDP+ estimate for an RF learner paired with fixed-X knockoffs is not
# validated (documented-invalid regime).
.warn_rf_fixedx_knockoff_nulls <- function() {
  warning(
    "base_learner = \"rf\" with artificial_type = \"knockoff\" (fixed-X): ",
    "fixed-X knockoffs are generated once on the full design, so their ",
    "exchangeability does not survive the bootstrap row subsampling; the ",
    "FDP+ estimate for this pairing is not validated. Use a model-X ",
    "knockoff type instead (artificial_type = \"knockoff_equi\" or ",
    "\"knockoff_mvr\"), whose knockoff rows are i.i.d. draws from a fitted ",
    "Gaussian model and remain valid under row subsampling.",
    call. = FALSE
  )
}

# ---- Internal batch adapter ---------------------------------------------------
# Signature matches the other batch adapters:
#   function(x, y, lambda_grid) -> logical matrix (n_features x n_lambdas)
# ONE ranger forest is fitted per bootstrap subsample; masks for the whole
# grid are extracted from that single importance vector, so the bootstrap
# loop costs n_bootstraps forests total (not n_bootstraps x n_lambdas).
.make_ranger_batch_adapter <- function(family = "gaussian", num_trees = 500L) {
  num_trees <- .validate_scalar_integer_like(num_trees, "num_trees", min = 1L)
  function(x, y, lambda_grid) {
    # ranger infers the task from `y`: factor -> classification.  Coerce
    # 0/1-numeric or character class labels so binomial/multinomial fits are
    # classification forests, mirroring the glmnet family semantics.
    if (identical(family, "binomial") || identical(family, "multinomial")) {
      if (!is.factor(y)) y <- factor(y)
    }
    n_lambdas  <- nrow(lambda_grid)
    n_features <- ncol(x)
    thresholds <- .rf_thresholds_from_grid(lambda_grid)

    # ranger runs on its own C++ RNG and does not consume R's stream.  Draw
    # the forest seed from the caller's (already locally seeded) RNG stream —
    # stabl_fit() wraps every adapter call in .with_local_seed(iter_seeds[[i]],
    # ...) — so each bootstrap iteration is reproducible under `random_state`
    # and bit-identical across sequential and parallel execution.
    ranger_seed <- sample.int(.Machine$integer.max, 1L)

    importance <- .rf_fit_permutation_importance(x, y, num_trees, ranger_seed)
    if (is.null(importance) || length(importance) != n_features ||
        anyNA(importance)) {
      stop(
        "ranger did not return usable permutation importances for this ",
        "bootstrap subsample.",
        call. = FALSE
      )
    }

    masks <- vapply(
      thresholds,
      function(threshold) importance >= threshold,
      logical(n_features)
    )
    dimnames(masks) <- NULL
    masks
  }
}
