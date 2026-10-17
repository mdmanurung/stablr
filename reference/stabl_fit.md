# Fit STABL (Stability-Penalized Feature Selection)

Performs the STABL bootstrap stability-selection procedure over a grid
of regularisation parameters, optionally injecting artificial features
for data-adaptive FDP+ threshold diagnostics. The fitting rule chooses
the FDP+ minimizer; it is not a universal false-discovery guarantee.
This is the R counterpart of `Stabl.fit()` in the Python STABL library.

## Usage

``` r
stabl_fit(
  x,
  y,
  lambda_grid,
  base_learner = "lasso",
  family = "gaussian",
  n_bootstraps = 1000L,
  artificial_type = "random_permutation",
  artificial_proportion = 1,
  sample_fraction = 0.5,
  replace = FALSE,
  hard_threshold = NULL,
  fdr_threshold_range = seq(0, 0.99, by = 0.01),
  explore = FALSE,
  n_explore = 5L,
  sample_groups = NULL,
  groups = NULL,
  stratify_bootstrap = FALSE,
  bootstrap_strata = NULL,
  n_lambda = 30L,
  l1_ratio = NULL,
  verbose = FALSE,
  workers = 1L,
  random_state = NULL,
  adaptive_gamma = 1,
  adaptive_epsilon = 1e-06,
  feature_groups = NULL,
  corr_group_threshold = NULL,
  knockoff_fallback = c("stop", "allow"),
  splsda_ncomp = 2L
)

# S3 method for class 'stabl_fit'
print(x, ...)
```

## Arguments

- x:

  A numeric matrix or `data.frame` (samples \\\times\\ features) with
  row names used as sample IDs.

- y:

  A named numeric/factor vector whose names are sample IDs, or a
  matrix-like outcome (for example
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html)) with
  row names as sample IDs.

- lambda_grid:

  A pre-expanded `data.frame` (each row = one parameter combination)
  with at least a `lambda` column. Pass `"auto"` to derive a data-driven
  grid via
  [`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md)
  (uses `family` and `n_lambda`).

- base_learner:

  Character; `"lasso"` (alpha = 1), `"elastic_net"` (reads `alpha` from
  the `lambda_grid` `alpha` column), or `"adaptive_lasso"`, or
  `"sparse_group_lasso"`, or `"rf"` (random forest via `ranger`;
  requires `lambda_grid` from
  [`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md)
  or `"auto"`), or `"oblique_rf"` (oblique random forest via `aorsf`;
  requires `lambda_grid` from
  [`auto_oblique_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_oblique_threshold_grid.md)
  or `"auto"`), or `"ordinal"` (elastic-net cumulative logit via
  `ordinalNet` for ordered-band outcomes; requires a factor `y` and
  [`auto_ordinal_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_ordinal_lambda_grid.md)
  or `"auto"`). Default: `"lasso"`. `"rf"`, `"oblique_rf"`, and
  `"ordinal"` are experimental extensions beyond Python STABL: FDP+
  calibration for (oblique) random forests is empirical, not guaranteed,
  and is validated separately (see the package's RF validation study);
  they must be paired with model-X knockoff artificial features. For
  `"oblique_rf"` the per-feature selection frequencies are noisier than
  for `"rf"` because importance is spread across the linear combinations
  used in each split, or `"splsda"` (sparse PLS discriminant analysis
  via `mixOmics`; classification only, requires `family = "binomial"`;
  sweeps a per-component `keepX` grid — use
  [`auto_splsda_keepx_grid()`](https://mdmanurung.github.io/stablr/reference/auto_splsda_keepx_grid.md)
  or `"auto"`, and control the number of latent components with
  `splsda_ncomp`). `"splsda"` is an experimental extension beyond Python
  STABL whose FDP+ calibration is empirical, not guaranteed. Unlike the
  RF learners, the random-permutation pairing is acceptable for
  `"splsda"`: on the package's multi-cohort validation benchmark,
  sPLS-DA + random permutation was the best-calibrated configuration
  (empirical FDP / reported FDP+ ratio 1.01), while the `knockoff_equi`
  pairing ran conservative there; validate knockoff pairings per dataset
  (see the package's MINT comparison study).

- family:

  Character; `glmnet` response family (`"gaussian"`, `"binomial"`,
  `"multinomial"`, or `"cox"`), or `"ordinal"` for
  `base_learner = "ordinal"`. `base_learner = "splsda"` requires
  `"binomial"`. Default: `"gaussian"`.

- n_bootstraps:

  Positive integer; bootstrap iterations per lambda. Default: `1000L`.

- artificial_type:

  Character or `NULL`; `"random_permutation"`, `"knockoff"`,
  `"knockoff_equi"`, `"knockoff_mvr"`, or `NULL` (no artificial features
  — requires `hard_threshold`). Default: `"random_permutation"`.

- artificial_proportion:

  Numeric in `(0, 1]`; fraction of original features to inject as
  artificial noise. Default: `1.0`.

- sample_fraction:

  Positive numeric; fraction of samples drawn per bootstrap. Default:
  `0.5`.

- replace:

  Logical; sample with replacement? Default: `FALSE`.

- hard_threshold:

  Numeric in `(0, 1]` or `NULL`. When supplied, FDP+ control is bypassed
  and this value is the stability-score cut-off. Default: `NULL`.

- fdr_threshold_range:

  Numeric vector swept when computing FDP+. Default:
  `seq(0, 0.99, by = 0.01)`, matching Python STABL's
  `np.arange(0., 1., .01)`.

- explore:

  Logical; if `TRUE` and no features pass the threshold, fall back to
  the top `n_explore` features. Default: `FALSE`.

- n_explore:

  Positive integer; fallback feature count. Default: `5L`.

- sample_groups:

  Named vector of group IDs (same names as `rownames(x)`) or `NULL`.
  When supplied,
  [`group_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/group_bootstrap_indices.md)
  is used instead of
  [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md)
  – use this when samples are clustered (e.g., repeated measures,
  multiple sites) and bootstrap resampling must keep clusters intact.
  Default: `NULL`.

- groups:

  Deprecated; use `sample_groups` instead.

- stratify_bootstrap:

  Logical; if `TRUE`, bootstrap subsamples are stratified by the outcome
  class. This is intended for classification tasks with small or
  imbalanced classes. Default `FALSE` preserves the original STABL
  sampling behavior.

- bootstrap_strata:

  Optional categorical stratification design for bootstrap sampling.
  Provide a named vector, matrix, `data.frame`, or list with one
  row/value per sample. Multiple columns are combined as a joint
  interaction stratum, for example outcome class by study group. When
  supplied, this overrides `stratify_bootstrap`.

- n_lambda:

  Integer; number of lambda values when `lambda_grid = "auto"`. Ignored
  otherwise. Default: `30L`.

- l1_ratio:

  Numeric scalar, numeric vector, or `NULL`. Passed to
  [`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md)
  when `lambda_grid = "auto"`. Use this for
  `base_learner = "elastic_net"` so the generated grid contains the
  elastic-net `alpha` column. Default `NULL` preserves the lasso-like
  auto grid used by existing workflows.

- verbose:

  Logical; emit progress messages. Default: `FALSE`.

- workers:

  Positive integer; parallel workers for
  [`furrr::future_map()`](https://furrr.futureverse.org/reference/future_map.html).
  Actual parallelism requires the caller to invoke
  `future::plan(multisession, workers = N)` first. Default: `1L`.

- random_state:

  Integer or `NULL`; top-level seed. Default: `NULL`.

- adaptive_gamma:

  Positive numeric scalar. Used only when
  `base_learner = "adaptive_lasso"`.

- adaptive_epsilon:

  Positive numeric scalar. Used only when
  `base_learner = "adaptive_lasso"`.

- feature_groups:

  Optional feature-group definition for
  `base_learner = "sparse_group_lasso"`. Provide a vector/factor of
  length `ncol(x)` assigning each original feature to a group.

- corr_group_threshold:

  Optional numeric percentile in `(0, 100]` used to derive sparse-group
  feature groups from pairwise correlations of original features. Used
  only when `base_learner = "sparse_group_lasso"`.

- knockoff_fallback:

  Character; policy when a knockoff constructor fails on a
  supported-dimensional input (audit A1). `"stop"` (the default) raises
  an error because a knockoff-to-random-permutation fallback silently
  changes the FDP+ calibration target; `"allow"` restores the pre-0.2.2
  warn-and-fall-back behavior (recorded in the artificial provenance's
  `fallback_history`).

- splsda_ncomp:

  Positive integer; number of latent components retained in every
  [`mixOmics::splsda`](https://rdrr.io/pkg/mixOmics/man/splsda.html)
  fit. Used only when `base_learner = "splsda"`. Default: `2L`.

- ...:

  Ignored; present for S3 `print` method compatibility.

## Value

An S3 object of class `"stabl_fit"` containing:

- `stabl_scores_`:

  Numeric matrix (features \\\times\\ lambdas) of per-bootstrap
  selection frequencies (stability scores).

- `stabl_scores_artificial_`:

  Analogous matrix for artificial features, or `NULL` when
  `artificial_type = NULL`.

- `fitted_lambda_grid`:

  The `data.frame` of lambda combinations actually used.

- `fdr_min_threshold_`:

  FDP+-optimal stability threshold, or `NULL`.

- `FDRs_`:

  Numeric vector of FDP+ estimates per threshold, or `NULL`.

- `min_fdr_`:

  Minimum FDP+ achieved, or `NULL`.

- `fdrs_table_`:

  Per-lambda FDP+ matrix, or `NULL`.

- `hard_threshold`:

  As supplied.

- `artificial_type`:

  As supplied.

- `artificial_provenance`:

  List returned by
  [`make_artificial_features()`](https://mdmanurung.github.io/stablr/reference/make_artificial_features.md)
  recording `actual_type`, all selected generator modes, and the
  complete ordered fallback history. This is the scientific source of
  truth for the generator used; or `NULL` when `artificial_type = NULL`.

- `artificial_proportion`:

  As supplied.

- `explore`:

  As supplied.

- `n_explore`:

  As supplied.

- `feature_names`:

  Character vector of original feature names.

- `n_features_in_`:

  Integer number of original features.

## Methods (by generic)

- `print(stabl_fit)`: Print a concise summary of a fitted `stabl_fit`
  object; invisibly returns `x`.

## Examples

``` r
set.seed(42)
n <- 50L
p <- 10L
x <- matrix(
  rnorm(n * p),
  nrow = n,
  dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p)))
)
y <- setNames(0.8 * x[, 1] - 0.5 * x[, 2] + rnorm(n, sd = 0.3), rownames(x))
lambda_grid <- data.frame(lambda = seq(0.02, 0.30, length.out = 4L))

# Default lasso path
fit_lasso <- stabl_fit(
  x = x,
  y = y,
  lambda_grid = lambda_grid,
  n_bootstraps = 6L,
  artificial_type = "random_permutation",
  random_state = 1L
)
get_feature_names(fit_lasso)
#> [1] "f1" "f2"

# Adaptive lasso path
fit_adaptive <- stabl_fit(
  x = x,
  y = y,
  lambda_grid = lambda_grid,
  base_learner = "adaptive_lasso",
  adaptive_gamma = 1.0,
  adaptive_epsilon = 1e-6,
  n_bootstraps = 6L,
  artificial_type = "random_permutation",
  random_state = 2L
)
get_importances(fit_adaptive)
#>        f1        f2        f3        f4        f5        f6        f7        f8 
#> 1.0000000 1.0000000 0.5000000 0.8333333 0.8333333 0.3333333 0.6666667 0.6666667 
#>        f9       f10 
#> 0.3333333 0.5000000 

# Sparse-group lasso path (optional dependency)
if (requireNamespace("sparsegl", quietly = TRUE)) {
  feature_groups <- rep(seq_len(5L), each = 2L)
  fit_sgl <- stabl_fit(
    x = x,
    y = y,
    lambda_grid = lambda_grid,
    base_learner = "sparse_group_lasso",
    feature_groups = feature_groups,
    n_bootstraps = 5L,
    artificial_type = "random_permutation",
    random_state = 3L
  )
  get_support(fit_sgl)
}
#>    f1    f2    f3    f4    f5    f6    f7    f8    f9   f10 
#>  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 

# Multinomial path with three classes
y_multi <- setNames(
  factor(sample(c("A", "B", "C"), n, replace = TRUE)),
  rownames(x)
)
fit_multi <- suppressWarnings(stabl_fit(
  x = x,
  y = y_multi,
  lambda_grid = lambda_grid,
  family = "multinomial",
  n_bootstraps = 5L,
  artificial_type = "random_permutation",
  random_state = 4L
))
get_feature_names(fit_multi)
#>  [1] "f1"  "f2"  "f3"  "f4"  "f5"  "f6"  "f7"  "f8"  "f9"  "f10"
```
