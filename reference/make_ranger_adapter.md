# Build a Random-Forest Learner Adapter (Calibrated, Scoped)

Returns a closure that fits ONE `ranger` random forest on a single
bootstrap subsample and returns a logical selection mask over the
features, mirroring the role of
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md)
for tree-based learners.

## Usage

``` r
make_ranger_adapter(num_trees = 500L, family = "gaussian")
```

## Arguments

- num_trees:

  Positive integer; number of trees per forest (`ranger`'s `num.trees`).
  Default `500L`.

- family:

  Character scalar; response family passed to the batch adapter. Default
  `"gaussian"`.

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)` for use
inside
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
where `lambda_val` is a 1-row `data.frame` from the threshold grid (see
[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md)).

## Details

This adapter factory is primarily an internal backend used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
(`base_learner = "rf"`). The response family is inferred by `ranger`
from `y` (numeric -\> regression, factor -\> classification, `Surv` -\>
survival);
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
coerces `y` for the binomial/multinomial families. The forest seed is
drawn from the caller's RNG stream so that results are reproducible
under
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)'s
`random_state` and identical across sequential and parallel execution.

## Selection rule

Instead of a penalty path, the RF adapter sweeps an
**importance-threshold path**. The grid (see
[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md))
contains quantiles of the pooled permutation-importance distribution;
given a 1-row grid `lambda_val`, a feature is selected when its
permutation importance is `>= ` the row's threshold. Lower thresholds
therefore give supersets of the selections at higher thresholds — the
same sparse-to-dense nesting the glmnet lambda path provides, so the
stability-accumulation and FDP+ logic in
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
work unchanged.

## Statistical guardrails (read before using)

- **Permutation importance, not Gini.** Gini importance is biased
  towards high-cardinality features; permutation importance is the
  better-studied choice for null calibration (Boruta likewise builds on
  mean-decrease-accuracy).

- **Pair with model-X knockoff artificial features.** A random forest is
  a multivariate learner: real null features correlated with true
  signals receive leaked importance, and permuted null features
  (`artificial_type = "random_permutation"`) do **not** reproduce that
  leakage, so the FDP+ estimate is anti-conservative. Model-X knockoff
  nulls (`"knockoff_equi"` / `"knockoff_mvr"`) preserve the correlation
  structure and are the only statistically sound pairing. Combining
  `base_learner = "rf"` with `artificial_type = "random_permutation"` in
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  emits a package warning to this effect.

- **Calibrated within a measured scope (WS-2).** With knockoff_equi
  decoys the mean empirical FDP is 0.025 (linear design) and 0.065
  (low-dimensional nonlinear design), meeting the predeclared criterion;
  it is NOT met at p = 50 with strong correlation (0.44) or for binary
  outcomes (0.26). Simulation evidence, not a guarantee.

## See also

[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md),
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(200L), 20L, 10L,
              dimnames = list(NULL, paste0("f", 1:10)))
  y <- x[, 1] * x[, 2] + rnorm(20L)
  grid  <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L,
                                  seed = 1L)
  adapter <- stablr:::make_ranger_adapter(num_trees = 100L)
  mask <- adapter(x, y, grid[1L, , drop = FALSE])
  cat("selected at sparsest threshold:", sum(mask), "of", ncol(x), "\n")
}
#> selected at sparsest threshold: 0 of 10 
```
