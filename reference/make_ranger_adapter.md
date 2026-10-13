# Build a Random-Forest Learner Adapter (Experimental)

Experimental learner adapter extending STABL to non-linear signals with
random forests (`ranger`). One forest is fit per bootstrap subsample;
the "penalty path" is an importance-threshold path built from quantiles
of the pooled permutation-importance distribution (see
[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md)),
so the selection set grows monotonically as the threshold decreases —
the same sparse-to-dense nesting the lambda path provides for glmnet
learners.

## Usage

``` r
make_ranger_adapter(num_trees = 500L)
```

## Arguments

- num_trees:

  Positive integer; number of trees per ranger forest.

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)`,
mirroring
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md):
the mask selects features whose permutation importance meets the
threshold carried by `lambda_val`.

## Details

**Experimental:** this is a stablr extension beyond Python STABL, and
FDP+ calibration for random-forest learners is empirical (validated
separately by the F5 benchmark), not a theoretical guarantee.

**Statistical guardrails:** permutation importance is used (not Gini).
Random forests must be paired with model-X knockoff artificial features
(`knockoff_equi` / `knockoff_mvr`): permuted nulls do not receive the
correlation leakage that real null features receive in a multivariate
learner, so `random_permutation` nulls make FDP+ anti-conservative —
stablr warns when this combination is requested.

## See also

[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md),
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
(`base_learner = "rf"`).

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(x[, 1] * x[, 2] * 0.5 + rnorm(60, sd = 0.3), rownames(x))
  adapter <- make_ranger_adapter(num_trees = 100L)
  grid <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  mask <- adapter(x, y, grid[1L, ])  # sparsest end of the path
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
#> selected: 0 of 6 features
```
