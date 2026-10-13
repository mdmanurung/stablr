# Build an Adaptive Lasso Learner Adapter

Returns a closure that computes feature-specific penalty weights from a
ridge initialisation on each bootstrap subsample, then fits a lasso
model with those `penalty.factor` values to obtain the selected-feature
mask.

## Usage

``` r
make_adaptive_lasso_adapter(
  family = "gaussian",
  gamma = 1,
  epsilon = 1e-06,
  bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
)
```

## Arguments

- family:

  Character; the `glmnet` response family, for example `"gaussian"`,
  `"binomial"`, `"multinomial"`, or `"cox"`.

- gamma:

  Positive numeric scalar; controls how sharply the weights
  down-penalise features with large ridge coefficients. Larger values
  make the penalty more selective. Default `1.0` matches the Python
  reference implementation.

- epsilon:

  Positive numeric scalar; added to the denominator of the weight to
  avoid division by zero. Default `1e-6`.

- bootstrap_threshold:

  Positive numeric; absolute-coefficient cutoff used to decide that a
  feature is selected in a given bootstrap. Features with
  `|coef| > bootstrap_threshold` are counted as selected.

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)`.

## Details

Adaptive lasso improves on standard lasso by assigning stronger
penalties to features with small initial coefficients (likely noise) and
weaker penalties to features with large initial coefficients (likely
signal). This asymmetric penalisation achieves the oracle property under
regularity conditions, selecting the true support more reliably than
plain lasso when signal features have moderate to large effect sizes.

Weights are defined as: \$\$w_j = 1 / (\|\hat\beta_j^{\mathrm{init}}\| +
\epsilon)^\gamma\$\$ where \\\hat\beta_j^{\mathrm{init}}\\ comes from a
ridge regression on the same bootstrap subsample. The `epsilon` floor
avoids division by zero for features with near-zero ridge coefficients.

This adapter factory is primarily an internal backend used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).
For end-to-end feature selection workflows, prefer
`base_learner = "adaptive_lasso"` in
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## See also

[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
[`make_sgl_adapter()`](https://mdmanurung.github.io/stablr/reference/make_sgl_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(2L)
  x <- matrix(rnorm(200L), 20L, 10L)
  y <- rnorm(20L)
  adapter <- make_adaptive_lasso_adapter(family = "gaussian", gamma = 1, epsilon = 1e-6)
  mask <- adapter(x, y, data.frame(lambda = 0.1))
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
#> selected: 7 of 10 features
```
