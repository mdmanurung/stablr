# Build a glmnet Learner Adapter

Returns a closure that fits a `glmnet` model on a single bootstrap
subsample and returns a logical selection mask over the features. The
closure mirrors the role of `fit_bootstrapped_sample()` in the Python
STABL library.

## Usage

``` r
make_glmnet_adapter(
  family = "gaussian",
  alpha_fixed = NULL,
  bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
)
```

## Arguments

- family:

  Character; the `glmnet` response family, for example `"gaussian"`,
  `"binomial"`, `"multinomial"`, or `"cox"`.

- alpha_fixed:

  Numeric scalar or `NULL`. When not `NULL`, this value overrides any
  `alpha` column in `lambda_val`.

- bootstrap_threshold:

  Positive numeric; absolute-coefficient cutoff used to decide that a
  feature is selected in a given bootstrap. Features with
  `|coef| > bootstrap_threshold` are counted as selected.

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)` for use
inside
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Details

Learner adapters decouple the modelling back-end from the STABL
bootstrap loop, making it easy to substitute different
penalised-regression solvers without changing the stability-accumulation
logic. This factory produces the standard lasso / elastic-net adapter
that is the default back-end for
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

The returned function accepts a pre-expanded lambda-grid row (a 1-row
`data.frame`) and applies `glmnet` at that exact penalty value. The
`alpha` elastic-net mixing parameter is resolved in priority order:
`alpha_fixed` argument \> `alpha` column in `lambda_val` \> default of 1
(pure lasso).

This adapter factory is primarily an internal backend used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).
Most users should configure `base_learner`, `family`, and `lambda_grid`
directly in
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
rather than calling this factory manually.

## See also

[`make_adaptive_lasso_adapter()`](https://mdmanurung.github.io/stablr/reference/make_adaptive_lasso_adapter.md),
[`make_sgl_adapter()`](https://mdmanurung.github.io/stablr/reference/make_sgl_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(200L), 20L, 10L)
  y <- rnorm(20L)
  adapter <- make_glmnet_adapter(family = "gaussian", bootstrap_threshold = 1e-5)
  mask <- adapter(x, y, data.frame(lambda = 0.1))
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
#> selected: 5 of 10 features
```
