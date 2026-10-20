# Build a Sparse Group Lasso Learner Adapter

Returns a closure that fits
[`sparsegl::sparsegl()`](https://dajmcdon.github.io/sparsegl/reference/sparsegl.html)
on a single bootstrap subsample and returns a logical selection mask
over the features.

## Usage

``` r
make_sgl_adapter(
  family = "gaussian",
  feature_groups,
  alpha_fixed = NULL,
  bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
)
```

## Arguments

- family:

  Character; response family (`"gaussian"`, `"binomial"`, or
  `"multinomial"`). Cox regression is not supported by `sparsegl`.

- feature_groups:

  Integer or factor vector of length `p` (number of features) assigning
  each feature to a group. Values are coerced via
  `as.integer(as.factor(...))`, so any type that has a natural ordering
  is accepted.

- alpha_fixed:

  Numeric scalar in `[0, 1]` or `NULL`. When not `NULL`, overrides any
  `alpha` column in `lambda_val`. Controls the lasso / group-lasso
  mixing weight (`asparse` in `sparsegl`).

- bootstrap_threshold:

  Positive numeric; absolute-coefficient cutoff used to decide that a
  feature is selected in a given bootstrap.

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)`.

## Details

Sparse-group lasso is useful when features have a known (or inferred)
block structure, for example gene pathways, omic layers, or correlated
feature clusters. It imposes simultaneous sparsity within and between
groups: the group-level penalty encourages whole groups to be zeroed
out, while the within-group lasso penalty allows groups to have only a
sparse subset of active features. This can substantially improve
stability in structured high-dimensional settings.

The `alpha` value (`asparse` in `sparsegl`) controls the balance between
the within-group lasso penalty and the group-level penalty: 0 = pure
group lasso; 1 = pure lasso (no group penalty). The default of `0.05`
used when no `alpha` column is present matches the Python STABL
reference.

This adapter factory is primarily an internal backend used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).
For end-to-end feature selection workflows, prefer
`base_learner = "sparse_group_lasso"` in
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
and provide `feature_groups` (or `corr_group_threshold`) there.

## See also

[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
[`make_adaptive_lasso_adapter()`](https://mdmanurung.github.io/stablr/reference/make_adaptive_lasso_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("sparsegl", quietly = TRUE)) {
  set.seed(3L)
  x <- matrix(rnorm(200L), 20L, 10L)
  y <- rnorm(20L)
  groups  <- rep(1:5, each = 2L)
  adapter <- stablr:::make_sgl_adapter(family = "gaussian", feature_groups = groups)
  mask <- adapter(x, y, data.frame(lambda = 0.1))
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
#> selected: 0 of 10 features
```
