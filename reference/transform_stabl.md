# Transform New Data to the Selected STABL Feature Set

Subsets a matrix or data frame to the features selected by a fitted
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
object. Columns are returned in fitted feature order, rows are
preserved, and the output remains two-dimensional even when no features
are selected.

## Usage

``` r
transform_stabl(object, x, new_hard_threshold = NULL)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- x:

  A matrix or `data.frame` with feature names in columns.

- new_hard_threshold:

  Numeric in `(0, 1]` or `NULL`. When supplied, overrides the fitted
  threshold for this transformation only.

## Value

A matrix when `x` is a matrix, or a `data.frame` when `x` is a
`data.frame`, containing selected columns in fitted feature order.

## See also

[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md),
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)

## Examples

``` r
set.seed(1L)
x <- matrix(
  rnorm(30 * 5), 30, 5,
  dimnames = list(paste0("s", 1:30), paste0("f", 1:5))
)
y <- setNames(rnorm(30), rownames(x))
fit <- stabl_fit(
  x, y, lambda_grid = data.frame(lambda = c(0.2, 0.1)),
  n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
  random_state = 1L
)
transformed <- transform_stabl(fit, x)
dim(transformed)
#> [1] 30  3
try(transform_stabl(fit, unname(x)))
#> Error : `x` must have non-empty column names.
```
