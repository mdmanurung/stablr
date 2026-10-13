# Support and Feature Names of a STABL Selection

Beginner-friendly generics for the selection: `support()` returns the
selected features' names or mask (delegating to
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md))
and `feature_names()` returns the selected feature names at the chosen
threshold (delegating to
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md)).
Methods are provided for `stabl`/`stabl_fit` objects.

## Usage

``` r
support(object, ...)

feature_names(object, ...)

# S3 method for class 'stabl_fit'
support(object, ...)

# S3 method for class 'stabl'
support(object, ...)

# S3 method for class 'stabl_fit'
feature_names(object, ...)

# S3 method for class 'stabl'
feature_names(object, ...)
```

## Arguments

- object:

  A fitted `stabl`/`stabl_fit` object.

- ...:

  Passed through (e.g. `new_hard_threshold`).

## Value

See
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
/
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md).

## See also

[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md),
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md),
[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)

## Examples

``` r
# \donttest{
set.seed(1L)
df <- data.frame(m1 = rnorm(40), m2 = rnorm(40), m3 = rnorm(40))
df$y <- 0.8 * df$m1 + rnorm(40, sd = 0.5)
fit <- stabl(df, y = "y", n_bootstraps = 50L, random_state = 1L)
feature_names(fit)
#> character(0)
# }
```
