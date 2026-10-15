# Apply a fitted preprocessing transform

Applies a transform fitted by
[`stabl_preprocess_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit.md)
to new data. The fitted parameters are used unchanged; assessment values
never influence them.

## Usage

``` r
stabl_preprocess_apply(object, x, ...)
```

## Arguments

- object:

  A `"stabl_preprocessing"` object.

- x:

  Numeric matrix (samples × features) with column names. Must contain
  every retained feature; extra columns are dropped.

- ...:

  Reserved for future use.

## Value

A numeric matrix with the retained features, transformed.

## Examples

``` r
x <- matrix(rnorm(40 * 5), 40, 5,
            dimnames = list(paste0("s", 1:40), paste0("f", 1:5)))
pp <- stabl_preprocess_fit(x, transforms = c("center", "scale"))
xt <- stabl_preprocess_apply(pp, x)
```
