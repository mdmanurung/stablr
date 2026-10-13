# Print a Fitted STABL Object (Beginner Style)

Prints a plain-language answer: what was fitted, which features were
selected, and how reliable the selection is estimated to be. The
underlying terse `stabl_fit` print remains available by calling
[`print.stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
explicitly or via
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
objects.

## Usage

``` r
# S3 method for class 'stabl'
print(x, ...)
```

## Arguments

- x:

  A `"stabl"` object returned by
  [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md).

- ...:

  Ignored; present for S3 `print` method compatibility.

## Value

The object, invisibly.

## See also

[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md),
[`summary.stabl()`](https://mdmanurung.github.io/stablr/reference/summary.stabl.md),
[`plot.stabl()`](https://mdmanurung.github.io/stablr/reference/plot.stabl.md)
