# Plot Selected-Feature Stability

Point plot of the stability score for each selected feature, ordered
along the y axis.

## Usage

``` r
# S3 method for class 'stabl_interpretation'
plot(x, ...)
```

## Arguments

- x:

  A `stabl_interpretation` object from
  [`stabl_interpret()`](https://mdmanurung.github.io/stablr/reference/stabl_interpret.md).

- ...:

  Ignored.

## Value

A `ggplot` object. Raises an error ("No selected features to plot.")
when the selection is empty.
