# Summarise a STABL Interpretation

Returns the evidence table restricted to the selected features, ordered
by decreasing stability score.

## Usage

``` r
# S3 method for class 'stabl_interpretation'
summary(object, ...)
```

## Arguments

- object:

  A `stabl_interpretation` object from
  [`stabl_interpret()`](https://mdmanurung.github.io/stablr/reference/stabl_interpret.md).

- ...:

  Ignored.

## Value

A `data.frame`: the rows of `object$evidence` with `selected == TRUE`,
ordered by decreasing `stability`. Zero-row when the selection is empty.
