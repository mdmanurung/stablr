# Summarize a Fitted STABL Object (Beginner Style)

Returns and prints a structured summary with the selection, the FDP+
estimate, and the chosen threshold, in plain language.

## Usage

``` r
# S3 method for class 'stabl'
summary(object, ...)
```

## Arguments

- object:

  A `"stabl"` object returned by
  [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md).

- ...:

  Ignored; present for S3 `summary` method compatibility.

## Value

A `"summary.stabl_fit"` object (see
[`summary.stabl_fit()`](https://mdmanurung.github.io/stablr/reference/summary.stabl_fit.md)).

## See also

[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md),
[`print.stabl()`](https://mdmanurung.github.io/stablr/reference/print.stabl.md)
