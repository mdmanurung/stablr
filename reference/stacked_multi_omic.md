# Deprecated: use [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)

`stacked_multi_omic()` is deprecated and will be removed in a future
release; use
[`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)
instead. (The Python parity fixtures keep the historical
`stacked_multi_omic` name on the Python side; the R entry point is
[`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md).)

## Usage

``` r
stacked_multi_omic(...)
```

## Arguments

- ...:

  See
  [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md).

## Examples

``` r
set.seed(1L)
y <- rep(c(0L, 1L), 10)
pred <- data.frame(
  omic_a = y + rnorm(20, sd = 0.5),
  omic_b = y + rnorm(20, sd = 1.5)
)
# warns; use stabl_multiomic_stacked()
old <- suppressWarnings(stacked_multi_omic(pred, y, task_type = "binary",
                                           n_iter = 200L, random_state = 42L))
new <- stabl_multiomic_stacked(pred, y, task_type = "binary",
                               n_iter = 200L, random_state = 42L)
stopifnot(identical(names(old), names(new)))
```
