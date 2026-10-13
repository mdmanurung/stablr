# Deprecated: use [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md)

`classic_bootstrap_indices()` is deprecated and will be removed in a
future release; use
[`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md)
instead.

## Usage

``` r
classic_bootstrap_indices(...)
```

## Arguments

- ...:

  See
  [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md).

## Examples

``` r
y <- c(rep(0, 10), rep(1, 10))
# warns; use sample_bootstrap_indices()
same <- suppressWarnings(identical(
  classic_bootstrap_indices(y, n_subsamples = 5L, random_state = 1L),
  sample_bootstrap_indices(y, n_subsamples = 5L, random_state = 1L)
))
stopifnot(same)
```
