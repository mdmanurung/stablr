# Deprecated: use [`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)

`pearson_similarity_measure()` is deprecated and will be removed in a
future release; use
[`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)
instead.

## Usage

``` r
pearson_similarity_measure(list_of_lists, d, stat = "median")
```

## Examples

``` r
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
# warns; use pearson_similarity_summary()
same <- suppressWarnings(identical(
  pearson_similarity_measure(sets, d = 10L),
  pearson_similarity_summary(sets, nb_total_elements = 10L)
))
stopifnot(same)
```
