# Deprecated: use [`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)

`adjusted_similarity_measure()` is deprecated and will be removed in a
future release; use
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
instead.

## Usage

``` r
adjusted_similarity_measure(list_of_lists, nb_total_elements, stat = "median")
```

## Arguments

- list_of_lists, nb_total_elements, stat:

  See
  [`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md).

## Examples

``` r
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
# warns; use adjusted_similarity_summary()
same <- suppressWarnings(identical(
  adjusted_similarity_measure(sets, nb_total_elements = 10L),
  adjusted_similarity_summary(sets, nb_total_elements = 10L)
))
stopifnot(same)
```
