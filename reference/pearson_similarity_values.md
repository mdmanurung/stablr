# Deprecated: use [`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)

`pearson_similarity_values()` is deprecated and will be removed in a
future release; use
[`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)
instead.

## Usage

``` r
pearson_similarity_values(list_of_lists, d)
```

## Examples

``` r
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
# warns; use pearson_similarity_pairwise()
same <- suppressWarnings(identical(
  pearson_similarity_values(sets, d = 10L),
  pearson_similarity_pairwise(sets, nb_total_elements = 10L)
))
stopifnot(same)
```
