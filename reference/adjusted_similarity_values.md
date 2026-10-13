# Deprecated: use [`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md)

`adjusted_similarity_values()` is deprecated and will be removed in a
future release; use
[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md)
instead.

## Usage

``` r
adjusted_similarity_values(list_of_lists, nb_total_elements)
```

## Examples

``` r
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
# warns; use adjusted_similarity_pairwise()
same <- suppressWarnings(identical(
  adjusted_similarity_values(sets, nb_total_elements = 10L),
  adjusted_similarity_pairwise(sets, nb_total_elements = 10L)
))
stopifnot(same)
```
