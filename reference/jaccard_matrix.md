# Deprecated: use [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)

`jaccard_matrix()` is deprecated and will be removed in a future
release; use
[`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)
instead.

## Usage

``` r
jaccard_matrix(list_of_lists, remove_diag = TRUE)
```

## Arguments

- list_of_lists:

  See
  [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md).

## Examples

``` r
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
# warns; use jaccard_similarity_matrix()
same <- suppressWarnings(identical(
  jaccard_matrix(sets),
  jaccard_similarity_matrix(sets)
))
stopifnot(same)
```
