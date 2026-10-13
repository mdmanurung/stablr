# Upper-Triangle Pearson Similarity Values

Computes all pairwise Pearson-corrected similarities and returns the
N\\\times\\(N-1)/2 upper-triangle values as a flat vector (same layout
as
[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md),
enabling direct comparison between the two metrics on the same data).

## Usage

``` r
pearson_similarity_pairwise(list_of_lists, nb_total_elements, d = NULL)
```

## Arguments

- list_of_lists:

  A list of character/integer vectors, one per run.

- nb_total_elements:

  Integer; total number of candidate features.

- d:

  Deprecated; use `nb_total_elements` instead.

## Value

Numeric vector of length N\*(N-1)/2.

## See also

[`pearson_similarity()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity.md)
for the pairwise function,
[`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)
for a one-number summary,
[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md)
for an alternative chance-correction.

## Examples

``` r
sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
pearson_similarity_pairwise(sets, nb_total_elements = 10L)
#> [1] 0.52380952 0.52380952 0.04761905
```
