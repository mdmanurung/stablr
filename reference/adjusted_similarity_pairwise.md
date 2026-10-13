# Upper-Triangle Adjusted Similarity Values

Computes the full pairwise adjusted-similarity matrix for N feature sets
and returns only the N\\\times\\(N-1)/2 off-diagonal upper-triangle
values as a flat vector. This format is convenient for computing summary
statistics (see
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md))
or for Wilcoxon/permutation tests comparing reproducibility across
methods.

## Usage

``` r
adjusted_similarity_pairwise(list_of_lists, nb_total_elements)
```

## Arguments

- list_of_lists:

  A list of character/integer vectors, one per STABL run.

- nb_total_elements:

  Integer; total number of candidate features.

## Value

Numeric vector of length N\*(N-1)/2 containing the pairwise
adjusted-similarity values for all unique pairs (row-major
upper-triangle order, matching the Python convention).

## See also

[`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md)
for the pairwise function,
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
for a one-number summary,
[`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)
for an alternative chance-correction.

## Examples

``` r
sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
adjusted_similarity_pairwise(sets, nb_total_elements = 10L)
#> [1] 0.36666667 0.36666667 0.03333333
```
