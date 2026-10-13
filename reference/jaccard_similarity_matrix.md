# Pairwise Jaccard Matrix from a List of Feature Sets

Computes all N\\\times\\N pairwise Jaccard similarities in one call.
Useful for visualising the reproducibility of STABL across multiple
cross-validation folds or repeated runs: a boxplot of the off-diagonal
values summarises how consistently the same features are selected.

## Usage

``` r
jaccard_similarity_matrix(list_of_lists, remove_diag = TRUE)
```

## Arguments

- list_of_lists:

  A list of character/integer vectors, one per STABL run (e.g. one per
  cross-validation fold).

- remove_diag:

  Logical; if `TRUE` (default) the diagonal column (self-similarity = 1)
  is removed from the output.

## Value

Numeric matrix of dimension N\\\times\\N (or N\\\times\\(N-1) when
`remove_diag = TRUE`). Row/column order matches `list_of_lists`.

## Details

By default the self-similarity diagonal (always 1) is removed, so the
returned matrix has \\N \times (N - 1)\\ columns and each row contains
the \\N - 1\\ similarities of run \\i\\ with every other run.

## See also

[`jaccard_similarity()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity.md),
[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md),
[`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)

## Examples

``` r
sets <- list(c("A","B","C"), c("B","C","D"), c("A","C","E"))
jaccard_similarity_matrix(sets)           # 3 x 2 (diagonal removed)
#>      [,1] [,2]
#> [1,]  0.5  0.5
#> [2,]  0.5  0.2
#> [3,]  0.5  0.2
jaccard_similarity_matrix(sets, remove_diag = FALSE)  # 3 x 3
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.5  0.5
#> [2,]  0.5  1.0  0.2
#> [3,]  0.5  0.2  1.0
```
