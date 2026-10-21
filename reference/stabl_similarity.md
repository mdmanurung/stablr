# Similarity Between Feature Sets

Computes the similarity between two sets of selected features, or
between all pairs of sets in a list, using one of three methods. Use
this to assess the reproducibility of STABL (or any other
feature-selection method) across cross-validation folds or repeated
runs.

## Usage

``` r
stabl_similarity(
  x,
  y = NULL,
  method = c("jaccard", "adjusted", "pearson"),
  n_universe = NULL,
  as_matrix = FALSE,
  remove_diag = FALSE
)
```

## Arguments

- x:

  Either a character/integer vector of selected feature identifiers
  (then `y` must also be supplied and the result is a single scalar), or
  a list of such vectors, one per STABL run (then `y` must be `NULL` and
  all pairwise similarities are computed).

- y:

  `NULL` (default) or a character/integer vector of selected feature
  identifiers. Supply `y` for scalar mode, or a list as `x` for pairwise
  mode — not both.

- method:

  Character; one of `"jaccard"` (default), `"adjusted"`, or `"pearson"`.

- n_universe:

  Integer; total number of candidate features in the universe (i.e. the
  number of columns in the original predictor matrix). Required for
  `method = "adjusted"` and `method = "pearson"`; ignored for
  `method = "jaccard"`.

- as_matrix:

  Logical; pairwise mode only. If `FALSE` (default), returns the \\N
  (N - 1) / 2\\ off-diagonal upper-triangle values as a flat vector in
  row-major order (matching the Python convention). If `TRUE`, returns
  the full symmetric \\N \times N\\ matrix.

- remove_diag:

  Logical; only used when `as_matrix = TRUE`. If `TRUE` the
  self-similarity diagonal is removed row-wise, yielding an \\N \times
  (N - 1)\\ matrix in which row \\i\\ contains the similarities of run
  \\i\\ with every other run.

## Value

In scalar mode, a numeric scalar. In pairwise mode with
`as_matrix = FALSE`, a numeric vector of length \\N (N - 1) / 2\\. In
pairwise mode with `as_matrix = TRUE`, a numeric matrix of dimension \\N
\times N\\ (or \\N \times (N - 1)\\ when `remove_diag = TRUE`).

The matrix diagonal holds the self-similarity of each set: exactly 1 for
`"jaccard"` and `"pearson"`, and the computed self-value for
`"adjusted"` (\\1 - k / d\\ for a set of size \\k \le d / 2\\), which is
below 1 by construction because the measure is chance-corrected.

## Details

Three methods are available:

- `"jaccard"` (default):

  Set-overlap fraction \\\|A \cap B\| / \|A \cup B\|\\. The most
  interpretable measure; it does not correct for the chance overlap
  expected between small or large sets. When both sets are empty the
  value is 0 (convention matching the Python implementation).

- `"adjusted"`:

  Chance-corrected overlap, analogous to Cohen's kappa for sets:
  \$\$S\_{\text{adj}}(A, B) = \frac{r - \mathbb{E}\[r\]}{\min(k_1,
  k_2) - \max(0, k_1 + k_2 - d)}\$\$ where \\r = \|A \cap B\|\\, \\k_i =
  \|A_i\|\\, \\d =\\ `n_universe`, and \\\mathbb{E}\[r\] = k_1 k_2 /
  d\\. Returns 0 when either set is empty or equals the full universe.

- `"pearson"`:

  Pearson-correlation-inspired chance correction (the phi coefficient):
  \$\$S\_{\text{Pearson}}(A, B) = \frac{r - k_i k_j / d}{d \cdot
  \upsilon_i \upsilon_j}\$\$ with \\\upsilon_i = \sqrt{\pi_i (1 -
  \pi_i)}\\ and \\\pi_i = k_i / d\\. Returns 1 when both sets are empty
  or both equal the universe; 0 when exactly one of them is.

## See also

[`stabl_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/stabl_similarity_summary.md)
for reducing pairwise values to a single location statistic.

## Examples

``` r
# Scalar mode: two sets
stabl_similarity(c("A", "B", "C"), c("B", "C", "D"))                # 0.5
#> [1] 0.5
stabl_similarity(c("f1", "f2", "f3"), c("f2", "f3", "f4"),
                 method = "adjusted", n_universe = 10L)
#> [1] 0.3666667

# Pairwise mode: a list of sets, one per run
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
stabl_similarity(sets)                                        # flat vector
#> [1] 0.5 0.5 0.2
stabl_similarity(sets, as_matrix = TRUE)                      # 3 x 3
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.5  0.5
#> [2,]  0.5  1.0  0.2
#> [3,]  0.5  0.2  1.0
stabl_similarity(sets, as_matrix = TRUE, remove_diag = TRUE)  # 3 x 2
#>      [,1] [,2]
#> [1,]  0.5  0.5
#> [2,]  0.5  0.2
#> [3,]  0.5  0.2
stabl_similarity(sets, method = "adjusted", n_universe = 10L)
#> [1] 0.36666667 0.36666667 0.03333333
```
