# Pearson-Corrected Similarity Between Two Feature Sets

Computes a Pearson-correlation-inspired similarity that corrects for
expected random intersection. This is a second chance-correction
approach (alongside
[`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md))
that normalises by the geometric mean of the within-set variances under
independent Bernoulli sampling.

## Usage

``` r
pearson_similarity(
  set1 = NULL,
  set2 = NULL,
  n_universe = NULL,
  list_i = NULL,
  list_j = NULL,
  d = NULL
)
```

## Arguments

- set1:

  Character or integer vector of selected feature identifiers.

- set2:

  Character or integer vector of selected feature identifiers.

- n_universe:

  Integer; total number of candidate features in the universe.

- list_i:

  Deprecated; use `set1` instead. If supplied by name, a deprecation
  warning is emitted and the value is used as `set1` when the latter is
  `NULL`. It will be removed in a future release.

- list_j:

  Deprecated; use `set2` instead. If supplied by name, a deprecation
  warning is emitted and the value is used as `set2` when the latter is
  `NULL`. It will be removed in a future release.

- d:

  Deprecated; use `n_universe` instead. If supplied by name, a
  deprecation warning is emitted and the value is used as `n_universe`
  when the latter is `NULL`. It will be removed in a future release.

## Value

Numeric scalar. Positive values indicate more overlap than chance; the
maximum is typically close to 1 for perfectly matching sets.

## Details

The formula is: \$\$S\_{\text{Pearson}}(A, B) = \frac{r - k_i k_j / d}{d
\cdot \upsilon_i \upsilon_j}\$\$ where \\r = \|A \cap B\|\\, \\k_i =
\|A_i\|\\, \\d\\ is the universe size, and \\\upsilon_i = \sqrt{\pi_i
(1 - \pi_i)}\\ with \\\pi_i = k_i / d\\.

Edge cases: returns 1 when both sets are empty or both equal the
universe; returns 0 when one set is empty or equals the universe.

## See also

[`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)
for computing all pairwise values at once,
[`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)
for a summary statistic,
[`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md)
for an alternative chance-correction.

## Examples

``` r
pearson_similarity(c("f1","f2","f3"), c("f2","f3","f4"), n_universe = 10L)
#> [1] 0.5238095
pearson_similarity(c("f1","f2"),      c("f1","f2"),      n_universe = 10L) # near 1
#> [1] 1
pearson_similarity(character(0),      c("f1","f2"),      n_universe = 10L) # 0
#> [1] 0
```
