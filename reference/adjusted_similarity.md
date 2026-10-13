# Adjusted Similarity Between Two Feature Sets

Computes a chance-corrected similarity between two feature sets. Unlike
[`jaccard_similarity()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity.md),
the adjusted measure accounts for the expected random overlap given the
sizes of both sets and the total feature universe, so it does not
systematically penalise methods that select many features.

## Usage

``` r
adjusted_similarity(
  set1 = NULL,
  set2 = NULL,
  n_universe = NULL,
  list1 = NULL,
  list2 = NULL,
  nb_total_elements = NULL
)
```

## Arguments

- set1:

  Character or integer vector of selected feature identifiers.

- set2:

  Character or integer vector of selected feature identifiers.

- n_universe:

  Integer; total number of candidate features in the universe (i.e. the
  number of columns in the original predictor matrix).

- list1:

  Deprecated; use `set1` instead. If supplied by name, a deprecation
  warning is emitted and the value is used as `set1` when the latter is
  `NULL`. It will be removed in a future release.

- list2:

  Deprecated; use `set2` instead. If supplied by name, a deprecation
  warning is emitted and the value is used as `set2` when the latter is
  `NULL`. It will be removed in a future release.

- nb_total_elements:

  Deprecated; use `n_universe` instead. If supplied by name, a
  deprecation warning is emitted and the value is used as `n_universe`
  when the latter is `NULL`. It will be removed in a future release.

## Value

Numeric scalar in \\(-1, 1\]\\. Values above 0 indicate more overlap
than expected by chance; 1 means perfect agreement; negative values
indicate less overlap than chance.

## Details

The formula is analogous to Cohen's kappa for sets:
\$\$S\_{\text{adj}}(A, B) = \frac{r - \mathbb{E}\[r\]}{\min(k_1, k_2) -
\max(0, k_1 + k_2 - d)}\$\$ where \\r = \|A \cap B\|\\, \\k_i =
\|A_i\|\\, \\d =\\ `n_universe`, and \\\mathbb{E}\[r\] = k_1 k_2 / d\\.

Returns 0 when either set is empty or equals the full universe (edge
cases where the correction denominator is zero).

## See also

[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md)
for computing all pairwise values at once,
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
for a summary statistic,
[`jaccard_similarity()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity.md)
for the non-chance-corrected alternative.

## Examples

``` r
# 10-feature universe; sets A and B share 2 out of 3 features each
adjusted_similarity(c("f1","f2","f3"), c("f2","f3","f4"), n_universe = 10L)
#> [1] 0.3666667
adjusted_similarity(c("f1","f2"),      c("f1","f2"),      n_universe = 10L) # 0.8
#> [1] 0.8
adjusted_similarity(character(0),      c("f1"),           n_universe = 10L) # 0
#> [1] 0
```
