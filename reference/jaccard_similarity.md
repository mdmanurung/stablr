# Jaccard Similarity Between Two Feature Sets

Measures the overlap between two sets of selected features as a simple
fraction of shared features out of all features that appeared in either
set. This is the most interpretable pairwise similarity when no prior
expectation of set sizes exists.

## Usage

``` r
jaccard_similarity(list1, list2)
```

## Arguments

- list1:

  Character or integer vector of selected feature identifiers.

- list2:

  Character or integer vector of selected feature identifiers.

## Value

Numeric scalar in \\\[0, 1\]\\. A value of 1 means the two sets are
identical; 0 means they share no features.

## Details

Use Jaccard when you want a quick, easy-to-explain measure. Prefer
[`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md)
or
[`pearson_similarity()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity.md)
when comparing across runs that select very different numbers of
features, because Jaccard does not correct for the chance overlap
expected between small or large sets.

When both sets are empty the function returns 0 (convention matching the
Python implementation).

## See also

[`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)
for computing all pairwise Jaccard similarities at once.

## Examples

``` r
jaccard_similarity(c("A", "B", "C"), c("B", "C", "D"))  # 0.5
#> [1] 0.5
jaccard_similarity(c("A", "B"),      c("A", "B"))        # 1
#> [1] 1
jaccard_similarity(character(0),     character(0))       # 0 (both empty)
#> [1] 0
```
