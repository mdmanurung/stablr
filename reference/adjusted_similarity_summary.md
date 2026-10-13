# Summary Statistic of Adjusted Similarity Values

Convenience wrapper that computes all pairwise adjusted similarities
(via
[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md))
and reduces them to a single location statistic with an associated
spread measure. Useful for reporting a single reproducibility number per
STABL configuration in benchmarking tables.

## Usage

``` r
adjusted_similarity_summary(list_of_lists, nb_total_elements, stat = "median")
```

## Arguments

- list_of_lists:

  A list of character/integer vectors, one per run.

- nb_total_elements:

  Integer; total number of candidate features.

- stat:

  Character; `"median"` (default, robust to outliers) or `"mean"`.

## Value

A named list with two elements:

- `statistic`:

  The median (or mean) of adjusted-similarity values.

- `err`:

  For `"median"`: the 25th and 75th percentile vector (IQR bounds). For
  `"mean"`: the root-mean-squared deviation (RMSD / population SD).

## See also

[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md),
[`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)

## Examples

``` r
sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
adjusted_similarity_summary(sets, nb_total_elements = 10L)
#> $statistic
#> [1] 0.3666667
#> 
#> $err
#> [1] 0.2000000 0.3666667
#> 
adjusted_similarity_summary(sets, nb_total_elements = 10L, stat = "mean")
#> $statistic
#> [1] 0.2555556
#> 
#> $err
#> [1] 0.1571348
#> 
```
