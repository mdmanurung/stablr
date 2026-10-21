# Summary Statistic of Pairwise Similarity Values

Computes all pairwise similarities between the feature sets in `x` (via
[`stabl_similarity()`](https://mdmanurung.github.io/stablr/reference/stabl_similarity.md))
and reduces them to a single location statistic with an associated
spread measure. Useful for reporting a single reproducibility number per
STABL configuration in benchmarking tables.

## Usage

``` r
stabl_similarity_summary(
  x,
  method = c("jaccard", "adjusted", "pearson"),
  n_universe = NULL,
  stat = c("median", "mean")
)
```

## Arguments

- x:

  A list of character/integer vectors, one per STABL run.

- method:

  Character; one of `"jaccard"` (default), `"adjusted"`, or `"pearson"`.
  See
  [`stabl_similarity()`](https://mdmanurung.github.io/stablr/reference/stabl_similarity.md)
  for method descriptions.

- n_universe:

  Integer; total number of candidate features in the universe. Required
  for `method = "adjusted"` and `method = "pearson"`; ignored for
  `method = "jaccard"`.

- stat:

  Character; `"median"` (default, robust to outliers) or `"mean"`.

## Value

A named list with two elements:

- `statistic`:

  The median (or mean) of the pairwise similarity values.

- `err`:

  For `"median"`: the 25th and 75th percentile vector (IQR bounds). For
  `"mean"`: the root-mean-squared deviation (RMSD / population SD).

## See also

[`stabl_similarity()`](https://mdmanurung.github.io/stablr/reference/stabl_similarity.md)

## Examples

``` r
sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
stabl_similarity_summary(sets)
#> $statistic
#> [1] 0.5
#> 
#> $err
#> [1] 0.35 0.50
#> 
stabl_similarity_summary(sets, method = "adjusted", n_universe = 10L)
#> $statistic
#> [1] 0.3666667
#> 
#> $err
#> [1] 0.2000000 0.3666667
#> 
stabl_similarity_summary(sets, method = "pearson", n_universe = 10L,
                         stat = "mean")
#> $statistic
#> [1] 0.3650794
#> 
#> $err
#> [1] 0.2244783
#> 
```
