# Summary Statistic of Pearson Similarity Values

Convenience wrapper analogous to
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
but using the Pearson-corrected similarity. See
[`pearson_similarity()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity.md)
for a description of the underlying metric and when to prefer it over
the adjusted similarity.

## Usage

``` r
pearson_similarity_summary(
  list_of_lists,
  nb_total_elements,
  stat = "median",
  d = NULL
)
```

## Arguments

- list_of_lists:

  A list of character/integer vectors, one per run.

- nb_total_elements:

  Integer; total number of candidate features.

- stat:

  Character; `"median"` (default) or `"mean"`.

- d:

  Deprecated; use `nb_total_elements` instead.

## Value

A named list with `statistic` and `err` (see
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
for the exact definitions).

## See also

[`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md),
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)

## Examples

``` r
sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
pearson_similarity_summary(sets, nb_total_elements = 10L)
#> $statistic
#> [1] 0.5238095
#> 
#> $err
#> [1] 0.2857143 0.5238095
#> 
pearson_similarity_summary(sets, nb_total_elements = 10L, stat = "mean")
#> $statistic
#> [1] 0.3650794
#> 
#> $err
#> [1] 0.2244783
#> 
```
