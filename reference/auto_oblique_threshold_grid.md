# Auto-Generate an Importance-Threshold Grid for Oblique-Random-Forest STABL (Experimental)

Builds the oblique-random-forest analog of
[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md):
one aorsf oblique forest with negation importance is fit on the full
data, and the grid is taken from quantiles of the pooled importance
distribution (`probs = seq(0.5, 0.999, length.out = n_grid)`).
Experimental extension beyond Python STABL.

## Usage

``` r
auto_oblique_threshold_grid(x, y, n_grid = 30L, num_trees = 500L, seed = NULL)
```

## Arguments

- x:

  A `matrix` or `data.frame` of features (samples in rows).

- y:

  Outcome vector.

- n_grid:

  Positive integer; number of grid points.

- num_trees:

  Positive integer; number of trees per aorsf oblique forest.

- seed:

  Optional integer seed for the forest fit.

## Value

A `data.frame` with columns `lambda` (the importance threshold) and
`threshold` (identical), sorted in *decreasing* order so the path runs
sparse-to-dense, matching the lambda-path convention used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Details

Negation importance is used instead of permutation importance because it
discriminates relevant from irrelevant predictors better in oblique
forests (Jaeger et al. 2023). Because importance is spread across the
linear combinations used in each split, the resulting per-feature
selection frequencies are noisier than for axis-aligned random forests —
prefer more bootstraps when using `base_learner = "oblique_rf"` in
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## See also

[`make_aorsf_adapter()`](https://mdmanurung.github.io/stablr/reference/make_aorsf_adapter.md),
[`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md),
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md).

## Examples

``` r
if (requireNamespace("aorsf", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(0.8 * x[, 1] + 0.8 * x[, 2] + 0.5 * x[, 1] * x[, 2] +
                rnorm(60, sd = 0.3), rownames(x))
  grid <- auto_oblique_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  grid
}
#>       lambda  threshold
#> 1 1.26845711 1.26845711
#> 2 1.13219414 1.13219414
#> 3 0.80114657 0.80114657
#> 4 0.18370775 0.18370775
#> 5 0.03284498 0.03284498
```
