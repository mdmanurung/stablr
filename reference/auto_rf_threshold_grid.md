# Auto-Generate an Importance-Threshold Grid for Random-Forest STABL (Experimental)

Builds the random-forest analog of
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md):
one ranger forest with permutation importance is fit on the full data,
and the grid is taken from quantiles of the pooled importance
distribution (`probs = seq(0.5, 0.999, length.out = n_grid)`).
Experimental extension beyond Python STABL.

## Usage

``` r
auto_rf_threshold_grid(x, y, n_grid = 30L, num_trees = 500L, seed = NULL)
```

## Arguments

- x:

  A `matrix` or `data.frame` of features (samples in rows).

- y:

  Outcome vector.

- n_grid:

  Positive integer; number of grid points.

- num_trees:

  Positive integer; number of trees per ranger forest.

- seed:

  Optional integer seed for the forest fit.

## Value

A `data.frame` with columns `lambda` (the importance threshold) and
`threshold` (identical), sorted in *decreasing* order so the path runs
sparse-to-dense, matching the lambda-path convention used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## See also

[`make_ranger_adapter()`](https://mdmanurung.github.io/stablr/reference/make_ranger_adapter.md),
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md).

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(x[, 1] * x[, 2] * 0.5 + rnorm(60, sd = 0.3), rownames(x))
  grid <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  grid
}
#>        lambda   threshold
#> 1 0.056898277 0.056898277
#> 2 0.049516610 0.049516610
#> 3 0.036206947 0.036206947
#> 4 0.014181369 0.014181369
#> 5 0.004090671 0.004090671
```
