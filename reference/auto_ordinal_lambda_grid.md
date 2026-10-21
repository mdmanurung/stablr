# Build a Data-Driven Lambda Grid for the Ordinal Learner (Experimental)

Builds the ordinal analog of
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md):
one `ordinalNet` cumulative-logit model with parallel terms is fit on
the full data, and the grid is taken from its decreasing lambda path.
The grid is computed once on the full data; each bootstrap subsample is
subsequently fit at these fixed penalties with its own model.

## Usage

``` r
auto_ordinal_lambda_grid(x, y, n_lambda = 30L, alpha = 1)
```

## Arguments

- x:

  A `matrix` or `data.frame` of features (samples in rows).

- y:

  Factor outcome; the level order is the band order. An ordered factor
  is recommended.

- n_lambda:

  Positive integer; number of lambda values. Default `30L`.

- alpha:

  Numeric in `(0, 1]`; the elastic-net mixing weight. Default `1` (pure
  lasso).

## Value

A `data.frame` with a single `lambda` column, sorted in decreasing order
so the path runs sparse-to-dense, matching the convention used by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## See also

[`make_ordinal_adapter()`](https://mdmanurung.github.io/stablr/reference/make_ordinal_adapter.md),
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("ordinalNet", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(400L), 40L, 10L,
              dimnames = list(NULL, paste0("f", 1:10)))
  eta <- x[, 1] - 0.5 * x[, 2]
  y   <- cut(eta + rnorm(40L, sd = 0.5), breaks = 4,
             labels = paste0("L", 1:4), ordered_result = TRUE)
  grid <- stablr:::auto_ordinal_lambda_grid(x, y, n_lambda = 5L)
  head(grid)
}
#>        lambda
#> 1 0.340240255
#> 2 0.107593416
#> 3 0.034024025
#> 4 0.010759342
#> 5 0.003402403
```
