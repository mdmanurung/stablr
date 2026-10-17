# Build an Automatic keepX Grid for the sPLS-DA Learner

Returns a pre-expanded lambda grid for `base_learner = "splsda"` whose
rows sweep sPLS-DA's per-component `keepX` (the number of features
retained per latent component). The grid is geometric in `p = ncol(x)`,
spanning from roughly `p / 50` up to `p / 3`, which keeps the sparse end
informative for FDP+ thresholding while still probing supports large
enough to matter for prediction.

## Usage

``` r
auto_splsda_keepx_grid(x, n_lambda = 10L)
```

## Arguments

- x:

  Numeric matrix; the real feature matrix (only `ncol(x)` is used).

- n_lambda:

  Positive integer; number of grid rows. Default `10L`. Duplicated
  rounded values are dropped, so the returned grid may be shorter than
  `n_lambda`.

## Value

A `data.frame` with columns `lambda` and `keepX` (identical positive
integer values, sorted increasingly).

## Details

The returned `data.frame` carries the keepX values in both the `lambda`
column (the column
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
requires) and a dedicated `keepX` column (which
[`make_splsda_adapter()`](https://mdmanurung.github.io/stablr/reference/make_splsda_adapter.md)
reads preferentially).

## See also

[`make_splsda_adapter()`](https://mdmanurung.github.io/stablr/reference/make_splsda_adapter.md),
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md)

## Examples

``` r
grid <- auto_splsda_keepx_grid(matrix(0, 10, 150), n_lambda = 8)
grid
#>   lambda keepX
#> 1      3     3
#> 2      4     4
#> 3      7     7
#> 4     10    10
#> 5     15    15
#> 6     22    22
#> 7     33    33
#> 8     50    50
```
