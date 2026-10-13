# Get the Fitted Lambda Grid

Returns the penalty grid that the fit consumed, after any `"auto"`
resolution. Useful for reproducing a fit exactly or for inspecting
auto-generated grids.

## Usage

``` r
get_lambda_grid(object)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Value

The `data.frame` penalty grid actually used by the fit (one row per
lambda, with columns `lambda` and, for elastic net, `alpha`).

## See also

[`get_fdrs()`](https://mdmanurung.github.io/stablr/reference/get_fdrs.md)
for the FDP+ diagnostics across the swept thresholds.

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(30 * 5), 30, 5,
            dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
y <- setNames(rnorm(30), rownames(x))
fit <- stabl_fit(x, y,
                 lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
                 n_bootstraps   = 4L,
                 hard_threshold = 0.3,
                 random_state   = 1L)
get_lambda_grid(fit)
#>   lambda
#> 1   0.20
#> 2   0.10
#> 3   0.05
```
