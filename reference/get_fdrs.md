# Get FDP+ Diagnostic Estimates Across the Swept Thresholds

Returns the vector of FDP+ diagnostic estimates computed by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
across the swept stability thresholds.

## Usage

``` r
get_fdrs(object)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Value

Named numeric vector with one FDP+ estimate per swept stability
threshold (names taken from the fit's `fdr_threshold_range`).

## Details

These are *diagnostic* estimates used to choose the stability threshold,
**not** a guarantee that false discovery proportion is controlled at any
level. Calibration depends on the artificial features acting as fair
null competitors. Fails with an error when the fit was produced without
artificial features (`artificial_type = NULL`).

## See also

[`get_stabl_scores()`](https://mdmanurung.github.io/stablr/reference/get_stabl_scores.md)
for the full selection-frequency path,
[`get_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/get_lambda_grid.md)
for the penalty grid used.

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(40 * 6), 40, 6,
            dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
y <- setNames(rnorm(40), rownames(x))
fit <- stabl_fit(x, y,
                 lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
                 n_bootstraps   = 4L,
                 artificial_type = "random_permutation",
                 random_state   = 1L)
fdrs <- get_fdrs(fit)  # named numeric, one FDP+ estimate per threshold
fdrs
#>        0     0.01     0.02     0.03     0.04     0.05     0.06     0.07 
#> 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 
#>     0.08     0.09      0.1     0.11     0.12     0.13     0.14     0.15 
#> 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 
#>     0.16     0.17     0.18     0.19      0.2     0.21     0.22     0.23 
#> 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 
#>     0.24     0.25     0.26     0.27     0.28     0.29      0.3     0.31 
#> 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 
#>     0.32     0.33     0.34     0.35     0.36     0.37     0.38     0.39 
#> 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 
#>      0.4     0.41     0.42     0.43     0.44     0.45     0.46     0.47 
#> 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 1.166667 
#>     0.48     0.49      0.5     0.51     0.52     0.53     0.54     0.55 
#> 1.166667 1.166667 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 
#>     0.56     0.57     0.58     0.59      0.6     0.61     0.62     0.63 
#> 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 
#>     0.64     0.65     0.66     0.67     0.68     0.69      0.7     0.71 
#> 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 
#>     0.72     0.73     0.74     0.75     0.76     0.77     0.78     0.79 
#> 1.000000 1.000000 1.000000 2.000000 2.000000 2.000000 2.000000 2.000000 
#>      0.8     0.81     0.82     0.83     0.84     0.85     0.86     0.87 
#> 2.000000 2.000000 2.000000 2.000000 2.000000 2.000000 2.000000 2.000000 
#>     0.88     0.89      0.9     0.91     0.92     0.93     0.94     0.95 
#> 2.000000 2.000000 2.000000 2.000000 2.000000 2.000000 2.000000 2.000000 
#>     0.96     0.97     0.98     0.99 
#> 2.000000 2.000000 2.000000 2.000000 
```
