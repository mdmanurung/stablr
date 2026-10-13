# Summarize a Fitted stabl_fit Object

Compact human-readable summary of a STABL fit: the selected feature set,
the stability threshold, the minimal FDP+ diagnostic, and the dimensions
of the fit. The `min_fdr` field is an FDP+ *diagnostic*, not an FDR
guarantee.

## Usage

``` r
# S3 method for class 'stabl_fit'
summary(object, ...)

# S3 method for class 'summary.stabl_fit'
print(x, ...)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- x:

  A `"summary.stabl_fit"` object as returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  Unused; present for S3 consistency.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns a list of
class `"summary.stabl_fit"` with elements `selected_features`,
`n_selected`, `stabl_threshold`, `min_fdr`, `lambda_grid_used`,
`n_features_in` and `n_bootstraps`.
[`print()`](https://rdrr.io/r/base/print.html) returns the object
invisibly.

## See also

[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
for the binary selection mask,
[`get_fdrs()`](https://mdmanurung.github.io/stablr/reference/get_fdrs.md)
for the full FDP+ diagnostic vector.

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
summary(fit)
#> <stabl_fit summary>
#>   Features in:       5 
#>   Features selected: 5 
#>   Threshold:         0.3 
#>   Min FDP+:          0.6 
#>   Lambda grid:      3 x 1
#>   Bootstraps:        4 
#>   Selected:          f1, f2, f3, f4, f5 
```
