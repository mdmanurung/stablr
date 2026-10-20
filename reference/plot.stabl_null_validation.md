# Plot Outcome-Null Validation Results

Plot Outcome-Null Validation Results

## Usage

``` r
# S3 method for class 'stabl_null_validation'
plot(
  x,
  which = c("panel_sizes", "feature_frequency", "reported_vs_realized", "saturation"),
  ...
)
```

## Arguments

- x:

  A `stabl_null_validation` object from
  [`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md).

- which:

  Character vector of panels: `"panel_sizes"` (null panel-size
  distribution), `"feature_frequency"` (per-feature null selection
  frequency), `"reported_vs_realized"` (reported FDP+ minima against the
  by-construction realized FDP of 1 on non-empty null panels), and
  `"saturation"` (real vs artificial at-selection proportions per
  replicate).

- ...:

  Ignored.

## Value

Invisibly, a named list of the drawn `ggplot` objects.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1)
  x <- matrix(rnorm(60 * 8), 60, 8,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
  y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
  nv <- stabl_null_check(
    x, y,
    fit_spec = list(n_bootstraps = 10,
                    artificial_type = "random_permutation"),
    n_reps = 3, seed = 7
  )
  plot(nv, which = "panel_sizes")
}
```
