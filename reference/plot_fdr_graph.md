# Plot the FDP+ Diagnostic Curve

Displays how the FDP+ diagnostic (a decoy-based estimate of the false
discovery proportion under the STABL procedure) changes across the full
range of candidate stability thresholds, and marks the threshold
selected by the first-minimum rule.

## Usage

``` r
plot_fdr_graph(object, title = "FDP+ diagnostic", reference_fdp_plus = 0.05)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- title:

  Character scalar; plot title. Default `"FDP+ diagnostic"`.

- reference_fdp_plus:

  Numeric scalar or `NULL`; a purely visual horizontal reference line
  only. It is not used for fitting or threshold selection and marks no
  control guarantee. Use `NULL` to omit the line. Default `0.05`.
  (Renamed from `fdr_target` in 1.0.0: the line is a visual reference,
  not an FDR-control target.)

## Value

A `ggplot` object. The curve shows the FDP+ diagnostic at each candidate
threshold; a vertical dashed line marks the selected threshold stored in
`object$fdr_min_threshold_`; a horizontal dashed line marks
`reference_fdp_plus` when supplied.

## Details

This diagnostic shows why a particular stability threshold was chosen
during fitting: the current rule selects the first minimum of the
displayed FDP+ curve. Curve shape is descriptive and does not by itself
establish FDP control, calibration, or signal/noise separation; the
minimized value is a data-dependent diagnostic, not a controlled or
nominal FDR. Interpretation depends on the artificial-feature
assumptions and the fitted data.

Requires that `object` was fitted with `artificial_type` set (not
`NULL`).

## See also

[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 8), 60, 8,
               dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
  y <- setNames(rnorm(60), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid        = data.frame(lambda = c(0.3, 0.1, 0.05)),
                   n_bootstraps       = 6L,
                   hard_threshold     = NULL,
                   artificial_type    = "random_permutation",
                   random_state       = 1L)
  plot_fdr_graph(fit)
}

# }
```
