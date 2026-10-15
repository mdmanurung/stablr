# Plot the STABL Stability Path

Produces a ggplot2 line chart showing the stability score (frequency of
selection) of each feature across the regularisation path. Selected
features are highlighted in red; noise features are shown in dark grey.
When artificial features were used during fitting, their stability path
is overlaid as thin dotted grey lines. The FDP+-optimal (or hard)
threshold is shown as a dashed horizontal line.

## Usage

``` r
plot_stabl_path(
  object,
  new_hard_threshold = NULL,
  title = "STABL Stability Path"
)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object from
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- new_hard_threshold:

  Numeric in `(0, 1]` or `NULL`. When supplied, overrides the threshold
  stored in `object`.

- title:

  Character; plot title (default `"STABL Stability Path"`).

## Value

A `ggplot` object.

## Details

When the lambda grid contains an `alpha` column (elastic-net mixed-alpha
path), the plot is automatically faceted by `alpha`.

The path is always drawn in decreasing-lambda order (high penalty on the
left, low penalty on the right) within each `alpha` level, independently
of the row order of the supplied lambda grid. The x-axis is the position
along that ordered path, so auto grids (log-spaced lambda) are evenly
spaced.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(40 * 6), 40, 6,
               dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
  y <- setNames(rnorm(40), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid = data.frame(lambda = c(0.3, 0.1, 0.05)),
                   n_bootstraps = 6L, hard_threshold = 0.3, random_state = 1L)
  plot_stabl_path(fit)
}
```
