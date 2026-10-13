# Plot a Fitted stabl_fit Object

Convenience S3 method so a fitted STABL model can be visualised with
`plot(fit)`; delegates to
[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md).

## Usage

``` r
# S3 method for class 'stabl_fit'
plot(x, ...)
```

## Arguments

- x:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- ...:

  Arguments forwarded to
  [`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md).

## Value

A `ggplot` object (invisibly), as returned by
[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md).

## See also

[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md)
for the full signature and styling options.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(30 * 5), 30, 5,
              dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
  y <- setNames(rnorm(30), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
                   n_bootstraps   = 4L,
                   hard_threshold = 0.3,
                   random_state   = 1L)
  plot(fit)
}
```
