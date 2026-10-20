# Scatterplots of Selected Features Against a Continuous Outcome

Produces a ggplot2 figure with one facet per selected feature, showing
the raw data points together with a LOESS smooth and 95% confidence
band. The smooth helps reveal non-linear relationships that might be
missed by reporting correlation coefficients alone. Use after STABL
feature selection for regression tasks.

## Usage

``` r
scatterplot_features(features, x, y, title = "Selected Features", ncol = 3L)
```

## Arguments

- features:

  Character vector of feature names to plot. Features not present as
  columns in `x` are silently skipped.

- x:

  Numeric matrix or data frame of predictors. Column names must include
  all elements of `features`.

- y:

  Numeric vector of continuous outcome values. Must have the same number
  of elements as `nrow(x)`.

- title:

  Character scalar; plot title. Default `"Selected Features"`.

- ncol:

  Positive integer; number of facet columns in the grid. Default 3.

## Value

A `ggplot` object with one facet panel per feature.

## See also

[`boxplot_features()`](https://mdmanurung.github.io/stablr/reference/boxplot_features.md)
for classification tasks,
[`save_stabl_results()`](https://mdmanurung.github.io/stablr/reference/save_stabl_results.md)
which calls this automatically.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(40 * 4), 40, 4,
               dimnames = list(paste0("s", 1:40), c("A", "B", "C", "D")))
  y <- x[, "A"] * 0.8 + rnorm(40, 0, 0.5)
  stablr:::scatterplot_features(c("A", "B"), x, y)
}
```
