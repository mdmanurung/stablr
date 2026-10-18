# Get the Names of Selected Features from a Fitted STABL Object

Convenience wrapper around
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
that returns only the names of the features that pass the stability
threshold, ready for use as column selectors in downstream modelling
(e.g. `x[, get_feature_names(fit)]`).

## Usage

``` r
get_feature_names(object, new_hard_threshold = NULL)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- new_hard_threshold:

  Numeric in `(0, 1]` or `NULL`. Forwarded to
  [`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md);
  see that function for the full threshold resolution order.

## Value

Character vector of selected feature names. An empty character vector is
returned when no feature passes the threshold.

## See also

[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
for the binary mask,
[`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md)
for ranked scores.

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
sel <- get_feature_names(fit)   # character vector
x_sel <- x[, sel, drop = FALSE]    # subset to selected columns
```
