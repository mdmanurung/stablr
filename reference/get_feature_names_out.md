# Deprecated: use [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md)

`get_feature_names_out()` is deprecated and will be removed in a future
release; use
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md)
instead.

## Usage

``` r
get_feature_names_out(object, new_hard_threshold = NULL)
```

## Arguments

- object, new_hard_threshold:

  See
  [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md).

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(30 * 4), 30, 4,
            dimnames = list(paste0("s", 1:30), paste0("f", 1:4)))
y <- setNames(x[, "f1"] + rnorm(30, sd = 0.5), rownames(x))
fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = 0.1),
                 n_bootstraps = 10L, artificial_type = NULL,
                 hard_threshold = 0.5)
# warns; use get_feature_names()
same <- suppressWarnings(identical(
  get_feature_names_out(fit),
  get_feature_names(fit)
))
stopifnot(same)
```
