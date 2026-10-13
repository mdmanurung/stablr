# Get Cooperative-Fusion Selected Features

Returns the feature names selected by the cooperative-fusion branch of a
fitted multi-omic workflow. This accessor gives downstream code a stable
public surface instead of reaching into `$cooperative_fusion` directly.

## Usage

``` r
get_cooperative_features(object, view = NULL)
```

## Arguments

- object:

  A fitted `"stabl_multiomic_fit"` object returned by
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  with `cooperative_fusion = TRUE`, or a `"stabl_multiomic_cv"` object
  returned by
  [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  with `cooperative_fusion = TRUE`.

- view:

  Optional character scalar naming one omic view. When supplied for a
  `"stabl_multiomic_fit"`, only that view's cooperative feature names
  are returned.

## Value

For `"stabl_multiomic_fit"`, a named list of character vectors, or a
character vector when `view` is supplied. For `"stabl_multiomic_cv"`, a
named list keyed by fold, where each element has the same structure as
the `"stabl_multiomic_fit"` method.

## See also

[`get_cooperative_diagnostics()`](https://mdmanurung.github.io/stablr/reference/get_cooperative_diagnostics.md),
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md),
[`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)

## Examples

``` r
# \donttest{
set.seed(1L)
n <- 40L
x_list <- list(
  omic1 = matrix(rnorm(n * 5), n, 5,
                 dimnames = list(paste0("s", 1:n), paste0("g", 1:5))),
  omic2 = matrix(rnorm(n * 4), n, 4,
                 dimnames = list(paste0("s", 1:n), paste0("p", 1:4)))
)
y <- setNames(rnorm(n), paste0("s", 1:n))
train <- 1:30
valid <- 31:40
x_train <- lapply(x_list, function(x) x[train, , drop = FALSE])
x_valid <- lapply(x_list, function(x) x[valid, , drop = FALSE])
lam <- data.frame(lambda = c(0.3, 0.1))
fit <- stabl_multiomic_train_validate(
  x_train, y[train], lambda_grid = lam,
  x_valid_list = x_valid, y_valid = y[valid],
  family = "gaussian",
  n_bootstraps = 4L,
  hard_threshold = 0.2, cooperative_fusion = TRUE,
  random_state = 1L
)
get_cooperative_features(fit)           # named list of feature vectors
#> $omic1
#> character(0)
#> 
#> $omic2
#> character(0)
#> 
get_cooperative_features(fit, "omic1") # one view only
#> character(0)
# }
```
