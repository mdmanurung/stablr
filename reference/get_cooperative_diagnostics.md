# Get Cooperative-Fusion Tuning Diagnostics

Returns the cooperative-fusion tuning diagnostics from a fitted
multi-omic workflow. For train/validation fits this is the per-candidate
tuning table; for outer cross-validation fits this is the fold
diagnostics table restricted to cooperative diagnostic columns.

## Usage

``` r
get_cooperative_diagnostics(object)
```

## Arguments

- object:

  A `"stabl_multiomic_fit"` or `"stabl_multiomic_cv"` object with
  cooperative fusion enabled.

## Value

A `data.frame` of cooperative tuning diagnostics.

## See also

[`get_cooperative_features()`](https://mdmanurung.github.io/stablr/reference/get_cooperative_features.md),
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
get_cooperative_diagnostics(fit)   # data.frame of tuning results
#>   rho    lambda metric_value selected
#> 1   0 0.3397982     1.092782     TRUE
# }
```
