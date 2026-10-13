# Get the Full Stability Score Matrix from a Fitted STABL Object

Returns the raw stability-score matrix accumulated over all bootstrap
iterations. Inspecting this matrix is useful for diagnostics: you can
check how stable features are across the regularisation path, identify
features that are consistently selected at many lambda values
(robustness), and spot features that peak only at one extreme of the
path (fragility).

## Usage

``` r
get_stabl_scores(object)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Value

Numeric matrix with one row per original feature and one column per
lambda in the fitted grid. Row names are the feature names from the
training matrix; column names are not set (use
`object$fitted_lambda_grid` to map column indices to lambda values).

## Details

The stability score for feature \\j\\ at regularisation strength
\\\lambda_k\\ is defined as the fraction of bootstrap subsamples in
which feature \\j\\ received a non-zero coefficient when the model was
fitted at \\\lambda_k\\. Values lie in \\\[0, 1\]\\.

## See also

[`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md)
for the per-feature maximum score (scalar summary),
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
for the binary selection mask.

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
scores <- get_stabl_scores(fit)  # features x lambdas matrix
dim(scores)
#> [1] 5 3
```
