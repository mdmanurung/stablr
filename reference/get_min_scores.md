# Get Minimum-Across-Environment Stability Scores

Accessor for the transportability summary of an environment-stability
analysis: the weakest-environment stability score per feature (the
ICP-intersection analog).

## Usage

``` r
get_min_scores(object)
```

## Arguments

- object:

  A `"stabl_env_stability"` object returned by
  [`stabl_environment_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_environment_stability.md).

## Value

Named numeric vector: for each feature, the minimum of its
max-over-lambda selection frequencies across environments.

## See also

[`stabl_environment_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_environment_stability.md).

## Examples

``` r
set.seed(1L)
n <- 90L
env <- factor(rep(c("batchA", "batchB", "batchC"), each = n / 3))
x <- matrix(rnorm(n * 6), n, 6,
            dimnames = list(paste0("s", seq_len(n)), paste0("f", 1:6)))
y <- setNames(1.2 * x[, "f1"] + 0.8 * x[, "f2"] + rnorm(n, sd = 0.5),
              rownames(x))
res <- stabl_environment_stability(
  x, y, env,
  lambda_grid   = data.frame(lambda = c(0.6, 0.4, 0.25)),
  n_bootstraps  = 10L,
  hard_threshold = 0.5,
  random_state  = 1L
)
get_min_scores(res)
#> f1 f2 f3 f4 f5 f6 
#>  1  1  0  0  0  0 
```
