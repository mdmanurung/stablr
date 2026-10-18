# Print an Environment-Stability Analysis

Compact summary of an environment-stability analysis: the environments
analyzed, the number of features, and the top features by
minimum-across-environment stability score.

## Usage

``` r
# S3 method for class 'stabl_env_stability'
print(x, ...)
```

## Arguments

- x:

  A `"stabl_env_stability"` object returned by
  [`stabl_environment_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_environment_stability.md).

- ...:

  Unused; included for S3 consistency.

## Value

The input object, invisibly.

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
print(res)
#> <stabl_env_stability>
#>   Environments:    3 (batchA, batchB, batchC)
#>   Features in:      6 
#>   Min score >= 0.5: 2 feature(s)
#>   Top by min score: f1, f2, f4, f3, f5 
```
