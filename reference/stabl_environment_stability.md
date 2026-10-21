# Environment-Stability Workflow: Reliable and Transportable Feature Sets

Fits STABL separately in each environment and aggregates per-feature
selection frequencies by their *minimum* across environments: a feature
is only as reliable as its weakest environment. This delivers
ICP-inspired stability-across-environments (Peters et al.) as a
*predictive* selection heuristic — without causal assumptions, without
an environment-variable requirement beyond labels, and without new
dependencies.

## Usage

``` r
stabl_environment_stability(x, y, environments, lambda_grid = "auto",
                            random_state = NULL, ...)
```

## Arguments

- x:

  A `matrix` or `data.frame` of features (samples in rows).

- y:

  Named outcome vector; sample IDs must match `rownames(x)`.

- environments:

  Factor or character vector with one label per sample (batch, cohort,
  site, ...). Must contain at least two non-empty levels and no `NA`.

- lambda_grid:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md);
  either `"auto"` (per-environment grids) or an explicit grid. Note that
  `"auto"` grids reach the saturated small-lambda regime, which makes
  max-over-lambda scores lenient; explicit top-heavy grids discriminate
  better.

- random_state:

  Optional integer seed; per-environment seeds are derived
  deterministically from it via the package's hashed multi-stream
  derivation (the same pattern as
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  and the multiomic workflows; changed from the pre-1.0 consecutive
  `random_state + i` derivation, F-019).

- ...:

  Additional arguments forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  in every environment.

## Value

A list of class `"stabl_env_stability"` with elements `env_fits` (named
list of per-environment `"stabl_fit"` objects), `env_scores` (features x
environments matrix of max-over-lambda selection frequencies),
`min_scores` (minimum across environments), `mean_scores`,
`environments` (levels used), `pairwise_adjusted_similarity`
(environment x environment similarity of selected sets), and `call`.

## Details

Environment labels may be batch, cohort, or site identifiers. Group
integrity within each per-environment fit is preserved by passing
environment labels through the usual grouped-bootstrap machinery. The
minimum-aggregation rule is deliberately conservative: features whose
signal is environment-specific receive a low minimum score. This is a
stability heuristic, **not** a causal-effect estimate and not an FDR
guarantee.

## See also

[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
[`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md),
[`stabl_similarity()`](https://mdmanurung.github.io/stablr/reference/stabl_similarity.md).

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
res$min_scores                      # transportability summary
#>  f1  f2  f3  f4  f5  f6 
#> 1.0 1.0 0.0 0.2 0.0 0.0 
get_min_scores(res)                 # same, via accessor
#>  f1  f2  f3  f4  f5  f6 
#> 1.0 1.0 0.0 0.2 0.0 0.0 
```
