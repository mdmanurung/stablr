# Build an Ordinal-Regression Learner Adapter (Experimental)

Returns a closure that fits ONE elastic-net cumulative-logit model
(`ordinalNet`) on a single bootstrap subsample and returns a logical
selection mask over the features, mirroring the role of
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md)
for ordered outcomes.

## Usage

``` r
make_ordinal_adapter(alpha = 1)
```

## Arguments

- alpha:

  Numeric in `(0, 1]`; the elastic-net mixing weight inside
  `ordinalNet`. Default `1` (pure lasso).

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)` for use
inside
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
where `lambda_val` is a 1-row `data.frame` from the lambda grid (see
[`auto_ordinal_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_ordinal_lambda_grid.md)).

## Selection rule

With `parallelTerms = TRUE` the proportional-odds model gives one slope
per feature. A feature is selected at a grid row when the absolute value
of its slope exceeds the same small coefficient cutoff the glmnet
adapters use. The grid is a lambda path (see
[`auto_ordinal_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_ordinal_lambda_grid.md));
lower penalties give supersets of the selections at higher penalties, so
the stability-accumulation and FDP+ logic in
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
work unchanged.

## Statistical guardrails (read before using)

- **Outcome type.** `y` must be a factor; its level order is the band
  order. An ordered factor is recommended. Numeric outcomes must be
  banded by the caller (for example with
  [`cut()`](https://rdrr.io/r/base/cut.html)).

- **Experimental.** This learner is an extension beyond Python STABL and
  has no Python-parity counterpart. FDP+ calibration for the ordinal
  learner is validated by the WS-9 simulation benchmark; consult the
  benchmark summary before treating FDP+ as calibrated.

- **Decoys.** Artificial features are outcome-agnostic; both
  `random_permutation` and model-X knockoffs apply unchanged.

## Calibration status (experimental)

Simulated calibration evidence (10 seeds x 2 proportional-odds designs)
shows knockoff_equi mean FDP+ of 0.00 (nominal 0.1) but mean TPR of only
0.37 on both designs – the learner is strongly conservative and did not
meet the predeclared graduation bar. It stays experimental; see
benchmarks/ws9_summary.md.

## See also

[`auto_ordinal_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_ordinal_lambda_grid.md),
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("ordinalNet", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(400L), 40L, 10L,
              dimnames = list(NULL, paste0("f", 1:10)))
  eta <- x[, 1] - 0.5 * x[, 2]
  y   <- cut(eta + rnorm(40L, sd = 0.5), breaks = 4,
             labels = paste0("L", 1:4), ordered_result = TRUE)
  grid   <- stablr:::auto_ordinal_lambda_grid(x, y, n_lambda = 5L)
  adapter <- stablr:::make_ordinal_adapter()
  mask <- adapter(x, y, grid[1L, , drop = FALSE])
  cat("selected at sparsest lambda:", sum(mask), "of", ncol(x), "\n")
}
#> selected at sparsest lambda: 0 of 10 
```
