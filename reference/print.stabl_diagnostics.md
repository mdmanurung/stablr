# Print STABL Selection Diagnostics

Print STABL Selection Diagnostics

## Usage

``` r
# S3 method for class 'stabl_diagnostics'
print(x, ...)
```

## Arguments

- x:

  A `stabl_diagnostics` object from
  [`stabl_diagnose()`](https://mdmanurung.github.io/stablr/reference/stabl_diagnose.md).

- ...:

  Ignored; present for S3 compatibility.

## Value

The input `x`, invisibly.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
fit <- stabl_fit(x, y, n_bootstraps = 20, lambda_grid = "auto",
                 artificial_type = "random_permutation")
print(stabl_diagnose(fit))
#> STABL selection diagnostics
#>   learner: lasso | family: gaussian | decoys: random_permutation
#>   features: 8 real / 8 artificial | lambdas: 23 | bootstraps: 20
#>   selection rule: fdp_plus_minimum | threshold: 1 | comparator: strictly greater (>)
#>   FDP+ minimum: 1.12 at threshold 0 (95 minimizing thresholds)
#>   evidence:
#>     - [empirically_concerning] lasso x random_permutation x gaussian (sparse independent, linear signal, low dimension)
#>     - [empirically_concerning] lasso x any x gaussian (pure null (no signal))
#>     - [validated_in_stablr] lasso x any x gaussian (block-correlated AR(0.6), weak signal, n=100 (Phase 6 D3 conditions))
#>   warning codes:
#>     - STABL_DIAG_TIED_MINIMA
#>     - STABL_DIAG_FLAT_MINIMUM
#>     - STABL_DIAG_BOUNDARY_MINIMUM
#>     - STABL_DIAG_HIGH_SCORE_DEPENDENCE
#>     - STABL_DIAG_EMPTY_SUPPORT
#>   assumption notes:
#>     - The first FDP+ minimum lies at the edge of fdr_threshold_range; consider widening the range before trusting the minimum.
#>     - Mean absolute score-profile correlation among real features is 0.79: features co-select across the regularization path. Interpret per-feature selection cautiously; this probe does not by itself invalidate the decoy design.
#>   calibration scope:
#>     - The FDP+ minimum is a stochastic-bound diagnostic, not a calibrated estimate of realized FDP (F-001). On pure-null designs the minimizing rule still selected 18-33 features with realized FDP 0.92-0.93 while reporting 0.42-0.66 (audit Phase 6, F-P6-1): an empty support is the only outcome a null dataset cannot produce by this rule.
#>     - Calibration is regime-dependent in both directions (audit Phase 6, F-P6-2): approximately calibrated for knockoff decoys on independent signal, anti-conservative for permutation decoys, conservative under block correlation. No single interpretation of the FDP+ minimum is safe across regimes.
```
