# Diagnose selection reliability of a fitted STABL object

Inspects a `stabl_fit` result and returns a structured
`stabl_diagnostics` object: score-distribution summaries,
decoy-generator provenance, the full threshold path of real/decoy
selection counts, FDP+ curve features (selected cut-off, ties, flat
minima, boundary solutions, low-decoy-count instability), a
null-saturation probe, a score-dependence assumption probe, and a
learner x decoy x outcome x dimensionality evidence profile with
benchmark provenance.

## Usage

``` r
stabl_diagnose(
  fit,
  ...,
  saturation_threshold = 0.05,
  saturation_proportion = 0.5,
  low_decoy_count = 5L
)
```

## Arguments

- fit:

  A fitted
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  object (including
  [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
  results, which inherit from `stabl_fit`).

- ...:

  Reserved for future use; passing arguments here is an error.

- saturation_threshold:

  Numeric in \[0, 1\]; thresholds at or below this value are considered
  "near-zero" for the null-saturation probe. Default `0.05`.

- saturation_proportion:

  Numeric in (0, 1\]; the probe flags null-like saturation when a
  near-zero threshold selects at least this proportion of real *and* of
  artificial features. Default `0.5`.

- low_decoy_count:

  Non-negative integer; the FDP+ minimum is flagged as unstable when the
  artificial-feature count at the selected threshold is below this
  value. Default `5L`.

## Value

An object of class `stabl_diagnostics` with components:

- score_summary:

  List of distribution summaries (quantiles, mean) for real and
  artificial per-feature maximum scores.

- decoy_provenance:

  The fit's `artificial_provenance` (requested vs actual generator,
  fallback counts), or `NULL`.

- threshold_path:

  `data.frame` with `threshold`, `n_real`, `n_artificial`, `prop_real`,
  `prop_artificial`, and `fdp_plus` along the full
  `fdr_threshold_range`.

- per_lambda:

  `data.frame` of per-lambda real/artificial selected counts at the
  selected threshold (regularization-path view).

- cutoff:

  List describing the selected threshold, the selection rule in effect
  (`"user_hard_threshold"`, `"fdp_plus_minimum"`, or `"none"`), the
  comparator, and the FDP+ denominator convention.

- fdp_curve:

  List with the FDP+ curve, all minimizing thresholds, and
  tie/flat/boundary flags.

- null_saturation:

  List with the null-saturation probe result.

- dependence_probe:

  List with score-profile correlation summaries (assumption probe only).

- evidence_profile:

  `data.frame` of matching evidence rows with `status` and `provenance`.

- warning_codes:

  Character vector of machine-readable codes.

- assumption_notes:

  Character vector of assumption statements.

- calibration_scope:

  Character vector delimiting what the FDP+ diagnostic does and does not
  establish.

## Details

The diagnostic is read-only: it never changes the fit or its selected
support. It distinguishes a user-supplied `hard_threshold` from the
minimized-FDP+ rule and reports the effective threshold comparator
(strict `>`) and the FDP+ denominator convention
`(n_artificial / artificial_proportion + 1) / pmax(1, n_real)`.

## See also

[`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md),
[`plot_fdr_graph()`](https://mdmanurung.github.io/stablr/reference/plot_fdr_graph.md),
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
fit <- stabl_fit(x, y, n_bootstraps = 20, lambda_grid = "auto",
                 artificial_type = "random_permutation")
d <- stabl_diagnose(fit)
d$warning_codes
#> [1] "STABL_DIAG_TIED_MINIMA"           "STABL_DIAG_FLAT_MINIMUM"         
#> [3] "STABL_DIAG_BOUNDARY_MINIMUM"      "STABL_DIAG_HIGH_SCORE_DEPENDENCE"
#> [5] "STABL_DIAG_EMPTY_SUPPORT"        
```
