# Reliability diagnostics, validation, and final models

## Contents

- 1 Introduction
- 2 A reference selection
- 3 Diagnosing a selection with
  [`stabl_diagnose()`](https://mdmanurung.github.io/stablr/reference/stabl_diagnose.md)
- 4 Outcome-null validation with
  [`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md)
- 5 Cohort replication with
  [`stabl_replication_profile()`](https://mdmanurung.github.io/stablr/reference/stabl_replication_profile.md)
- 6 Validated development and a portable final model
- 7 Putting it together

## 1 Introduction

A stability-selection run answers “which features?”. The reliability
toolchain answers the four questions that come next:

1.  **Is this selection itself well-behaved?**
    [`stabl_diagnose()`](https://mdmanurung.github.io/stablr/reference/stabl_diagnose.md)
    inspects the fitted object: score distributions, decoy provenance,
    the FDP+ curve along the full threshold path, and saturation or
    boundary behaviour that weakens the error estimate.
2.  **Would the procedure select features when there is nothing to
    find?**
    [`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md)
    re-runs the same fit specification under an explicit outcome-null
    exchangeability design and reports how often anything is selected.
3.  **Does the selection replicate across cohorts?**
    [`stabl_replication_profile()`](https://mdmanurung.github.io/stablr/reference/stabl_replication_profile.md)
    partitions the data by study label and reports pooled stability,
    leave-one-study-out training replication, and optional held-out
    performance.
4.  **How do I ship the final model?**
    [`stabl_finalize()`](https://mdmanurung.github.io/stablr/reference/stabl_finalize.md)
    packs a development result into a portable `"stabl_predictor"` —
    with preprocessing stored, no training observations serialized — and
    [`predict()`](https://rdrr.io/r/stats/predict.html) applies it to
    new data under a strict schema.

All four are read-only with respect to the fitted object: they never
refit, rethreshold, or otherwise modify the selection they inspect.

## 2 A reference selection

A small synthetic problem with three true signals and five noise
features, fit with the expert interface:

``` r

library(stablr)

set.seed(1)
ids <- paste0("s", 1:60)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(ids, paste0("f", 1:8)))
y <- setNames(2 * x[, 1] - x[, 2] + 0.5 * x[, 3] + rnorm(60), ids)

fit <- stabl_fit(x, y, n_bootstraps = 50, lambda_grid = "auto",
                 artificial_type = "random_permutation",
                 random_state = 1)
fit
#> <stabl_fit>
#>   Features in:      8 
#>   Features selected: 0 
#>   Min FDP+:         1.125 
#>   FDP+ threshold:   1 
#>   Artificial:       random_permutation
```

## 3 Diagnosing a selection with `stabl_diagnose()`

[`stabl_diagnose()`](https://mdmanurung.github.io/stablr/reference/stabl_diagnose.md)
takes a fitted `"stabl_fit"` and returns a machine-readable reliability
report:

``` r

d <- stabl_diagnose(fit)
d$warning_codes
#> [1] "STABL_DIAG_TIED_MINIMA"           "STABL_DIAG_FLAT_MINIMUM"         
#> [3] "STABL_DIAG_BOUNDARY_MINIMUM"      "STABL_DIAG_HIGH_SCORE_DEPENDENCE"
#> [5] "STABL_DIAG_EMPTY_SUPPORT"
```

`warning_codes` is the headline field: each entry is a stable,
machine-readable string naming a specific concern (for example a flat or
boundary-hugging FDP+ curve, tie-heavy score distributions, or a low
effective decoy count), so downstream code can branch on it without
parsing prose. An empty vector means no concern fired — it is not a
guarantee of validity.

Two further fields carry the detail:

``` r

d$cutoff
#> $selected_threshold
#> [1] 1
#> 
#> $selection_rule
#> [1] "fdp_plus_minimum"
#> 
#> $hard_threshold
#> NULL
#> 
#> $fdr_min_threshold
#> [1] 1
#> 
#> $comparator
#> [1] "strictly greater (>)"
#> 
#> $denominator
#> $denominator$artificial_proportion
#> [1] 1
#> 
#> $denominator$formula
#> [1] "(n_artificial / artificial_proportion + 1) / pmax(1, n_real)"
#> 
#> $denominator$plus_one_offset
#> [1] TRUE
#> 
#> $denominator$zero_real_guard
#> [1] "pmax(1, n_real)"
head(d$threshold_path)
#>   threshold n_real n_artificial prop_real prop_artificial fdp_plus
#> 1      0.00      8            8         1               1    1.125
#> 2      0.01      8            8         1               1    1.125
#> 3      0.02      8            8         1               1    1.125
#> 4      0.03      8            8         1               1    1.125
#> 5      0.04      8            8         1               1    1.125
#> 6      0.05      8            8         1               1    1.125
d$calibration_scope
#> [1] "The FDP+ minimum is a stochastic-bound diagnostic, not a calibrated estimate of realized FDP (F-001). On pure-null designs the minimizing rule still selected 18-33 features with realized FDP 0.92-0.93 while reporting 0.42-0.66 (audit Phase 6, F-P6-1): an empty support is the only outcome a null dataset cannot produce by this rule."
#> [2] "Calibration is regime-dependent in both directions (audit Phase 6, F-P6-2): approximately calibrated for knockoff decoys on independent signal, anti-conservative for permutation decoys, conservative under block correlation. No single interpretation of the FDP+ minimum is safe across regimes."
```

The `cutoff` block discloses the conventions the FDP+ estimate relies on
— the selection rule (`fdp_plus_minimum`), the comparator (strictly
greater), and the exact denominator formula with its plus-one offset and
zero-real guard; `threshold_path` walks the full threshold grid with
per-threshold selected counts and FDP+ estimates, which is what the
saturation and boundary probes in `warning_codes` summarize. The
evidence profile (`d$evidence_profile`) places the fit’s learner x decoy
x outcome x dimensionality regime against the calibration benchmarks —
it is context for interpretation, not a binary safe/unsafe verdict. Read
`calibration_scope` before trusting any single number: the FDP+ minimum
is a stochastic-bound diagnostic, not a calibrated estimate of realized
FDP, and a low minimum on the curve is not “controlled FDR” (the
plot-reading article develops this distinction).

## 4 Outcome-null validation with `stabl_null_check()`

A selection procedure should select nothing when the outcome carries no
signal.
[`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md)
tests exactly that by refitting the same specification under an explicit
exchangeability design:

``` r

nv <- stabl_null_check(
  x, y,
  design = "independent",
  fit_spec = list(n_bootstraps = 10,
                  artificial_type = "random_permutation"),
  n_reps = 5, seed = 7
)
nv$p_any_selection
#> [1] 0
```

`p_any_selection` is the fraction of null replicates in which anything
was selected at all; on a problem with real signal fitted under the
independent design it should be near zero, because permuting the outcome
destroys the signal. The four designs — `independent`, `within_block`,
`whole_cluster`, and `simulated` — encode different exchangeability
assumptions, and the function refuses fail-closed combinations rather
than silently running an invalid null. Asking for a global permutation
while supplying block structure is an error:

``` r

stabl_null_check(
  x, y,
  design = "independent",
  blocks = rep(1:3, each = 20),
  fit_spec = list(n_bootstraps = 10),
  n_reps = 2, seed = 7
)
#> Error:
#> ! design = "independent" applies a GLOBAL permutation of the outcome, which is not exchangeable under block/cluster structure. Use design = "within_block" or design = "whole_cluster", or drop the structure only if you can defend sample independence.
```

Every replicate seed is derived replayably from `seed`, so a null check
is exactly reproducible.

## 5 Cohort replication with `stabl_replication_profile()`

When samples come from several studies,
[`stabl_replication_profile()`](https://mdmanurung.github.io/stablr/reference/stabl_replication_profile.md)
measures whether the selection is a property of the pooled data or of
individual cohorts:

``` r

studies <- setNames(rep(paste0("study", 1:3), each = 20L), ids)
rp <- stabl_replication_profile(
  x, y, studies,
  fit_spec = list(n_bootstraps = 50,
                  artificial_type = "random_permutation"),
  seed = 5
)
head(rp$features)
#>   feature n_measured_cohorts n_eligible_cohorts within_study_frequency
#> 1      f1                  3                  3              1.0000000
#> 2      f2                  3                  3              1.0000000
#> 3      f3                  3                  3              0.3333333
#> 4      f4                  3                  3              0.0000000
#> 5      f5                  3                  3              0.3333333
#> 6      f6                  3                  3              0.0000000
#>   training_replication fold_frequency in_pooled_support loso_retention
#> 1            1.0000000      1.0000000             FALSE             NA
#> 2            0.6666667      0.6666667             FALSE             NA
#> 3            1.0000000      1.0000000             FALSE             NA
#> 4            1.0000000      1.0000000             FALSE             NA
#> 5            0.6666667      0.6666667             FALSE             NA
#> 6            1.0000000      1.0000000             FALSE             NA
#>   dropped_on_removal_of
#> 1                  <NA>
#> 2                  <NA>
#> 3                  <NA>
#> 4                  <NA>
#> 5                  <NA>
#> 6                  <NA>
```

The `features` table places each feature’s pooled support
(`in_pooled_support`) alongside its leave-one-study-out
`training_replication` and `loso_retention`; `dropped_on_removal_of`
names the cohort whose removal costs the feature its support, so a
selection carried by a single cohort is visible directly. Preprocessing
is recomputed inside each training partition — nothing fitted on
held-out cohorts leaks into selection, and no held-out cohort chooses
thresholds. With `effects =` supplied, the profile additionally pools
user-provided per-study effect estimates (inverse-variance pooling with
Cochran’s Q); with `heldout_performance = TRUE`, a downstream
[`stats::glm`](https://rdrr.io/r/stats/glm.html) refit is evaluated on
each held-out study.

## 6 Validated development and a portable final model

For a held-out performance estimate during development,
`stabl_analyze(mode = "validated")` wraps selection in nested
cross-validation: selection runs inside every outer training fold, and
the outer test folds produce held-out predictions. The gaussian path is
fully supported alongside binomial:

``` r

val <- stabl_analyze(
  x, y,
  family = "gaussian", mode = "validated",
  outer_v = 4L, inner_v = 3L,
  n_bootstraps = 30L,
  random_state = 42,
  final_refit = TRUE
)
val$performance
#> $pooled
#>    n      mae     rmse        r2
#> 1 60 1.465443 2.072181 0.2382401
#> 
#> $by_repeat
#>    n      mae     rmse        r2
#> 1 60 1.465443 2.072181 0.2382401
```

`final_refit = TRUE` refits on the full data after validation, which is
what makes the result eligible for finalization.
[`stabl_finalize()`](https://mdmanurung.github.io/stablr/reference/stabl_finalize.md)
then packs the selected support, the stored preprocessing, and the
refitted model into a portable `"stabl_predictor"`:

``` r

final <- stabl_finalize(val, x, y, preprocessing = c("center", "scale"),
                        random_state = 2)
final
#> <stabl_predictor: gaussian >
#>   engine:      cvglmnet
#>   panel:       7 of 8 features
#>   selected:    f1, f2, f3, f5, f6, ...
#>   requires:    8 input column(s)
#>   preprocessing: center, scale
#>   development: stabl_analysis_validated (full-data refit; not cross-validated performance)
#>   stablr:      0.4.0
predict(final, x[1:3, ], type = "response")
#> [1] -3.5231401  0.7149864 -2.1361474
```

No training observations are serialized into the object — one-feature
models are packed as bare coefficients and preprocessing templates are
stripped — so the predictor can be shipped without disclosing the
development cohort. [`predict()`](https://rdrr.io/r/stats/predict.html)
validates new data against the stored schema: missing, renamed, or
duplicated columns always error, and extra columns error unless
explicitly dropped with `extra_columns = "drop"`.

## 7 Putting it together

The intended order of operations:

1.  Fit and **diagnose**
    ([`stabl_diagnose()`](https://mdmanurung.github.io/stablr/reference/stabl_diagnose.md)):
    resolve any `warning_codes` before trusting the selection.
2.  **Null-check**
    ([`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md))
    under the exchangeability design that matches the sampling
    structure.
3.  **Replicate**
    ([`stabl_replication_profile()`](https://mdmanurung.github.io/stablr/reference/stabl_replication_profile.md))
    when multiple cohorts exist.
4.  **Validate and finalize** (`stabl_analyze(mode = "validated")` with
    `final_refit = TRUE`, then
    [`stabl_finalize()`](https://mdmanurung.github.io/stablr/reference/stabl_finalize.md))
    to obtain a portable model with a held-out performance estimate.

Each step is read-only over the previous one’s output, so the chain can
be scripted, rerun, and audited end to end.
