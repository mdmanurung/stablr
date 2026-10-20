# Cohort Replication and Heterogeneity Profile

Characterizes cross-cohort replication of feature selection: pooled
selection stability, training-cohort replication across
leave-one-study-out fits, leave-one-study-out sensitivity of the pooled
panel, optional direction/effect replication from user-supplied
per-study effect estimates, and optional held-out cohort performance via
a documented downstream-refit design.

## Usage

``` r
stabl_replication_profile(
  x,
  y,
  study_labels,
  fit_spec = list(),
  pooled_fit = NULL,
  effects = NULL,
  heldout_performance = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  Numeric predictor matrix with non-empty row names (sample ids) and
  column names (feature ids), as for
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- y:

  Outcome vector, named to match `rownames(x)`.

- study_labels:

  Named vector (factor-coercible) assigning each sample to a
  study/cohort; names must match `rownames(x)` exactly. At least three
  studies are required.

- fit_spec:

  Named list of further arguments for
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
  e.g. `list(n_bootstraps = 50, artificial_type = "knockoff_equi")`. The
  per-fit `random_state` is always set by `stabl_replication_profile`;
  `lambda_grid` defaults to `"auto"` when not supplied.

- pooled_fit:

  Optional pre-computed pooled
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  result on `(x, y)`; validated for feature-name agreement and refitted
  if absent.

- effects:

  Optional data frame with columns `feature`, `study`, `estimate`, `se`
  giving per-study effect estimates for direction replication. Estimates
  must come from a documented per-study estimation design; they are used
  as supplied. Duplicated (feature, study) rows are an error.

- heldout_performance:

  Logical; if `TRUE`, evaluate each LOSO fold's selected panel on the
  held-out study via a downstream
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) refit (gaussian:
  MSE; binomial: balanced error rate). Empty panels contribute `NA` with
  reason `"empty_support"`.

- seed:

  Optional base seed for the per-fit derived seeds (recorded; fully
  replayable). If `NULL`, drawn from the current RNG state and recorded.

- ...:

  Reserved; unused arguments are an error.

## Value

A `stabl_replication_profile` object: a list with

- features:

  Per-feature data frame: `n_measured_cohorts`, `n_eligible_cohorts`,
  `within_study_frequency`, `training_replication`, `fold_frequency`,
  `in_pooled_support`, `loso_retention`, `dropped_on_removal_of`, and
  (when `effects` is supplied) `n_estimates`, `sign_concordance`,
  `pooled_estimate`, `pooled_se`, `pooled_ci_lo`, `pooled_ci_hi`,
  `het_Q`, `het_df`, `het_p`.

- estimands:

  Named list describing the four estimands and where each lives in the
  output.

- pooled_support:

  Character vector of the pooled selected panel.

- heldout:

  When `heldout_performance = TRUE`, a per-study data frame of held-out
  metrics with `reason` for `NA`s.

- eligibility:

  Per feature-by-cohort eligibility with `reason` for ineligibility
  (`"all_missing"`, `"constant_in_cohort"`).

- study_labels, fit_spec, seed, fit_seeds, call:

  Provenance.

## Interpretation caveat

The four estimands answer different questions. Selection agreement
across cohorts is not evidence of causal or directional agreement;
direction/effect replication is only as valid as the supplied per-study
estimation design. Proportions are cohort-level observations over
dependent folds – no binomial confidence intervals are computed.

## Examples

``` r
set.seed(1)
ids <- paste0("s", 1:60)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(ids, paste0("f", 1:8)))
y <- setNames(x[, 1] * 3 + rnorm(60, sd = 0.8), ids)
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
#> 3      f3                  3                  3              1.0000000
#> 4      f4                  3                  3              0.3333333
#> 5      f5                  3                  3              1.0000000
#> 6      f6                  3                  3              0.3333333
#>   training_replication fold_frequency in_pooled_support loso_retention
#> 1            1.0000000      1.0000000              TRUE      1.0000000
#> 2            0.6666667      0.6666667              TRUE      0.6666667
#> 3            1.0000000      1.0000000              TRUE      1.0000000
#> 4            0.6666667      0.6666667             FALSE             NA
#> 5            0.6666667      0.6666667              TRUE      0.6666667
#> 6            0.3333333      0.3333333              TRUE      0.3333333
#>   dropped_on_removal_of
#> 1                  <NA>
#> 2                study1
#> 3                  <NA>
#> 4                  <NA>
#> 5                study3
#> 6         study1,study3
```
