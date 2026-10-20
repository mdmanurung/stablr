# Outcome-Null Calibration Check for STABL

Repeatedly reruns the complete
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
pipeline on outcome-null versions of the data (permuted or explicitly
simulated) and summarizes null selection behavior: the
selected-panel-size distribution, the probability of selecting anything
at all, per-feature selection frequencies, reported-FDP+ versus
by-construction realized FDP, and null-saturation proportions.

## Usage

``` r
stabl_null_check(
  x,
  y,
  design = c("independent", "within_block", "whole_cluster", "simulated"),
  fit_spec = list(),
  n_reps = 20L,
  blocks = NULL,
  clusters = NULL,
  null_simulator = NULL,
  decoy_reps = 1L,
  seed = NULL,
  max_failures = NULL,
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

- design:

  Exchangeability unit prescribing the null: one of `"independent"`
  (global permutation; refused when `blocks` or `clusters` are
  supplied), `"within_block"` (permute within `blocks`),
  `"whole_cluster"` (permute whole cluster outcome vectors across
  `clusters`), or `"simulated"` (use `null_simulator`).

- fit_spec:

  Named list of further arguments for
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
  e.g. `list(n_bootstraps = 50, artificial_type = "knockoff_equi")`. The
  per-replicate `random_state` is always set by `stabl_null_check`;
  `lambda_grid` defaults to `"auto"` when not supplied.

- n_reps:

  Number of null replicates (permutations or simulated nulls). Total
  fits = `n_reps * decoy_reps`.

- blocks:

  Optional factor/vector of length `nrow(x)` giving block membership;
  required for `design = "within_block"`.

- clusters:

  Optional factor/vector of length `nrow(x)` giving cluster membership;
  required for `design = "whole_cluster"`.

- null_simulator:

  Required for `design = "simulated"`: a function
  `function(y, x, rep_id, seed)` returning a null outcome vector (named
  like `y`). Under this design the truth is synthetic: every selected
  feature is a false positive by construction.

- decoy_reps:

  Number of fit replicates per null replicate with distinct
  decoy-generator seeds; values \> 1 separate sensitivity to the random
  decoy realization from sensitivity to the null permutation.

- seed:

  Optional base seed. Per-replicate permutation and fit seeds are
  derived from it by the package's hashed multi-stream derivation and
  stored in the result, so any replicate is exactly replayable. If NULL,
  a seed is drawn from the current RNG state and recorded.

- max_failures:

  Bound on tolerated fit failures. Every failure is recorded with its
  replicate id, seeds, and error message; if failures exceed
  `max_failures` the run aborts (fail-closed) rather than summarizing a
  mostly-failed design.

- ...:

  Reserved; unused arguments are an error.

## Value

A `stabl_null_validation` object with elements:

- design, truth_status:

  The exchangeability unit and whether the null truth is
  `"synthetic_global_null"` (design `"simulated"`) or
  `"induced_global_null_permutation"` (permutation designs).

- n_reps, decoy_reps, n_failures, failures:

  Replication counts and a data frame of bounded failure records.

- replicates:

  Per-replicate results: seeds, panel size, selected features, reported
  FDP+ minimum and threshold, at-selection real/artificial proportions,
  and decoy-generator provenance.

- panel_sizes:

  Integer vector of selected-panel sizes.

- p_any_selection:

  Probability the null pipeline selects at least one feature.

- feature_frequency:

  Named, decreasing vector of per-feature null selection frequencies.

- reported_fdp_plus:

  Reported FDP+ minima across replicates; under the global null the
  realized FDP is 1 whenever the panel is non-empty, by construction.

- fit_spec, call:

  The executed fit specification and call.

## Interpretation caveat

Null resampling quantifies empirical selection under a globally null
DGP. It does not validate FDP control in partly non-null observational
data, and label permutation is not a substitute for modeled confounders.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
nv <- stabl_null_check(
  x, y,
  design = "independent",
  fit_spec = list(n_bootstraps = 10,
                  artificial_type = "random_permutation"),
  n_reps = 3, seed = 7
)
nv$p_any_selection
#> [1] 0
```
