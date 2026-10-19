# attic/ — archived experimental code (not part of the built package)

This directory retains source and tests removed from the stablr core API so
the record of what was evaluated stays inspectable. Nothing here is built,
exported, or checked (`attic/` is listed in `.Rbuildignore`); code here may
not run against the current package without modification.

## Contents

- `learner_oblique_rf.R` (formerly `R/learner_oblique_rf.R`) — the
  experimental oblique random-forest learner adapter (`make_aorsf_adapter()`,
  `auto_oblique_threshold_grid()`, and the internal `.make_aorsf_batch_adapter()`
  / `.aorsf_fit_negation_importance()` helpers), backed by the `aorsf`
  package.
- `tests/test-oblique-rf-adapter.R` (formerly
  `tests/testthat/test-oblique-rf-adapter.R`) — its testthat suite.

## Why the oblique random forest was retired (1.0.0, Phase 13 PR-G)

The WS-3 benchmark (`benchmarks/ws3_summary.md`,
`benchmarks/ws3_diagnosis.csv`) evaluated the oblique-RF adapter against the
predeclared FDP+ calibration bar; no predeclared remedy met the criterion on
the tested regimes (knockoff decoy importance is inflated for oblique
forests, so FDP+ is anti-conservative). Per the 1.0.0 specification (§5.5),
the algorithm is retired from the core rather than shipped with a
calibration claim it did not earn. `stabl_fit(base_learner = "oblique_rf")`
now fails closed with a pointer here; use `base_learner = "rf"` for the
calibrated-in-scope axis-aligned forest. Re-opening requires predeclared
positive new calibration evidence.

Historical benchmark scripts under `benchmarks/` that reference
`oblique_rf` are preserved as the record of the experiments that were run;
they describe the pre-retirement package.
