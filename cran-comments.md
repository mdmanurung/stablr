# CRAN submission comments — DRAFT

This file is intentionally provisional. Do not submit it until the exact
`stablr_0.3.0.tar.gz` release candidate has completed current R release/devel,
Windows, macOS, Win-builder, and sanitizer checks. Final error, warning, note,
test, and skip counts will be copied from those immutable-candidate results.
Preliminary sandbox results (single Linux host, conda toolchain): 0 errors,
0 warnings, 3 environmental NOTEs; full test suite green. These are not the
final numbers.

This will be the first CRAN submission of `stablr`.

## What changed in 0.3.0

First release-candidate version. Relative to the development history in the
public repository, 0.3.0 adds: a two-tier user interface (beginner `stabl()`
wrapper over the expert `stabl_fit()`), unambiguous renames with warn-once
deprecation shims, study-disjoint nested cross-validation
(`study_labels` on `stabl_multiomic_nested_cv()`), per-study diagnostics and
study-disjoint train/validation split helpers, per-study performance
reporting, a sample-i.i.d. vs study-disjoint scheme comparison helper,
selection-stability summaries across study-disjoint folds, fail-closed
guardrails for degenerate grouped/classification folds, and a
public-interface audit with deterministic contract tests.

## Vendored cooperative S3 registrations

The package vendors a parity-locked subset of the cooperative-learning engine.
Seven cooperative S3 registrations are maintained manually in `NAMESPACE`
because the vendored source files (`R/cooperative-*.R`) cannot carry roxygen
`@export` tags without breaking parity provenance: `coef.multiview`,
`coef.cv.multiview`, `family.multiview`, `family.cv.multiview`,
`plot.multiview`, `predict.multiview`, and `predict.cv.multiview`.

## Tests skipped on CRAN (3)

`tests/testthat/test-fdp-calibration.R`,
`tests/testthat/test-methodology-validation-runner.R`, and
`tests/testthat/test-signal-recovery.R` carry `skip_on_cran()` guards. The
long-form calibration, methodology-validation, and signal-recovery sweeps run
in CI (GitHub Actions R-CMD-check on Windows, macOS, and Ubuntu
release/devel/oldrel) and are skipped on CRAN to keep check time within
policy limits.

## Compiled code

The package vendors a GPL-2 cooperative-learning implementation and glmnetpp
headers under `src/`, declares `SystemRequirements: C++17`, and therefore is
R-native without requiring a Python runtime rather than pure R.
