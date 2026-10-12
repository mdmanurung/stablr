# Reproducible multiview-parity fixtures

These fixtures lock the vendored cooperative-learning backend
(`R/cooperative-*.R`, a fork of CRAN `multiview` 1.0) against reference
values produced by **CRAN multiview 1.0 itself** on fixed input data, and
pin the three deliberate divergences documented in
`inst/PATCHES.cooperative.md` (D4/B1, D5/B2, D6/B3). They extend the
SHA-256-manifest pattern of `tests/testthat/fixtures/python_parity/` and
make the B4 divergence log testable. The reference engine is needed only
to **regenerate** the fixtures; package checks consume only the committed
RDS/CSV/JSON files and never load `multiview`.

## Reproduce

```sh
R_LIBS=<lib> Rscript tests/testthat/fixtures/multiview_parity/generate.R
(cd tests/testthat/fixtures/multiview_parity && sha256sum -c manifest.sha256)
```

Requirements at generation time: CRAN `multiview` 1.0, `glmnet`,
`pkgload` (loads the vendored backend including its compiled C++
routines), and `jsonlite`. The generator runs from the package root,
loads the vendored backend via `pkgload::load_all()`, regenerates every
input from seed 2026, and rewrites all outputs plus `manifest.sha256`.

## Contract map

- `inputs.rds` is the single frozen input bundle: seed 2026, n = 80,
  two views (`a`, `b`) of p = 6 features each, gaussian and binomial
  outcomes, unit and heavy-head/light-tail case weights, zero / mild
  non-zero / aggressive separating offsets, a fixed 5-fold `foldid`,
  user-specified lambda grids, and rho = 0.3.
- `reference_parity_gaussian.csv`, `reference_parity_binomial.csv`, and
  `reference_parity_weighted.csv` freeze coefficients and selection
  schedules per lambda for both engines (`engine` = `vendored` or
  `cran_1_0`) on the paths where the fork must match upstream exactly
  (no offset, unit or non-unit weights). The tests recompute the
  vendored fit live and require it to reproduce both the frozen vendored
  rows (fixture lock) and the frozen CRAN rows (parity lock).
- `reference_b1_predmat_offset.csv` and `reference_b1_cv_offset.csv`
  freeze the vendored `build_predmat()` fold predictions and the
  `cv.multiview` curve for a fit **with an offset**. Upstream multiview
  1.0 self-references `offset_sub` before it exists
  (`R/build_predmat.R:26`) and errors on the same input; the vendored
  fix subsets rows correctly (D4/B1).
- `reference_b2_weighted_offset.csv` freezes weighted fits **with a
  non-zero offset** — the only regime where the halved-step warm-start
  residual weighting (D5/B2) can diverge — including an aggressive
  separating-offset variant when both engines complete it.
- `reference_b3_lower_limits.csv` freezes a gaussian fit with a binding
  `lower.limits` entry (feature `b2` clamped at -0.2; its unbounded
  coefficient is near -0.76), pinning the `lower.limits` plumbing
  (D6/B3).
- `provenance.json` is the machine-readable record of the reference
  engine, versions, seed, and the observed upstream-vs-vendored behavior
  per divergence (the B4 divergence log input). `manifest.sha256` is the
  integrity manifest covering every other file in this directory.

## Deliberate limitation

Coefficient parity is only claimed for the frozen inputs and the pinned
engine versions (R, glmnet, multiview 1.0); the tests allow a small
numerical tolerance (1e-6) for cross-platform BLAS differences while the
selection schedules must match exactly. The D5/B2 divergence is only
observable when step-halving actually triggers; where it cannot be
forced deterministically the fixtures pin current vendored behavior
instead of claiming an observed difference.
