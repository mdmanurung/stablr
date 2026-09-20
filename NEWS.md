# stablr 0.1.1

## Release correctness

- Late-fusion weights now default to leakage-safe out-of-fold training. Each
  stacking fold reruns STABL selection and its downstream learner using only
  fold-training samples; `late_fusion_training = "python_legacy"` preserves
  the historical Python-compatible in-sample algorithm.
- Named outcomes, groups, bootstrap strata, and validation outcomes are
  canonically aligned to omic sample identifiers before fusion.
- OOF provenance records named fold assignments, seeds, selected features,
  artificial-feature metadata, warnings, fallbacks, and full-refit diagnostics.
- Scalar stacking now rejects malformed outcomes and infinite observed
  predictions instead of allowing invalid metric inputs.
- Binary event levels and multiclass probability columns are now validated
  explicitly; probability rows must be complete, finite, and normalized, and
  metric calculations no longer silently discard malformed observations.

## Architecture, profiling, and calibration hardening

- **Model-X equicorrelated knockoffs:** `artificial_type = "knockoff_equi"`
  now uses a scikit-learn 1.7.1-compatible Ledoit-Wolf covariance estimator
  and knockpy 1.3.5-compatible equicorrelated S-matrix. Provenance records the
  estimator, compatibility versions, shrinkage, and S-matrix settings. Results
  from this generator can change relative to earlier 0.1.1 candidates.
  `"knockoff"` retains fixed-X semantics, while `"knockoff_mvr"` retains its
  existing covariance and solver semantics.
- **Fixed-X knockoff correctness:** fixed-X construction now rejects row-
  augmented knockoffs and uses the documented random-permutation fallback.
  Truncating only the knockoff matrix to the original row count violated the
  fixed-X Gram identities when `p < n < 2p`; fallback provenance records the
  reason.
- **FDP+ performance and contracts:** threshold counts now sort once per score
  vector and use interval lookup, preserving strict `>` tie semantics and
  first-minimum threshold selection while avoiding feature-by-threshold logical
  allocations. Public inputs now reject malformed, non-finite, out-of-range,
  or column-incompatible score matrices.
- **Learner-module depth:** public single-lambda glmnet and adaptive-lasso
  adapters now delegate to the same path-fitting implementations used by the
  batched bootstrap interface, removing duplicate fitting and coefficient-
  extraction logic without changing seeded fitted outputs.
- **Methodology runner:** artificial-feature strategies are compared on paired
  simulated data and fit seeds. Replicate artifacts record both seeds, and
  summaries report SD and Monte Carlo SE for selection count, empirical FDP,
  and TPR, plus Monte Carlo SE for FDP-exceedance and fallback rates.
- **ROC/PRC diagnostics:** `plot_roc()` and `plot_prc()` reject ambiguous or
  incomplete labels and probabilities, and aggregate tied scores so metrics
  are row-order invariant.
- **Python parity evidence:** reproducible Python 3.11 fixtures pin the public
  STABL reference commit and include frozen matrices, bootstrap schedules,
  artificial features, lambda grids, masks, FDP+ curves, stacking candidates,
  solver-ranking contracts, provenance, and a SHA-256 manifest.
- **MVR knockoffs:** solver checks now cover PSD feasibility, objective
  behaviour, determinism, and caller RNG isolation. Fallbacks and approximate
  high-dimensional chunking carry explicit provenance.
- **Release infrastructure:** expanded the R CMD check platform matrix,
  installed-package methodology runner resolution, always-uploaded check
  artifacts, one-pass Cobertura/Codecov reporting, release-gated pkgdown
  deployment, the canonical manual R-hub v2 workflow for sanitizer and
  specialist-platform checks, GPL and vendored-code provenance, and
  reproducible bundled OOL data preparation.
- **Metadata and documentation:** set the release version to 0.1.1, refreshed
  the citation and method DOIs, credited the STABL and multiview contributors,
  recorded the bundled OOL data source and CC0 terms, and replaced "Pure-R"
  with the accurate "R-native without a Python runtime" description.
  Documentation now describes FDP+ as a data-dependent diagnostic minimizer,
  states its assumptions and fallback limitations, and distinguishes frozen-
  score parity from end-to-end solver ranking parity.
- **Release evidence:** added deterministic tests for sample-order invariance,
  grouped OOF separation, validation-outcome independence, fold perturbations,
  sequential/parallel agreement, legacy stacking parity, MVR feasibility, and
  locked statistical release gates. Validation commands, versions, hashes,
  warnings, profiling limitations, and platform results are recorded in
  `docs/release/0.1.1-validation.md`.
  Locked long-running validators now require a clean Git source tree, record
  both start and end provenance, and reject evidence if the tracked tree moves
  during computation.
- **Release stop condition:** SCI-04 remains independently unaccepted. The
  review packet identifies that the draft 100-replicate Hoeffding bound cannot
  satisfy its proposed 0.10 threshold even when no features are selected. The
  package contract therefore remains non-promotable pending a reviewed,
  jointly feasible statistical decision.
- **Pre-review release controls:** added versioned SCI-04 and OOL v3 review
  bundles, detached-signature reviewer authentication, exact 1,317-feature OOL
  score alignment, independent packet verification, artifact-bound candidate
  platform receipts, and an external evidence-only promotion validator. These
  controls do not accept either scientific decision or make v0.1.1 releasable.

# stablr 0.1.0

## Publication-readiness blocker remediation

- **License:** corrected `License:` field from `MIT + file LICENSE` to `GPL (>= 2)`
  and removed the MIT `LICENSE` file. The vendored cooperative-learning engine
  under `src/glmnetpp/` is derived from glmnet/multiview (GPL-2); the compiled
  combined work must therefore be GPL. `inst/COPYING.cooperative` already
  documented the GPL-2 component.
- **Repository URLs:** updated `URL:` and `BugReports:` in DESCRIPTION and the
  package URL in `inst/CITATION` from the upstream Python repo
  (`gregbellan/Stabl`) to the canonical R package repo (`mdmanurung/stablr`).
- **Contributors:** added the original STABL method authors (Hédou, Marić, Bellan,
  Gfeller — Nat Biotechnol 2024) and the cooperative-learning/multiview authors
  (Ding, Li, Narasimhan, Tibshirani — PNAS 2022) as `ctb` in `Authors@R`,
  consistent with the already-correct `inst/CITATION` entries.
- **Cox in late fusion:** `late_fusion = TRUE` with `family = "cox"` now throws
  a clear error instead of silently scoring with R-squared and OLS predictors
  (both invalid for censored data). Early-fusion and per-omic cox remain valid.
  Mirrors the existing cooperative-fusion cox guard.
- **README:** documented the C++17 compiler / toolchain requirement (Rtools on
  Windows, Xcode CLT on macOS, gcc/clang on Linux) and added a canonical
  `remotes::install_github()` install path.
- **Vignettes:** added reader-facing `eval = FALSE` callout blocks to the three
  heavy vignettes (`stablr-multiomic`, `stablr-tcga`, `stablr-tcga-nestedcv`)
  explaining that code is illustrative and how to run it.
- **MCC metric:** added Matthews Correlation Coefficient (Gorodkin 2004 multiclass
  formula) to the list returned by `.classification_metrics()`. MCC is a standard
  metric for imbalanced classification, complementing the existing accuracy, BER,
  and macro-F1. Returns 0 for degenerate (all-one-class) predictions.
- **Elastic-net parity disclosure:** the Python-parity fixtures for elastic-net
  models assert *signal ranking* (which features are selected), not bit-identical
  coefficients. R and Python elastic-net implementations use different coordinate
  descent schedules; ranking parity is the appropriate criterion and is explicitly
  documented in the test file.
- **Verification:** `mixOmics` (Suggests) is used only in TCGA vignettes and is
  already guarded by `requireNamespace()` + global `eval = FALSE`; no change
  needed.
- **Release metadata:** aligned pkgdown URLs, article listings, and the public
  reference index with the current `mdmanurung/stablr` repository. Removed the
  stale Python-parity article link; the Python-to-R mapping now lives in
  `docs/PYTHON_TO_R_MAPPING.md`.
- **Parallel tests:** the parallel RNG determinism test now probes
  `future::multisession`/PSOCK availability and skips with a clear reason in
  socket-restricted check sandboxes.
- **Release hardening:** added shared scalar validators, duplicate sample/feature
  name guards, effective FDP+ decoy scaling, `transform_stabl()` for selected
  feature subsetting, `matrixStats`-backed row maxima, bundled Python-parity
  fixture checks, fixed-seed methodology/stacking scripts, validator rollout
  across workflow/plot/export helpers, and parity-preserving chunked
  `stacked_multi_omic()` evaluation for binary/regression tasks. The multiclass
  stacking path remains scalar by design to preserve probability-normalization
  parity.
## Quality and publication-readiness sweep

- **Infrastructure:** corrected `License:` field; lowered `Depends: R (>=
  4.1.0)` from 4.4.0; moved package-level doc to `R/stablr-package.R`; guarded
  `save_stabl_results()` against a missing ggplot2 Suggests dep; removed stray
  `Rplots.pdf` build artifacts.
- **Documentation:** documented all four `print.*` S3 methods; added `@examples`
  to all 54 exported symbols; replaced duplicated `@param` blocks with
  `@inheritParams stabl_fit`; added `@seealso` cross-links across the metrics
  family and `compute_fdp_plus`.
- **Correctness:** guarded cooperative lambda-index lookup against float
  mismatch; fixed integer overflow in `.derive_nested_seed` for large fold
  counts; added explicit `stop()` default to the cooperation `type_measure`
  switch; replaced `seq(start+1, start+n)` with `seq_len(n)` idiom at all index
  ranges; added tie-policy comments at every `which.min`/`which.max` site.
- **Performance:** centralised `rowMaxs` helper with edge-case fast paths for
  0-row and 1-column matrices; eliminated quadratic `c(accumulator, chunk)`
  vector growth in three bootstrap helpers; extracted
  `.BOOTSTRAP_COEF_THRESHOLD` named constant; halved pairwise-similarity call
  count in `jaccard_matrix`, `adjusted_similarity_values`, and
  `pearson_similarity_values` via upper-triangle loops.
- **Internal quality:** extracted `.require_ggplot2()` to replace six identical
  inline guards; collapsed `.cooperative_selected_features` two-pass loop to a
  single `lapply`; vectorised the finite-check in `.multiclass_log_loss`;
  extracted `.scores_to_long_df` to deduplicate the real- and
  artificial-feature long-data-frame loops in `plot_stabl_path`.
- **Tests:** added tests for `get_stabl_scores` and `load_ool_data`; modernised
  `skip()` to `skip_if()` idiom. Current local no-vignette release check:
  FAIL 0 | WARN 0 | SKIP 1 | PASS 1611.
  No-manual/no-vignette `R CMD check` is `Status: OK`; all six source
  vignettes rendered to `/tmp/stablr-vignette-review`; pkgdown built to
  `/tmp/stablr-pkgdown`.
- **CI:** added `.github/workflows/R-CMD-check.yaml` (ubuntu release+devel,
  macOS release) with a separate covr/Codecov coverage job; gated heavy
  vignettes (`stablr-multiomic`, `stablr-tcga`, `stablr-tcga-nestedcv`) with
  `eval = FALSE` in their setup chunks.

## Native cooperative learning (multiview engine)

- Vendored CRAN `multiview` v1.0 cooperative-learning engine (GPL-2 submodule in
  `inst/COPYING.cooperative`) for gaussian and binomial cooperative fusion.
- No runtime dependency on the external `multiview` package for cooperative
  workflows.
- Compiled C++ core (`wls_exp`, glmnetpp headers) with `NeedsCompilation: yes`.
- New `vignette("stablr-advanced")` documents Cox, knockoffs, grouped bootstrap,
  reproducibility metrics, export helpers, and outer CV.

## Application-note release candidate

- R-native STABL implementation without a Python runtime, with glmnet backends (lasso, elastic net,
  adaptive lasso, optional sparse group lasso).
- FDP+ thresholding with random-permutation and knockoff artificial features.
- Multi-omic train/validation, outer CV, and nested CV workflows.
- Early, late, and cooperative fusion integration paths.
- Bundled Onset of Labor (OOL) tutorial subset in `inst/extdata/`.
- Python parity fixtures and vignettes for cross-language validation.
- Publication reproduction script: `inst/analysis/generate_publication_parity_table.R`.

## Naming

The package is **`stablr`**. Cooperative multi-omic fusion is a workflow inside
`stablr`, not a separate package (`costablr` is not used).
