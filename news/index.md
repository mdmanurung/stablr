# Changelog

## stablr 0.3.0

### Cross-study reliability: study-disjoint resampling (B1–B3)

- **New `study_labels` argument on
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)**
  declares the sampling unit of the cross-validation. When supplied (a
  named vector covering exactly the sample IDs), whole studies become
  the CV unit: outer folds assign whole studies (`outer_v` equal to the
  number of studies gives leave-one-study-out), inner
  candidate-selection folds are study-disjoint within each outer
  training split, and `stratified = TRUE` is reinterpreted at the study
  level through a deterministic largest-first greedy pass that minimises
  class-proportion imbalance. The default `study_labels = NULL`
  preserves the sample-i.i.d. scheme byte-identically
  (`fold_scheme = "sample_iid"`). Study disjointness concerns the CV
  folds only — `bootstrap_strata` (sampling design) and SGL
  `feature_groups` (features) remain orthogonal.
- **Fail-closed validation** for the study-disjoint scheme: label names
  must cover exactly the sample IDs; at least `max(outer_v, 2)` studies;
  each outer-training split needs at least `inner_v + 1` studies; every
  validation fold (outer and inner) must hold both outcome classes; and
  1-sample studies are rejected.
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
  also gained a majority-class fallback when a fold’s selection leaves
  fewer than two features (previously a glmnet error).
- **Provenance and reporting:** the nested-CV object records
  `fold_scheme` (`"sample_iid"` / `"study_disjoint"`) and the validated
  `study_labels`; under the study-disjoint scheme each fold carries
  `study_ids` and held-out `class_counts`, and
  [`print()`](https://rdrr.io/r/base/print.html) shows the scheme.
- **New
  [`stabl_multiomic_study_diagnostics()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_study_diagnostics.md)**
  (inferred scope, flagged): per-study outcome composition and
  per-feature absolute standardized mean differences versus the pooled
  rest, quantifying the study-level shift that motivates study-disjoint
  folds.
- **New
  [`stabl_multiomic_study_split()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_study_split.md)**
  (inferred scope, flagged): a deterministic study-disjoint
  train/validation split (greedy, seeded, class-balanced) feeding
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).
  [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  already supports group-disjoint assessment folds through its `groups`
  argument and now records `fold_scheme` = `"group_disjoint"` for
  provenance.
- **Benchmark evidence:** `inst/benchmarks/b1_study_disjoint_cv.R`
  compares held-out BER under sample-i.i.d. versus study-disjoint outer
  folds on multi-study simulations, isolating two batch mechanisms:
  class-independent location shifts on null features (a deliberate null
  probe — the gap is approximately zero there, as expected) and
  outcome-prevalence drift across studies (gap ≈ +0.04 BER over 10
  seeds, positive in 9 of 10 seeds — the sample-i.i.d. estimate is
  optimistic precisely when study-level outcome composition differs).
  Per-fold selection sets are persisted (`b1_selected_features.csv`) for
  the Round C stability-of-stability analysis.

### Cross-study reliability: reporting, comparison, and guardrails (C1–C5)

All Round C items are inferred scope (the approved plan flags them for
audit reconciliation):

- **Per-study performance (C1):**
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
  now returns `performance_by_study` when `study_labels` is supplied —
  per-study outer held-out accuracy and balanced error rate pooled
  across folds and repeats (`NA` balanced error rate for studies whose
  pooled held-out predictions contain a single outcome class). `NULL`
  under the sample-i.i.d. scheme.
- **Scheme comparison helper (C2):** new
  [`stabl_multiomic_scheme_comparison()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_scheme_comparison.md)
  runs the same nested CV under both sampling units (sample-i.i.d. and
  study-disjoint) with matched arguments and returns both objects plus a
  comparison table and the BER gap (positive = the sample-i.i.d.
  estimate was optimistic).
- **Harmonization guardrail (C3, documentation):** the CV docs now warn
  that batch harmonization (e.g. ComBat) fitted on the full dataset
  before cross-validation leaks study-level information across folds;
  fit harmonization within training folds or interpret metrics as
  conditional on that step. Study-disjoint folds do not repair leakage
  from pooled preprocessing.
- **Stability of stability (C4):** new
  [`stabl_multiomic_selection_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_selection_stability.md)
  computes fold x fold Jaccard similarity of the selected feature sets
  across the outer folds of a study-disjoint nested CV (the deferred
  selection-stability heatmap analysis; render with
  [`stats::heatmap()`](https://rdrr.io/r/stats/heatmap.html)).
- **Grouped-CV guardrail (C5):**
  [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  now fails closed when a grouped assessment fold holds a single outcome
  class under a classification family (mirroring the study-disjoint
  validation); regression is unaffected.

### Two-tier API: beginner `stabl()` + expert `stabl_fit()`

- **New
  [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md):**
  a one-stop-shop entry point for first-time and application users. Pass
  a `data.frame` (or a matrix plus an outcome vector; the outcome may
  also be named as a column of the `data.frame`) and it runs STABL with
  statistically safe defaults — model-X knockoff artificial features
  (`artificial_type = "knockoff_equi"`), FDP+ calibration, strict
  knockoff fallback (`knockoff_fallback = "stop"`), and an automatic
  lambda grid. Roughly eight arguments cover the common cases;
  everything else is forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  through `...`.
- **New S3 methods for the `stabl` class:**
  [`print.stabl()`](https://mdmanurung.github.io/stablr/reference/print.stabl.md)
  (plain-language answer: what was fitted, which features were selected,
  the achieved FDP+ estimate),
  [`summary.stabl()`](https://mdmanurung.github.io/stablr/reference/summary.stabl.md),
  and
  [`plot.stabl()`](https://mdmanurung.github.io/stablr/reference/plot.stabl.md)
  (stability paths plus the FDP+ curve). `stabl` objects also remain
  full `stabl_fit` objects, so every accessor and workflow applies
  unchanged.
- **New beginner generics
  [`support()`](https://mdmanurung.github.io/stablr/reference/support.md)
  and
  [`feature_names()`](https://mdmanurung.github.io/stablr/reference/support.md)**
  with methods for `stabl` and `stabl_fit`; they delegate to
  [`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
  and
  [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md).

### Renamed functions and arguments (old names kept as deprecation shims)

Each old name still works: it warns once per session and delegates to
the new name. Positional calls are unaffected.

| Old name | New name | Why |
|----|----|----|
| `stabl_fit(groups =)` | `stabl_fit(sample_groups =)` | `groups` collided with `feature_groups`; sample- vs feature-grouping is now unambiguous |
| [`get_feature_names_out()`](https://mdmanurung.github.io/stablr/reference/get_feature_names_out.md) | [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md) | drops the sklearn-ism; reads naturally in R |
| [`stacked_multi_omic()`](https://mdmanurung.github.io/stablr/reference/stacked_multi_omic.md) | [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md) | joins the `stabl_multiomic_*` family and fixes the multi_omic/multiomic casing split |
| [`jaccard_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_matrix.md) | [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md) | states what the matrix contains |
| [`adjusted_similarity_values()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_values.md) | [`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md) | the triplet now reads pair / all-pairs matrix / summary consistently |
| [`adjusted_similarity_measure()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_measure.md) | [`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md) | as above |
| [`pearson_similarity_values()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_values.md) | [`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md) | as above |
| [`pearson_similarity_measure()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_measure.md) | [`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md) | as above |
| [`classic_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/classic_bootstrap_indices.md) | [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md) | states what is resampled (samples, not features or groups) |
| `pearson_similarity_pairwise(d =)` / `_summary(d =)` | `nb_total_elements =` | aligns the element-count argument with the adjusted-similarity family (`d =` still accepted with a deprecation warning) |

### Documentation

- pkgdown reference index reorganized into tiers: Getting started, Core,
  Workflows, Expert, Metrics, Plotting, Data, and a dedicated Deprecated
  section.
- New vignette “Getting started with stabl” — a beginner walkthrough
  from a `data.frame` to a plain-language selection.
- README gains a two-tier quick start (beginner
  [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
  first, expert
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  second).
- All help pages regenerated; every shim’s help page points to its
  replacement.

## stablr 0.2.2

### Scientific hardening (audit A1-A3, E1-E3)

- **Strict fallback policy (flagged behavior change):** when a knockoff
  constructor fails on a supported-dimensional input, the previous
  behavior was to warn and silently fall back to random-permutation
  features — which do not support FDP+ calibration, so the fit’s FDP+
  estimates silently changed target.
  [`make_artificial_features()`](https://mdmanurung.github.io/stablr/reference/make_artificial_features.md)
  and
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  now take `knockoff_fallback = "stop"` (the default): the fallback is
  an error with an actionable message. `knockoff_fallback = "allow"`
  restores the warn-and-fall-back behavior (recorded in the provenance
  `fallback_history`). The benign `knockoff_mvr` -\> `knockoff_equi`
  fallback is unaffected (both are model-X constructions with valid FDP+
  calibration) and remains a warning under either policy. Related fix:
  the linear-algebra infeasibility classifier now recognizes the
  hyphenated “not positive-definite” and “cannot solve SDP” messages
  emitted by the knockoff package’s SDP solver, so those failures are
  typed and handled consistently instead of escaping as raw errors.
- **RF/oblique-RF confirmatory guard extended:** the runtime guardrail
  for multivariate RF learners now also fires for the fixed-X `knockoff`
  pairing (fixed-X exchangeability does not survive bootstrap row
  subsampling), not only for `random_permutation`. The model-X pairings
  (`knockoff_equi`, `knockoff_mvr`) remain warning-free.
- **Covariance terminology:** the fixed-X knockoff documentation now
  states the actual contract (second-moment/Gram structure on the input
  scale) rather than “covariance structure”.
- **SCI-04 feasibility clause (E1):** the draft
  `SCI-04-null-selected-fraction-v2` contract now documents that the
  proposed release criterion `U <= 0.10` (Hoeffding, `alpha = 0.05`) is
  vacuous at the previously considered n = 100 replicates (margin 0.1224
  \> 0.10) and predeclares `n >= 150` as a design-time fail-closed
  feasibility requirement, with deterministic reference tests pinning
  the arithmetic. The gate remains `proposed_unaccepted` pending
  independent review; no locked-run data was used.
- **F5 benchmark expansion (E2/E3):** the RF validation suite now
  includes the oblique-RF adapter (`aorsf`) when installed, giving the
  adapter its first F5 FDP evidence, and the start-of-run knockoff probe
  verifies the 0.2.2 shipped-pair scale contract on a
  heterogeneous-scale design (column sds 0.01-100) instead of a
  standardized one.

### Audit hotfixes (source audit of 0.2.1, PR [\#3](https://github.com/mdmanurung/stablr/issues/3) follow-up)

- **Fixed-X shipped-pair scale contract completed (flagged behavior
  change):** 0.2.1 rescaled the fixed-X knockoffs to the *standardized*
  design but still shipped them next to the caller’s *raw* matrix — a
  mismatched pair whenever input columns were heterogeneous. The
  construction now runs on the raw design (create.fixed normalizes
  internally, scale only) and the knockoffs are back-transformed by the
  input column norms, so the shipped pair `(x, xk)` satisfies the
  fixed-X Gram condition on the input’s own scale: `t(xk) xk = t(x) x`
  exactly and `t(x) xk = t(x) x - diag(s)` with `s_j >= 0`. The real
  block is untouched; the pre-construction centering/scaling step is
  removed. Gram/cross-Gram identities and paired-block scale matching
  are now asserted on the actual `x_augmented` matrix under
  heterogeneous scales (column sds 0.01-100). Provenance field
  `knockoff_rescaled_to_standardized_scale` is replaced by
  `knockoff_rescaled_to_input_scale`.
- **Windows CI fixtures:** `.gitattributes` now enforces LF (and marks
  `.rds` binary) for the `multiview_parity` fixtures; the Windows
  `R-CMD-check` manifest failures were checkout line-ending conversion,
  not hash drift.
- **pkgdown index:** the `stablr-environments` vignette and the
  `print.stabl_env_stability` topic are indexed; the site builds clean.
- **Oblique-RF tests:** a stale probe (measuring `create.fixed`’s
  internal unit-norm scale against a sd ~ 1 expectation, and gating a
  *model-X* assertion on a *fixed-X* build check) was removed; the
  `knockoff_equi` signal-recovery assertion now runs wherever the
  knockoff package is installed.
- **0.2.1 claim corrections:** the F5 measurement of empirical FDP 0.665
  was for lasso + *model-X equicorrelated* knockoffs (not fixed-X), and
  the “invisible to the learner” causal statement is moderated: for the
  lasso backend, glmnet’s internal standardization moderates the
  selection-level impact of a scale mismatch. The contract violation
  remains the reason for the fix; the 0.2.1 notes above are corrected in
  place.

## stablr 0.2.1

### Fixed-X knockoff scale and alignment fix (flagged behavior change)

- **Scale contract:**
  [`knockoff::create.fixed()`](https://rdrr.io/pkg/knockoff/man/create.fixed.html)
  column-normalizes its input and returns knockoffs on that internal
  scale (entries ~ 1/√n). The fixed-X path previously paired these
  knockoffs with the sd-1 standardized design, violating the fixed-X
  construction contract (the paired blocks lived on different scales).
  Knockoffs are now rescaled so that `t(Xk) Xk = t(X) X` holds on the
  paired scale. (0.2.2 correction: the F5 measurement of empirical FDP
  0.665 cited here was for lasso + *model-X equicorrelated* knockoffs,
  not fixed-X; and for the lasso backend the selection-level impact of a
  scale mismatch is moderated by `glmnet`’s internal standardization —
  the contract violation, not a measured FDP figure, is the reason for
  this fix.) Fits using fixed-X knockoffs
  (`artificial_type = "knockoff"`) may select different features, and
  FDP+ estimates on such fits change.
- **Column alignment:** artificial features are now shipped in aligned
  order (injected column k is the artificial counterpart of feature k’s
  selected subset). The previous random column permutation broke the
  knockoff swap-symmetry in the bootstrap design, because stability
  scores are accumulated positionally. Applies to the fixed-X and
  model-X knockoff types; random-permutation nulls are unaffected.
- **Probes and tests:** the knockoff sanity probe (tests and F5
  benchmark) now probes the shipped contract (rescaled scale + Gram
  identity) instead of the internal normalized scale; the formerly
  CI-gated knockoff quality and calibration assertions now run wherever
  the knockoff package is available. The end-to-end calibration test was
  redesigned to the bootstrap-valid model-X regime (`knockoff_equi`)
  with an explicit lambda grid and powered signals; fixed-X knockoffs
  under row subsampling remain documented as an invalid combination and
  are not asserted.

### New experimental oblique random-forest learner (plan item D2)

- `base_learner = "oblique_rf"` via
  [`make_aorsf_adapter()`](https://mdmanurung.github.io/stablr/reference/make_aorsf_adapter.md)
  /
  [`auto_oblique_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_oblique_threshold_grid.md)
  (aorsf backend, negation importance), mirroring the random-forest
  adapter: importance-threshold penalty path, per-tree seeds drawn from
  the caller’s seeded RNG stream, and the same model-X-knockoff-pairing
  guardrail. As with `rf`, selection frequencies are noisier than the
  glmnet family’s; the adapter is experimental and FDP+ calibration is
  empirical (F5-style benchmarks). aorsf is a Suggests dependency; the
  feature is guarded by `.require_pkg()`.

### Multiview-parity fixtures (plan item F4)

- SHA-256-pinned fixtures (`tests/testthat/fixtures/multiview_parity/`)
  lock the vendored cooperative-learning backend’s outputs to CRAN
  multiview 1.0 on gaussian, binomial, and weighted fits, plus the
  documented B1/B2/B3 divergence sites (offset handling in
  `build_predmat`, warm-start residual weighting, `lower.limits`
  clamping). Upstream `build_predmat` errors with an offset (divergence
  D4 confirmed live); the vendored fix is pinned by fixture. Makes the
  PATCHES.cooperative.md divergence log testable.

## stablr 0.2.0

Breaking renames ship with deprecation shims that warn for one minor
release before removal.

### Scientific correctness

- **Behavior change (flagged):** the fixed-X knockoff path now
  standardizes the input matrix before calling
  [`knockoff::create.fixed()`](https://rdrr.io/pkg/knockoff/man/create.fixed.html),
  matching that function’s documented contract. Fits on heterogeneous
  feature scales may select different features. Standardization details
  (means, sds, zero- variance handling) are recorded in the fit
  provenance, and the fixed-X documentation now states the
  bootstrap-reuse caveat.
- **FDP+ communication:** new “What FDP+ is and is not” vignette section
  and documented interpretation caveats on
  [`get_fdrs()`](https://mdmanurung.github.io/stablr/reference/get_fdrs.md).

### Vendored cooperative-backend fixes

- **Offset subsetting in `build_predmat`:** the cross-validation
  prediction path now subsets the `offset` argument per fold instead of
  referencing an undefined symbol. The path was previously unreachable
  because offsets are not exposed through the public API.
- **Halved-step residual weighting:** after step-halving, the warm-start
  residual now weights the offset consistently with the non-halved
  branch.
- **`lower.limits` symmetry:** the cooperative path now passes
  `lower.limits` to `elnet.fit()` alongside `upper.limits`.
- **Divergence log:** `inst/PATCHES.cooperative.md` records every
  deviation from CRAN multiview 1.0.

All backend fixes are inherited from CRAN multiview 1.0 and are
behavior-neutral under unit weights, the only weighting stablr exercises
today.

### API consistency

- **`seed` deprecated:**
  [`classic_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/classic_bootstrap_indices.md)
  and
  [`group_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/group_bootstrap_indices.md)
  now take `random_state`; `seed` still works but warns, and will be
  removed in the following minor release. The helpers no longer clobber
  the global RNG state (`withr` promoted from Suggests to Imports).
- **Similarity-metric argument renames:**
  [`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md)
  and
  [`pearson_similarity()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity.md)
  (and their `_values`/`_measure` variants) rename their arguments to
  `set1`/`set2`/`n_universe`; old names keep working through deprecation
  shims.
- **New accessors:**
  [`get_fdrs()`](https://mdmanurung.github.io/stablr/reference/get_fdrs.md)
  and
  [`get_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/get_lambda_grid.md)
  expose fitted FDP+ estimates and the lambda grid, with documented
  interpretation caveats.
- **New S3 methods:**
  [`plot.stabl_fit()`](https://mdmanurung.github.io/stablr/reference/plot.stabl_fit.md)
  and
  [`summary.stabl_fit()`](https://mdmanurung.github.io/stablr/reference/summary.stabl_fit.md)
  (with
  [`print.summary.stabl_fit()`](https://mdmanurung.github.io/stablr/reference/summary.stabl_fit.md))
  for fitted STABL objects.
- **Internal:** late-fusion internal formals are dot-prefixed so user
  arguments passed through `...` cannot collide with them; no
  user-visible impact.

### New capabilities (extensions beyond Python STABL, labeled experimental)

- **Random-forest base learner:** `base_learner = "rf"` via
  [`make_ranger_adapter()`](https://mdmanurung.github.io/stablr/reference/make_ranger_adapter.md)
  and
  [`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md)
  — permutation- importance threshold paths so the full STABL machinery
  runs unchanged on non-linear signals. Random forests must be paired
  with model-X knockoff artificial features; the anti-conservative
  random-permutation pairing now warns. Calibration is validated
  empirically by the F5 benchmark
  (`inst/benchmarks/f5_rf_validation.R`).
- **Environment-stability workflow:**
  [`stabl_environment_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_environment_stability.md)
  fits STABL per environment (batch, cohort, site) and aggregates
  selection frequencies by their minimum across environments — an
  ICP-inspired transportability heuristic without causal assumptions or
  new dependencies — with
  [`get_min_scores()`](https://mdmanurung.github.io/stablr/reference/get_min_scores.md),
  a print method, and a new vignette.

### Tests and documentation

- Golden-path rendering tests for all four plotting functions (plot-data
  correctness: rows, mappings, threshold geometries, FDP+ curve values).
- Lambda-grid leakage guidance for cross-validated multi-omic workflows,
  a note on the classification-bootstrap subsample deviation, dead-code
  removal in the MVR knockoff solver, and the release-evidence glossary
  moved from the repository root to `docs/CONTEXT.md`.

## stablr 0.1.1

### Release correctness

- Late-fusion weights now default to leakage-safe out-of-fold training.
  Each stacking fold reruns STABL selection and its downstream learner
  using only fold-training samples;
  `late_fusion_training = "python_legacy"` preserves the historical
  Python-compatible in-sample algorithm.
- Named outcomes, groups, bootstrap strata, and validation outcomes are
  canonically aligned to omic sample identifiers before fusion.
- OOF provenance records named fold assignments, seeds, selected
  features, artificial-feature metadata, warnings, fallbacks, and
  full-refit diagnostics.
- Scalar stacking now rejects malformed outcomes and infinite observed
  predictions instead of allowing invalid metric inputs.
- Binary event levels and multiclass probability columns are now
  validated explicitly; probability rows must be complete, finite, and
  normalized, and metric calculations no longer silently discard
  malformed observations.

### Architecture, profiling, and calibration hardening

- **Fixed-X knockoff correctness:** fixed-X construction now rejects
  row- augmented knockoffs and uses the documented random-permutation
  fallback. Truncating only the knockoff matrix to the original row
  count violated the fixed-X Gram identities when `p < n < 2p`; fallback
  provenance records the reason.
- **FDP+ performance and contracts:** threshold counts now sort once per
  score vector and use interval lookup, preserving strict `>` tie
  semantics and first-minimum threshold selection while avoiding
  feature-by-threshold logical allocations. Public inputs now reject
  malformed, non-finite, out-of-range, or column-incompatible score
  matrices.
- **Learner-module depth:** public single-lambda glmnet and
  adaptive-lasso adapters now delegate to the same path-fitting
  implementations used by the batched bootstrap interface, removing
  duplicate fitting and coefficient- extraction logic without changing
  seeded fitted outputs.
- **Methodology runner:** artificial-feature strategies are compared on
  paired simulated data and fit seeds. Replicate artifacts record both
  seeds, and summaries report SD and Monte Carlo SE for selection count,
  empirical FDP, and TPR, plus Monte Carlo SE for FDP-exceedance and
  fallback rates.
- **ROC/PRC diagnostics:**
  [`plot_roc()`](https://mdmanurung.github.io/stablr/reference/plot_roc.md)
  and
  [`plot_prc()`](https://mdmanurung.github.io/stablr/reference/plot_prc.md)
  reject ambiguous or incomplete labels and probabilities, and aggregate
  tied scores so metrics are row-order invariant.
- **Python parity evidence:** reproducible Python 3.11 fixtures pin the
  public STABL reference commit and include frozen matrices, bootstrap
  schedules, artificial features, lambda grids, masks, FDP+ curves,
  stacking candidates, solver-ranking contracts, provenance, and a
  SHA-256 manifest.
- **MVR knockoffs:** solver checks now cover PSD feasibility, objective
  behaviour, determinism, and caller RNG isolation. Fallbacks and
  approximate high-dimensional chunking carry explicit provenance.
- **Release infrastructure:** expanded the R CMD check platform matrix,
  installed-package methodology runner resolution, always-uploaded check
  artifacts, one-pass Cobertura/Codecov reporting, release-gated pkgdown
  deployment, the canonical manual R-hub v2 workflow for sanitizer and
  specialist-platform checks, GPL and vendored-code provenance, and
  reproducible bundled OOL data preparation.
- **Metadata and documentation:** set the release version to 0.1.1,
  refreshed the citation and method DOIs, credited the STABL and
  multiview contributors, recorded the bundled OOL data source and CC0
  terms, and replaced “Pure-R” with the accurate “R-native without a
  Python runtime” description. Documentation now describes FDP+ as a
  data-dependent diagnostic minimizer, states its assumptions and
  fallback limitations, and distinguishes frozen- score parity from
  end-to-end solver ranking parity.
- **Release evidence:** added deterministic tests for sample-order
  invariance, grouped OOF separation, validation-outcome independence,
  fold perturbations, sequential/parallel agreement, legacy stacking
  parity, MVR feasibility, and locked statistical release gates.
  Validation commands, versions, hashes, warnings, profiling
  limitations, and platform results are recorded in
  `docs/release/0.1.1-validation.md`. Locked long-running validators now
  require a clean Git source tree, record both start and end provenance,
  and reject evidence if the tracked tree moves during computation.

## stablr 0.1.0

### Publication-readiness blocker remediation

- **License:** corrected `License:` field from `MIT + file LICENSE` to
  `GPL (>= 2)` and removed the MIT `LICENSE` file. The vendored
  cooperative-learning engine under `src/glmnetpp/` is derived from
  glmnet/multiview (GPL-2); the compiled combined work must therefore be
  GPL. `inst/COPYING.cooperative` already documented the GPL-2
  component.

- **Repository URLs:** updated `URL:` and `BugReports:` in DESCRIPTION
  and the package URL in `inst/CITATION` from the upstream Python repo
  (`gregbellan/Stabl`) to the canonical R package repo
  (`mdmanurung/stablr`).

- **Contributors:** added the original STABL method authors (Hédou,
  Marić, Bellan, Gfeller — Nat Biotechnol 2024) and the
  cooperative-learning/multiview authors (Ding, Li, Narasimhan,
  Tibshirani — PNAS 2022) as `ctb` in `Authors@R`, consistent with the
  already-correct `inst/CITATION` entries.

- **Cox in late fusion:** `late_fusion = TRUE` with `family = "cox"` now
  throws a clear error instead of silently scoring with R-squared and
  OLS predictors (both invalid for censored data). Early-fusion and
  per-omic cox remain valid. Mirrors the existing cooperative-fusion cox
  guard.

- **README:** documented the C++17 compiler / toolchain requirement
  (Rtools on Windows, Xcode CLT on macOS, gcc/clang on Linux) and added
  a canonical `remotes::install_github()` install path.

- **Vignettes:** added reader-facing `eval = FALSE` callout blocks to
  the three heavy vignettes (`stablr-multiomic`, `stablr-tcga`,
  `stablr-tcga-nestedcv`) explaining that code is illustrative and how
  to run it.

- **MCC metric:** added Matthews Correlation Coefficient (Gorodkin 2004
  multiclass formula) to the list returned by
  `.classification_metrics()`. MCC is a standard metric for imbalanced
  classification, complementing the existing accuracy, BER, and
  macro-F1. Returns 0 for degenerate (all-one-class) predictions.

- **Elastic-net parity disclosure:** the Python-parity fixtures for
  elastic-net models assert *signal ranking* (which features are
  selected), not bit-identical coefficients. R and Python elastic-net
  implementations use different coordinate descent schedules; ranking
  parity is the appropriate criterion and is explicitly documented in
  the test file.

- **Verification:** `mixOmics` (Suggests) is used only in TCGA vignettes
  and is already guarded by
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) + global
  `eval = FALSE`; no change needed.

- **Release metadata:** aligned pkgdown URLs, article listings, and the
  public reference index with the current `mdmanurung/stablr`
  repository. Removed the stale Python-parity article link; the
  Python-to-R mapping now lives in `docs/PYTHON_TO_R_MAPPING.md`.

- **Parallel tests:** the parallel RNG determinism test now probes
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html)/PSOCK
  availability and skips with a clear reason in socket-restricted check
  sandboxes.

- **Release hardening:** added shared scalar validators, duplicate
  sample/feature name guards, effective FDP+ decoy scaling,
  [`transform_stabl()`](https://mdmanurung.github.io/stablr/reference/transform_stabl.md)
  for selected feature subsetting, `matrixStats`-backed row maxima,
  bundled Python-parity fixture checks, fixed-seed methodology/stacking
  scripts, validator rollout across workflow/plot/export helpers, and
  parity-preserving chunked
  [`stacked_multi_omic()`](https://mdmanurung.github.io/stablr/reference/stacked_multi_omic.md)
  evaluation for binary/regression tasks. The multiclass stacking path
  remains scalar by design to preserve probability-normalization parity.
  \## Quality and publication-readiness sweep

- **Infrastructure:** corrected `License:` field; lowered
  `Depends: R (>= 4.1.0)` from 4.4.0; moved package-level doc to
  `R/stablr-package.R`; guarded
  [`save_stabl_results()`](https://mdmanurung.github.io/stablr/reference/save_stabl_results.md)
  against a missing ggplot2 Suggests dep; removed stray `Rplots.pdf`
  build artifacts.

- **Documentation:** documented all four `print.*` S3 methods; added
  `@examples` to all 54 exported symbols; replaced duplicated `@param`
  blocks with `@inheritParams stabl_fit`; added `@seealso` cross-links
  across the metrics family and `compute_fdp_plus`.

- **Correctness:** guarded cooperative lambda-index lookup against float
  mismatch; fixed integer overflow in `.derive_nested_seed` for large
  fold counts; added explicit
  [`stop()`](https://rdrr.io/r/base/stop.html) default to the
  cooperation `type_measure` switch; replaced `seq(start+1, start+n)`
  with `seq_len(n)` idiom at all index ranges; added tie-policy comments
  at every `which.min`/`which.max` site.

- **Performance:** centralised `rowMaxs` helper with edge-case fast
  paths for 0-row and 1-column matrices; eliminated quadratic
  `c(accumulator, chunk)` vector growth in three bootstrap helpers;
  extracted `.BOOTSTRAP_COEF_THRESHOLD` named constant; halved
  pairwise-similarity call count in `jaccard_matrix`,
  `adjusted_similarity_values`, and `pearson_similarity_values` via
  upper-triangle loops.

- **Internal quality:** extracted `.require_ggplot2()` to replace six
  identical inline guards; collapsed `.cooperative_selected_features`
  two-pass loop to a single `lapply`; vectorised the finite-check in
  `.multiclass_log_loss`; extracted `.scores_to_long_df` to deduplicate
  the real- and artificial-feature long-data-frame loops in
  `plot_stabl_path`.

- **Tests:** added tests for `get_stabl_scores` and `load_ool_data`;
  modernised `skip()` to `skip_if()` idiom. Current local no-vignette
  release check: FAIL 0 \| WARN 0 \| SKIP 1 \| PASS 1611.
  No-manual/no-vignette `R CMD check` is `Status: OK`; all six source
  vignettes rendered to `/tmp/stablr-vignette-review`; pkgdown built to
  `/tmp/stablr-pkgdown`.

- **CI:** added `.github/workflows/R-CMD-check.yaml` (ubuntu
  release+devel, macOS release) with a separate covr/Codecov coverage
  job; gated heavy vignettes (`stablr-multiomic`, `stablr-tcga`,
  `stablr-tcga-nestedcv`) with `eval = FALSE` in their setup chunks.

### Native cooperative learning (multiview engine)

- Vendored CRAN `multiview` v1.0 cooperative-learning engine (GPL-2
  submodule in `inst/COPYING.cooperative`) for gaussian and binomial
  cooperative fusion.
- No runtime dependency on the external `multiview` package for
  cooperative workflows.
- Compiled C++ core (`wls_exp`, glmnetpp headers) with
  `NeedsCompilation: yes`.
- New
  [`vignette("stablr-advanced")`](https://mdmanurung.github.io/stablr/articles/stablr-advanced.md)
  documents Cox, knockoffs, grouped bootstrap, reproducibility metrics,
  export helpers, and outer CV.

### Application-note release candidate

- R-native STABL implementation without a Python runtime, with glmnet
  backends (lasso, elastic net, adaptive lasso, optional sparse group
  lasso).
- FDP+ thresholding with random-permutation and knockoff artificial
  features.
- Multi-omic train/validation, outer CV, and nested CV workflows.
- Early, late, and cooperative fusion integration paths.
- Bundled Onset of Labor (OOL) tutorial subset in `inst/extdata/`.
- Python parity fixtures and vignettes for cross-language validation.
- Publication reproduction script:
  `inst/analysis/generate_publication_parity_table.R`.

### Naming

The package is **`stablr`**. Cooperative multi-omic fusion is a workflow
inside `stablr`, not a separate package (`costablr` is not used).
