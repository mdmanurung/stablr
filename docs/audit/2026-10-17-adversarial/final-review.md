# Final pre-publication per-function review (2026-10-18)

Scope: **every exported function of stablr** (74 exports resolved from the
hand-maintained `NAMESPACE` against `asNamespace("stablr")`), audited against
four questions: (1) does it do what the documentation says, (2) is it
implemented well, (3) is it scientifically correct, (4) is it robust; plus
argument/defaults calibration and ensemble (cross-function) coherence. This
is the last review-and-revision pass before publication.

Method: four delegated per-function slices (core engine/FDP/bootstrap/
accessors; learner adapters/grids/preprocessing; multiomic workflows/
validation/interpretation/data; environment stability/similarity/plots/
export), each doing a full source read against the man pages with runtime
probes; the two remaining exports (`get_cooperative_features`,
`get_cooperative_diagnostics`) were audited inline with a runtime probe.
All findings below are NEW relative to the existing ledger (F-001..F-026,
F-CV-01..08, P4-F24..29, F-P6-1..3, F-P7-1..3, F-P10-1..2, A1-A6); where a
finding is a known defect class at a new location this is stated. Slice
working notes (per-function verdict tables, probe scripts incl.
`probe_D.R`) are preserved in the session audit workspace; the probe-backed
claims are summarised here.

## Per-function verdict rollup (74 exports)

| Slice | Functions | OK | OK-with-notes | issue |
|---|---|---|---|---|
| A — core engine, FDP+, bootstrap, accessors | 18 | 11 | 5 | 2 |
| B — learner adapters, grids, preprocessing | 14 | 1 | 9 | 4 |
| C — multiomic workflows, validation, interpretation, data | 14 | 4 | 8 | 2 |
| D — env stability, similarity family, plots, export | 26 | 11 | 11 | 4 |
| inline — cooperative-fusion accessors | 2 | 1 | 1 | 0 |
| **total** | **74** | **28** | **34** | **12** |

"issue" = at least one defect requiring code or doc correction; "OK-with-
notes" = correct but carries a P3/P4 note or a known ledger finding.
Scientific-correctness spot-verifications that PASSED: `compute_fdp_plus`
(formula, strict-`>` counting, argmin tie-break, `fdrs_table` dims),
`make_artificial_features` (provenance/fallback), adjusted-similarity
boundedness (−1, 1], `pearson_similarity` = phi coefficient, F-beta/TPR
denominators, ROC/PRC conventions, all 5 deprecation shims (warn + delegate
identically), leakage-safe OOF late fusion, honest evidence labeling in
`stabl_analyze`/`stabl_interpret`, fail-closed study splits and cooperative
accessors.

## New findings and dispositions

Severities: P1 = default-path failure; P3 = misleading doc or silent
wrong-output risk; P4 = hygiene/robustness.

### Fixed in this commit (code + tests)

| ID | Sev | Where | Defect | Fix |
|---|---|---|---|---|
| B-01 | P1 | `R/learner_rf.R` | `auto_rf_threshold_grid()` default `seed=NULL` always errored: `.rf_fit_permutation_importance` computed `seed + (attempt-1L)*104729L` → `numeric(0)` → `ranger(seed=numeric(0))` "Expecting a single value: [extent=0]"; propagated to `stabl_fit(base_learner="rf", lambda_grid="auto")` without `random_state` (stabl_fit.R:309). Tests always passed explicit seeds, so untested. | Resolve `seed <- seed %||% sample.int(.Machine$integer.max, 1L)` in `.rf_fit_permutation_importance` (honours the documented "seeded from the caller's RNG stream" contract; reproducible under `random_state`); regression test with default seed. |
| B-02 | P3 | `R/learner_rf.R` | Retry-failure returned NULL importance → `.rf_threshold_grid_from_importance` silently produced a 0-row grid (batch adapter stops; grid path did not). | Fail-closed `stop()` in `auto_rf_threshold_grid()` when importance is NULL. |
| N-02 | P3 | `R/stabl_accessors.R` | `get_support()`/`get_feature_names()`/`transform_stabl()` never range-validated `new_hard_threshold` despite documented (0, 1]: `1.5` silently selected nothing, `-1` everything (probe-verified). | `.validate_proportion(new_hard_threshold, ...)` at the `get_support.stabl_fit` choke point (covers the accessor family); boundary tests. |
| D-1 | P3 | `R/visualization.R` | `boxplot_features`/`scatterplot_features` did not enforce documented `length(y) == nrow(x)`; divisible mismatches silently recycled labels into a wrong figure (probe: 40-row x with 20 labels plotted silently). | Length check with informative error in both functions; regression tests. |
| C-3 | P3 | `R/multiomic_workflows.R` | `stabl_multiomic_train_validate(y_valid=...)` without `x_valid_list` was accepted and silently ignored (probe-verified). | Fail-closed `stop()`; regression test. |
| C-5 | P3 | `R/multiomic_workflows.R` | `.cooperative_prediction_matrix` silently recycled/truncated wrong-length predictions (`matrix(out, nrow=n_rows)`) feeding AUC/mse tuning. | Fail-closed `stop()` on row mismatch; regression test. |
| D-3 | P3 | `vignettes/stablr-advanced.Rmd` | The Jaccard heatmap fed the default `remove_diag=TRUE` N×(N−1) matrix through `as.table()`, labelling a fake second axis "random_state" — the deployed figure's y-axis was meaningless. | Heatmap chunk now uses `remove_diag=FALSE` with explicit seed dimnames on both axes; site re-renders on push. |

### Fixed in this commit (documentation)

| ID | Sev | Where | Defect | Fix |
|---|---|---|---|---|
| N-01 | P3 | `man/stabl.Rd` | Documented "Automatic lambda grid when none is supplied" — no such option exists; `lambda_grid` is hard-coded `"auto"` and passing it via `...` errors ("matched by multiple actual arguments", probe-verified). | Reworded: the grid is always automatic; manual grids go through `stabl_fit()`. |
| B-03 | P3 | `man/auto_ordinal_lambda_grid.Rd` | alpha documented as `[0, 1]`; code rejects 0 (`min_open=TRUE`). | Corrected to `(0, 1]`. |
| B-05 | P4 | `man/stabl_preprocess_fit.Rd` | "transforms to fit, in order" was false — application order is fixed (zv → impute → center → scale). | Documented the fixed order. |
| C-1 | P3 | `man/stabl_multiomic_scheme_comparison.Rd` | Documented a nonexistent `stabl_workers` argument (actual formal: `workers`; `stabl_workers` is a returned element). | Corrected to `workers`. |
| C-2 | P3 | `man/load_ool_data.Rd` | Claimed "up to 27 validation samples"; the bundled valid split is 21×100 (probe-verified). | Corrected to 21 (after ID intersection). |
| C-6 | P3 | `man/stabl_multiomic_train_validate.Rd` | `oof_score` optimism (stacking weights tuned on the same OOF matrix, max over `n_iter_lf` draws) was disclosed only in a source comment (late_fusion.R:2-9). | Caveat + base-model-OOF definition added to the @return section, with a pointer to nested CV for honest estimates. |
| D-9 | P4 | `man/pearson_similarity.Rd` | "maximum is typically close to 1" understated: the statistic is the phi coefficient, exactly 1 for identical non-degenerate sets. | Corrected. |
| D-11 | P4 | `man/export_stabl_to_csv.Rd` | Doc stated the artificial-files condition as `artificial_type != NULL`; code tests `stabl_scores_artificial_`. | Reworded to the implemented condition. |

### Deferred to Phase 13 (remediation matrix additions)

Rationale for deferral: API/behaviour design decisions (user-facing
signature or default changes, convention changes that alter benchmark
outputs) must go through the remediation program, not a pre-publication
commit.

| ID | Sev | Summary |
|---|---|---|
| N-03 | P4 | Exported bootstrap samplers skip integer/weight validation (`n_subsamples=10.5` silently truncated; `0` → misleading 1000-retry error; negative `class_weights` → cryptic error). |
| N-04 | P4 | `stabl_fit` returns 5 undocumented fields; `stratify_bootstrap` field reports *effective* stratification (argument/field semantic mismatch). |
| N-05 | P4 | `auto_lambda_grid(l1_ratio = numeric(0))` silently returns NULL, violating the documented return contract. |
| N-06 | P4 | `stabl()` y-alignment doc omits the unnamed-x + named-y case (fails loudly, but the documented rule does not cover it). |
| B-04 | P3 | `make_ranger_adapter`/`make_aorsf_adapter` drop `family`; binomial y-coercion unreachable via public factories; 0/1-numeric y silently fits a regression forest (distinct from F-026). |
| B-06 | P4 | `auto_splsda_keepx_grid` with p=1 emits keepX=2 > p (mixOmics errors without stablr context). |
| B-07 | P4 | sPLS-DA batch adapter silently emits all-FALSE masks for bootstraps with constant columns (source comment only). |
| B-08 | P4 | Adaptive-lasso ridge init takes the smallest-path lambda (≈OLS in p≫n, forfeiting shrinkage); penalty factors up to 1/ε applied to a grid calibrated without `penalty.factor`. |
| C-4 | P3 | `stabl_multiomic_cv` folds never outcome-stratified; single-class guardrail fires only for grouped folds. |
| C-7/C-7b | P4 | `stabl_interpret` silently overrides `loss` for gaussian+log_loss; NA in `x` propagates silently into importance NAs. |
| C-8 | P4 | `stabl_interpret`/`stabl_analyze` permutation-importance helpers clobber the global RNG without restoring (F-CV-02 class at new locations). |
| C-9 | P4 | Defaults drift: `n_bootstraps` 1000 (core) vs 100 (all wrappers); `family`/`lambda_grid` defaults inconsistent across wrappers. |
| C-10..C-14 | P4 | Doc gaps and silent degenerate outputs in multiomic helpers (cooperative diagnostics columns, single-study diagnostics, `-Inf` stacked score, numeric-y coercion, block-pooling semantics). |
| D-2 | P3 | `save_stabl_results` default `task_type="binary"` renders a degenerate per-unique-y-value boxplot for regression fits; no cross-check vs `object$family`. |
| D-4 | P3 | `.pairwise_env_similarity` hard-codes diagonal 1 although `adjusted_similarity` is not reflexive (self-similarity = 1−k/d); convention differs from `jaccard_similarity_matrix`. |
| D-5..D-8, D-10, D-12, D-13 | P4 | Edge cases and hygiene: `jaccard_similarity_matrix(list())` obscure error; unbounded `fdr_target` rescales plot axis; unvalidated `fscore` beta; NA/duplicate tokens in similarity metrics; lexicographic level ordering; internal `rowMaxs` shadowing the matrixStats import; late `figure_fmt` validation leaves partial export dirs. |

## Ensemble-coherence notes (cross-function)

1. **Threshold semantics split:** glmnet/sgl/ordinal/splsda use strict `>`;
   rf/oblique use `>=`. Under exactly-zero importances (null data) `>=`
   selects everything — a contributing mechanism to F-P7-1. Unifying is a
   behaviour change → Phase 13.
2. **Grid conventions differ three ways:** rf/oblique prefer a `threshold`
   column, splsda prefers `keepX`, everything else is `lambda`-only.
   Documented per function; a unified convention is a breaking change.
3. **Factory signature inconsistency:** the sPLS-DA factory exposes
   `family`/`threshold`; ranger/aorsf expose neither (B-04); `num_trees` is
   hard-coded 500L in `stabl_fit`'s auto paths (stabl_fit.R:896-903) and
   ordinal alpha hard-coded 1 — not user-reachable through `stabl_fit`.
4. **Defaults drift:** `n_bootstraps` 1000 (core) vs 100 (all multiomic
   wrappers); `family` gaussian vs multinomial; `lambda_grid` required vs
   `"auto"` (C-9). Each man page shows its signature default, so no
   doc/code mismatch — but the 10× resolution drop between core and
   wrappers deserves a vignette note.
5. **Grouping-argument naming:** `groups_train` / `groups` /
   `study_labels` / `sample_groups` name two distinct roles across the
   family; forwarding is correct but easy to misroute.
6. **Empty-set conventions disagree across similarity families**
   (`adjusted(∅,∅)=0` vs `pearson(∅,∅)=1`); each is documented locally, the
   cross-family picture is not.
7. **Diagonal conventions disagree** (D-4): Jaccard matrix computes its
   diagonal; env-stability hard-codes 1 for a non-reflexive metric.
8. **Return-shape asymmetry:** `selected_features` is a named list in
   train_validate/cv but a long data.frame in nested CV;
   `*_summary()$err` is length-2 (IQR) for median and length-1 (RMSD) for
   mean.

## Defaults-calibration verdict

- `stabl()`'s opinionated defaults (knockoff_equi, FDP+ on, strict knockoff
  fallback, `n_bootstraps = 1000`) are consistent with the package's own
  benchmark evidence (`benchmarks/learner_comparison/guide_table.csv`;
  the vignette defaults table matches the actual formals).
- `n_lambda = 30`, `sample_fraction = 0.5` match the Python reference and
  benchmark configurations.
- RF grid defaults (500 trees, quantile probabilities `seq(0.5, 0.999)`)
  match the man pages and WS-2 calibration scope; the F-P7-1 caveat (fixed
  full-data quantiles do not calibrate refit noise; max-over-lambda
  scoring saturates under the null) stands and is a documented limitation,
  not a defaults bug.
- The `n_bootstraps = 100` wrapper default is the one calibration concern
  (C-9): documented, but a silent 10× resolution drop relative to the
  core — Phase 13 candidate for alignment or a doc note.
- Ordinal stays experimental (WS-9), consistent with its man-page framing.

## Verification after revisions

- Full test suite (`testthat::test_local`, stop_on_failure = FALSE): 54
  test files, **0 failures**. The only failures in the first run were two
  pre-existing assertions in `test-public-interface-audit.R` (a hard-coded
  72-export list and registry CSV that the sPLS-DA commit a5d1470 outgrew
  without updating the audit test) — corrected in this commit (registry
  rows added for `make_splsda_adapter`/`auto_splsda_keepx_grid`, count
  72 → 74) and re-verified: 168/168 expectations pass. 1 pre-existing skip
  (knockoff_equi recovery asserted on CI only); all warnings pre-existing
  (deprecation shims, MVR fallback, phase7 warning-path exercises).
- `R CMD check` (tarball built with vignettes): **0 errors, 0 warnings,
  2 NOTEs** — (1) installed size 6.3 Mb (doc 3.6 Mb from the nine built
  vignettes), (2) `vignettes/figure/` "looks like a leftover from
  'knitr'" (committed rendered PNGs; pre-existing repo hygiene — added to
  the Phase 13 list). Tests inside the check: FAIL 0 | WARN 19 | SKIP 20 |
  PASS 3067 (warnings are the pre-existing deprecation-shim / MVR-fallback
  set).
- Vignette render verification: R CMD check's "re-building of vignette
  outputs ... OK" re-rendered all 9 vignettes including the modified
  `stablr-advanced.Rmd` (D-3 fix renders cleanly).
- Regression tests added: `tests/testthat/test-final-review-regressions.R`
  (B-01, B-02, N-02, D-1, C-3, C-5).

## Publication-readiness verdict

After this revision pass there are **no known P1 defects** and every P3
finding is either fixed or explicitly deferred with a rationale in the
remediation matrix. The remaining deferred items are robustness/hygiene
(P4) and four behaviour-design P3s (B-04, C-4, D-2, D-4) that require
user-facing API or convention decisions. From the AI-audit standpoint the
package is publication-ready **pending**: (1) SCI-04 — independent human
confirmation, which no AI audit can substitute; (2) Phases 12-15
(synthesis, remediation execution, referee report).
