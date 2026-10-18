# Remediation priority matrix (Phase 13 draft — pre-Phase 7)

Consolidated from findings F-001..F-026, F-CV-01..08, P4-F24..F29, F-P6-1..3,
F-P10-1..2. Priorities: P0 = blocks any release; P1 = scientific-integrity
fixes before release; P2 = fix or explicitly fence before release; P3 =
hygiene/robustness, batch at leisure. Every remediation must preserve the
Python-parity contract or document the deliberate divergence (no silent
workarounds).

## Phase 13 progress

- **PR-1 (2026-10-18) — DONE, blocking set cleared:** F-011, F-012, F-013
  fixed in `R/cooperative-multiview.path.R` + `R/cooperative-get_start.R`
  (shared `.jerr_code_msg` mapper), with compiled-boundary regression tests
  in `tests/testthat/test-compiled-boundary-regressions.R` (9 blocks / 25
  expectations). Full suite after the fix: 56 files, 593 tests, 0 failed,
  0 errors, 3416 passed, 1 pre-existing skip. Numerical behavior on valid
  inputs is unchanged (validation and warnings only — parity contract
  intact). Next: PR-2 (RNG-hygiene cluster).

- **PR-2 (2026-10-18) — DONE, RNG-hygiene cluster cleared:** F-CV-01,
  F-CV-02, F-CV-03, P4-F28 fixed across `R/stabl_fit.R` (hermetic
  `.Random.seed` save/restore), `R/nested_cv.R` (pre-fork per-fold seeds;
  seeded downstream `cv.glmnet` foldid), `R/multiomic_workflows.R` +
  `R/late_fusion.R` (per-consumer stream derivation via new
  `.derive_stream_seeds`; offsets 4000/41000/42000). Two additional
  reachable unseeded `cv.glmnet` foldid sites found and fixed during
  implementation (corrects F-CV-07's "unreachable" judgment — see
  findings.md). New regression file `tests/testthat/test-rng-hygiene.R`
  (9 blocks). Full suite after the fix: 57 files, 602 tests, 0 failed,
  0 errors, 1 pre-existing skip, 52 warnings (identical to pre-PR-2).
  R-internal multiomic characterization golden values (char-F1, char-F2-cv)
  deliberately updated with dated rationale — old values pinned the
  seed-reuse behavior the audit flagged; per-omic selections unchanged;
  Python-parity fixtures (68 checks) and multiview-parity (63) pass clean,
  so the parity contract is intact. CRAN-comments draft updated for the
  new `skip_on_cran()` guard (4 files). Next: PR-3 (documentation/fence
  cluster).

- **RNG sweep (2026-10-18, PR-2 follow-up):** repo-wide enumeration of all
  RNG primitives + implicit consumers across R/ and src/, classified
  against the PR-2 graph, with empirical probes (6/6 as expected). Core
  fit/CV paths confirmed complete. Three NEW P3 findings outside the core
  paths — F-CV-09 (`make_artificial_features`), F-CV-10 (`stabl_interpret`
  permutation importance), F-CV-11 (`stabl_analyze` permutation
  importance) — all bare-`set.seed` clobber sites, all empirically
  confirmed; queued for the PR-4 P3 batch. One vendored residual
  documented (cv.multiview internal foldid, unreachable in-package,
  parity-locked).

## P0 — blocks any release

| ID | Defect | Fix | Est. |
|---|---|---|---|
| ~~F-011~~ | Sparse path OOB read (SIGABRT) when `attr(x,"xm"/"xs")` missing; silent heap OOB under -DNDEBUG | **FIXED (PR-1, 2026-10-18):** length + finite/positive validation in `elnet.fit` before `.Call`; regression tests | DONE |

## P1 — scientific integrity, fix before release

| ID | Defect | Fix | Est. |
|---|---|---|---|
| ~~F-012~~ | NA in `x` silently drops columns (`jerr=0`, no warning) | **FIXED (PR-1, 2026-10-18):** `anyNA` checks on x (dense+sparse)/y/weights in `elnet.fit`; regression tests | DONE |
| ~~F-013~~ | maxit non-convergence returns bare `list(jerr=-1)`, no warning; `jerr.multiview` conflates codes | **FIXED (PR-1, 2026-10-18):** boundary warning in `elnet.fit`; shared `.jerr_code_msg` with distinct messages per code range in both mappers | DONE |
| F-P6-1 | `min_fdr_` anti-conservative 1.5-2.2× on pure nulls (realized FDP 0.92-0.93 vs reported 0.42-0.66; 18-33 features selected on no-signal data) | Documentation-first: quantify null behavior in `?stabl_fit`/`?get_fdrs`/vignette; report the F-001 bound alongside `min_fdr_`; consider a null-calibration diagnostic. Do NOT silently change the threshold rule (parity contract) | M |

## P2 — fix or explicitly fence before release

| ID | Defect | Fix | Est. |
|---|---|---|---|
| F-005/L1 | RP-under-dependence warning fires only for rf/oblique_rf; default lasso+RP path un-warned | Extend warning to all learners (or document in `artificial_type` help); keep RF-specific wording | S |
| F-002 | Environment-aware min/mean aggregation is a novel unvalidated extension | Either a validation study (like WS-2/WS-5) or explicit "experimental, unvalidated" labeling in help + pkgdown | M |
| F-CV-01 | `random_state=NULL` + `cv_workers>1`: forked workers share `.Random.seed` → correlated folds, worker-count-dependent results | Derive per-worker seeds even when NULL (or warn loudly) | **FIXED (PR-2, 2026-10-18)** — per-fold seeds drawn pre-fork via `.derive_stream_seeds`; regression test asserts worker-count independence |
| F-CV-02 | `stabl_fit` permanently clobbers global RNG | Save/restore `.Random.seed` on exit (pattern already exists in `.with_local_seed`) | **FIXED (PR-2, 2026-10-18)** — `on.exit()` save/restore both branches; hermeticity regression tests |
| F-CV-03 / P4-F28 | One `random_state` reused across folds/omics/cooperative foldid | Derive independent streams per use (stabl_fit's own 3-seed pattern) | **FIXED (PR-2, 2026-10-18)** — `.derive_stream_seeds` per-consumer streams at all fan-out sites (offsets 4000/41000/42000); R-internal multiomic golden values deliberately updated with dated rationale; Python-parity contract intact |
| F-CV-09 | `make_artificial_features()` (exported): bare `set.seed(random_state)`, never restored — clobbers global RNG (F-CV-02 class) | Same `on.exit()` save/restore as PR-2 | OPEN (P3; found in 2026-10-18 repo-wide RNG sweep; empirically confirmed) |
| F-CV-10 | `stabl_interpret()` → `.interpret_permutation_importance()`: bare `set.seed`, never restored | Same `on.exit()` save/restore as PR-2 | OPEN (P3; found in 2026-10-18 repo-wide RNG sweep; empirically confirmed) |
| F-CV-11 | `stabl_analyze()` → `.permutation_importance_sets()`: bare `set.seed`, never restored | Same `on.exit()` save/restore as PR-2 | OPEN (P3; found in 2026-10-18 repo-wide RNG sweep; empirically confirmed) |
| F-CV-04 / P4-F24 | `cooperation_selection="validation"` tunes on the assessment fold; optimistic metrics, no guard | Default-safe already; add explicit warning + rename option to make bias legible | S |
| F-P6-2 | Calibration regime-dependent in both directions (anti-conservative independent-null; conservative correlated) | Document regime dependence; consider a regime diagnostic (e.g. reported selection-rate ratio real:artificial) | M |
| F-P10-1 | Locked release contract is 0.1.1-era vs DESCRIPTION 0.3.0 | Re-ratify contract at 0.3.0 (owner action) or amend contract version | M |

## P3 — hygiene / robustness (batch)

- F-014: copy-guard all warm-start buffers passed to C++ (aliasing hazard).
- F-015: fix or rename `elnet.fit$dev.ratio` (not a deviance ratio; units mismatch cold/warm).
- F-016 (corrected): `fit$warm$a` → `fit$warm_fit$a` (works today only via `$` partial matching; one added `warm*` element from breakage).
- F-017: zero-guard `xs_(k)` in `spwls_exp` + constant-column detection.
- F-018: drop dangling `-Iglmnetpp/src`; pin vendored glmnetpp version/commit; name coxdev in COPYING.cooperative; prune ~90 dead vendored files; keep exception `what()` in jerr=10001 message.
- F-019: derive per-environment seeds via hashed streams, not `random_state + i`.
- F-020: replace message-regex error classification with condition classes where possible; keep regex as fallback.
- F-021: message when `workers > 1` silently degrades to sequential.
- F-022: document group-sampler undershoot possibility (all groups exhausted).
- F-023: document SGL noise-column group inheritance.
- F-024: tolerance-guard the ordinal path-trim fallback (reuse `.near_lambda_match`).
- F-025: validate ordinal outcome once per `stabl_fit`, not per bootstrap.
- F-026: validate declared family vs `y` type for RF/oblique (factor y + gaussian → error or warn).
- F-CV-05: fix nested-seed offset collisions.
- F-CV-06: upstream cvstats `nlams` discard (cvsd at path end) — note, vendor-tracked.
- F-CV-07: seed the latent internal foldid in `.mv_cv_multiview` (unreachable today).
- F-CV-08: warn when nested CV lacks a repeated-measures fold unit.
- P4-F25: harmonize binomial mse scale across cooperative branches.
- P4-F26: derive a fused-model lambda grid (or document first-omic reuse).
- P4-F27: per-view `hard_threshold` support.
- F-P10-2: comment in DESCRIPTION for Suggests entries used only by inst//vignettes.

## Release gates (unchanged, human tasks)

1. `fresh-release-run-pending` — fresh n=150 release run after remediation.
2. `sci04-independent-confirmation-pending` — independent human statistician
   confirmation; the AI audit (this program) does NOT satisfy SCI-04.
3. Phase 7 (aging case-study adversarial audit) — blocked on the v2
   bit-for-bit reproduction run (in flight).

## Sequencing proposal

1. P0 + P1 code fixes (F-011..F-013) + regression tests — one PR, small.
2. P2 RNG-hygiene cluster (F-CV-01..03, P4-F28) — one PR, mechanical.
3. P2 documentation cluster (F-P6-1, F-P6-2, F-005 warning scope, F-002
   labeling, F-CV-04 warning) — one PR.
4. P3 batch — one or two PRs.
5. Re-ratify release contract (F-P10-1) → fresh release run → SCI-04 human
   confirmation → referee report (Phase 15) finalizes.

## Final-review additions (2026-10-18 — Phase 13 inputs)

From the final pre-publication per-function review (`final-review.md`);
fixed-there items (B-01, B-02, N-02, D-1, C-3, C-5, D-3 and the doc set
N-01, B-03, B-05, C-1, C-2, C-6, D-9, D-11) are not repeated here.

| ID | Sev | Defect | Proposed fix | Est. |
|---|---|---|---|---|
| B-04 | P3 | ranger/aorsf factories drop `family`; 0/1-numeric y silently fits regression forest | Add `family` param to both factories; coerce y like the internal batch adapters | S |
| C-4 | P3 | `stabl_multiomic_cv` folds never outcome-stratified; single-class guard only fires for grouped folds | Pass outcome-derived strata (or add `stratified`); extend guard to ungrouped folds | M |
| D-2 | P3 | `save_stabl_results` default `task_type="binary"` renders degenerate boxplots for regression fits | Cross-check `task_type` vs `object$family` (warn or derive default) | S |
| D-4 | P3 | `.pairwise_env_similarity` hard-codes diagonal 1 for a non-reflexive metric | Compute the diagonal with `adjusted_similarity(S_i, S_i, d)` or document the convention (changes outputs → needs decision) | S |
| N-03 | P4 | Exported bootstrap samplers skip integer/weight validation | `.validate_scalar_integer_like(n_subsamples, ..., min = 1L)`; reject non-positive `class_weights` | S |
| N-04 | P4 | 5 undocumented `stabl_fit` return fields; `stratify_bootstrap` field semantics | Document fields; consider `effective_stratification` rename | S |
| N-05 | P4 | `auto_lambda_grid(l1_ratio = numeric(0))` silently returns NULL | Validate `length(l1_ratio) >= 1` | S |
| N-06 | P4 | `stabl()` y-alignment doc omits unnamed-x + named-y case | One sentence in `@param y` | S |
| B-06 | P4 | `auto_splsda_keepx_grid` p=1 emits keepX=2 > p | Cap `hi` at p; pre-validate in adapter | S |
| B-07 | P4 | sPLS-DA batch adapter silently emits all-FALSE masks on constant columns | Document; tolerance-based constancy detection | S |
| B-08 | P4 | Adaptive-lasso ridge init takes smallest-path lambda (≈OLS in p≫n); penalty factors vs unpenalized grid calibration | Document init choice; consider mid-path init and penalized-grid recalibration | M |
| C-7/C-7b | P4 | `stabl_interpret` silent loss override; NA propagation into importances | Validate `loss` vs `family`; finite-value check or `na_reason` path | S |
| C-8 | P4 | `stabl_interpret`/`stabl_analyze` importance helpers clobber global RNG | Wrap in `.with_local_seed()` | S |
| C-9 | P4 | Defaults drift: `n_bootstraps` 1000 core vs 100 wrappers; `family`/`lambda_grid` defaults | Align or add "defaults differ" notes to wrapper Rds | S |
| C-10..C-14 | P4 | Multiomic doc gaps + silent degenerate outputs (cooperative diagnostic cols, single-study diagnostics, `-Inf` stacked score, numeric-y coercion, block pooling) | Batch doc + guard fixes | M |
| D-5 | P4 | `jaccard_similarity_matrix(list())` obscure error; n=1 + remove_diag → 1×0 | Early returns / informative error | S |
| D-6 | P4 | `plot_fdr_graph` unbounded `fdr_target` rescales axis | `max = 1` in validator | S |
| D-7 | P4 | `fscore_similarity` unvalidated `beta` | Positive-scalar check | S |
| D-8 | P4 | NA/duplicate tokens in similarity metrics | Validate or document set semantics | S |
| D-10 | P4 | Lexicographic level ordering; factor silently coerced to codes in scatterplot | Numeric-aware sort; reject factors with pointer | S |
| D-12 | P4 | Internal `rowMaxs` shadows dead `importFrom(matrixStats, rowMaxs)` | Rename to `.row_maxs` or drop the import | S |
| D-13 | P4 | Late `figure_fmt` validation leaves partial export dirs | Validate up front | S |
| FR-hygiene | P4 | `vignettes/figure/` committed knitr leftovers trigger an R CMD check NOTE | Verify nothing references them, then remove from the repo | S |

Ensemble-coherence decisions also owed to Phase 13 (see
`final-review.md` §Ensemble-coherence notes): threshold-semantics split
(strict `>` vs `>=`), three-way grid-convention divergence, factory
signature inconsistency, `n_bootstraps` 1000-vs-100 alignment, grouping
argument naming, empty-set and diagonal convention documentation.
