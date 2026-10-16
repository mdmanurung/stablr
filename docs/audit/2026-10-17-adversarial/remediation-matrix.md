# Remediation priority matrix (Phase 13 draft — pre-Phase 7)

Consolidated from findings F-001..F-026, F-CV-01..08, P4-F24..F29, F-P6-1..3,
F-P10-1..2. Priorities: P0 = blocks any release; P1 = scientific-integrity
fixes before release; P2 = fix or explicitly fence before release; P3 =
hygiene/robustness, batch at leisure. Every remediation must preserve the
Python-parity contract or document the deliberate divergence (no silent
workarounds).

## P0 — blocks any release

| ID | Defect | Fix | Est. |
|---|---|---|---|
| F-011 | Sparse path OOB read (SIGABRT) when `attr(x,"xm"/"xs")` missing; silent heap OOB under -DNDEBUG | Length-validate `xm`/`xs` vs `ni` and `xs > 0` in `elnet.fit` before `.Call` (or C++ entry → `jerr=10001`); add regression test | S |

## P1 — scientific integrity, fix before release

| ID | Defect | Fix | Est. |
|---|---|---|---|
| F-012 | NA in `x` silently drops columns (`jerr=0`, no warning) | `anyNA` check on x/y/weights in `elnet.fit`; regression test | S |
| F-013 | maxit non-convergence returns bare `list(jerr=-1)`, no warning; `jerr.multiview` conflates codes | Warn at `elnet.fit` boundary; distinct messages per code range | S |
| F-P6-1 | `min_fdr_` anti-conservative 1.5-2.2× on pure nulls (realized FDP 0.92-0.93 vs reported 0.42-0.66; 18-33 features selected on no-signal data) | Documentation-first: quantify null behavior in `?stabl_fit`/`?get_fdrs`/vignette; report the F-001 bound alongside `min_fdr_`; consider a null-calibration diagnostic. Do NOT silently change the threshold rule (parity contract) | M |

## P2 — fix or explicitly fence before release

| ID | Defect | Fix | Est. |
|---|---|---|---|
| F-005/L1 | RP-under-dependence warning fires only for rf/oblique_rf; default lasso+RP path un-warned | Extend warning to all learners (or document in `artificial_type` help); keep RF-specific wording | S |
| F-002 | Environment-aware min/mean aggregation is a novel unvalidated extension | Either a validation study (like WS-2/WS-5) or explicit "experimental, unvalidated" labeling in help + pkgdown | M |
| F-CV-01 | `random_state=NULL` + `cv_workers>1`: forked workers share `.Random.seed` → correlated folds, worker-count-dependent results | Derive per-worker seeds even when NULL (or warn loudly) | S |
| F-CV-02 | `stabl_fit` permanently clobbers global RNG | Save/restore `.Random.seed` on exit (pattern already exists in `.with_local_seed`) | S |
| F-CV-03 / P4-F28 | One `random_state` reused across folds/omics/cooperative foldid | Derive independent streams per use (stabl_fit's own 3-seed pattern) | S |
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
