# Adversarial Audit Findings — living report

Format per Phase 14: Finding ID, severity (P0/P1/P2/P3/Deferred), location,
reproduction, evidence, scientific consequence, root cause, proposed change,
regression test, independent validation, residual limitations.

Status vocabulary: OPEN / CONFIRMED / REFUTED / ALREADY-RESOLVED / INCONCLUSIVE.

<!-- Findings are appended below as audit phases complete. IDs are assigned
     in discovery order and never reused. Reconciliation (Phase 13 Pass F)
     may reclassify but not delete entries. -->

## Findings

(none confirmed yet — Phase 2 first pass below records audit leads)

### Phase 1 findings (independent methodologist pass — integrated and spot-verified)

Full report: `specification-matrix.md` (20-row traceability matrix, primary
references with DOIs, formal specifications). Lead-auditor spot-checks
confirmed the key claims against source (strict `>` at
R/stabl_accessors.R:58; RP default at R/stabl_fit.R:203; RF guardrails at
R/stabl_fit.R:362-377; no quantile aggregation in package R/).

- **F-001 (P1, interpretive) — `min_fdr_` alpha semantics (G1).** Stabl has
  no user-facing FDR target: theta = argmin_t FDP+(t), q+ = min_t FDP+(t)
  is random, only stochastically relatively bounded
  (P(FDP(theta) >= (1+D)q+) <= 1/(1+D), Hédou et al. 2024 Eq. 7).
  Implementation is faithful and docstrings disclaim universal control,
  but `min_fdr_` invites a "controlled FDR" misreading downstream.
  Proposed: documentation + accessor naming guidance; simulation evidence
  of the (q+, realized FDP) joint distribution (Phase 6).
- **F-002 (P1, unvalidated extension) — environment-aware min/mean
  aggregation (G2).** No literature estimand; per-environment FDP+
  thresholds not propagated to the aggregated score; no error control at
  the aggregated level. The package's own docs acknowledge leniency.
  Proposed: keep, but label explicitly as uncalibrated research feature in
  all user-facing docs; Phase 6 simulation (environment-label permutation,
  environment-shift DGPs) to quantify calibration.
- **F-003 (P2, spec deviation, trivially fixable) — strict `>` support rule
  vs paper's `>=` (G6).** Discrete stability scores (multiples of 1/B) make
  exact ties at grid thresholds systematic; features at f_j = theta are
  excluded — anti-conservative in |S(theta)| and FDP+ counts. Both stablr
  AND Python STABL use strict `>` (parity-preserving), so fixing requires a
  documented divergence from the reference implementation or an upstream
  fix. Decision needed: parity vs paper fidelity.
- **F-004 (P2, documented approximation) — second-order Gaussian knockoffs
  (G4).** Exact model-X only under (approximate) multivariate Gaussianity;
  single knockoff realization reused across all B bootstraps (Python-parity)
  leaves unaveraged Monte Carlo noise in FDP+. Documentation + Phase 6
  calibration study under non-Gaussian margins.
- **F-005 (P2, silent-risk default) — RP nulls under dependence for glmnet
  learners (G5).** `stabl_fit` defaults to `random_permutation` (line 203);
  permuted features are exchangeable with nulls only under (near-)
  independence. RF learners warn (guardrails at 362-377); glmnet learners
  do not and cannot know. Proposed: extend the guardrail warning to all
  learners when `artificial_type = "random_permutation"` and correlation
  screening suggests dependence, or document the scope limit at the
  user level.
- **F-006 (P2, empirical-only calibration) — RF/oblique-RF/ordinal learners
  (G7).** Outside R1; package warns and restricts to model-X pairings.
  Keep restricted; Phase 6 calibration curves per learner x artificial-type.
- **F-007 (P3, resampling-distribution deviation) — classification
  feasibility repair (G8).** Redraw/stratified fallback when a subsample
  has <2 per class conditions the resampling scheme (deviation from i.i.d.
  subsampling of R1/R2). Benign for balanced problems; validate on
  imbalanced DGPs. Failed fits themselves fail closed (no denominator
  distortion) — the safe choice, confirmed.
- **F-008 (P3, numeric parity) — rounding and proportion semantics (G10).**
  `round` vs `int` for n_injected; effective vs nominal artificial
  proportion in 1/pi scaling. Identical at pi=1; sub-1% count differences
  otherwise. Document; cross-language test at pi in {0.3, 0.5, 0.7}.
- **F-009 (Deferred, validation task) — cooperative engine objective (G9).**
  Faithful for gaussian/binomial per the methodologist pass; numerically
  verify against the direct augmented-matrix solve (R3 Eqs. 23-24) in
  Phase 6; poisson/cox dead code documented as unsupported.
- **F-010 (P3, documentation) — lambda-grid provenance (G10).** stablr's
  auto lambda grids are glmnet-derived, not Python's linspace on the
  inverse-penalty C — comparable only via the stability-path/FDP+ layer.
  Document in PYTHON_TO_R_MAPPING.md if not already there.

**Summary judgment (Phase 1, verified):** core Stabl machinery is a faithful
Python-parity translation; the two scientifically consequential unverified
extensions are environment-aware aggregation (F-002) and RF-family/ordinal
learners (F-006); the most dangerous interpretive gap is `min_fdr_`
semantics (F-001).

### Phase 2 first-pass leads (status: OPEN leads, not confirmed defects)

- **L1 (RP exchangeability under correlated predictors)** — `make_rp_features`
  (R/artificial_features.R:16-45) shuffles copies of randomly chosen real
  columns: artificial features are independent of *everything*, while null
  real features correlated with signal features are selected by proxy. The
  FDP+ estimator's exchangeability premise (artificial features exchangeable
  with *null* real features) therefore weakens under correlated designs.
  Verify what Hédou et al. claim for the permutation variant and what the
  package documents; the dispatcher doc already disclaims universal control.
  Chase in Phase 1/6.
  **Update (baseline evidence):** the test suite itself states this
  limitation for the RF learner specifically (test-rf-adapter.R:178: real
  null features correlated with signals receive leaked importance while
  permuted nulls do not — FDP+ anti-conservative; package recommends
  model-X modes). Open question narrowed to: (a) is the warning raised
  everywhere it applies (all learners, all entry points), (b) does the
  paper's permutation variant carry the same caveat, (c) is the
  documentation at user level (not just test level).
- **L2 (PD repair changes the fitted model)** — `.estimate_pd_sigma`
  (R/artificial_features.R:396-410) applies corpcor shrinkage or a diagonal
  ridge when cov(x) is not PD; model-X knockoffs are then drawn from an
  approximate Gaussian model. Approximation is standard practice but the
  equi/MVR docs do not state it as explicitly as the fixed-X caveat does.
  Documentation-level unless simulations show calibration impact (Phase 6).
- **L3 (novel MVR solver needs independent validation)** — `solve_mvr()` is
  an R-native coordinate-descent port of `knockpy.mrc._solve_mvr_ungrouped`
  (novel to stablr, no Python-STABL counterpart). Requires independent
  validation against knockpy or against S-matrix feasibility properties
  (0 < s_j <= 2*lambda_min). Phase 6 task.
  **Update (2026-10-17, independent numerical verification — PARTIALLY
  RESOLVED):** property tests executed at this commit
  (`/workspace/audit_mvr_checks.R`, all passed):
  (i) `.mvr_cholupdate` matches direct `chol(V ± xx')` to 3.6e-15 over 200
  random trials for both update and downdate — the nonstandard
  `c = r/Rkk, s = xk/Rkk` Givens form is the same rotation as the standard
  `c = Rkk/r, s = xk/r` (algebraically verified: `(R_tail + s*x_tail)/c`
  with the code's c,s equals `c_std*R_tail + s_std*x_tail`);
  (ii) 20/20 random PD Sigmas: reproducible under fixed seed, S > 0,
  2*Sigma - S PSD at tol, MVR loss strictly decreased from the
  min-eig*I initialization (e.g. 122.01 -> 90.61);
  (iii) 2x2 case: solver S (0.7394, 0.4692) matches brute-force grid
  argmin (0.74, 0.47), identical loss 6.2196 — coordinate descent reaches
  the global optimum there.
  First-order condition independently re-derived: the code's quadratic is
  the negation of the derived stationary equation, same roots; both
  branches retained and minimized by explicit evaluation (correct).
  Residual: comparison against the actual knockpy Python implementation
  (not installed in this environment — blocked, optional) and downstream
  FDP+ calibration impact (Phase 6 simulations).
- **L4 (max-over-lambda + threshold minimization optimism)** —
  `compute_fdp_plus` takes rowMaxs over the lambda grid and then minimizes
  FDP+ over 100 thresholds. Whether the paper's procedure includes the same
  double optimization (and any correction) must be established in Phase 1;
  the realized FDP at the chosen threshold is a post-selection quantity.
- **L5 (REFUTED lead: FDP+ test coverage)** — ledger token-matching showed
  `fdp_control.R` with no name-matched tests; direct grep found
  `test-fdp-plus-invariants.R`, `test-python-parity-fixtures.R` (Python
  reference fixtures `fdp_plus_curve.csv`, `fdp_plus_by_lambda.csv`), and
  coverage in `test-stabl-fit.R`. Ledger corrected manually; matcher
  limitation noted (tokens < 5 chars filtered).

### Phase 2 first-pass positive observations (evidence of correctness)

- `compute_fdp_plus` strict-exceedance counting via sort + `findInterval`
  complement is exactly `#(score > tau)` including tie behavior; tie-breaking
  in `which.min` matches Python `argmin`; min_fdr > 1 caps the threshold at
  1.0 (empty selection — fail-closed).
- Input contracts on both core functions are thorough (finite, [0,1],
  dimension agreement, positive proportion).
- Artificial-feature provenance (requested vs actual type, chunk records,
  ordered fallback history) and the strict `knockoff_fallback = "stop"`
  default are scientifically motivated: permutation features do not support
  FDP+ calibration, so a silent fallback would change the estimand.
- Fixed-X knockoff bootstrap-reuse limitation is explicitly documented
  (exchangeability does not survive row subsampling) with a preference
  ordering toward model-X modes.
- MVR fallback chain (MVR -> equi -> RP) records each hop in
  `fallback_history`; the mvr->equi hop stays a warning under both policies
  (both model-X).

## Audit leads carried forward (to verify against current commit, not assume)

| Lead | Prior session evidence | Re-verification plan |
|---|---|---|
| Empirical FDP calibration failures (learner-specific) | S1/S2 benchmarks: M4/M5 saturate 59/59 under all truth regimes; empirical FDP 0.83-1.00 | Phase 2D/6: independent FDP+ derivation + synthetic DGP calibration check |
| Null-selection saturation in aging benchmark | S0 null: 50/50 features, CP 0.942 | Phase 7: re-derive from committed artifacts; attribute algorithmic cause |
| Historical cross-language parity failures | python_parity fixtures; IV CI smoke | Phase 2/9: fixture provenance + current parity status |
| Multiview parity warnings | CRAN multiview deprecations in parity tests | Phase 4: vendored engine vs reference outputs |
| Ordinal graduation criteria | WS-9 incomplete runs; experimental status | Phase 5: check predeclared criteria vs current evidence |
| Release contract version skew | locked contract is 0.1.1-era vs DESCRIPTION 0.3.0 | Phase 11/15: release-readiness assessment |
| SCI-04 independent confirmation pending | provisional owner acceptance only | Phase 12/15: must remain a human task; AI audit must not self-satisfy |

## Phase 8C — Compiled-code & numerical-integrity audit (delegated subagent; lead-auditor verification 2026-10-17)

Full report: `specification-matrix.md` companion — compiled-layer report delivered
separately; runtime confirmations executed against a package copy installed to
`/tmp/rlib` (repo untouched). Lead auditor independently reproduced the P0 and
both P1 runtime claims and statically confirmed P2-3 before integration.

- **F-011 (P0, CONFIRMED — independently reproduced):** Sparse-path
  out-of-bounds read when `attr(x,"xm"/"xs")` are missing. `elnet.fit` builds
  `xm`/`xs` via `as.double(attr(x,"xm"))` (`R/cooperative-multiview.path.R:1071-1072`);
  a plain `dgCMatrix` yields `numeric(0)`, and `spwls_exp`'s `Eigen::Map` then
  reads out of bounds in `compute_grad`/`compute_xv` (`sp_gaussian_wls.hpp`).
  Reproduced: R process aborts with Eigen `index < size()` assertion (SIGABRT,
  exit 134). Under `-DNDEBUG` the same call is a silent heap OOB read. The
  shipped `multiview.fit` path always sets the attributes, so current workflows
  are safe; the boundary is one refactor away from memory unsafety. Fix: length
  validation (`length(xm)==ni`, `length(xs)==ni`, `xs>0`) at the R wrapper or
  C++ entry.
- **F-012 (P1, CONFIRMED — independently reproduced):** NA in `x` is silently
  accepted. NaN gradients fail all strong-rule/KKT comparisons
  (`check_kkt`: `g > l1_regul*penalty` is FALSE for NaN), so the affected
  column is silently dropped with `jerr=0` and no warning. Verified: NA column
  returns all-zero coefficients, `jerr=0`, no diagnostic.
- **F-013 (P1, CONFIRMED — independently reproduced):** maxit non-convergence
  returns a bare, class-less `list(jerr=-1)` with no warning at the
  `elnet.fit` boundary (visible only as path-truncation warning in
  `multiview.path`). `jerr.multiview` conflates all negative codes as
  "convergence not reached" and all positive codes < 7777 as "memory
  allocation error" (misleading for `jerr=10001`, any `std::exception`).
- **F-014 (P2, CONFIRMED static):** C++ mutates caller-owned warm-start R
  memory in place via non-const `Eigen::Map`s; only `a` is copy-guarded on the
  R side (`:1063-1064`). Benign in the current call graph (warm state is
  recomputed before each call) but bypasses R copy-on-write; latent aliasing
  hazard.
- **F-015 (P2, CONFIRMED runtime):** `elnet.fit$dev.ratio` is not a deviance
  ratio — cold-start returns the raw unnormalized `rsqc` accumulator while the
  warm branch seeds a fraction (units mismatch). Both current consumers
  overwrite it before use, so no current impact; misleading if consumed raw.
- **F-016 (P2, CONFIRMED static; reachability SUSPECT for gaussian-only):**
  `fit$warm$a` at `:657` reads a non-existent element (result carries
  `warm_fit`, never `warm`) → "non-conformable" error on any IRLS iteration
  ≥ 2 (non-gaussian families). Gaussian fits converging in one pass never
  execute the line. Upstream multiview bug carried over. Fix: `warm_fit$a`.
- **F-017 (P2, SUSPECT):** `spwls_exp` divides by `xs_(k)` with no zero guard;
  safe today only because the caller hard-codes `xs=1` and the wrapper comment
  admits "assume that there are no constant variables" — constant columns are
  not detected, unlike glmnet.
- **F-018 (P3, CONFIRMED):** Build/vendoring hygiene: dangling
  `-Iglmnetpp/src` include; no `-DNDEBUG` (assertions enabled — turned P0-1
  into a clean abort, but an explicit bounds check would be strictly better);
  coxdev + elnet_driver/elnet_path trees are dead vendored code (~90 of 101
  files unused, enlarging the GPL-2 surface); vendored glmnetpp version/commit
  not pinned anywhere (only "multiview 1.0" in `inst/COPYING.cooperative`);
  coxdev sub-library not named in COPYING.cooperative; `driver.h` discards
  exception `what()` text; unused `no`/`ni` parameters.
- **Positive (clean bill for the core):** coordinate descent, soft-thresholding,
  strong rules, KKT re-check loop, and 1-index bookkeeping are faithful and
  internally consistent; no RNG, no threading, no `-ffast-math`, no POSIX-only
  headers; `Makevars.win` byte-identical to `Makevars`; RcppExports argument
  order matches all R callers; attribution via `inst/COPYING.cooperative` is
  adequate for the combined work (GPL-2 consistent with DESCRIPTION).

## Phase 2 completion — core machinery audit (lead auditor, 2026-10-17)

Files audited this pass: `stabl_fit.R`, `environment_stability.R`,
`learner_adapters.R`, `conditions.R`, `stabl_accessors.R`,
`bootstrap_helpers.R` (completing fdp_control.R, artificial_features.R,
knockoff_mvr.R from the first pass). All nine core-machinery files read in
full.

**L1 scope question (a) ANSWERED — warning-scope gap confirmed:** the
RP-under-dependence warning (`.warn_rf_random_permutation_nulls()`) fires
only when `base_learner` is `rf`/`oblique_rf` (`stabl_fit.R:362-377`). The
package-default configuration — `base_learner = "lasso"` with
`artificial_type = "random_permutation"` (also the default) — never warns,
although F-005/L1's Assumption-1 concern (permuted nulls break the null
distribution under feature dependence) applies to every learner. The
default path is the un-warned one. Remediation candidate: extend the
warning to all learners (or document the limitation in `stabl_fit()`'s
`artificial_type` help entry, which currently does not mention it).

**Positive observations (machinery):**
- RNG architecture is sound: three independent seed streams (artificial
  features, bootstrap indices, per-iteration adapter seeds) derived from one
  seeded draw; per-iteration seeds make sequential and `furrr` parallel
  paths bit-identical; `.with_local_seed` restores `.Random.seed` correctly
  including the no-prior-seed case.
- Effective (post-rounding) artificial proportion — not the requested value —
  is forwarded to `compute_fdp_plus`, so the FDP+ weighting matches the
  actually injected count.
- Threshold resolution order (new_hard_threshold > hard_threshold >
  fdr_min_threshold_) is fail-closed; `get_fdrs()` errors when no artificial
  features exist; `transform_stabl()` rejects missing/duplicate/empty column
  names and preserves 2-dimensionality.
- The explore fallback's tie-over-selection bug is fixed with a documented
  rationale (top-n by `order()`, not `sort() - 0.01`).
- Bootstrap samplers: whole-group invariant never trims partial groups;
  degenerate-bootstrap guard retries up to 1000 then fails closed;
  stratified counts use largest-remainder allocation with minimum
  enforcement and a final replace=FALSE feasibility check; strata interaction
  uses `sep = "\r"` to avoid delimiter collisions.
- Batch adapters reduce the loop to one glmnet path call per bootstrap (per
  unique alpha); the on-grid fast path is tolerance-guarded
  (`100*.Machine$double.eps`) with per-lambda fallback, and on-grid
  bit-identity is covered by `test-coef-batch-ongrid.R`.
- `conditions.R` typed-error design fails closed: unmatched errors propagate
  as hard errors; numerical-infeasibility classification is the only
  conversion.

**New P3 observations (recorded, non-blocking):**
- **F-019 (P3):** `stabl_environment_stability()` derives per-environment
  seeds as consecutive integers (`random_state + i`), unlike `stabl_fit()`'s
  hashed multi-stream derivation. Practically safe under R's MT seed
  scrambler but internally inconsistent and weaker than the package's own
  best practice.
- **F-020 (P3):** `.is_glmnet_numerical_error()` /
  `.is_linear_algebra_numerical_error()` classify by English message regex —
  brittle across glmnet/linalg version changes. Failure direction is safe
  (hard error instead of fold-skip), but a glmnet message change would
  silently convert skip-fold behavior into aborts.
- **F-021 (P3):** `workers > 1` silently falls back to sequential execution
  when furrr/future are not installed (no message).
- **F-022 (P3):** `group_bootstrap_indices()` documentation guarantees
  realised count >= n_subsamples, but when all groups are exhausted before
  the target the unstratified path can return fewer (undershoot possible
  with few groups and large sample_fraction); the stratified-group path
  errors instead. Documentation gap only.
- **F-023 (P3):** `.append_noise_groups()` assigns each injected noise
  column to its source feature's SGL group (new group only for
  non-derived columns). Group-level penalties can therefore be influenced
  by noise columns sharing a group with real features; SGL is experimental
  and this is a design choice worth documenting.

## Phase 5 — outcome-family audit (lead auditor, 2026-10-17)

Scope: family dispatch across glmnet adapters (`auto_lambda_grid`,
`.feature_abs_coefs{,_batch}`), Cox/Surv flow through `stabl_fit`, ordinal
(`learner_ordinal.R`), RF/oblique family coercion (`learner_rf.R`,
`learner_oblique_rf.R`).

- **Clean:** intercept handling is family-correct in both fast and slow
  coefficient paths (multinomial: per-class (p+1)-row matrices, intercept
  dropped per class, max-abs across classes; Cox: no-intercept p-row layout,
  `row_sel = TRUE` correct; gaussian/binomial: `-1L` drop). Fast-path
  (`fit$beta`) and slow-path (`coef.glmnet`) layouts agree.
- **Clean (runtime-verified):** Cox/Surv outcomes flow correctly through the
  bootstrap sampler — `length(Surv)` returns the number of subjects
  (survival ≥ 3.x defines `length.Surv`), so `sample.int(n = length(y))`
  indexes subjects, not columns. End-to-end `stabl_fit(family = "cox")` with
  a `Surv` outcome completes and selects coherently. Cox ties are pinned to
  Efron via the package constant `.COX_TIES` (parity with Python default).
- **Clean:** the ordinal learner is honestly labeled experimental with a
  predeclared WS-9 graduation criterion and preserved adverse evidence
  (knockoff_equi mean FDP+ 0.00 — nominal 0.1 — but mean TPR 0.37 on both
  simulated designs; strongly conservative; did not graduate). Fail-closed
  outcome validation (factor required, NA rejected, unordered-factor
  warning).
- **F-024 (P3):** the ordinal path-trim fallback matches trimmed lambdas by
  `which.min(abs(fit$lambdaVals - s))` with no tolerance guard (unlike
  `.near_lambda_match`), silently mapping grid rows to distant path points
  when ordinalNet trims the path.
- **F-025 (P3):** `.validate_ordinal_outcome()` re-warns on every bootstrap
  iteration for unordered factors (up to n_bootstraps identical warnings).
- **F-026 (P3):** RF/oblique adapters infer the task from `y` (ranger/aorsf
  inference), so a factor `y` with `family = "gaussian"` silently fits a
  classification forest while `stabl_fit`'s bootstrap sampling treats it as
  regression (`.modeled_class_levels` returns NULL). Declared-family/task
  mismatch is not validated.

## Phase 9 — test-suite audit (lead auditor, 2026-10-17)

Inventory: 59 test files, 562 `test_that` blocks, 13,935 lines, 3,325
passing assertions / 0 failures / 1 skip (fresh baseline).

- **All 7 ledger "no tests" flags are token-matcher false leads.** Every
  flagged file has direct coverage: `fdp_control.R` (Phase 2), `conditions.R`
  (test-nested-cv.R, test-late-fusion-oof.R, test-artificial-feature-provenance.R),
  `exports.R` (test-phase7.R, test-public-interface-audit.R), `metrics.R`
  (test-environment-stability.R, test-stabl-wrapper.R), `s3_methods.R`
  (test-accessors-s3.R), `RcppExports.R` (Phase 8C static verification),
  `stablr-package.R` (docs only, n/a). Ledger corrected.
- **Parity fixtures are strong:** fixture files are hash-verified
  (`expect_equal` on SHA hashes — bit-level provenance), and score
  comparisons use `tolerance = 1e-12` against Python-derived references.
- **Determinism claims are directly tested:** sequential reproducibility and
  parallel/sequential bit-identity (`test-rng-determinism.R`, furrr-gated)
  cover the M-2/V-7 seed architecture verified statically in Phase 2.
- **Guardrails are tested as guardrails:** the RF+RP test asserts that the
  anti-conservative-FDP+ warning fires (documents the limitation) rather
  than canonizing the anti-conservative result as expected behavior.
- **Gap (links to F-011..F-013):** no tests cover the compiled-boundary
  failure modes (missing sparse attributes, NA inputs, maxit
  non-convergence) — the P0/P1 behaviors were found by adversarial probing,
  not by the suite. Remediation: add regression tests alongside the boundary
  validation fixes.
- **Gap (minor):** the single skip (test-rf-adapter.R) means RF tests
  silently vanish in environments without ranger; acceptable for an
  optional-dependency learner, worth noting for coverage accounting.

## Phases 3, 4, 6 integration + F-016 correction (lead auditor, 2026-10-17)

### CORRECTION — F-016 downgraded from P2-crash to P3-fragility (lead-auditor verification)

The Phase 8C claim that `fit$warm$a` (`cooperative-multiview.path.R:657`)
"crashes any IRLS iteration ≥ 2" is **REFUTED by runtime experiment**. R's
`$` list operator performs partial matching: `fit$warm` resolves to
`fit$warm_fit` (the unique element whose name starts with "warm"), so the
coefficient read succeeds. Instrumented runtime trace (binomial cooperative
fit, rho = 0.5, 90 `multiview.fit` calls, every call warm-started with a
`glmnetfit`): completes without error. **P4-F29 (binomial cooperative fails
on realistic signal) is likewise refuted** — the same probe completed on a
strong-signal binomial design. The fragility is real and worth fixing
(`fit$warm$a` → `fit$warm_fit$a`): any future list element whose name starts
with `warm` and sorts before `warm_fit` silently breaks the read, and
`fit[["warm", exact = TRUE]]` semantics would crash. No current impact.

### Phase 3 — CV leakage audit (delegated independent pass; report: audit_cv_leakage.md)

**Clean bill on leakage integrity:** feature selection, artificial-feature
generation, FDP+ threshold resolution, auto lambda grids, fold preprocessing,
and held-out prediction are strictly fold-training-local in both nested CV
and `stabl_multiomic_cv`; folds are disjoint, stratified, seed-reproducible;
the vendored `build_predmat` predicts fold-i-rows with the fold-i model with
no view misalignment. No path exists for a fold to inherit a full-dataset
threshold (degenerate FDP+ caps at 1.0 — empty selection, conservative).

- **F-CV-01 (P2, CONFIRMED runtime):** with `random_state = NULL`, nested-CV
  results depend on `cv_workers` (1 vs 2 give different predictions); under
  `mclapply` every forked worker inherits the identical `.Random.seed`, so
  outer folds share correlated RNG streams.
- **F-CV-02 (P2, CONFIRMED runtime):** `stabl_fit` permanently clobbers the
  global RNG (bare `set.seed()` calls, never restored) — later unseeded
  randomness becomes a deterministic function of the last `random_state`.
- **F-CV-03 (P2, CONFIRMED static):** `stabl_multiomic_cv` reuses one
  `random_state` for fold construction, every fold's `stabl_fit`, and the
  cooperative internal foldid — correlated fold noise (not held-out leakage).
- **F-CV-04 (P2, CONFIRMED static):** `cooperation_selection = "validation"`
  inside `stabl_multiomic_cv` tunes rho/lambda on the outer assessment fold
  and reports the tuned predictions as held-out — optimistic metrics, no
  guard (default `"cv"` is safe). Same defect as P4-F24.
- **F-CV-05..08 (P3):** nested-seed offset collisions; vendored `cvstats`
  discards per-fold `nlams` (slightly wrong cvsd at path end, upstream-
  inherited); latent unseeded internal foldid (unreachable in-package);
  no repeated-measures fold unit in nested CV except `study_labels`
  (documented footgun, no warning).

### Phase 4 — multi-omic audit (delegated independent pass; report: audit_multiomic.md)

**Clean bill:** strict name+order alignment across views (runtime-confirmed
hard error on reversed order); fail-closed duplicate/NULL-threshold handling;
CV fold integrity across all views simultaneously; no full-data fits in
evaluation; leakage-safe OOF late fusion; deterministic reruns; correct
original-scale coefficient restoration in the vendored engine (gaussian and
binomial rho>0 paths).

- **P4-F24 (P2, CONFIRMED):** cooperative `cooperation_selection="validation"`
  tunes rho AND lambda on the validation set and reports that metric as
  `score` — selection-biased (documented but optimistic). = F-CV-04.
- **P4-F25 (P3, CONFIRMED):** binomial mse is link-scale in the cooperative
  validation branch vs response-scale in the CV branch.
- **P4-F26 (P3, CONFIRMED):** early fusion silently reuses the first omic's
  lambda grid for the fused model.
- **P4-F27 (P3, CONFIRMED):** `hard_threshold` is global-only; no per-view
  thresholds.
- **P4-F28 (P3, CONFIRMED):** seed reuse across folds/omics/cooperative
  foldid (hygiene; = F-CV-03).

### Phase 6 — independent synthetic-DGP validation (lead auditor; artifacts: /workspace/phase6_out/)

Designs: D1 sparse independent (10/150 signal), D2 pure null, D3
block-AR(0.6) correlated; 20 seeds each × {random_permutation, knockoff_equi};
n=100, 100 bootstraps, auto lambda grid (20 values). D4: environment-label
permutation null (10 seeds).

**Headline result — min_fdr_ is NOT a calibrated estimate of realized FDP in
several regimes** (mean realized FDP of the FDP+-selected support vs mean
reported min_fdr_):

| design | artificial | realized FDP | reported min_fdr_ | mean n_selected |
|---|---|---|---|---|
| D1_sparse | knockoff_equi | 0.458 | 0.473 | 26.4 |
| D1_sparse | rand_perm | 0.210 | 0.149 | 9.1 |
| D2_null | knockoff_equi | 0.918 | 0.423 | 18.2 |
| D2_null | rand_perm | 0.934 | 0.662 | 32.6 |
| D3_correlated | knockoff_equi | 0.354 | 0.522 | 27.6 |
| D3_correlated | rand_perm | 0.031 | 0.098 | 10.2 |

- **F-P6-1 (P1, headline):** on a PURE NULL (no signal), the machinery still
  selects 18-33 features with realized FDP 0.92-0.93 while reporting
  min_fdr_ 0.42-0.66 — the estimate is anti-conservative by 1.5-2.2× exactly
  where a user is most likely to trust an empty result. This is consistent
  with F-001's semantics (min over a random q+ curve; only the probabilistic
  bound holds) but quantifies the practical consequence: **an empty support
  is the only outcome a null dataset cannot produce by this rule**. The
  package's "diagnostic, not a guarantee" framing is accurate but the degree
  of null anti-conservatism is material and undocumented.
- **F-P6-2 (P2):** calibration is regime-dependent in BOTH directions:
  anti-conservative on independent nulls (both decoy types), approximately
  calibrated for knockoff_equi on independent signal (0.473 vs 0.458),
  anti-conservative for rand_perm on independent signal (0.149 vs 0.210),
  and conservative under block correlation (both types). No single
  interpretation of min_fdr_ is safe across regimes; the RP-vs-knockoff
  ordering predicted by F-005/L1 was NOT observed in the correlated design
  at this signal level (RP was the more conservative there) — the F-005
  mechanism is assumption-level, and its empirical signature is
  design-dependent.
- **F-P6-3 (positive):** D4 environment-label permutation null: max-over-
  features min_score = 0.168 (sd 0.071) — far below the 0.5 heuristic;
  `stabl_environment_stability` does not manufacture cross-environment
  stability on nulls.
- Caveats: weak-signal regime (per-feature effect 0.25, n=100), 100
  bootstraps, 20 seeds; realized-FDP-on-null conflates estimate quality with
  selection behavior; results characterize the DEFAULT threshold rule, not
  the bound's validity.

## Phase 10-11 — documentation & distribution audit (lead auditor, 2026-10-17)

- **Clean:** dependency hygiene is correct — Imports holds only true hard
  deps (glmnet, Matrix, matrixStats, methods, RColorBrewer, stats, utils,
  withr); all optional learners/generators (knockoff, ranger, aorsf,
  ordinalNet, sparsegl, ggplot2) are Suggests with `.require_pkg()` guards
  and informative messages. `License: GPL (>= 2)` is consistent with the
  vendored GPL-2 engine + `inst/COPYING.cooperative` attribution.
- **Clean:** NEWS.md preserves adverse results verbatim (ordinal graduation
  bar NOT met; oblique-RF fenced; redundancy co-selection gate NOT met;
  leaky full-data selection quantified) — exemplary honest changelog
  practice.
- **Clean:** 9 vignettes including the aging case study; man/ has 82 pages.
- **F-P10-1 (P2, carried):** the locked release contract is 0.1.1-era while
  DESCRIPTION says 0.3.0 (release-contract version skew, from the prior
  session record) — a distribution-governance defect to resolve before any
  release (Phase 12-15).
- **F-P10-2 (P3):** mixOmics/jsonlite/rlang appear in Suggests but are used
  only by `inst/analysis/release_evidence.R` and vignettes, not by R/.
  Acceptable practice; worth a comment in DESCRIPTION.

### Phase 6 addendum — winner's-curse mechanism + F-001 bound check (2026-10-17)

Mechanism demo (200 pure-null reps, exchangeable real/artificial scores,
package's own formula): the per-threshold estimator is approximately unbiased
(mean 1.022 at a fixed threshold), but the min over the 100-threshold curve
collapses to mean 0.674 with 35 features selected — the argmin hunts the
luckiest real-vs-artificial count crossing. The Phase 6 null anti-
conservatism (reported 0.42-0.66, realized 0.92-0.93) is this winner's curse
plus score-distribution effects, not a biased per-threshold estimator.

F-001 bound check on the actual D2-null runs (true FDP = 1.0, 40 seed-runs):
empirical P(FDP >= (1+Δ)·q+) = 0.60 / 0.47 / 0.17 / 0.03 at Δ = 0.5 / 1 / 2 /
5, against bounds 0.67 / 0.50 / 0.33 / 0.17 — **within the bound at every Δ,
and nearly tight at small Δ**. The observed anti-conservatism is therefore
consistent with the package's theory: the F-001 guarantee is a weak
Markov-type conditional tail bound, not a calibration statement, and it
permits almost exactly the observed behavior. Remediation implication for
F-P6-1: report the bound's practical meaning alongside min_fdr_ (e.g. "the
true FDP exceeds twice this value with probability up to ~50%") rather than
changing the threshold rule.

## Phase 7 — aging case-study adversarial audit (lead auditor, 2026-10-17)

Scope: `vignettes/stablr-aging-case-study.Rmd` (v1 historical case study +
v2 summary section) and the ws8_v2 saturation chain (claims-ledger V8, V12,
V13). Mandate from the leads table: re-derive the S0 null saturation from
committed artifacts and attribute the algorithmic cause. Evidence:
`phase7_validation/` (rederive.R + log; mechanism probe + log).

### Re-derivation — every checked number reproduces exactly

All values recomputed from the committed CSVs in the compendium
(`results/ws8_v2/`), not copied from the ledger:

- **V8 (S0 null, 8 cohorts):** 50/50 permutations select 59/59; max
  min-score median 1.000, range [0.980, 1.000]; same-design real-data
  reference 35/59 native. Matches the ledger and vignette §10.4 exactly.
- **V12 (S2a SPAC-excluded null):** 50/50 select 59/59; median 1.000,
  range [0.995, 1.000]; real reference 30/59. Matches.
- **V13 (S1 synthetic):** M4/M5 select 59/59 native features under all four
  truth regimes (homogeneous, heterogeneous, null, sparse-strong); empirical
  FDP of native panels 0.83–1.00; M1 selects 0 features under the null.
  Matches.
- **Vignette §10.3:** macro R@10 per selector (0.395/0.386/0.321/0.286/
  0.142/0.118) and all five M4 contrasts with CIs and sign-flip p-values
  reproduce exactly from `run/heldout_replication.csv` and
  `run/method_comparisons.csv`.
- **Vignette §10.5:** top-10 panel macro/micro MAE 12.19/13.56 vs all-59
  12.35/14.35 vs median baseline 20.11/22.74; predicted-on-observed slopes
  0.317–0.675 ("0.32 to 0.68"). Matches.
- **Vignette §4 (v1):** 24/60 selected; the eight predeclared min-agg
  scores (1.000, 1.000, 0.995, 0.755, 0.740, 0.535, 0.505, 0.445); 7/8
  predeclared selected; all eight sign-consistent 8/8; 11/24 consistent
  features significant. All reproduce from `benchmarks/ws8_feature_stats.csv`.

The v1 defects (A1–A4) remain flagged in place in the vignette; no
unflagged v1 claim was found to contradict the committed artifacts.

### F-P7-1 (P2, mechanism CONFIRMED by probe): the null saturation is caused by the rf threshold-path floor + max-over-lambda scoring + a decision rule that never consults FDP+

Algorithmic cause chain, verified in code and by a controlled probe
(`phase7_validation/phase7_mechanism_probe.R`):

1. **Grid floor at the median.** `stabl_fit(lambda_grid = "auto",
   base_learner = "rf")` resolves to `auto_rf_threshold_grid()`, whose path
   is `probs = seq(0.5, 0.999, length.out = n_grid)` — quantiles of the
   full-data permutation-importance distribution of the real features. The
   densest path point is therefore the **median** importance. Under a pure
   null this bar is numerically zero (probe: median importance −0.0003,
   max 0.036).
2. **Max-over-lambda scoring.** The per-cohort stability score is the
   max-over-thresholds of per-bootstrap selection frequencies
   (`environment_stability.R` aggregates max-over-lambda scores; the min
   across cohorts inherits them). At the densest thresholds the probe's
   mean selection frequency is 0.94–1.00 for reals AND decoys, so
   max-over-lambda scores are ≈ 1 for everything: reals median 1.000
   (min 0.980), decoys median 1.000. The probe reproduces the committed
   null signature exactly (min 0.980 = the committed range floor).
3. **The native panel rule ignores the decoy contrast.** The panel rule is
   a fixed bar, `min_score >= 0.5` — the FDP+ machinery is computed but
   never consulted. The decoy exchangeability itself is intact: per-lambda
   real/decoy mean frequencies agree within 0.02 at every threshold, which
   is why per-cohort FDP+ ≈ 1 (vignette Table 1: decoys ≥ 0.5 in every
   cohort). An FDP+-gated rule would select ~nothing under the null —
   exactly what M1 (per-cohort voting) does in S1 (0 features under s1c).
4. **Deeper layer — fixed thresholds do not calibrate refit noise.** A
   sparse-only path (quantiles 0.9–0.999) still yields max-over-lambda
   ≈ 0.76 under the null: a fixed threshold set at the full-data 90th
   percentile is exceeded in ~74% of bootstrap refits. Fixed full-data
   importance quantiles are not calibrated for bootstrap-refit importance
   scale. The min-across-cohorts aggregation is innocent (min of ≈1 is ≈1);
   the leniency is created per cohort.

Ruled out as mechanisms: SPAC's PD-projected knockoff covariance (V12:
saturation without SPAC) and the equicorrelated knockoff construction
(V15/S2c for R@10). The vignette's §10.4 statement that "the projection is
not the mechanism" is confirmed; the actual mechanism (items 1–3) was
previously unstated anywhere in the repo.

Remediation implications (for Phase 12, not silent fixes): (a) the
env-stability native-panel rule should either gate on the per-cohort FDP+
thresholds or document that it is a descriptive filter (the vignette now
says this, but the CODE presents `min_score >= 0.5` as a selection rule);
(b) `auto_rf_threshold_grid`'s median floor plus max-over-lambda scoring
should carry an explicit leniency warning analogous to the existing
random-permutation warning; (c) any future rf-based env-stability
calibration study must use FDP+-gated selection, not the fixed 0.5 bar.

### F-P7-2 (P3, CONFIRMED): inconsistent Clopper–Pearson conventions between V8 and V12

The ledger and vignette §10.4 quote two "Clopper–Pearson 95% lower bounds"
for the same 50/50 saturation statistic with different values: V8 → 0.929,
V12 → 0.942. The difference is purely conventional: 0.929 = 0.025^(1/50)
(two-sided convention, lower tail α/2) and 0.942 = 0.05^(1/50) (one-sided,
`binom.test(alternative = "greater")`, `s2a_null_merge.R:43`). V8's bound
is not computed in `null_probe.R` at all (provenance outside the committed
scripts). Neither the ledger nor the vignette states the convention.
Remediation: state the convention (or unify on one) in both places; the
saturation conclusion is unaffected (both bounds ≥ 0.92).

### F-P7-3 (P3, CONFIRMED): the WS-5 min-agg calibration does not cover the rf pairing used by the case study

Vignette §3.1 invokes WS-5's FDP+ 0.083–0.142 as the source of the min-agg
rule. WS-5 (`ws5_transportability.R`) measured that calibration with the
**lasso** learner, **random_permutation** artificials, binomial family,
p = 30, E = 3, 100 bootstraps. The case study applies the same aggregation
with **rf + knockoff_equi**, gaussian, p = 60, E = 8 (200 bootstraps) — a
learner/artificial-type/family/shape combination for which no min-agg
calibration existed before WS-8 v2's own S1 (which found saturation). The
vignette honestly reports non-transfer (§4.1, §8 limit 1) but does not note
that the cited calibration was measured on a different configuration.
Remediation: one clarifying sentence in §3.1 (Phase 12 batch).

### Scope boundaries

- S1 rests on the 5-replicate pilot of a planned 20 (flagged in the ledger
  and §10.8); extending to 20 reps is the already-designed weak-signal/
  p≫n sweep, deliberately not run inside Phase 7.
- The v1 pipeline's donor-level derivation script was never archived
  (vignette §2.3 discloses this); Phase 7 cannot audit that step. The v2
  rebuild regenerates features from the raw file with a committed script.
- The AI audit does not satisfy SCI-04; the independent human confirmation
  remains outstanding.

## Final pre-publication per-function review (2026-10-18)

All 74 exported functions audited (four delegated slices + two inline);
verdict rollup, the 42 new findings (N-01..N-06, B-01..B-08, C-1..C-14,
D-1..D-13), ensemble-coherence notes, defaults-calibration verdict, and
dispositions are recorded in `final-review.md` (not duplicated here).
Headline: one P1 found and fixed in the same commit (B-01 —
`auto_rf_threshold_grid` default `seed=NULL` always failed, propagating to
`stabl_fit(base_learner = "rf", lambda_grid = "auto")` without
`random_state`); five fail-closed robustness fixes with regression tests
(B-02, N-02, D-1, C-3, C-5); the advanced-vignette Jaccard-heatmap axis
defect fixed (D-3); nine documentation corrections; 27 findings deferred
to Phase 13 with rationale (remediation-matrix.md). The pre-existing
`test-public-interface-audit` failure (hard-coded 72-export list vs the
74 exports introduced by the sPLS-DA commit a5d1470) was also corrected
here. SCI-04 human confirmation remains outstanding.
