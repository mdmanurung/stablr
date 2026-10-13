# Feature-maturity roadmap (0.3.0 → mature)

Date: 2026-10-14 · Branch: `docs/pkgdown-0.3.0` · Status: PROPOSED

## Principle

Maturity is earned by evidence, not declared. A component is **mature** when
it has (a) validation or parity evidence on predeclared designs, (b) full
test coverage, (c) documented interpretation limits, and (d) a governed
release-evidence run. FDP+ calibration is a property of the
**(learner × artificial-decoy × design)** triple, not of a feature in
isolation — the roadmap therefore always names the design.

Current evidence base: testthat 2997/0/1; `R CMD check --as-cran` 0/0/3;
F5 RF-validation simulation (non-linear signals, n=120, p=10, 10 seeds,
100 bootstraps, nominal FDP+ 0.1); B1 study-disjoint benchmark (6 studies,
10 seeds); sealed 0.1.1 full-candidate OOL parity run (negative); SCI-04
gate `proposed_unaccepted`.

## Maturity matrix (current state)

| # | Component | Status | Evidence today | Gap to mature |
|---|-----------|--------|----------------|---------------|
| 1 | Core engine (`stabl_fit`, `auto_lambda_grid`, `compute_fdp_plus`) | **Mature** (linear-signal regimes) | Python parity (ranking/support contracts); full suite | None beyond WS-6 documentation |
| 2 | Lasso / elastic net / adaptive lasso adapters | **Mature** (linear regimes) | Parity + suite; F5: lasso+knockoff_equi FDP 0.578 on the non-linear design | WS-6: document the learner-regime boundary |
| 3 | Random-permutation decoys | **Mature, with documented anti-conservatism** | F5: permutation nulls lack correlation leakage → inflated FDP (rf 0.322); guard rationale confirmed | Fold per-learner table into WS-6 |
| 4 | RF adapter (`ranger`) | **Calibrated-experimental** | F5 full: rf+knockoff_equi FDP 0.065, TPR 1.000 (10/10 seeds) | WS-2: broaden designs, codify permutation-path guard, promote |
| 5 | Oblique RF adapter (`aorsf`) | **Experimental — not calibrated** | F5 full: knockoff_equi 0.242, permutation 0.387 vs nominal 0.10 | WS-3: root-cause, recalibrate or fence scope |
| 6 | Sparse group lasso | **Mature** (linear regimes) | Suite + parity scope | Include in WS-2 calibration matrix |
| 7 | Multi-omic early/late fusion | **Mature** (software); **evidence-gated** (claims) | OOF leakage-safe defaults; suite; OOL tutorial | WS-1: full-candidate parity run has never passed |
| 8 | Cooperative fusion (gaussian/binomial) | **Shipped; benchmark gap** | Suite; vendored multiview engine; vignette | WS-4: fusion-regime benchmark; Cox deferred |
| 9 | Cox cooperative fusion | **Deferred** | — | WS-4: implement or record permanent boundary |
| 10 | Environment-stability workflow | **Experimental** | Suite only | WS-5: transportability simulation |
| 11 | Study-disjoint nested CV (B1) | **Evidence-backed** | B1: null probe gap ≈ 0 (−0.0096); prevalence-drift gap +0.0395 | Fold into WS-1 release evidence |
| 12 | Knockoffs > 3,000 predictors | **Fail-closed boundary** | Documented; no globally valid chunked construction | WS-7: research or permanent boundary |
| 13 | SCI-04 release gate | **`proposed_unaccepted`** | Packet merged (PR #7); reference tests green | WS-0: external acceptance + n ≥ 150 run |

## Workstreams

### WS-0 — Release gate (blocking; governance)
External reviewer completes the SCI-04 packet (6 assumptions + response
template) → immutable accepted decision record → contract flip
(`promotable`, replicates 100 → ≥ 150; margin(150) = 0.0999 ≤ 0.10) →
`release_evidence.R` check swap → fresh locked release run at accepted n →
win-builder/sanitizer → finalize cran-comments → CRAN submit.
**Exit:** gate accepted; CRAN submission sent. **Owner:** maintainer +
external reviewer (cannot be self-accepted). **Wall-clock:** external.

### WS-1 — Full-candidate OOL parity (must-pass for release claims)
The frozen protocol failed on the 0.1.1 candidate: support Jaccard 0.778
(≥ 0.80 required), score Spearman 0.785 (≥ 0.90), FDP-threshold delta 0.08
(≤ 0.05); the actual-generator gate failed via the typed
`random_permutation` fallback (knockoff infeasible, n > p required).
Tasks: (a) root-cause the score-correlation gap — lambda-grid alignment,
solver tolerance, preprocessing order; (b) decide per-gate remediation;
(c) re-run the frozen protocol on the 0.3.0 candidate from a hashed tarball.
**Exit:** all gates pass on a sealed run — or gates re-predeclared with
justification through the governance path. **Effort:** 1–2 diagnosis
sessions + reruns (45–90 min each, sandbox-scale).

### WS-2 — RF adapter graduation (plan item D1)
Broaden F5 beyond the single non-linear design: n × p grid, correlated
nulls, gaussian + binomial + Cox, both decoy types; predeclare the design
matrix before running. Codify the permutation-path guard (default-on,
documented rationale: permuted nulls do not receive correlation leakage).
Promote `make_ranger_adapter` from experimental to calibrated with an
explicit design scope in man pages and pkgdown.
**Exit:** predeclared calibration study committed; labels updated.
**Effort:** 2–3 benchmark sessions.

### WS-3 — Oblique RF adapter: fix or fence
Root-cause FDP 0.242 (knockoff_equi) / 0.387 (permutation). Hypotheses:
noisier per-feature selection frequencies under heterogeneous predictor
scales; threshold instability at small n_selected. Candidate remedies:
scale standardization before frequency aggregation, wider bootstrap counts,
threshold-pool smoothing — each predeclared before the next F5-style run.
If recalibration fails, fence: restrict documented scope and mark
not-recommended for the affected designs.
**Exit:** calibrated on predeclared designs, or fenced with explicit scope.
**Effort:** 2–4 sessions; genuine research risk — may honestly end at
"fenced".

### WS-4 — Cooperative fusion completion
(a) Benchmark cooperative vs early vs late fusion across view-correlation
regimes (complementary vs redundant views), extending the F5 pattern to
multi-omic designs. (b) Cox cooperative fusion: implement on the vendored
engine or record a permanent deferral decision with rationale.
**Exit:** benchmark committed; Cox shipped or boundary recorded.
**Effort:** 2–3 sessions.

### WS-5 — Environment-stability graduation
Transportability simulation using B1's isolated-mechanism pattern: location
shifts, prevalence drift, and mixed mechanisms across environments;
evaluate whether min-aggregated selection frequencies control FDP+ on
transported cohorts. Document when the workflow is calibrated.
**Exit:** benchmark + docs; promote from experimental.
**Effort:** 1–2 sessions.

### WS-6 — Learner-regime documentation
Characterize the lasso+knockoff_equi 0.578 row: a main-effects linear
learner on non-linear signals leaves structure for decoys to absorb, so
FDP+ is anti-conservative — this is a learner-regime property, not a bug.
Ship a "which learner for which signal" guidance table (docs + vignette
cross-links).
**Exit:** guidance shipped. **Effort:** 1 session.

### WS-7 — Knockoff high-dimension boundary
Survey chunked / ghost-knockoff constructions for p > 3,000; implement
behind fail-closed validation if a globally valid construction can be
established, otherwise record the boundary as permanent in docs.
**Exit:** decision record. **Effort:** open-ended (research).

## Sequencing

- **Now, in parallel:** WS-0 (external review, wall-clock-bound) and WS-1
  diagnosis (engineering-bound).
- **Next:** WS-2, WS-3 (benchmark compute; WS-3 may inform WS-2's guard
  policy), WS-6 (cheap, unblocks honest messaging).
- **Then:** WS-4, WS-5.
- **Last:** WS-7 (research-y; the fail-closed boundary is safe to keep).

## What maturity does NOT mean

FDP+ is a diagnostic minimizer, not a universal false-discovery guarantee;
calibration results are design-relative and must be cited with their design.
These limits stay in the README, man pages, and vignettes regardless of
workstream outcomes.
