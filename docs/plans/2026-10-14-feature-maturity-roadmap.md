# Feature-maturity roadmap (0.3.0 → mature)

Date: 2026-10-14 · Branch: `docs/pkgdown-0.3.0` · Status: EXECUTED (2026-10-14; outcomes below; WS-0 external, WS-1 environment-blocked)

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
| 4 | RF adapter (`ranger`) | **Calibrated (scoped)** — WS-2 done | WS-2 (60 runs): rf+knockoff_equi FDP 0.025 (linear) / 0.065 (nonlinear p=10); NOT calibrated at p=50 (0.44) or binomial (0.26); scope in man pages | Scope restriction documented; no further work |
| 5 | Oblique RF adapter (`aorsf`) | **Fenced** — WS-3 done | WS-3: best remedy (boot 200) FDP 0.184/0.392 vs predeclared 0.15; decoy importance inflated (mean 0.26 vs rf 0.23) | Not recommended with knockoff decoys; use rf |
| 6 | Sparse group lasso | **Mature** (linear regimes) | Suite + parity scope | Include in WS-2 calibration matrix |
| 7 | Multi-omic early/late fusion | **Mature** (software); **evidence-gated** (claims) | OOF leakage-safe defaults; suite; OOL tutorial | WS-1: full-candidate parity run has never passed |
| 8 | Cooperative fusion (gaussian/binomial) | **Shipped; benchmarked** — WS-4a done | WS-4 (60 runs): no mode dominates; cooperative competitive (0.24–0.26 BER); late-fusion gaussian weakness documented (1.09 MSE) | Late-gaussian caveat in docs |
| 9 | Cox cooperative fusion | **Permanently deferred** — WS-4b done | No CV-loss precedent for Cox cooperative learning; fail-closed stops at 3 layers | Cox users keep per-omic + early fusion paths |
| 10 | Environment-stability workflow | **Promoted (scoped)** — WS-5 done | WS-5 (60 runs): min-agg FDP+ 0.08–0.14 across shift/drift designs; pooled drifts to 0.16; per-env rows anti-conservative (0.28–0.54) | Min-agg is the calibrated read; documented |
| 11 | Ordinal regression learner (WS-9, added 2026-10-14) | **Experimental (benchmarked; bar not met)** | WS-9 (40 runs, proportional-odds ground truth): knockoff_equi mean FDP+ 0.00 (nominal 0.1) but mean TPR 0.37 on both designs — strongly conservative; predeclared graduation bar (FDP+ <= 0.15 AND TPR >= 0.8) NOT met | Stays experimental; scope documented in ?make_ordinal_adapter and benchmarks/ws9_summary.md |
| 11 | Study-disjoint nested CV (B1) | **Evidence-backed** | B1: null probe gap ≈ 0 (−0.0096); prevalence-drift gap +0.0395 | Fold into WS-1 release evidence |
| 12 | Knockoffs > 3,000 predictors | **Permanent boundary** — WS-7 done | Decision record: docs/plans/2026-10-14-ws7-knockoff-boundary-decision.md | Reopen only with new literature-backed construction |
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
**OUTCOME (2026-10-14): ENVIRONMENT-BLOCKED.** The sealed Python reference
outputs live in the separate stablr-experiments workspace with its own HPC
runtime; they are not present in this repository. The bundled
python_metrics_parity gate passes. No speculative root-cause is recorded
without the sealed artifacts. Release claims stay gated on WS-0/WS-1.

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
**OUTCOME (2026-10-14): NOT MET GLOBALLY — SCOPED PROMOTION.** 60 seed-runs,
zero errors (benchmarks/ws2_summary.md). rf + knockoff_equi meets the
predeclared criterion (mean FDP <= 0.15, TPR >= 0.8) on A_linear (0.025) and
B_nonlinear (0.065) but fails C_highdim_nonlinear (0.44; 3/10 seeds select
~all features) and D_nonlinear_binomial (0.26). Disposition: promote rf with
an explicit design scope (linear and low-dimensional nonlinear, gaussian);
document not-recommended regimes; no removal (literature-backed, competitive
where calibrated).

Broaden F5 beyond the single non-linear design: n × p grid, correlated
nulls, gaussian + binomial + Cox, both decoy types; predeclare the design
matrix before running. Codify the permutation-path guard (default-on,
documented rationale: permuted nulls do not receive correlation leakage).
Promote `make_ranger_adapter` from experimental to calibrated with an
explicit design scope in man pages and pkgdown.
**Exit:** predeclared calibration study committed; labels updated.
**Effort:** 2–3 benchmark sessions.

### WS-3 — Oblique RF adapter: fix or fence
**OUTCOME (2026-10-14): FENCED.** All predeclared remedies failed the
criterion (benchmarks/ws3_summary.md): baseline 0.242/0.450,
R1 standardize (no change), R2 boot-200 0.184/0.392, R3 both 0.184/0.392 —
none <= 0.15 on both designs. Diagnosis: knockoff decoy selection-frequency
noise is inflated for oblique forests (decoy mean importance 0.26 vs rf
0.23; decoy max 0.60-0.66). Disposition: oblique_rf is fenced — not
recommended with knockoff decoys; rf is the recommended knockoff-calibrated
learner.

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
**OUTCOME (2026-10-14): BENCHMARK DONE; COX PERMANENTLY DEFERRED.**
(a) 60-run fusion benchmark (benchmarks/ws4_summary.md): no mode
dominates. Complementary views: early 0.234 / late 0.264 / cooperative 0.243
BER. Redundant views: cooperative best on gaussian (0.258) and competitive on
binomial (0.309). Late-fusion gaussian is the weakness (1.09 MSE, 4x worse).
Cooperative fusion is KEPT (competitive, literature-backed). (b) Cox
cooperative fusion: the vendored engine fits Cox paths, but no cv.coxnet-style
CV loss exists; explicit fail-closed stops guard 3 layers. Implementing one
is novel method development without literature precedent — recorded as a
permanent deferral; Cox users keep per-omic and early-fusion paths.

(a) Benchmark cooperative vs early vs late fusion across view-correlation
regimes (complementary vs redundant views), extending the F5 pattern to
multi-omic designs. (b) Cox cooperative fusion: implement on the vendored
engine or record a permanent deferral decision with rationale.
**Exit:** benchmark committed; Cox shipped or boundary recorded.
**Effort:** 2–3 sessions.

### WS-5 — Environment-stability graduation
**OUTCOME (2026-10-14): PROMOTED WITH DOCUMENTED SCOPE.** 60-run
transportability simulation (benchmarks/ws5_summary.md), corrected
selection contract (min score >= 0.5). Min-aggregation controls FDP+ across
null-shift (0.083-0.142), signal-shift (0.075-0.163), and prevalence-drift
(0.086-0.142) designs; pooled aggregation drifts to 0.163. Per-environment
permutation rows are anti-conservative (0.28-0.54), confirming the documented
caveat. TPR cost of min-aggregation (0.53-0.73 vs pooled 0.90-0.93) is the
transportability price. Disposition: promote stabl_environment_stability from
experimental; min-agg is the calibrated read.

Transportability simulation using B1's isolated-mechanism pattern: location
shifts, prevalence drift, and mixed mechanisms across environments;
evaluate whether min-aggregated selection frequencies control FDP+ on
transported cohorts. Document when the workflow is calibrated.
**Exit:** benchmark + docs; promote from experimental.
**Effort:** 1–2 sessions.

### WS-6 — Learner-regime documentation
**OUTCOME (2026-10-14): SHIPPED.** Guidance table added to the intro
vignette, section 7.3, grounded in F5 + WS-2 + WS-3 measured numbers.

Characterize the lasso+knockoff_equi 0.578 row: a main-effects linear
learner on non-linear signals leaves structure for decoys to absorb, so
FDP+ is anti-conservative — this is a learner-regime property, not a bug.
Ship a "which learner for which signal" guidance table (docs + vignette
cross-links).
**Exit:** guidance shipped. **Effort:** 1 session.

### WS-7 — Knockoff high-dimension boundary
**OUTCOME (2026-10-14): DECIDED — PERMANENT BOUNDARY.** Decision record:
docs/plans/2026-10-14-ws7-knockoff-boundary-decision.md (literature
grounding: Candes 2018; Bates 2020; Barber-Candes-Samworth 2020; Blain 2024;
Spector-Janson 2020).

Survey chunked / ghost-knockoff constructions for p > 3,000; implement
behind fail-closed validation if a globally valid construction can be
established, otherwise record the boundary as permanent in docs.
**Exit:** decision record. **Effort:** open-ended (research).

## Sequencing

- **Done (2026-10-14):** WS-2 (scoped promotion), WS-3 (fenced), WS-4
  (benchmark + Cox deferral), WS-5 (promoted), WS-6 (shipped), WS-7
  (permanent boundary).
- **Open, external:** WS-0 (SCI-04 reviewer; wall-clock-bound).
- **Open, environment-blocked:** WS-1 (sealed reference outputs live outside
  this repository).

## What maturity does NOT mean

FDP+ is a diagnostic minimizer, not a universal false-discovery guarantee;
calibration results are design-relative and must be cited with their design.
These limits stay in the README, man pages, and vignettes regardless of
workstream outcomes.
