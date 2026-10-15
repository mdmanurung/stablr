# Predeclared design — interpretation & validated-workflow benchmark (IV)

Date: 2026-10-15. Written BEFORE any run of this benchmark, per the
release-governance pattern (predeclare, then run, then report — including
adverse outcomes). This benchmark supports Gate B (binary MVP): it exercises
the PR1 interpretation layer and the PR3 `stabl_analyze()` exploratory and
validated modes on synthetic data with known ground truth.

Scope and honesty caveats (fixed):

- Simulation-only evidence; no application-data claim.
- FDP+ remains a learner x decoy x design diagnostic. Nothing here asserts
  universal finite-sample FDR control.
- The leaky-vs-nested contrast (IV-5) is directional on average across
  seeds. We do NOT demand that every random replicate displays a fixed
  optimism gap.
- Empirical FDP/TPR, where reported, use the known simulated truth set and
  are distinct from the fitted FDP+ estimate.

## Common settings

- Binary outcomes throughout: `eta = x %*% beta`, centred so mean
  prevalence is 0.5; `y ~ Bernoulli(plogis(eta))`, levels A/B with B the
  event class.
- Seeds 1..10 per design; `n_bootstraps = 50`; `artificial_type =
  "random_permutation"`; `hard_threshold = NULL`; `lambda_grid =
  data.frame(lambda = c(0.3, 0.2, 0.1))`.
- Validated arm: `stabl_analyze(mode = "validated")` with `outer_v = 3`,
  `inner_v = 3`, `preprocessing = c("center", "scale")`, `n_repeats = 10`
  for held-out importance. IV-4 uses study-disjoint folds (`outer_v = 4`
  leave-one-study-out, `inner_v = 2`).
- Exploratory arm: `stabl_analyze(mode = "exploratory")` on the full data.
- Arm coverage by design: IV-0..IV-3 run both the exploratory and the
  validated arm; IV-4 runs the validated arm only (study-disjoint); IV-5
  runs the leaky and nested arms only. Planned runs: 4 x 2 x 10 + 10 +
  20 = 110.

## Design matrix

| Design | n | p | Truth | Purpose |
|---|---|---|-------|---------|
| IV-0 null | 120 | 10 | none (y independent) | null behaviour of validated metrics |
| IV-1 sparse | 120 | 10 | beta = (1.5, -1.2, 1.0) on f1..f3 | effect-sign recovery; signal vs noise importance |
| IV-2 substitutes | 120 | 10 | f1a (1.5) with correlated copies f1b/f1c (rho 0.85); f2 (-1.2) | correlation grouping; group vs individual importance |
| IV-3 single | 80 | 5 | f1 only (1.5) | one-feature downstream path in validated mode |
| IV-4 study shift | 160 | 10 | f1 (1.5), f2 (-1.2); study intercept shifts (-1, -0.3, 0.3, 1); study 4 feature scale x3 | study-disjoint validated mode under prevalence/scale drift |
| IV-5 leaky vs nested | IV-1 data | as IV-1 | leaky arm: selection on FULL data, then fold-wise downstream fit on the leaked panel over the SAME outer folds; nested arm: `stabl_analyze` validated | directional optimism of leaky selection |

The IV-5 leaky arm reproduces the nested arm's exact outer folds by calling
the package's internal fold generator with the same seed
(`.make_repeated_cv_folds(y, v = 3, repeats = 1, stratified = TRUE,
random_state = seed)`), so the comparison isolates selection leakage only.

## Recorded per run

Design, seed, arm (exploratory / validated / leaky); selected feature
identities and count; per-feature stability; exploratory coefficients and
effect directions; correlation-group assignments; validated pooled and
per-repeat accuracy, BER, macro-F1, MCC, Brier, clipped log-loss, AUROC;
held-out importance summary (mean/sd per feature, folds selected); error
class/message for failed runs; runtime; package version; git SHA.

## Gate criteria (fixed before running)

- G-IV-1 completeness: every planned design x seed x arm run is accounted
  for in `iv_results.csv` or `iv_failures.csv`; the benchmark exits
  non-zero if any planned run is missing.
- G-IV-2 effect-direction recovery (IV-1): among selected signal features,
  the exploratory `effect_direction` matches the true sign in >= 90% of
  cases pooled across seeds.
- G-IV-3 redundancy grouping (IV-2): f1a/f1b/f1c share one correlation
  group in >= 80% of exploratory runs.
- G-IV-4 leakage direction (IV-5): leaky-arm mean AUROC >= nested-arm mean
  AUROC across seeds (pooled directional check, not per-replicate).
- G-IV-5 null sanity (IV-0): no run errors; median validated AUROC within
  [0.35, 0.65] (a loose sanity band, not a calibration claim).
- G-IV-6 one-feature path (IV-3): all seeds complete; whenever exactly one
  feature is selected in a fold, that fold's event probabilities are
  non-degenerate (not all equal to the training prevalence).

Outcome reporting: each criterion is reported as MET / NOT MET with the
actual numbers. Adverse results are committed, not suppressed. Graduation
language for Gate B is reserved for the Gate B report, which combines this
benchmark with the full test-suite evidence and the WS-8/WS-9 repairs.

## Addendum 2026-10-15 (post-run repairs; design, arms, and gates unchanged)

An initial execution of this exact predeclared design exposed two
implementation defects. Neither the design matrix, the arms, the seeds,
nor the gate criteria were changed; the benchmark was rerun end-to-end
after the repairs and only the rerun outputs are committed.

1. Package bug (repaired in PR0d, commit 2eb6929): folds with an empty
   selection aborted `.analyze_fold_importance()` ("arguments imply
   differing number of rows: 1, 0"), failing 3 of 110 runs (iv0_null
   validated, seeds 6/9/10). The empty-selection fold row is now a true
   zero-row data.frame; a regression test reproduces the exact IV-0
   seed-6 settings. Rerun completed 110/110.
2. Summary-table defect (benchmark script only): the per-design x arm
   mean table silently dropped arms that record only a subset of metrics
   (`aggregate`'s default `na.omit`). Fixed with `na.action = na.pass`;
   gate values are computed from the raw per-run results and were
   unaffected. Outputs were regenerated from an identical deterministic
   rerun.
