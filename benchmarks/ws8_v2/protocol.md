# WS-8 v2 — Locked Protocol

**Status:** Locked at M0 (2026-10-15). Prospective revision of the exploratory WS-8
case study. Not a retrospective preregistration. Changes after this lock require a
numbered amendment appended to `ws8_showcase_predesign.md` and a note here.

**Source snapshot:** repository `mdmanurung/stablr`, branch `ws8-v2`, base commit
`1ff3cd33c4b084201e481a88177031cc2c2ad432` (main after PR #14). Working tree clean
except untracked `benchmarks/ws8_loco_silo_objects.rds` (committed under v2 for
provenance).

## 1. Aims and hypotheses

Primary aim: identify immune-cell composition features with reproducible associations
with chronological age across eight cohorts, and determine whether environment-aware
STABL improves independent cross-cohort replication and parsimony versus simpler
selectors.

- **H1 (primary):** a reproducible immune-aging core exists across the eight cohorts.
- **H2 (primary):** environment-aware STABL improves held-out association replication
  at a fixed feature budget versus simpler selectors.
- **H3 (secondary):** a sparse STABL panel retains age prediction in unseen cohorts.
- **H4 (secondary):** a reproducible core coexists with cohort-dependent components.
- **H5 (validation):** observed selection is not mostly produced by permissive tuning
  paths, correlations, or small-cohort instability.

Decision hierarchy: H1/H2 decisive; H3/sparsity/biology secondary. No method
superiority from a favorable secondary endpoint when the primary comparison fails.

## 2. Data

- Source: Terekhova et al. (2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012.
  Public portal file `all_cells_metadata.csv.gz`, MD5
  `50651968855bf6ca3efcdf199e2a9fd5`, retrieved 2026-10-13.
- Unit of analysis: donor. 2,609 donors, 8 cohorts (ABF300 166, AIDA 625, CIMA 421,
  KCL 188, OneK1K 981, SoundLife 96, SPAC 33, UCSF 99).
- Feature construction: per-donor per-cell-type proportion `(count + 0.5)/(total +
  0.5)`, logit transform. The `Bad quality` annotation class is excluded from both
  numerator and denominator (now verifiable from the raw file; the v1 derivation
  script was not archived — v2 rebuilds it).
- **Primary feature universe: 59 harmonized cell-type proportions only.** Sex, CMV,
  and ancestry are NOT discovery biomarkers; they enter only as secondary
  adjustment/stratification. Consequence (recorded assumption): v2 primary results are
  NOT directly comparable to v1's 60-feature (incl. sex) 24-feature panel.
- Complete-case exclusion uses only variables required for the specific analysis;
  donors are not dropped from the composition-age analysis for missing CMV.
- Outcome: donor age in years (gaussian). Environments: the 8 cohorts.

## 3. Primary endpoint (locked)

**Cohort-balanced held-out directional association replication at matched budget
K = 10.**

For each held-out cohort `e`:
1. Select/rank features using only the other seven cohorts.
2. Direction `d_j ∈ {-1,+1}` from the training-only association meta-estimate.
   Indeterminate (|rho| < 1e-6 or exact tie) → zero replication contribution, flagged.
3. Held-out Spearman `rho_{j,e}` between feature and age.
4. `R_e(K) = mean_{j in S_e(K)}(d_j * rho_{j,e})`.
5. Macro estimate `R(K) = mean_e R_e(K)` over the 8 cohorts.

Locked sub-rules:
- **K:** primary 10; sensitivity {5, 15, 20}.
- **M4 ranking:** min_score descending, tie-break mean_score descending, then feature
  name (resolves protocol gap E2). Native M4 panel = min_score >= 0.5 (unchanged).
- **M1 vote:** significant (BH-adjusted) in >= 5 of 7 training cohorts for native
  selection; rank by vote count then meta-evidence.
- **Multiplicity:** the primary confirmatory contrasts are M4 vs M0 (meta-analysis)
  and M4 vs M3 (pooled STABL) at K=10. All other K values and methods are exploratory.
- Report all eight per-cohort values, sign-agreement proportion, and rank magnitudes.
  Report negative/zero R faithfully.

## 4. Selector matrix (locked)

| ID | Selector | Training-only rule |
|---|---|---|
| M0 | Pooled age-association meta-analysis | Signed random-effects meta estimate across 7 training cohorts |
| M1 | Original-style cohort significance vote | Per-cohort test, BH, >= 5/7 vote; rank by count + meta |
| M2 | Pooled sparse multivariable (lasso/EN stability) | Same resampling budget as STABL where feasible |
| M3 | Pooled STABL | STABL ignoring cohort partition |
| M4 | Environment-aware STABL (PRIMARY) | Per-environment STABL, min-across-cohorts aggregation |
| M5 | Env-aware STABL, alternative aggregation (secondary) | Majority/quantile aggregation |

Controlled pair for attribution: M3 vs M4 share base learner (rf), artificial features
(knockoff_equi), bootstrap budget (200), and candidate K; only cohort aggregation
differs. M0/M1/M2 are scientific comparators, not isolation controls.

## 5. Predictive benchmark (secondary, H3)

- Ridge (`glmnet`, alpha = 0) on selected features; regularization tuned by
  training-only, cohort-respecting inner CV (leave-one-training-cohort-out preferred;
  documented fallback to random donor folds if infeasible).
- Baselines: training-median-age predictor; all-59-feature ridge; source-paper panel
  only if its full feature list is documented (see §8).
- Metrics: per-cohort and macro (primary) + micro (companion) MAE, RMSE, mean bias,
  rank correlation, and BOTH calibration slopes (observed-on-predicted and
  predicted-on-observed), explicitly named.
- Paired macro-MAE differences (STABL minus comparator) with cohort-level uncertainty
  as primary; donor-level stratified bootstrap as companion. No inferring paired
  significance from overlapping marginal CIs.
- Zero-shot target prediction is the primary predictive setting. Any few-shot
  recalibration uses a predefined target-cohort calibration subset and a disjoint
  untouched evaluation set, separately labeled.
- Empty native panel in a fold → training-median-age baseline for the native-panel
  comparison, `n_selected = 0` recorded, no silent top-K substitution.

## 6. Leakage containment (locked)

Eight fixed outer LOCO folds; identical donor IDs across methods for paired
comparison. All data-adaptive preprocessing (imputation, scaling, filtering,
ancestry-category collapsing, supervised filtering, covariance/knockoff fit) is fit on
the seven training cohorts only. No test age labels enter selection, tuning, sign
assignment, or model choice. Held-out predictions are written once and immutable.
Per-fold manifest records train/test IDs (or hashes), preprocessing fit IDs, method
config, selected features, seeds, software versions. Automated provenance check fails
if a test ID appears in any upstream fit.

## 7. Inference rules (locked)

- The 8 cohorts are the generalization units; cohort-level paired estimates are the
  primary uncertainty basis; show all eight.
- CIs on mean paired cohort differences state assumptions and the n=8 imprecision;
  exact paired sign-flip/permutation tests as sensitivity; no robust population-wide
  generalization claim.
- Clopper–Pearson intervals over 200 bootstrap selections are conditional resampling
  summaries, not CIs for independent-cohort reproducibility.

## 8. Source-paper comparator (locked)

- Exact reproduction of the 20/55 clock panels requires Table S1, which is paywalled
  and was NOT retrievable at M0. **Exact source-paper comparator is unavailable.**
- v2 uses a **conceptual** comparator: the paper's per-cohort association rule
  (linear model `% ~ age + sex + ancestry`, BH-adjusted, "consistent" = significant in
  >= 5 of 8 cohorts), adapted to 7 training cohorts as >= 5/7. Not called an exact
  replication.
- The paper's 8.2 y (full) and 8.9 y (reduced) test MAEs are context only, protocols
  clearly different; no head-to-head superiority assertion.
- The eight canonical populations (CD8 naive, CD8 RTE, CD4 RTE, CD4 Temra, MAIT, CD8
  Temra, CD8 Tem GZMK+, CD8 Tem GZMB+) are an a priori biological concordance set, not
  the paper's complete clock.

## 9. Null probes and simulations (locked)

- **S0:** within-cohort age permutation (preserves features, cohort structure, each
  cohort's age distribution), full selector refit per permutation, independent seeds.
  Smoke run of 10 permutations first to measure per-permutation runtime, then scale to
  the largest count (target 50–100) fitting a ~4 h chunked background cap; report
  Monte Carlo uncertainty. Descriptive + a predefined empirical null statistic; kept
  separate from calibrated FDR guarantees.
- **S1:** known-signal DGPs — universal linear, cohort-specific, correlated
  surrogates, compositional (Dirichlet-multinomial), 8-cohort size imbalance incl.
  SPAC p>n. Readouts: true/false selection, group recovery, empirical FDP and power,
  prediction degradation under shift. Diagnostic FDP+ labeled separately from true
  simulation FDP.
- **S2:** RF grid density/range; thresholds 0.4/0.5/0.6/0.7; leave-one-cohort-out
  selection sensitivity; several seeds; knockoff covariance diagnostics per cohort.
- No promotion of a new threshold/learner from S0–S2 into the primary analysis without
  an explicit amendment and separate validation.

## 10. Outputs and claims gate

All v2 outputs under `benchmarks/ws8_v2/` with a versioned `results/` subfolder.
Required tables: `data_manifest`, `folds`, `feature_scores`, `heldout_replication`,
`predictions`, `prediction_metrics`, `method_comparisons`, `null_results`,
`simulation_results`, `run_manifest.json`, `claims_ledger.md`. No claim of
"validated"/"replicated"/"FDR controlled"/"out-of-cohort improvement" without a linked
executed result and defensible design. A credible negative finding is a valid outcome.

## Amendments

- **Amendment 1 (2026-10-15)** — code conformance and M7 review remediation:
  M0/M1 implementations brought into conformance with locked §4/§8; S0 null
  extended to the locked 50-permutation target with a same-design real-data
  reference fit; §6/§10 persistence gaps closed; §7 interval labeling fixed;
  §5 empty-native fallback corrected to the median-age baseline. Full text
  appended to `ws8_showcase_predesign.md`; dispositions in
  `review_decisions.md`.
