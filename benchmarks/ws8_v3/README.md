# WS-8 v3 — Candidate mechanisms for cross-cohort reproducible biomarker discovery

Benchmark of three candidate selection mechanisms (plus the mandatory
batch-correction baseline) against the locked WS-8 v2 method matrix (M0–M5),
using the Terekhova et al. (2026, *Immunity*) multi-cohort immune-aging atlas
as the application case study.

**Application case study — NOT calibration evidence.**

## Arms

| ID | Selector | Role |
|----|----------|------|
| M0 | RE-meta ranking (DL random-effects on per-cohort Spearman, Fisher-z) | v2 baseline (strongest v2 arm) |
| M1 | Cohort significance vote (adjusted LM, >=5/7 BH) | v2 baseline |
| M2 | Pooled lasso stability | v2 baseline |
| M3 | Pooled STABL (rf, knockoff_equi) | v2 baseline |
| M4 | Env-aware STABL, min aggregation | v2 primary |
| M5 | Env-aware STABL, majority vote | v2 secondary |
| M6 | RE-meta + homogeneity filter (BH<0.05 AND Cochran's Q p>=0.10) | Candidate 1 |
| M7 | MINT-STABL (`mixOmics::mint.spls` regression inside STABL bootstrap + knockoff FDP+) | Candidate 2 |
| M7c | Tuned MINT, no stability | M7 control |
| M8 | Contrast-penalized lasso (common + cohort-deviation terms, omega=2) in STABL | Candidate 3 |
| M9 | ComBat(~age) harmonization + pooled STABL; held-out transform `eb_newbatch` | mandatory baseline |
| M9s | M9 with strict `standardize_only` held-out transform | transform sensitivity |
| M9b | M9 with `mod = NULL` (no age protection) | post-hoc sensitivity motivated by the M9 prediction failure mechanism |

## Endpoints

- **Primary**: held-out signed-rho replication R(K), macro over 8 LOCO folds,
  K=10 (sensitivity 5/15/20). Signs from training cohorts only.
- **Co-primary**: recovery of the paper's 20-population reduced clock panel
  (authors' repo notebook, `Age_clock_model/clean_notebook_for_aging_models.ipynb`,
  cell 32 `selected_by_marina`; all 20 map onto the harmonized 59 features).
  Concordance target, not ground truth.
- **Secondary**: cohort-respecting ridge prediction (macro MAE), reference
  population sign concordance, GZMK+/GZMB+ Tem components, null probes
  (M7/M8, 10 within-cohort age permutations, full n_boot=200).

## Package version policy (important)

All v3 numbers are produced with stablr pinned at commit
`a8f3a58d74f17658fd56467e0355aade1ad667ab`. The historical v2 results were
produced at earlier commits (M4: `618e8e3`, M3: `6d19afb`, ...); the package
fitting path changed materially between those commits and the pin, so v2
STABL-arm numbers do NOT reproduce under the pin (M4 macro R@10 0.288 ->
0.327). Deterministic arms (M0/M1/M2) reproduce bit-exactly. The v3 matrix
therefore reruns M3/M4/M5 like-for-like on the pin; v2 numbers are reported
only as historical reference (`version_drift_reference.csv`).

## ComBat implementation note

`selectors_v3.R` reimplements parametric ComBat (Johnson et al. 2007) with the
fitted parameters exposed, matching `sva::ComBat` numerics to <1.5e-6 on the
unit test (per-batch EB hyperpriors, iterative `it.sol` MAP updates, sva's
exact stopping rule). Held-out application is label-free: the new batch's own
moments and per-batch priors are estimated from its feature values
(`eb_newbatch`); only the grand mean and pooled variance are frozen from
training. With a covariate-free fit this reproduces the training adjustment to
<1e-8 (test T2c). With `mod = ~age` the covariate add-back term cannot be
applied to unlabeled new data — the mechanism behind M9's held-out prediction
failure (see report).

## Files

- `selectors_v3.R` — M6–M9b selectors (v2 contract) + ComBat fit/apply.
- `run_v3.R` — LOCO driver (per-(method, fold) checkpoints, extras persistence,
  `test_transform` hook for M9-family prediction).
- `null_probe_v3.R` — within-cohort age-permutation nulls (M7/M8) + real
  full-atlas fits (perm 0).
- `full_atlas_v3.R` — full-atlas fits for the paper-comparison co-primary
  (M5 derived from M4's fit in-session; M9s shares M9's selection).
- `assemble_v3.R` — primary/contrast/recovery/prediction/null tables.
- `figures_v3.R`, `heterogeneity_map_v3.R` — figures (PNG).
- `test_selectors_v3.R` — unit tests T1–T6 (ComBat-vs-sva equivalence, frozen
  apply exactness, M6==M0 ranking, M8 omega behavior, M9 contract, determinism).

## Data provenance

`donor_tables.rds` (2,609 donors x 59 harmonized cell types, 8 cohorts) is the
M1 artifact of WS-8 v2 (`benchmarks/ws8_v2/results/m1/donor_tables.rds`),
lineage-verified (MD5 `50651968855bf6ca3efcdf199e2a9fd5`) to the portal file
`all_cells_metadata.csv.gz` (13,137,775 cells). Design matrix:
`logit((ct+0.5)/(total+0.5))`, `CT:` prefix; per-fold zero-variance feature
filtering fit on training cohorts only. Seed 20261013 throughout.

## Results

Full report: `results/report_ws8_v3_candidate_benchmark.md`; assembled tables
in `results/tables/`; figures in `results/figures/`.

Headline (macro held-out signed-rho replication R@10, 8 LOCO folds, pinned
package): M0/M6 0.395 > M1 0.386 > M5 0.354 > M4 0.327 > M7c 0.234 >
M7 0.219 > M3 0.142 > M2 0.118 > M8 0.109 > M9/M9s 0.060 > M9b 0.056.
No candidate beats the plain RE-meta ranking; M6-M4 = +0.068 (p = 0.047);
M7, M8, M9 all significantly worse than M0 (p <= 0.016). M7 leads full-atlas
P20 clock-panel recovery (0.85, GZMK/GZMB ranked #1/#1) but does not
transport. M7/M8 native panels are not distinguishable from permutation
null. ComBat-merge-then-select is an anti-pattern for this use case
(replication collapse; held-out prediction at the median baseline).
