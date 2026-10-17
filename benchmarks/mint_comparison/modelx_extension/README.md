# Model-X knockoff extension: does the knockoff constructor repair the KE miscalibration?

Follow-up to the MINT comparison (parent directory). Question: the
multi-cohort benchmark showed `knockoff_equi` (KE) FDP+ estimates are
miscalibrated under cohort heterogeneity (lasso anti-conservative,
sPLS-DA conservative). Is that an artifact of the *equi-correlated*
covariance shortcut, or does it persist for a higher-quality model-X
constructor (`knockoff_mvr`)?

## Design

- Same DGP as the parent benchmark: K = 3 cohorts x 80, p = 150, 10 shared
  signals (beta = 1.2/sqrt(10)), batch location shifts on 20 nuisance
  features, prevalence offsets.
- Full factorial, native `stabl_fit()` for every arm:
  `base_learner` in {lasso, splsda} x `artificial_type` in
  {random_permutation, knockoff_equi, knockoff_mvr}, 20 seeds,
  `n_bootstraps = 100`, `lambda_grid = "auto"`, `random_state = seed`.
- Fixed-X knockoffs (`artificial_type = "knockoff"`) are not constructible
  on this DGP (n = 240 < 2p = 300; `knockoff::create.fixed` augments rows)
  and the strict `knockoff_fallback = "stop"` correctly refuses.
- Calibration = realized FDP / reported `min_fdr_` per seed; medians below
  (means are inflated by seeds with `min_fdr_` near 0).

## Results (20 seeds, median per-seed calibration)

| arm        | TPR  | realized FDP | n_sel | reported FDP+ | calibration | Jaccard | secs |
|------------|------|--------------|-------|---------------|-------------|---------|------|
| lasso+RP   | 0.47 | 0.33         | 7.6   | 0.21          | **1.00**    | 0.172   | 0.9  |
| lasso+equi | 0.49 | 0.46         | 15.5  | 0.25          | 1.40        | 0.124   | 1.4  |
| lasso+MVR  | 0.49 | 0.43         | 13.1  | 0.26          | 1.33        | 0.127   | 3.7  |
| sPLS-DA+RP | 0.46 | 0.33         | 7.6   | 0.21          | **1.00**    | 0.178   | 29.4 |
| sPLS-DA+equi | 0.45 | 0.27       | 9.8   | 0.50          | 0.42        | 0.139   | 29.3 |
| sPLS-DA+MVR  | 0.45 | 0.25       | 9.4   | 0.50          | 0.38        | 0.148   | 31.1 |

## Reconciling with the parent benchmark's lasso+KE severity

The parent MINT comparison reported lasso + knockoff_equi calibration of
~8.0; this extension reports 1.40 for the same nominal arm. Both are
correct within their own protocols, which differ in two ways: the parent
swept a 20-value `auto_lambda_grid` and aggregated as
mean(realized FDP) / mean(reported FDP+) over seeds, while this extension
uses the default 30-value auto grid and the median of per-seed ratios.
Coarser grids give the argmin FDP+ rule fewer exit thresholds, which
amplifies the winner's-curse anti-conservatism (realized FDP 0.83 in the
parent vs 0.46 here), and the ratio-of-means is dominated by seeds with
near-zero reported FDP+. The qualitative conclusion is identical under
both protocols — equi and MVR track each other, random permutation is
calibrated — but cross-benchmark severity comparisons must be made within
one protocol; all six arms in this extension share a single protocol.

## Verdict

**No.** Switching the constructor from equi-correlated to model-X MVR does
not repair the miscalibration: it tracks equi closely in both direction and
magnitude (lasso 1.33 vs 1.40 anti-conservative; sPLS-DA 0.38 vs 0.42
conservative), while random permutation is exactly calibrated (1.00) for
both learners. MVR buys a marginal power gain at matched TPR (FDP 0.43 vs
0.46 for lasso) at ~2.7x the compute, with no calibration benefit.

Interpretation: the failure is not the covariance-solve shortcut — it is
the model-X construction itself. Both constructors fit a single multivariate
Gaussian to the *pooled* batch-heterogeneous matrix, so knockoffs are not
exchangeable with real features under cohort mixture (consistent with the
Phase 6 mechanism probe: real null features systematically outscore their
knockoffs). Under known batch structure, use random permutation for these
learners; model-X knockoffs require a per-batch (or covariate-adjusted)
construction, which stablr does not implement.

Reproducibility: `Rscript knockoff_mx_test.R` (checkpoints one RDS per
arm-seed; `N_SEEDS` env var; `PILOT=1` for a smoke test). Raw per-seed
records in `mx_per_seed.csv`, aggregated in `mx_summary.csv`.
