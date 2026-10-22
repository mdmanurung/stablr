# M10 addendum — Invariant Causal Prediction (ICP) arm (2026-10-22)

Post-hoc addendum to the WS-8 v3 candidate benchmark. ICP (Peters, Bühlmann
& Meinshausen 2016) inspired the package's environment-stability workflow
(M4/M5) but was not among the benchmarked arms; this addendum provides
direct evidence. Same locked LOCO harness, same primary endpoint (held-out
signed-ρ replication R@10), same seed policy. Application case study — NOT
calibration evidence.

## Design

- `InvariantCausalPrediction::ICP()` (v0.8), the 7 training cohorts as
  environments (`ExpInd`), α = 0.05, package-default variable selection
  ("boosting" for p > 8), `stopIfEmpty = FALSE`.
- Ranking: ICP p-value descending, tie-broken by |maximin coefficient| desc,
  then feature name. `native` = accepted set (p ≥ α).
- 8-fold LOCO, identical fold construction, feature filtering, and endpoint
  code as the v3 run (`run_v3.R` engine).

## Result: ICP is degenerate on this atlas

| | M10 (ICP) | M0 (RE meta) |
|---|---|---|
| macro R@10 | **0.055** | 0.395 |
| accepted set size | 59/59 in every fold | — |

M0 − M10 = +0.340 [+0.249, +0.430], exact sign-flip p = 0.0078 (all eight
fold differences positive). M10 is the worst arm in the extended matrix,
marginally below the ComBat controls (M9b 0.056).

## Diagnosis: the ICP premise fails on this atlas

The degeneracy is not a tuning artifact — four configurations (default
boosting; lasso; lasso with α = 0.01; stability selection; gof relaxed to
0.01; maxNoObs raised to 1,000) all produce `modelReject = 1.00` and
p-values identically 1.00 for all 59 features. The mechanism:

1. ICP requires at least one candidate set S whose residuals are invariant
   across environments. The intercept-only model fails (cohort age
   distributions differ), and **no tested set (≤ 8 features) restores
   invariance** — every candidate model is rejected, so ICP has no valid
   model and returns the trivial accept-everything output
   (`stopIfEmpty = FALSE`).
2. Direct diagnostic: even the **full 59-feature OLS** has strongly
   cohort-dependent residuals (Kruskal–Wallis p = 2e-13; Levene p = 1.4e-17;
   per-cohort mean residuals −4.75 to +1.28 years). There is no invariant
   linear model in this feature space for ICP to find.
3. Larger candidate sets (maxNoVariables = 20) exceed 40 minutes per fold —
   computationally infeasible for the 8-fold harness — and would still need
   an invariant model to exist.

## Interpretation

ICP's premise — some subset of features yields an outcome mechanism that is
invariant across environments — is empirically false on this atlas, for the
same reason the other invariance-flavored arms failed: cohort identity
explains the largest variance share, effect heterogeneity is pervasive
(median I² = 0.81), and even the full-feature model's residuals are
cohort-dependent (consistent with the calibration slopes < 1 in every
cohort of the case study). The degenerate acceptance (59/59) mirrors the
M4/M8 null-saturation finding and the v2 S1 synthetic result that
saturation is structural in this feature/outcome regime.

Practical reading for `stablr` users: ICP-style invariance testing needs
environments that differ in mechanism while sharing an invariant outcome
model. Multi-cohort harmonized-composition atlases with cohort-shifted
outcome distributions do not satisfy this premise; the meta-analytic
workflow (per-cohort association + RE ranking + heterogeneity audit) is the
appropriate tool there.

## Files

- `heldout_replication_m10.csv`, `per_feature_m10.csv`,
  `fold_manifest_m10.csv`, `extras_m10.csv` — harness outputs (8 folds).
- Selector: `../m10_icp_v3.R` (reuses the locked v2 engine verbatim).
