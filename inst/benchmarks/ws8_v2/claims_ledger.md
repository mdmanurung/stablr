# WS-8 v2 — Claims Ledger

Every narrative claim, its source table/artifact, the exact comparison, status
(`supported` / `tentative` / `unsupported`), and limitations. Updated at each
milestone. Initialized at M0 (2026-10-15); finalized at M7 (2026-10-16) after the
adversarial review and the Amendment 1 re-run.

## M0 audit findings (source-code and historical-artifact claims)

| # | Claim | Evidence | Status | Limitation |
|---|---|---|---|---|
| A1 | v1 `mean_score` is entirely missing | `ws8_results.csv`: `mean_score` NA for 60/60 rows; `colMeans(env_scores)` returns 8 cohort means, name-indexed by 60 features → all NA | supported (artifact) | none |
| A2 | v1 `selected_lasso` is entirely missing | `ws8_results.csv`: `selected_lasso` NA for 60/60 rows; `as.numeric()` strips names before name-indexing | supported (artifact) | v1 vignette "lasso selects 6, only CD8+ T Naive overlaps" is NOT verifiable from the committed CSV |
| A3 | v1 conflates macro and micro MAE | script `mean(loco_all$mae)` = 14.03 (micro, donor-weighted); vignette 12.14 = macro (unweighted cohort mean) | supported (recomputed) | two estimands presented as one |
| A4 | SPAC knockoff covariance is a PD projection | `.estimate_pd_sigma()` uses `corpcor::make.positive.definite`; SPAC n=33 < p=60 rank-deficient | supported (code) | PD repair ≠ knockoff validity |
| A5 | Paper design matches | Terekhova 2026 Immunity: 2,609 donors, 8 cohorts, 59 cell types, ~12M cells; clock = 55 full / 20 robust; association rule = per-cohort `% ~ age + sex + ancestry`, BH, consistent ≥5/8 | supported [1,2] | exact 20/55 lists in Table S1, paywalled, not retrieved |
| A6 | Protocol §3.2 "5/8 frequency vote → 59/59" | inconsistent with `env_scores` (many features < 0.5 in individual cohorts, e.g. KCL 0.265) | unsupported (do not propagate) | the 59/59 rule is unspecified; treated as unverified historical claim |

## v2 outcome claims

Primary endpoint: macro signed-rho replication R(K) over 8 LOCO folds. Artifacts:
`results/run/heldout_replication.csv`, `method_comparisons.csv`,
`per_cohort_replication_K10.csv`. Selectors M0–M5, matched top-K panels.
CIs on paired differences are paired t-intervals over the 8 folds (Amendment 1,
review finding F3); p-values are sign-flip permutation p-values. M0 and M1 are the
Amendment 1 conformant versions: M0 = DerSimonian–Laird random-effects meta-analysis
of per-cohort Fisher-z correlations ranked by |z| magnitude; M1 = per-cohort
`lm(feature ~ age + sex + ancestry)` adjusted vote.

| # | Claim | Evidence | Status | Limitation |
|---|---|---|---|---|
| V1 | Environment-aware aggregation adds value over pooled STABL | M4 R@10 = 0.286 vs M3 = 0.142; paired diff +0.145 [+0.083, +0.206], sign-flip p = 0.0078 (confirmatory); positive in all 8 cohorts | supported | 8 cohorts; imprecise; association replication, not prediction. SPAC fold uses a PD-projected knockoff covariance (A4) |
| V2 | Env-aware STABL does NOT beat the adjusted per-cohort association vote | M4 0.286 vs M1 0.386; paired diff −0.100 [−0.163, −0.036], p = 0.023 (exploratory) | supported (negative for H2 superiority) | M1 is a comparator, not the controlled pair |
| V3 | Env-aware STABL LOSES to random-effects meta-analysis | M4 0.286 vs M0 0.395; paired diff −0.108 [−0.174, −0.043], p = 0.016 (confirmatory) | supported (negative for H2) | pre-amendment M0 (signed-desc ranking) ranked strong negative associations last and scored 0.203; the conformant magnitude ranking is much stronger |
| V4 | Env-aware STABL beats lasso stability | M4 0.286 vs M2 0.118; paired diff +0.168 [+0.067, +0.269], p = 0.0078 (exploratory) | supported | exploratory contrast |
| V5 | Majority aggregation (M5) recovers replication vs strict min (M4) | M5 0.321 vs M4 0.286; paired diff M4−M5 −0.034 [−0.068, −0.0003], p = 0.023 (exploratory) | tentative | strict-min is conservative; majority trades off |
| V6 | A sparse STABL panel retains age prediction in unseen cohorts (H3) | M4 top-10 ridge macro MAE 12.19 y vs all-59 12.35 y and median-age 20.11 y; paired top-10 − all-59 −0.165 [−2.844, +2.515], p = 0.93; top-10 − median −7.92 [−11.13, −4.70], p = 0.0078 | supported | prediction is secondary; parity with all-59 is driven by CIMA (the largest held-out cohort); not a superiority claim over all-59 |
| V7 | A reproducible core exists (H1) | M4 sign agreement 0.70–1.00 across folds at K=10 (mean 0.85); per-cohort R@10 all positive (0.183–0.414) | supported | magnitude modest; see V8 null. (Range corrected at M7 review: an earlier draft stated 0.80–0.95, which appears nowhere in the committed CSV.) |
| V8 | The min-agg selection is NOT non-null at the feature-count level (H5 confirmed) | S0 within-cohort age permutation, 50 permutations: M4 selects 59/59 features in 50/50 (median max min-score 1.000, range 0.980–1.000; Clopper–Pearson 95% lower bound for the saturation probability 0.929). Same-design real-data comparator: 35/59 native features. `results/null/null_results.csv`, `null_scores_perm*.csv`, `real_reference_scores.csv` | supported (adverse) | selection saturation means the panel COUNT carries no signal; the case rests on replication/prediction RANKINGS, not on selection being non-null. Consistent with v1 FDP+ ≈ 1. SPAC PD projection (A4) applies to null and real runs alike |
| V9 | A reproducible core coexists with cohort-dependent components (H4) | All 8 a priori reference populations are sign-consistent in 8/8 cohorts (mean rho −0.63 to +0.29); between-cohort sd of rho 0.10–0.21; the GZMK+/GZMB+ CD8 Tem ratio is NOT sign-consistent (6/8 cohorts negative, range −0.58 to +0.18). `results/biology/` | supported | descriptive, all-donor associations (not held-out); ratio is a secondary predefined readout; no causal claim |
| V10 | The benchmark is seed-level reproducible | Amendment 1 re-run: selectors whose code was unchanged (M2, M3, M4, M5) reproduced their committed per-cohort R, sign-agreement, and prediction-MAE values bit-exactly (max abs diff = 0); only the deliberately changed M0/M1 values moved | supported | same seed (20261013), same data (MD5 50651968855bf6ca3efcdf199e2a9fd5), persisted fold manifest (`results/folds.csv`) |

**Headline (H2): mixed/negative.** Environment-aware STABL min-aggregation improves
held-out association replication over pooled STABL — the controlled aggregation
contrast, +0.145 [+0.083, +0.206] — and over lasso, but LOSES to both simple
association selectors at matched budget: random-effects meta-analysis (−0.108
[−0.174, −0.043]) and the adjusted association vote (−0.100 [−0.163, −0.036]).
Selection count is null-saturated (50/50 permutations select all 59 features vs
35/59 on real data), so rankings — not panel membership — carry the signal. No
claim of method superiority is made.

## References

[1] Terekhova et al. (2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012
    (cell.com fulltext S1074-7613(26)00387-0).
[2] KCL Pure record, https://kclpure.kcl.ac.uk/portal/en/publications/single-cell-analyses-of-global-cohorts-outline-determinants-of-im/
