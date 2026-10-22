# WS-8 v3 — Benchmark of candidate mechanisms for cross-cohort reproducible biomarker discovery

**Date:** 2026-10-21 · **Package:** stablr @ `a8f3a58` (pinned) · **Data:** Terekhova et al. (2026, *Immunity*) multi-cohort immune-aging atlas [4] · **Seed:** 20261013 throughout

*Application case study — NOT calibration evidence.*

---

## 1. Summary

Three candidate mechanisms for transportable biomarker discovery were benchmarked against the locked WS-8 v2 matrix (M0–M5) and a mandatory ComBat-merge baseline, using 8-fold leave-one-cohort-out (LOCO) evaluation on 2,609 donors × 59 harmonized immune cell-type features [4]:

- **Candidate 1 — M6**, heterogeneity-filtered random-effects meta-analysis (MetaIntegrator-style [63, 64] with an explicit Cochran's-Q homogeneity gate).
- **Candidate 2 — M7**, MINT (`mixOmics::mint.spls` [13]) inside the STABL bootstrap + knockoff FDP+ machinery; **M7c** = tuned MINT without stability (control).
- **Candidate 3 — M8**, contrast-penalized lasso (common + cohort-deviation terms, ω = 2; Rashid/GMCP family [15, 56]) inside STABL.
- **Baseline — M9**, frozen-hyperparameter ComBat harmonization + pooled STABL [17]; **M9s** (strict transform) and **M9b** (no age protection; post-hoc, added after the M9 prediction failure was diagnosed) as sensitivities.

**Headline results (primary endpoint: macro held-out signed-ρ replication R@10):**

1. **No candidate beats the simplest comparator.** The plain RE-meta ranking (M0, 0.395) remains the strongest arm; M6 ties it exactly (0.395, by construction — the homogeneity filter changes the native panel, not the ranking). The v2 primary arm M4 (env-STABL) reruns at 0.327 on the pinned package and is significantly **worse** than M0 (M4−M0 = −0.068, 95% CI [−0.122, −0.015], sign-flip p = 0.047) — the v2 confirmatory finding replicates, attenuated.
2. **M7 (MINT-STABL) = 0.219** — significantly worse than both M0 (−0.176, p = 0.008) and M4 (−0.108, p = 0.016). Stability selection adds nothing over tuned MINT (M7−M7c = −0.015, p = 0.52).
3. **M8 (contrast-STABL) = 0.109** — significantly worse than M0 (p = 0.008) and M4 (p = 0.008), and its native-panel machinery is degenerate under permutation null (§7).
4. **M9 (ComBat+STABL) = 0.060** — the worst arm; harmonization actively destroys cross-cohort replication (M9−M3 = −0.082 vs pooled STABL on raw data) and its held-out prediction collapses to near the training-median baseline (§8).
5. **Co-primary (recovery of the paper's 20-population clock panel):** M0/M6 again lead on LOCO top-10 recovery (0.42); M7 leads on full-atlas top-20 recovery (0.85) and uniquely ranks both GZMK+/GZMB+ Tem components #1 — the paper's headline ratio [4] — but this does not translate into held-out replication.
6. **Heterogeneity is the dominant feature of this data:** 33/59 cell types are significantly age-associated (BH < 0.05) but 29/33 are significantly heterogeneous across cohorts (Q p < 0.10; median I² = 0.81; 37/59 sign-discordant). Only 4 features pass M6's homogeneity gate (§6).

**Recommendation for stablr:** ship Candidate 1 (per-cohort association + RE-meta ranking with an explicit heterogeneity report) as the transportability workflow; do not pursue MINT-STABL or contrast-STABL as implemented; document the ComBat-merge-then-select strategy as an anti-pattern for this use case.

---

## 2. Design

**Data.** 2,609 donors × 59 harmonized cell-type proportions from 8 PBMC scRNA-seq cohorts (ABF300 n=166, AIDA 625, CIMA 421, KCL 188, OneK1K 981, SoundLife 96, SPAC 33, UCSF 99) [4]. Design matrix `logit((ct+0.5)/(total+0.5))`, `CT:` prefix; per-fold zero-variance feature filtering fit on training cohorts only (leakage-proof engine reused verbatim from WS-8 v2). Outcome: chronological age.

**Evaluation.** 8 LOCO folds. Primary: signed-ρ replication R(K) — mean over the top-K panel of train-sign × held-out within-cohort Spearman ρ — macro-averaged over folds; K=10 primary, {5,15,20} sensitivity. Co-primary: recovery of the paper's 20-population reduced clock panel (authors' repository notebook `Age_clock_model/clean_notebook_for_aging_models.ipynb`, cell 32 `selected_by_marina`; all 20 populations verified to map onto the 59 harmonized features). This is a **concordance target, not ground truth**. Secondary: cohort-respecting ridge prediction (macro MAE), reference-population sign concordance (8 curated populations, sign-consistent across all 8 cohorts), GZMK/GZMB components, and within-cohort age-permutation null probes for the new STABL-machinery arms (10 permutations, full n_boot=200).

**Statistics.** Contrasts are paired across the 8 folds (the generalization units): mean difference with a t interval (7 df) and an exact two-sided sign-flip permutation p — identical methodology to v2. Predeclared confirmatory contrasts: each candidate vs M0 and vs M4; plus the v2 confirmatory pair (M4−M0, M4−M3) rerun on the pinned package. All other contrasts are exploratory.

**Package version policy.** All v3 numbers use stablr @ `a8f3a58`. The historical v2 runs used earlier commits whose fitting code differs materially; M0–M3 reproduce bit-exactly under the pin (M3's pooled-STABL path is unchanged), while the env-stability wrapper shifted (M4 0.288→0.327, M5 0.321→0.354). The v3 matrix therefore reruns M3/M4/M5 like-for-like; v2 numbers appear only as historical reference (`tables/version_drift_reference.csv`).

---

## 3. Primary endpoint: held-out signed-ρ replication

Macro R(K) across 8 LOCO folds (full per-fold table: `tables/per_cohort_replication_K10_v3.csv`; Fig 1–2):

| Arm | K=5 | **K=10** | K=15 | K=20 |
|---|---|---|---|---|
| **M0** RE-meta ranking | 0.480 | **0.395** | 0.347 | 0.285 |
| **M6** het-filtered RE-meta | 0.480 | **0.395** | 0.347 | 0.285 |
| M1 cohort vote | 0.442 | 0.386 | 0.336 | 0.285 |
| M5 env-STABL (vote) | 0.446 | 0.354 | 0.287 | 0.241 |
| M4 env-STABL (min) | 0.455 | 0.327 | 0.254 | 0.230 |
| M7c MINT tuned | 0.361 | 0.234 | 0.238 | 0.217 |
| M7 MINT-STABL | 0.219 | 0.219 | 0.275 | 0.273 |
| M3 pooled STABL | 0.216 | 0.142 | 0.212 | 0.214 |
| M2 pooled lasso | 0.096 | 0.118 | 0.170 | 0.168 |
| M8 contrast-STABL | 0.028 | 0.109 | 0.114 | 0.130 |
| M9 ComBat+STABL | 0.019 | 0.060 | 0.100 | 0.107 |
| M9s ComBat std-only | 0.019 | 0.060 | 0.100 | 0.107 |
| M9b ComBat(no-age)+STABL | 0.020 | 0.056 | 0.103 | 0.111 |

Confirmatory contrasts at K=10 (paired, 8 folds; `tables/method_comparisons_v3.csv`):

| Contrast | Δ macro R@10 | 95% CI | sign-flip p |
|---|---|---|---|
| M6 − M0 | 0.000 | [0, 0] | 1.0 |
| **M6 − M4** | **+0.068** | [+0.015, +0.122] | 0.047 |
| M7 − M0 | −0.176 | [−0.238, −0.114] | 0.008 |
| M7 − M4 | −0.108 | [−0.152, −0.065] | 0.016 |
| M8 − M0 | −0.286 | [−0.364, −0.208] | 0.008 |
| M8 − M4 | −0.218 | [−0.293, −0.143] | 0.008 |
| M9 − M0 | −0.335 | [−0.424, −0.245] | 0.008 |
| M9 − M4 | −0.266 | [−0.352, −0.181] | 0.008 |
| M4 − M0 (v2 confirmatory, rerun) | −0.068 | [−0.122, −0.015] | 0.047 |
| M4 − M3 (v2 confirmatory, rerun) | +0.185 | [+0.095, +0.274] | 0.008 |

Exploratory: M7−M7c = −0.015 (p = 0.52); M9−M9s = 0 (p = 1.0 — identical selection, as constructed); M9−M3 = −0.082 (p = 0.13); M9b−M9 = −0.005 (p = 1.0); M7c−M0 = −0.161 (p = 0.008).

**Reading.** The ranking-only arms (M0/M6/M1) dominate. Stability-based arms span the middle (M5/M4) to the bottom (M3), and every mechanism that *jointly* models or *merges* cohorts (M7, M8, M9) degrades replication — consistent with the extreme effect-size heterogeneity in this data (§6) and with the documented failure of invariant-learning approaches on tabular biomedical data [67, 76].

---

## 4. Co-primary: recovery of the paper's 20-population clock panel

LOCO top-10 panels, macro over folds (`tables/panel_recovery_macro_v3.csv`; Fig 3):

| Arm | recovery@10 | precision@10 |
|---|---|---|
| **M0 / M6** | **0.419** | 0.838 |
| M1 | 0.400 | 0.800 |
| M7 | 0.331 | 0.663 |
| M5 | 0.313 | 0.625 |
| M7c | 0.288 | 0.575 |
| M4 | 0.275 | 0.550 |
| M3 | 0.194 | 0.388 |
| M2 | 0.144 | 0.288 |
| M8 | 0.088 | 0.175 |
| M9 / M9s | 0.013 | 0.025 |
| M9b | 0.000 | 0.000 |

Full-atlas fits (all 8 cohorts; `tables/full_atlas_summary_v3.csv`): top-20 P20 recovery — **M7 0.85**, M0/M6 0.80, M8 0.75, M1 0.70, M3 0.60, M7c 0.55, M4/M5 0.50, M2 0.35, M9/M9b 0.20. Reference-population recovery in top-20: M0/M6/M7 7-of-8 (all sign-concordant); M8/M3 6-of-8; M1/M7c/M4/M5/M2 5-of-8; M9/M9b 2-of-8.

**GZMK⁺/GZMB⁺ Tem components** (the paper's headline ratio [4]): M7 ranks both features **#1** on the full atlas; M8 ranks GZMK⁺ #1 but GZMB⁺ #49; M0/M6 rank them 17/16; M4 ranks them 36/28. MINT's joint integrative objective does home in on the paper's signature biology — but its LOCO replication (0.219) shows the panel does not transport as an absolute-effect ranking. The ratio components are themselves strongly age-associated but heterogeneous (GZMK⁺: ρ_RE = +0.20, I² = 0.91; GZMB⁺: ρ_RE = +0.22, I² = 0.80), so no method that requires effect-size consistency will elevate them.

**Interpretation.** M0/M6 dominate the transportable co-primary (LOCO recovery); M7's full-atlas concordance is real but non-transportable — exactly the distinction the benchmark was designed to expose.

---

## 5. Secondary: held-out prediction

Macro MAE in years, top-10 panel, cohort-respecting ridge (`tables/prediction_summary_v3.csv`; Fig 5). Training-median baseline: 20.11 y.

| Arm | MAE | rank ρ |
|---|---|---|
| M4 | **11.87** | 0.666 |
| M3 | 12.87 | 0.613 |
| M7c | 13.02 | 0.669 |
| M5 | 13.11 | 0.659 |
| M1 | 13.19 | 0.674 |
| M7 | 14.05 | 0.622 |
| M2 | 14.29 | 0.636 |
| M0 / M6 | 14.32 | 0.665 |
| M8 | 14.65 | 0.585 |
| M9s | 17.61 | 0.212 |
| M9b | 18.50 | 0.175 |
| M9 | 18.92 | 0.201 |

Env-STABL (M4) is the best *predictor* despite being mid-table on replication — the stability–performance tradeoff [31] in action. All ComBat variants collapse to rank ρ ≈ 0.2 and MAE within ~2 y of the median baseline.

---

## 6. Heterogeneity map (why M6's native panel has only 4 features)

Full-atlas per-feature statistics (`tables/heterogeneity_m6_fullatlas.csv`; `figures/fig_m6_heterogeneity_map.png`):

- 33/59 features significant at BH < 0.05; **29/33 fail the homogeneity gate** (Cochran's Q p < 0.10). Median I² = 0.81; 37/59 features are sign-discordant across cohorts.
- The strongest associations are the least transportable: CD8⁺ T Naive (ρ_RE = −0.63, I² = 0.90, per-cohort range −0.83 to −0.42), CD8⁺ RTE (−0.60, I² = 0.96, −0.82 to −0.15), CD4⁺ RTE (−0.57, I² = 0.92). Cochran's Q tests *magnitude* homogeneity, not sign — the 8 curated reference populations are sign-consistent everywhere but magnitude-heterogeneous, so they are filtered.
- The 4 native features: **MAIT** (ρ_RE = −0.30, I² = 0.09 — also a reference population), **B Plasma cells** (−0.09, I² = 0), **CD4⁺ Treg cytotoxic** (+0.08, I² = 0), **CD4⁺ Proliferative** (+0.07, I² = 0): small but effect-size-transportable. None are in the paper's P20 panel.
- Caveat: with 8 cohorts, Q has limited power; Q p ≥ 0.10 is absence of *detectable* heterogeneity, not proof of homogeneity.

---

## 7. Null probes: M7/M8 native panels are not interpretable

10 within-cohort age permutations, full-atlas, locked n_boot = 200 (`tables/null_summary_v3.csv`; Fig 4):

- **M7:** real native panel 40 vs null median 10.5 (range 0–38; empirical p ≈ 0.09). More damning, the real maximum stability score (0.125) sits *inside* the null score distribution (null max 0.125): MINT loadings are so unstable across bootstrap resamples that no feature exceeds ~12% selection frequency, and the FDP+ argmin threshold (0.08–0.13) labels noise-level selections as native.
- **M8:** degenerate — null median native panel 59/59 (range 54–59) vs real 48. Max-over-λ stability scores saturate at 1.0 for nearly every feature on the 944-column expanded design, so the FDP+ curve is flat and the threshold admits everything, on real and permuted data alike.

These nulls calibrate the *native-panel* machinery only; they do not bear on the top-K rankings behind the primary endpoint. Framing: as implemented, these mechanisms do not calibrate — an implementation-level negative result, not a verdict on the underlying ideas [13, 15, 56].

---

## 8. M9 prediction failure mechanism (and the M9b verdict)

M9 fits ComBat with `mod = ~age` on the 7 training cohorts; the adjusted training data includes the covariate add-back `age·β` (a large offset at raw age scale). The label-free held-out transform cannot include that term, and its per-batch location estimate additionally absorbs the held-out cohort's mean age shift (label-free ComBat cannot distinguish technical offset from "this cohort is older"). Quantified on the KCL fold: train MAE 11.4 y vs held-out 33.8 y; mean prediction 27.5 y vs true mean 61.2 y — the −22.5 y collapse ≈ the entire mean-age add-back, amplified by ridge coefficients.

**M9b** (`mod = NULL`, consistent train/test transform) fixes the fold-specific catastrophe (KCL 33.8 → 20.5 y) but leaves macro MAE unchanged (18.5 vs 18.9) with rank ρ still ≈ 0.2: the dominant failure is the *irreducible* label-free removal of between-cohort age signal, not the covariate add-back mismatch. M9s (no new-batch moments at all) sits between (17.6). Conclusion: correct-then-predict cannot transport absolute age calibration across cohorts with different age distributions, under any tested transform variant. The primary replication endpoint is immune to this by design (unadjusted held-out matrix; within-cohort Spearman is invariant to monotone per-feature transforms).

---

## 9. Limitations

- Single dataset, single outcome (age), 59 compositional features; 8 folds give coarse contrast resolution (minimum achievable sign-flip p = 0.0078) and the t CIs are approximate at n = 8.
- The P20 panel is a concordance target recovered from the authors' repository notebook, not ground truth; the paywalled Table S1 may differ.
- M8's degeneracy is partly an interaction between the contrast-expanded design and max-over-λ scoring; a redesigned path (e.g., fixed-λ stability) might behave differently — not pursued.
- Null probes used 10 permutations (empirical p granularity ≈ 0.09); sufficient for the qualitative conclusions drawn.
- FDP+ is a diagnostic, not FDR control; all arms are evaluated on held-out replication, not on their internal error estimates.
- Runtimes (`tables/runtime_v3.csv`) are wall-clock under parallel contention and should be read as approximate.

---

## 10. Reproducibility

- Code module: `benchmarks/ws8_v3/` on branch `ws8-v3-candidates` (selectors, driver, null probe, full-atlas, assembly, figures, unit tests T1–T6 incl. ComBat-vs-`sva` numerical equivalence < 1.5e-6 and frozen-apply exactness < 1e-8).
- Data: `benchmarks/ws8_v2/results/m1/donor_tables.rds` (lineage MD5 `50651968855bf6ca3efcdf199e2a9fd5` → portal `all_cells_metadata.csv.gz`, 13,137,775 cells).
- Environment: R 4.4.2; stablr @ `a8f3a58`; mixOmics 6.30.0; sva 3.54.0; glmnet 5.0; ranger 0.18.0; knockoff 0.3.6. Seed 20261013; n_boot = 200; per-(method, fold) checkpoints.
- Session deliverables: `ws8_v3/tables/` (12 CSVs) and `ws8_v3/figures/` (6 PNGs).

## References

Inline citations refer to the session's literature search: MINT [13]; penalized mixed-model / contrast family [15, 56]; MetaIntegrator [63, 64]; ComBat-merge baseline evidence [17]; IRM-family failures on tabular data [67, 76]; stability–performance tradeoff [31]; immune-aging reproducibility context [42, 43]; blood aging-signature positive control [41]; the benchmarked paper [4].
