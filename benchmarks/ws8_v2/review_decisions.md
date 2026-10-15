# WS-8 v2 — M7 Adversarial Review: Decisions

Independent adversarial review of the v2 benchmark package (report held with
this ledger). Verdict of the review: no leakage path; primary endpoint matches
the locked protocol exactly; all headline numbers and sign-flip p-values
independently reproduced; but 7 major findings required remediation. This file
records the disposition of every finding. Remediation code landed in commit
`b85decb` ("Amendment 1"); results were regenerated in the amended re-run
(`results/run_m012b`, `run_m3b`, `run_m4b`, `run_m5b`).

## Major findings

1. **V7 sign-agreement range wrong (claimed 0.80–0.95; actual 0.70–1.00).**
   FIXED. Ledger V7 corrected and annotated with the erratum.
2. **V8's real-data comparator "24 features" was a hard-coded literal; null
   (8 cohorts) vs real (7-cohort LOCO) designs differed.** FIXED.
   `real_reference.R` runs the full-atlas real M4 fit through the identical
   code path as the null (same transform, same 8-cohort design) and commits
   `results/null/real_reference_scores.csv` (per-feature scores + native
   flags) and `real_reference_env_scores.csv`. `null_probe.R` and `figures.R`
   now read the comparator from that artifact. The LOCO native panels are
   additionally auditable per fold in `training_scores.csv`.
3. **"CIs" were spread quantiles of the 8 paired differences, not CIs for the
   mean.** FIXED. `assemble_primary.R` now reports paired t intervals for the
   mean cohort difference (7 df; approximate-symmetry assumption and n = 8
   imprecision stated in the header) alongside the exact sign-flip p; the
   2.5–97.5% quantiles of the differences are kept under the honest name
   `spread_low/spread_high`. Ledger intervals updated.
4. **M1 deviated from locked §8 (unadjusted Spearman vote instead of
   `% ~ age + sex + ancestry` LM).** FIXED by conformance: M1 now implements
   the locked adjusted-LM vote (per-cohort LM, BH, ≥ 5/7; constant covariate
   terms dropped per cohort — ancestry is constant in 4 of 8 cohorts; sex and
   ethnicity are complete). Filed as Amendment 1; M1 re-run. Outcome: M1
   macro R@10 = 0.386 (was 0.381 under the deviating variant); the headline
   negative result (M4 < M1) is unchanged.
5. **M0 deviated from locked §4 (equal-weight mean instead of random-effects
   meta).** FIXED by conformance: DerSimonian–Laird random-effects meta on
   Fisher-z rhos. The top-K ranking is by |estimate| (direction-agnostic,
   consistent with the signed endpoint and the other selectors); the previous
   signed-descending ranking ranked strong negative associations last, which
   weakened M0 as a comparator. Filed as Amendment 1; M0 re-run. Outcome: M0
   macro R@10 = 0.395 (was 0.203); the confirmatory contrast M4 vs M0 turns
   negative — see the ledger headline.
6. **S0 null stopped at 10 permutations vs locked 50–100 target; post-hoc
   "point mass" rationalization.** FIXED. Extended to 50 permutations
   (chunked, independent seeds; `WS8_PERM_START`/`WS8_N_PERM`), reported with
   Monte Carlo uncertainty. The incorrect "point mass / MC error irrelevant"
   sentence was removed from V8; the correct statement is that the empirical
   null of the count statistic sat at the ceiling in every permutation run.
7. **§6/§10 conformance gaps (no per-fold manifest, no folds table, training
   scores not persisted).** FIXED. `folds_manifest.R` writes
   `results/folds.csv`; `run.R` writes `fold_manifest.csv` (train/test ID
   hashes, drop sets, panel sizes, seed), `training_scores.csv` (all 59
   features, native flags), and `training_env_detail.csv` (per-cohort audit
   values: selection frequencies for M4/M5, Spearman rhos for M0, adjusted
   LM p-values for M1) per run.

## Minor findings

8. **§5 predictive inference incomplete; V6 framing.** FIXED.
   `assemble_prediction.R` computes micro MAE companions, paired per-cohort
   MAE contrasts (top-K vs all-59, top-K vs median, native vs all-59) with
   t intervals and sign-flip p, and calibration ranges. V6's limitation now
   states that the macro edge is driven by CIMA and that the top-10 panel is
   worse than all-59 in 5 of 8 cohorts.
9. **Null did not reproduce the real pipeline (no fold transform).** FIXED.
   Null and real reference both apply `fit_fold_transform`; the transform
   depends only on X, which permutation does not change, so the drop set is
   identical across permutations by construction (and empty on this atlas).
10. **SPAC p > n knockoff PD projection unacknowledged for v2.** ACKNOWLEDGED
    in the ledger (V1/V8 limitations) and the vignette limits: every
    per-environment fit containing SPAC uses a PD-projected covariance, and
    min-aggregation gives the smallest, most-projected cohort veto power.
    Per-cohort knockoff diagnostics remain in the deferred S2 scope.
11. **M0 native panel tautology (`abs(sc) >= 0` selects all).** FIXED. M0 has
    no native panel (`character(0)`); the empty-panel median-age fallback of
    §5 applies.
12. **M3 rows provenance (uncommitted log; K-sensitivity absent).** FIXED.
    The M3 regeneration was committed (`run_m3`, all K, real sign agreement)
    and the combined table re-assembled from executed artifacts; the amended
    re-run (`run_m3b`) adds score persistence. Smoke directories removed.

## Nits

13. **`loco.R` comment contradicted code ("donor-weighted").** FIXED
    ("cohort-balanced (unweighted)").
14. **Dead code in M1.** FIXED (removed in the M1 rewrite).
15. **`data_manifest` field vs script.** VERIFIED NOT A BUG: exactly one bad
    annotation level exists, whose literal name is "Bad quality Bad quality";
    `paste(bad, collapse = ";")` on a length-1 vector reproduces the manifest
    value byte-for-byte.
16. **cwd assumption undocumented.** FIXED: every script header states "run
    from the repository root"; README documents the run order and env vars.
17. **NaN in console macro table for M3 at K ≠ 10.** MOOT: `run_m3`/`run_m3b`
    provide all K values.

## Reproducibility verification (Amendment 1, item 7)

The amended re-run re-executed all six selectors. Implementations that did not
change (M2, M3, M4, M5) were required to reproduce their committed per-cohort
R@10 values exactly; M0/M1 changed by design. Outcome: M2 reproduced
bit-exactly (max |diff| = 0 across folds); M3/M4/M5 verification recorded in
`claims_ledger.md` (claim V10) once their re-runs landed.
