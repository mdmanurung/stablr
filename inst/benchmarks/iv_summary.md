# IV — Interpretation & validated-workflow benchmark summary

- Predeclared design: interpretation_validation_predesign.md (2026-10-15); 10 seeds.
- Git SHA: 2eb6929; stablr version: 0.3.0.
- Run completeness: 110 planned, 110 completed, 0 failed (iv_failures.csv).
- Simulation-only evidence; no application-data claim. FDP+ remains a
  learner x decoy x design diagnostic.

## Mean metrics per design x arm

            design         arm n_selected  accuracy       ber     brier
1         iv0_null exploratory        3.9        NA        NA        NA
2       iv1_sparse exploratory        4.5        NA        NA        NA
3  iv2_substitutes exploratory        4.3        NA        NA        NA
4       iv3_single exploratory        3.0        NA        NA        NA
5        iv5_leaky       leaky        4.5        NA        NA 0.1656023
6        iv5_leaky      nested        3.3 0.7441667 0.2578162 0.1731275
7         iv0_null   validated        4.1 0.5150000 0.5201398 0.2575107
8       iv1_sparse   validated        3.3 0.7441667 0.2578162 0.1731275
9  iv2_substitutes   validated        4.0 0.7458333 0.2580031 0.1730890
10      iv3_single   validated        3.0 0.7137500 0.2912628 0.1962912
11 iv4_study_shift   validated        2.1 0.7562500 0.2498244 0.1709588
    log_loss     auroc
1         NA        NA
2         NA        NA
3         NA        NA
4         NA        NA
5  0.5004355 0.8298268
6  0.5173541 0.8149938
7  0.7099139 0.4676965
8  0.5173541 0.8149938
9  0.5222720 0.8204933
10 0.5774612 0.7741191
11 0.5103026 0.8283751

## Gate criteria

- G-IV-2 effect-direction recovery >= 0.90 (IV-1): value 1 -> MET
- G-IV-3 substitutes in one group >= 0.80 (IV-2): value 0.6 -> NOT MET
- G-IV-4 leaky AUROC >= nested AUROC (IV-5): value 0.0148 -> MET
- G-IV-5 null median AUROC in [0.35, 0.65] (IV-0): value 0.473 -> MET
- G-IV-6 one-feature probabilities non-degenerate (IV-3): value 4 -> MET

Notes:
- IV-5 contrasts leaky full-data selection with correctly nested
  selection over identical outer folds; the criterion is directional
  (pooled across seeds), not per-replicate.
- Empirical selection truth uses the known simulated signal sets;
  these are distinct from the fitted FDP+ estimate.

