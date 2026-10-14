# WS-9 — Ordinal learner graduation summary (corrected)

- Predeclared design: MATURITY_PREDESIGN.md (WS-9); 10 seeds x 2 designs x 2 decoy types x 100 bootstraps.
- Nominal FDP+ level: 0.1; graduation: knockoff_equi mean FDP <= 0.15 AND TPR >= 0.8 on BOTH designs.
- Truth set: nonzero coefficients of beta (3 of 10 predictors). The
  historical script used all 10 indices; see the dated erratum in
  MATURITY_PREDESIGN.md.
- Run completeness: 40 planned, 38 completed, 2 failed (ws9_failures.csv).

## Mean empirical FDP / TPR per cell

           design    artificial_type n_selected       fdp tpr
1 O1_linear_indep      knockoff_equi   3.700000 0.1400000   1
2   O2_equicorr03      knockoff_equi   3.666667 0.1555556   1
3 O1_linear_indep random_permutation   4.000000 0.2150000   1
4   O2_equicorr03 random_permutation   4.000000 0.2000000   1

Outcome: NOT MET — ordinal stays experimental with documented scope

- Empirical FDP/TPR are over the corrected 3-predictor truth set; they are
  distinct from the fitted FDP+ estimate.
- Permutation rows are expected anti-conservative (characterization only).
- Simulation-only evidence; no application-data claim.

