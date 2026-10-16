# WS-9 — Ordinal learner graduation summary

- Predeclared design: MATURITY_PREDESIGN.md (WS-9); 10 seeds x 2 designs x 2 decoy types x 100 bootstraps.
- Nominal FDP+ level: 0.1; graduation: knockoff_equi mean FDP <= 0.15 AND TPR >= 0.8 on BOTH designs.

## Mean empirical FDP / TPR per cell

           design    artificial_type n_selected fdp       tpr
1 O1_linear_indep      knockoff_equi   3.700000   0 0.3700000
2   O2_equicorr03      knockoff_equi   3.666667   0 0.3666667
3 O1_linear_indep random_permutation   4.000000   0 0.4000000
4   O2_equicorr03 random_permutation   4.000000   0 0.4000000

Outcome: NOT MET — ordinal stays experimental with documented scope

- Permutation rows are expected anti-conservative (characterization only).
- Simulation-only evidence; no application-data claim.

