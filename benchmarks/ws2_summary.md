# WS-2 — RF adapter graduation summary

- Predeclared design: MATURITY_PREDESIGN.md
 - Knockoff sanity probe: ok 
 - Nominal FDP+ level: 0.1 

## Mean empirical FDP / TPR per cell

                design    artificial_type n_selected       fdp       tpr
1             A_linear      knockoff_equi        2.8 0.0250000 0.9000000
2          B_nonlinear      knockoff_equi        3.3 0.0650000 1.0000000
3  C_highdim_nonlinear      knockoff_equi       17.5 0.4376899 1.0000000
4 D_nonlinear_binomial      knockoff_equi        5.0 0.2646429 1.0000000
5             A_linear random_permutation        6.1 0.4433333 1.0000000
6          B_nonlinear random_permutation        4.9 0.3221429 1.0000000
7  C_highdim_nonlinear random_permutation        6.4 0.3630201 1.0000000
8 D_nonlinear_binomial random_permutation        5.6 0.3340368 0.9666667

## Predeclared graduation criterion (rf + knockoff_equi: mean FDP <= 0.15 AND mean TPR >= 0.8 on ALL designs)

Outcome: NOT MET — adapter keeps a design-scope restriction 

## Caveats

- FDP+ is a diagnostic minimizer, not a universal FDR guarantee; rf calibration is purely empirical (plan D1).
 - Permutation rows are expected anti-conservative (permuted nulls lack correlation leakage); they characterize the guard, not graduation.
 - Simulation-only evidence; no application-data claim.
