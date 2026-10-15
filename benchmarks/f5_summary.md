# F5 — RF validation simulation summary

- Mode: **full** (smoke = execution proof only, not evidence)
- Design: non-linear signals (x1*x2 interaction, x3^2 quadratic, main effects) with correlated nulls; n=120, p=10.
- Nominal FDP+ level: 0.1
- Knockoff sanity probe (0.2.2 shipped-pair contract on a heterogeneous-scale design, Gram identity + nonnegative cross-Gram shift): **ok**
- Boruta baseline: not installed - skipped gracefully

## Mean empirical FDP / TPR per configuration (at the FDP+-chosen threshold)

| learner | artificial_type | mean empirical FDP | mean TPR | mean n selected | nominal FDP |
|---|---|---|---|---|---|
| lasso | knockoff_equi | 0.578 | 0.867 | 7.6 | 0.1 |
| lasso | random_permutation | 0.398 | 0.600 | 5.4 | 0.1 |
| oblique_rf | knockoff_equi | 0.242 | 0.767 | 3.4 | 0.1 |
| oblique_rf | random_permutation | 0.387 | 0.900 | 5.0 | 0.1 |
| rf | knockoff_equi | 0.065 | 1.000 | 3.3 | 0.1 |
| rf | random_permutation | 0.322 | 1.000 | 4.9 | 0.1 |

## Caveats

- Knockoff build passed the sanity probe (shipped-pair contract on a heterogeneous-scale design: Gram identity Xk'Xk = X'X exact to 1e-6 relative, nonnegative cross-Gram diagonal shift).
- FDP+ is a diagnostic minimizer, not a universal FDR guarantee; for `rf` its calibration is purely empirical (plan D1) and is validated by this suite.
- With `random_permutation` nulls, RF FDP+ is expected to be anti-conservative (permuted nulls do not receive the correlation leakage real nulls receive); the rf+permutation rows should show inflated empirical FDP relative to the knockoff rows on a sane knockoff build.

