# WS-5 — Environment-stability transportability summary

- Predeclared: MATURITY_PREDESIGN.md; E = 3, n_env = 80, p = 30 (3 signal); binomial; random_permutation; 10 seeds; 100 bootstraps.

## Mean FDP / TPR by mechanism, method, environment

           mechanism      method environment n_selected        fdp       tpr
       M1_null_shift  env_minagg         all        2.1 0.14166667 0.5666667
     M2_signal_shift  env_minagg         all        2.5 0.09166667 0.7333333
 M3_prevalence_drift  env_minagg         all        1.8 0.08333333 0.5333333
       M1_null_shift      pooled         all        3.5 0.16250000 0.9000000
     M2_signal_shift      pooled         all        3.1 0.07500000 0.9333333
 M3_prevalence_drift      pooled         all        3.4 0.13000000 0.9333333
       M1_null_shift env_per_env        env1        3.5 0.54166667 0.5000000
     M2_signal_shift env_per_env        env1        5.5 0.28035714 0.6666667
 M3_prevalence_drift env_per_env        env1        6.6 0.53500000 0.5666667
       M1_null_shift env_per_env        env2        5.0 0.33880952 0.8333333
     M2_signal_shift env_per_env        env2        5.0 0.33880952 0.8333333
 M3_prevalence_drift env_per_env        env2        5.0 0.33880952 0.8333333
       M1_null_shift env_per_env        env3        3.7 0.44000000 0.6000000
     M2_signal_shift env_per_env        env3        5.4 0.53168831 0.6333333
 M3_prevalence_drift env_per_env        env3        2.3 0.33333333 0.5000000

## Predeclared reading

- M1 (null shifts) is the calibrated null probe: min-aggregation should keep FDP at or below the pooled fit.
 - M2 (signal shifts) is the honest cost probe: min-aggregation trades TPR for transportability.
 - M3 (prevalence drift) tests the prior-shift regime.
 - FDP+ is a diagnostic minimizer, not a universal FDR guarantee; simulation-only.
