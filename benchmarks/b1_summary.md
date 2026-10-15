# B1 — study-disjoint nested CV benchmark

Configuration: 6 studies x 25 samples, p = 50 (3 signal); outer_v = 3, inner_v = 2, n_bootstraps = 100; 10 seeds; lasso base learner; artificial_type = "random_permutation".

Two batch configurations (isolated mechanisms): `location_shift` =
per-study N(0, 0.6) location shift on 10 null features, balanced prevalence;
`prevalence_drift` = per-study outcome prevalence ~ Uniform(0.3, 0.7), no feature shifts.

| Configuration | sample-iid BER (SD) | study-disjoint BER (SD) | gap |
|---|---|---|---|
| location_shift | 0.2997 (0.0643) | 0.2901 (0.0681) | -0.0096 |
| prevalence_drift | 0.3002 (0.0508) | 0.3396 (0.0703) | +0.0395 |

Interpretation: the gap is sample-iid minus study-disjoint mean held-out
BER. A positive gap means the sample-i.i.d. folds understate held-out
error when studies carry the corresponding structure. The
location-shift-only row is a deliberate null probe: class-independent
location shifts on null features under balanced prevalence give
sample-iid folds nothing to leak, so an approximately zero gap is the
calibrated outcome there, and a positive gap is expected only under
outcome-prevalence drift (prior shift). See the honesty caveats in the
script header (simulation-only; two batch models; reduced effective
training size in the study-disjoint arm; no reharmonization claim).

