# G0 baseline inventory — unified roadmap

**Date:** 2026-10-15
**Roadmap baseline:** `main` at `15fad15fb7cca5b90e8ee2054f49d1f4d76de3a5` (package v0.3.0, 2026-10-14)
**Working branch:** `docs/pkgdown-0.3.0`, post-merge SHA `108465b` (merge of `origin/main` into the branch; main-side delta was DESCRIPTION ordinalNet registration only, already mirrored on the branch)

## Environment

- R 4.4.2 (2024-10-31), conda toolchain, `R_LIBS=/workspace/Rlib-shared`
- glmnet 5.0, knockoff 0.3.6, ordinalNet 2.14, ranger 0.18.0, aorsf 0.1.6,
  ggplot2 4.0.3, testthat 3.3.2, recipes 1.4.0, pkgdown 2.2.1
- sgl NOT installed (sparse-group adapter tests skip by design)
- Package version 0.3.0; `Config/testthat/edition: 3`; roxygen2 8.1.1

## Test baseline (exact commands and results)

Command pattern: `R_LIBS=/workspace/Rlib-shared Rscript -e 'testthat::test_local(".", filter=<f>, ...)'`

| Suite | fail | pass | skip |
|---|---|---|---|
| nested-cv | 0 | 150 | 0 |
| late-fusion-oof | 0 | 59 | 0 |
| learner-ordinal | 0 | 20 | 0 |
| multiomic-workflows | 0 | 185 | 0 |
| environment-stability | 0 | 41 | 0 |
| **full suite** | **0** | **3040** | **1** |

## Selection / cutoff / decoy contracts (recorded, not altered)

- `R/stabl_fit.R`: bootstrap learner contract; per-bootstrap fit on the full
  lambda grid; selection frequency = fraction of bootstrap resamples in which
  a feature has a nonzero coefficient at a given lambda; per-feature score =
  max over the grid. `hard_threshold` selects directly; otherwise the
  FDP+-optimal threshold `fdr_min_threshold_` is used.
- `R/stabl_accessors.R`: `get_support()`, `get_feature_names_out()`,
  `get_importances()` (max-over-lambda stability score, sorted),
  `get_stabl_scores()` (raw feature × lambda matrix), `get_fdrs()`,
  `.resolve_stabl_threshold()` (hard → FDP+-optimal → override).
- `R/fdp_control.R`: FDP+ = (1 + E[decoys above threshold]) / (selected
  count), capped at 1; decoys from `artificial_type` (random_permutation or
  knockoff variants). Parity fixtures sealed — do not alter.
- `R/nested_cv.R`: nested classification CV with study-disjoint folds;
  `.predict_selected_multinomial()` fallback chain (see PR0a).
- `R/multiomic_workflows.R` / `R/late_fusion.R`: per-omic fits, early/late
  fusion, OOF stacking weights; zero-feature binary fallback = 0.5 (see PR0a).

## Verified defects carried into the roadmap (from roadmap §2)

1. `.predict_selected_multinomial()` majority-class fallback for <2 selected
   features discards exactly-one-feature signal.
2. `.late_fusion_fit_omic()` zero-feature binary prediction fixed at 0.5
   regardless of training prevalence.
3. `late_fusion.R` `oof_predictions` is base-model OOF, not full-stack OOF.
4. WS-9 benchmark truth set uses all 10 indices though only 3/10 coefficients
   are nonzero; 38/40 runs saved (seed 8 missing in both correlated arms).
5. WS-9 generator is Gaussian-latent with data-derived `cut()` thresholds —
   not an exact cumulative-logit DGP.
6. WS-8 aging showcase selects on all eight cohorts before LOCO prediction;
   selected-panel MAE 11.72 y is post-selection/apparent.
7. `learner_ordinal.R` sorts requested lambdas and can substitute a nearest
   fitted lambda; ordinal tests wrongly require nested supports.
8. Multiview validation used for cooperation strength is not an untouched
   test set.

## Release restrictions in force

- SCI-04 null-selected-fraction gate: fail-closed at `proposed_unaccepted`
  pending external review (`inst/analysis/review/SCI-04-independent-review-packet.md`).
- Sealed 0.1.1 negative parity artifact retained; mutable-checkout runs do
  not qualify as release evidence.
- RF/knockoff scope (WS-2), oblique-RF fence (WS-3), >3,000-predictor
  restriction, and cooperative-fusion caveats (WS-4) remain as documented.

**Gate G0: PASSED** — reproducible baseline inventory recorded; full suite green.
