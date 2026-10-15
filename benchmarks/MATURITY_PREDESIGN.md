# Predeclared designs — maturity workstreams WS-2/3/4/5

Date: 2026-10-14. Written BEFORE any benchmark run of this cycle, per the
release-governance pattern (predeclare, then run, then report — including
adverse outcomes). Reuses the F5 harness conventions: paired data/seed per
configuration, empirical FDP = |selected nulls| / |selected| at the FDP+
threshold, nominal level 0.1, honesty caveats carried in script headers and
summaries. All runs are simulation benchmarks; none constitute application
evidence.

## WS-2 — RF adapter graduation (primary learner: `rf`/ranger)

Predeclared design matrix (all cells run; adverse cells reported as-is):

| Design | n | p | Signal | Family | Seeds | Bootstraps |
|---|---|---|--------|--------|-------|------------|
| A linear | 120 | 10 | 3 main effects, AR(1) rho=0.6 nulls | gaussian | 10 | 100 |
| B non-linear (F5 original) | 120 | 10 | x1*x2, x3^2, main effects | gaussian | 10 | 100 |
| C higher-dim non-linear | 200 | 50 | x1*x2, x3^2, main effects | gaussian | 10 | 100 |
| D non-linear binary | 150 | 30 | x1*x2, x3^2, main effects | binomial | 10 | 100 |

Decoys: `knockoff_equi` and `random_permutation` in every cell.
Learners: `rf` in every cell; `oblique_rf` included only as the WS-3
vehicle. Graduation criterion (predeclared): mean empirical FDP <= 0.15
(1.5x nominal) with TPR >= 0.8 on designs A–D for rf+knockoff_equi; the
permutation rows are reported to characterize the guard, not to graduate.
If any cell fails, `make_ranger_adapter` keeps an explicit design-scope
restriction in its man page instead of full promotion.

## WS-3 — Oblique RF adapter: fix-or-fence

Diagnosis first (no package changes): quantify per-feature selection-
frequency noise — SD of each feature's selection frequency across
bootstrap replicates, real vs decoy, `rf` vs `oblique_rf`, on design B.
Predeclared remedies, tested at benchmark level only (input-level, no
runtime edits): R1 = column-standardize x (unit variance) before
`stabl_fit`; R2 = n_bootstraps 200; R3 = R1 + R2. Each remedy runs on
design B and on a heterogeneous-scale variant of B (column SDs
log-uniform on [0.1, 10]). Success criterion: mean empirical FDP <= 0.15
on both designs for oblique_rf+knockoff_equi. If no remedy meets the
criterion, the adapter is FENCED: man page and pkgdown gain an explicit
not-recommended scope, and the roadmap records the adverse result.

## WS-4 — Cooperative vs early vs late fusion

Simulation: 2 views (x1 informative, x2 informative in the complementary
arm / pure noise in the redundant arm), n = 150, p_view = 25, binomial and
gaussian, 5 seeds, 100 bootstraps. Modes: early fusion, late fusion,
cooperative fusion (built-in engine). Report held-out BER (validation
split), selection-set Jaccard across seeds, and view-wise selection
shares. Cox cooperative fusion: investigate the vendored cv.multiview
layer's cox path; record ship-or-defer with rationale in the roadmap —
no implementation in this cycle.

Evaluation rule (fixed before any run): held-out BER/MSE is computed from
native `valid_predictions` for late and cooperative fusion; for early
fusion, from a downstream `stats::glm` refit on the concatenated selected
features (`early_fusion$selected_train`), because `stabl_fit()` stores no
reusable model. Selection Jaccard is computed on the mode-appropriate
selection set (early: concatenated early-fusion selection; late: per-omic
stabl selections; cooperative: multiview per-view selection). Disposition
ladder for all workstreams: promote on met criteria; scope-fence (explicit
not-recommended scope in man/pkgdown) when criteria fail but the feature is
literature-backed; remove the feature when criteria fail AND it lacks
literature backing or actively worsens default behavior.

## WS-5 — Environment-stability transportability

B1-pattern isolated mechanisms across E = 3 environments, n_env = 80,
p = 30 (3 signal): M1 location shift on 10 null features; M2 location
shift on 2 signal features; M3 outcome prevalence drift. Compare pooled
single-fit STABL vs `stabl_environment_stability()` (min-aggregation):
per-environment empirical FDP / TPR, and cross-environment selection
Jaccard. Predeclared question: does min-aggregation keep empirical FDP
at or below the pooled fit's on shifted environments (the transportability
claim under test), and what TPR does it cost?

## WS-9 — Ordinal learner graduation (predeclared before any run)

New `base_learner = "ordinal"` (elastic-net cumulative logit via
`ordinalNet`, parallel terms). Simulations: proportional-odds ground truth —
latent score `eta = X %*% beta`, outcome bands are cuts of `eta + N(0, 0.6)`.
Two designs: O1 independent features (n = 150, p = 10, 3 signals,
beta = (1.2, 1.0, 0.8, 0, ...)); O2 equicorrelated features (rho = 0.3,
same beta). Both decoy types (knockoff_equi, random_permutation), 10 seeds
per cell, 100 bootstraps, `family = "ordinal"`. Graduation criterion
(predeclared): ordinal + knockoff_equi mean FDP <= 0.15 AND mean TPR >= 0.8
on BOTH designs. Permutation rows are characterization only (expected
anti-conservative). Outcome scale: band recovery, not band-MAE. If the
criterion fails, the disposition ladder applies (scope-fence; the learner is
literature-backed but ships experimental).

## Reporting

Each script writes `benchmarks/wsX_results.csv` + `wsX_summary.md`
with the honesty caveats of F5 carried forward. Adverse outcomes are
committed, not suppressed. Runtime behavior changes (label promotions,
fences, guards) happen only AFTER results land, in separate commits that
cite the summary files.

---

## Erratum / addendum — WS-9 (dated 2026-10-15)

This addendum corrects the WS-9 record. The historical preregistration above
is preserved unchanged; the corrections below apply to the executed benchmark
and its outputs, not to the predeclared design intent.

1. **Truth-set error.** The executed script defined the signal set as
   `paste0("f", seq_along(beta))` — all 10 predictors — although only 3 of
   the 10 coefficients are nonzero (beta = (1.2, 1.0, 0.8, 0, ...)). All
   empirical FDP/TPR values in the pre-correction outputs are therefore
   invalid. The corrected script uses `paste0("f", which(beta != 0))`
   (predictors f1, f2, f3).

2. **Incomplete run grid.** The saved pre-correction results contain 38 of
   the 40 planned runs: seed 8 is missing in both decoy arms of the
   correlated design (O2_equicorr03). Root cause: seed 8's data-derived
   `cut()` thresholds produce a singleton L1 class (count 1), and the
   ordinal learner requires at least 2 observations per class per bootstrap,
   so the fit is deterministically infeasible. The historical script dropped
   the failure silently via `try(...); next`. The corrected script records
   every failure with seed, error class, and message in `ws9_failures.csv`
   and fails the run (nonzero exit) if any planned combination is missing.

3. **Non-exact proportional-odds generator.** The WS-9 DGP uses Gaussian
   latent noise and data-derived `cut()` thresholds. It is therefore not an
   exact cumulative-logit / proportional-odds data-generating process. The
   historical experiment is preserved as-is; a genuinely cumulative-logit
   design with fixed ordered thresholds and logistic latent noise is a
   separate, separately predeclared benchmark (`ws9_revalidation.R`,
   PR6-S), not a rerun of this design.

The pre-correction outputs are archived under `benchmarks/archive/` and
must not be cited as valid evidence. No graduation inference from the
pre-correction run survives.
