# WS-8 — Application showcase predesign: stabl on the human immune aging atlas

Written BEFORE any showcase analysis. Data source: Terekhova et al.,
Immunity (2026), "Single-cell analyses of global cohorts outline determinants
of immune aging and link the ratio of GZMK+ to GZMB+ Tem cells to health
trajectories". Public portal file `all_cells_metadata.csv.gz`
(13,137,775 cell rows, 2,609 donors, 8 cohorts), downloaded 2026-10-13,
md5 50651968855bf6ca3efcdf199e2a9fd5.

## Role of this showcase

This is an APPLICATION CASE STUDY. It demonstrates the end-to-end stablr
workflow on real published data. It is NOT calibration evidence: simulation
benchmarks (WS-2/WS-4/WS-5, F5) remain the only calibration evidence. Any
FDP+ number reported here is a diagnostic computed on one real dataset.

## Data preparation (predeclared)

- Unit of analysis: donor. Cell-level metadata is aggregated to one row per
  donor (2,609 donors).
- Features: per-donor proportions of the 59 harmonized cell types
  (`annotation` column; the `Bad quality` class is excluded from both
  numerator and denominator). Proportions get add-half smoothing
  ((count + 0.5) / (total + 0.5)) and a logit transform, so every feature is
  an unconstrained numeric column.
- Covariate features enter the design matrix alongside cell types: sex
  (female = 1), CMV (positive = 1), and ancestry indicators from the
  `ethnicity` column (levels with fewer than 25 donors are pooled into
  `Other`). Selected covariate features are reported separately from cell
  types, never dropped silently.
- Outcome: donor age in years (gaussian).
- Environments: the 8 cohorts (`Dataset`): ABF300, AIDA, CIMA, KCL, OneK1K,
  SPAC, SoundLife, UCSF.

## Models (predeclared)

- Primary: `stabl_environment_stability` with `base_learner = "rf"`,
  `artificial_type = "knockoff_equi"`, 200 bootstraps per environment,
  `random_state = 20261013`. RF + equicorrelated knockoffs is the lead
  learner because it is the only pairing with graduated calibration evidence
  on linear and low-dimensional nonlinear designs (WS-2, F5). Known scope
  limits from WS-2 (high-dimensional C design, binomial D design) are
  restated in the report; this dataset is gaussian, p = 59 + covariates
  (< 100), and per-cohort n between 188 and 981.
- Sensitivity: the same call with `base_learner = "lasso"`. The F5 benchmark
  showed lasso + knockoff_equi FDP+ can be anti-conservative in the
  learner-regime benchmark (mean FDP 0.578 at p = 10); lasso results are
  reported as sensitivity, not headline.
- Consistency criterion (lead): min-aggregation — a feature is
  age-consistent when its min-over-cohorts max-over-lambda stability score
  is >= 0.5. This is stablr's documented contract and is stricter than the
  paper's "significant in >= 5 of 8 cohorts" rule. The paper's rule is
  reported as a comparison column, not as the headline.

## Age-prediction calibration (predeclared)

- Task: predict donor age from the selected features, evaluated
  leave-one-cohort-out (LOCO) across the 8 cohorts.
- Models: (a) ridge regression on all 59 cell-type features (unselected
  reference), (b) linear model on the min-agg selected set.
- Metrics: per-held-out-cohort MAE, pooled LOCO MAE, calibration slope and
  intercept of predicted vs observed age. Comparator context: the paper
  reports MAE 8.2 y (test, full 55-feature clock) and 8.9 y (reduced
  20-feature clock). Our clock uses the same kind of input (cell-type
  proportions) but a different feature screen; the claim is comparability,
  not superiority.

## Success criteria (predeclared)

1. The pipeline runs end to end on the real portal file with zero errors.
2. FDP+ diagnostics are computable for every per-environment fit.
3. The min-agg selected set is non-empty and each selected feature gets a
   directional read (sign of the per-cohort age association).
4. LOCO age MAE and calibration curves are reported for both models.
5. A comparison table maps the selected set against the paper's consistent
   findings (MAIT down, CD4+/CD8+ ratio up, CD8 Temra up, RTE down,
   CD8 naive down, Tem GZMK+/GZMB+ changes).

## Outputs

- `inst/benchmarks/ws8_results.csv` — per-feature per-cohort stability
  scores, selection flags, per-cohort age coefficients.
- `inst/benchmarks/ws8_summary.md` — narrative summary with the tables above.
- Figures (PNG, in the case-study vignette): stability path, environment
  min-score ranking, LOCO calibration curves, consistency heatmap.

---

## Addendum — LOCO evidence label (dated 2026-10-15)

This addendum documents a departure from the predeclared LOCO intent and the
planned correction. The historical preregistration above is preserved
unchanged.

**Departure.** The executed `ws8_aging_showcase.R` selects the 24-feature
panel using all eight cohorts *before* the leave-one-cohort-out (LOCO) ridge
prediction. The selected-panel LOCO MAE (historically 11.72 y) is therefore
a **post-selection** estimate: it is apparent for the selection procedure,
not an unbiased held-out error for the complete workflow. The all-feature
LOCO column (12.14 y) is unaffected and remains a genuine held-out
evaluation, because the held-out cohort plays no role in that feature set.

**Intended correction (PR5, not yet executed).** A selection-isolated LOCO
will re-run the full pipeline — filtering, `stabl_environment_stability()`
selection, ridge tuning — inside each seven-cohort training partition, and
predict the eighth cohort with the fold-specific panel. This is a
prospective correction; it is not claimed here as prespecified, and no
corrected numbers are inserted until the rerun completes.

**Preservation.** All historical MAE values, figures, CSVs, and the dataset
are preserved unchanged. The vignette and `ws8_summary.md` carry updated
evidence labels only; no historical result was edited or regenerated.

---

## Addendum 2 — execution deviations and post-review analyses (dated 2026-10-15)

An adversarial review of the vignette and artifacts found deviations from
this predesign and reporting gaps. This addendum documents them, records the
follow-up analyses executed on 2026-10-15, and predeclares the
selection-isolated LOCO protocol before its run. The historical
preregistration above is preserved unchanged.

### Execution deviations (review finding F7)

1. **Cohort sizes.** This predesign says "per-cohort n between 188 and
   981". The executed data has SPAC 33, SoundLife 96, and UCSF 99. The
   claim was wrong when written.
2. **Model (b).** Predeclared as "a linear model". Executed as ridge
   (`glmnet`, alpha = 0, 5-fold internal CV, `lambda.min`).
3. **Model (a) feature count.** Predeclared as 59 cell-type features.
   Executed with 60: the 59 cell types plus `COV:sex`.
4. **Stability-path figure.** Promised under Outputs; not delivered.
5. **Reduced-clock comparator.** The 8.9 y reduced-clock MAE was promised
   as comparator context; the vignette cited only the 8.2 y full-clock
   value. The vignette now states both are not directly comparable to this
   protocol.
6. **Covariate reporting.** The predesign promises covariates are "never
   dropped silently". CMV and the two ancestry indicators were dropped
   (zero variance inside at least one cohort) and are disclosed in the
   vignette. Sex entered the all-feature model; whether it survived
   selection was not stated. It did not: the 24-feature panel contains no
   sex feature. The vignette now states this.

### Reproducibility gap (review finding F9)

The donor-level count table (`donor_celltype_counts.tsv`,
`donor_meta.tsv`) was derived from `all_cells_metadata.csv.gz` before
`ws8_aging_showcase.R` runs. The derivation script is not archived in the
repository, so the predeclared exclusion of the `Bad quality` class cannot
be verified from the package source. The vignette states this gap instead
of asserting the exclusion.

### Review correction

The review stated SPAC is the binding cohort for 14 of the 24 selected
features. The artifact-backed count is **10 of 24 by strict argmin** (12 of
24 counting ties; two features have a selection frequency of 1.000 in all
eight cohorts, so all eight cohorts tie at the minimum). The vignette uses
the corrected counts.

### Criterion 2 closure: per-environment FDP+ diagnostics (executed)

The FDP+ diagnostics were computed inside every per-environment `stabl_fit`
(field `min_fdr_`) but never extracted or reported.
`ws8_review_fdp.R` extracts them (`ws8_fdp_diagnostics.csv`). Result:
**adverse**. FDP+ at the 0.5 bar is 1.03–1.05 in every cohort (decoy
features pass the bar as often as real features); the minimum over the
threshold grid is 0.89–1.02. The WS-5 calibration (0.083–0.142) does not
transfer to this atlas. The vignette reports the min-agg rule as a
descriptive reproducibility filter, not calibrated FDR control.

### Min-agg sensitivity without SPAC (executed)

`ws8_review_sensitivity_nospac.R` recomputes the min-agg rule over the
seven cohorts excluding SPAC (`ws8_minagg_sensitivity_nospac.csv`,
`ws8_binding_cohorts.csv`). Result: the panel grows from 24 to 33
features; 9 features join, none leave; CD8+ Tem GZMK+ rises 0.505 → 0.660;
the other seven predeclared scores are unchanged. Given the FDP+ result,
neither panel size is calibrated; the sensitivity shows panel membership
depends on the smallest cohort.

### Paired per-donor MAE differences (executed)

`ws8_review_paired_mae.R` computes the per-donor difference of absolute
errors (all-feature minus selected panel) with 2,000-resample bootstrap
intervals, seed 20261013 (`ws8_paired_mae_diff.csv`). Result: CIMA +7.98
[7.30, 8.67], UCSF +2.59 [1.91, 3.27], OneK1K −2.18 [−2.47, −1.90],
SoundLife −2.81 [−3.92, −1.65]; ABF300, KCL, SPAC inconclusive; pooled
+0.50 [0.29, 0.70]. The paired intervals resolve the marginal-interval
overlap at UCSF. The covariate asymmetry (sex in the all-feature model
only) is disclosed as a confounder of the comparison.

### Selection-isolated LOCO (PR5) — protocol predeclared before the run

For each held-out cohort h, in fixed alphabetical order:

1. Partition: train on the other seven cohorts; hold out h.
2. Selection: `stabl_environment_stability()` on the training partition
   with the primary settings (rf, `knockoff_equi`, 200 bootstraps,
   `sample_fraction = 1`, `replace = TRUE`, `random_state = 20261013`).
3. Panel: min-agg rule at >= 0.5 over the seven **training** cohorts.
4. Prediction: ridge (`glmnet`, alpha = 0, `nfolds = 5`, `lambda.min`) on
   the within-fold panel; predict h.
5. Primary endpoint: per-cohort MAE, compared with the all-feature and
   post-selection panel columns. Secondary: per-fold panel size, overlap
   with the full-data panel, calibration slope.
6. Fallback: if a fold selects no features, the fold uses the all-feature
   model and is flagged in the output.

Script: `ws8_loco_selection_isolated.R`. Output:
`ws8_loco_silo.csv` (per-cohort results) and
`ws8_loco_silo_features.csv` (per-fold panels). No corrected number enters
the vignette until this run completes.
