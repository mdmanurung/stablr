# stablr API Reference

This file is a human-readable map of the exported `stablr` API. The canonical
function documentation is generated from roxygen comments in `R/` and lives in
`man/`; the documentation website groups the same functions through
`_pkgdown.yml`.

## Core STABL Engine

### `stabl()`

The high-level entry point for new users: a data frame in, a plain-language
selection out. Runs STABL with statistically safe defaults — model-X
knockoff artificial features, FDP+ calibration, strict fallback policy,
automatic lambda grid — and prints a plain-language answer. Use
`get_feature_names()` / `get_support()` for the selection and
`get_fdrs()` for the FDP+ curve; `plot()` shows the stability path and
FDP+ diagnostic.

### `stabl_fit()`

Fits the core bootstrap stability-selection procedure for one feature matrix.

Core arguments:

- `x`: numeric matrix or data frame with samples in rows and features in columns.
  Non-empty row names are sample IDs; column names are feature names.
- `y`: named outcome vector, factor, or matrix-like outcome such as
  `survival::Surv`. Names or row names must identify the same samples as
  `rownames(x)`.
- `lambda_grid`: data frame, usually from `auto_lambda_grid()`. Each row is one
  learner setting evaluated during bootstrapping.
- `family`: current examples use `"gaussian"`, `"binomial"`,
  `"multinomial"`, and `"cox"`.
- `n_bootstraps`: finite integer-like scalar; number of bootstrap iterations.
- `artificial_type`: `"random_permutation"`, `"knockoff"`,
  `"knockoff_equi"`, `"knockoff_mvr"`, or `NULL`.
- `groups`: optional named vector whose names match `rownames(x)`; grouped
  bootstrapping keeps whole groups together.
- `random_state`: optional finite integer-like scalar seed.

Returned `stabl_fit` objects support `print()`, accessors, diagnostic plots,
and disk export helpers.

### `auto_lambda_grid()`

Builds data-driven `glmnet` lambda grids. Use it before `stabl_fit()` and pass
the returned data frame as `lambda_grid`.

### `compute_fdp_plus()`

Computes FDP+ threshold diagnostics from original and artificial feature
stability scores.

## Learner Selection

Learners are selected by name via the `base_learner` argument of
`stabl_fit()` (`"lasso"`, `"elastic_net"`, `"adaptive_lasso"`,
`"sparse_group_lasso"`, `"rf"`, `"ordinal"`, `"splsda"`). The fitting path
uses internal batch adapters; the former single-row factory functions were
internalized in 1.0.0. `auto_lambda_grid()` is the single exported grid
helper; the learner-specific grids (`auto_rf_threshold_grid()`,
`auto_ordinal_lambda_grid()`, `auto_splsda_keepx_grid()`) are internal
machinery behind `lambda_grid = "auto"`. The experimental oblique random
forest (`"oblique_rf"`) was retired in 1.0.0 (see
`docs/MIGRATION_0.3_TO_1.0.md`).

## Artificial Features

- Artificial-feature generation is internal machinery behind the
  `artificial_type` argument of the fit functions
  (`make_artificial_features()` and the per-type constructors
  `make_rp_features()`, `make_knockoff_features()`,
  `make_knockoff_equi_features()`, `make_knockoff_mvr_features()` are not
  exported); select the construction by string (table below) rather than
  calling factories directly.

Supported `artificial_type` values in fit functions are:

| Value | Meaning |
|---|---|
| `"random_permutation"` | Column-permuted decoys; no optional dependency. |
| `"knockoff"` | Fixed-X knockoffs. |
| `"knockoff_equi"` | Model-X equicorrelated Gaussian knockoffs. |
| `"knockoff_mvr"` | Model-X MVR Gaussian knockoffs. |
| `NULL` | No artificial features; FDP diagnostics are unavailable. |

## Input Validation and Bootstrapping

Input validation and bootstrap index generation are internal machinery
invoked by the fit functions (`validate_sample_alignment()`,
`validate_multiomic_inputs()`, `sample_bootstrap_indices()`,
`group_bootstrap_indices()` are not exported). Users control resampling
through the `stabl_fit()` arguments `n_bootstraps`, `stratify_bootstrap`,
`bootstrap_strata`, and `sample_groups`.

Current validation rules are strict and name based:

- `x` row names, `y` names or row names, and optional `groups` names must refer
  to the same unique sample IDs.
- Multi-omic matrices must have the same row names in the same order.
- Feature names, when present, must be non-empty and unique.
- Numeric counts and seeds are finite integer-like scalars.
- Plot titles, file paths, and file-format arguments are non-empty character
  scalars.

## Accessors

- `get_support(object, new_hard_threshold = NULL)` returns a named logical
  support mask. Use it for logical subsetting or counts such as
  `sum(get_support(fit))`.
- `get_feature_names(object, new_hard_threshold = NULL)` returns selected
  feature names as a character vector. Use it when displaying feature lists or
  selecting columns by name.
- `transform_stabl(object, x, new_hard_threshold = NULL)` subsets new data to
  selected features in fitted feature order and preserves a two-dimensional
  matrix or data frame.
- `get_stabl_scores()` returns the feature-by-lambda stability-score matrix.
- `get_importances()` returns maximum-over-lambda feature scores.
- `get_cooperative_features()` returns cooperative-fusion selected features.
- `get_cooperative_diagnostics()` returns cooperative-fusion tuning diagnostics.

## Multi-Omic Workflows

### `stabl_multiomic_train_validate()`

Runs per-omic STABL plus optional early, late, and cooperative fusion on
train/validation splits.

Important arguments:

- `x_train_list` and `x_valid_list`: named lists of aligned omic matrices.
- `y_train` and `y_valid`: named outcomes aligned to the omic sample IDs.
- `lambda_grid`: named list of per-omic lambda grids, or a shared grid where
  supported.
- `early_fusion`: concatenate omics and fit a joint STABL model.
- `late_fusion`: learn a weighted stack over per-omic downstream predictions.
- `late_fusion_training`: only `"oof"` (default) is supported; it reruns
  selection and downstream fitting inside each stacking fold. The historical
  in-sample `"python_legacy"` mode was removed from the public API in 1.0.0
  (retained unexported for parity tests).
- `late_fusion_nfolds`: number of leakage-safe stacking folds (default 5).
- `n_iter_lf`: finite integer-like scalar number of late-fusion weight draws.
- `cooperative_fusion`: run the native cooperative multiview branch.

`print.stabl_multiomic_fit()` reports branch presence and selected-feature
counts. When late fusion is present, it reports the late-fusion score rather than validation-metric labels from older drafts.

Binary and regression late fusion use parity-preserving batched stacking. OOF
mode costs approximately `nfolds + 1` complete fits because it adds a final
full-training refit after weights are frozen. Fold IDs, seeds, selected
features, artificial-feature provenance, warnings, and fallbacks are returned
under `$late_fusion$provenance`.
Multiclass stacking remains scalar by design to preserve probability
normalization and tie behavior.

If a fold selects no usable features or its downstream learner cannot be fit,
OOF mode predicts from that fold's training data only: the mean for regression,
the event prior for binary outcomes, and class priors for multiclass outcomes.
The reason is recorded rather than silently treated as a fitted learner.

### `stabl_multiomic_cv()`

Runs outer cross-validation over named multi-omic inputs. Optional early, late,
and cooperative fusion settings are propagated into each outer fold.

### `stabl_multiomic_nested_cv()`

Runs nested multi-omic cross-validation for benchmark workflows. Cooperative
fusion is not forwarded by this nested-CV helper. Nested CV estimates predictive
performance for the specified workflow; it does not establish dataset-wide FDP
control, and small outer-training folds can legitimately trigger the documented
late-fusion prior fallbacks.

### `stabl_multiomic_stacked()`

Performs random-search late-fusion weight selection for `task_type` values `"binary"`, `"regression"`, or `"multiclass"`. (The pre-0.3.0 name `stacked_multi_omic()` was removed in 1.0.0.)

The pinned cross-language solver evidence compares feature rankings and support,
not bit-identical coefficients. R and Python use different optimisation paths;
exact contracts are reserved for accumulation, FDP+, masks, metrics, and the
committed stacking-candidate matrix.

### `load_ool_data()`

Loads the bundled OOL example subset with aligned `cytof` and `proteomics`
matrices plus a named DOS outcome vector.

## Analysis Workflow

### `stabl_analyze(x, y, family, mode, ...)`

Explicit evidence-level analysis. The `mode` argument is required:

- `mode = "exploratory"` fits STABL on the full data and returns a
  canonical interpretation table labelled `post_selection_descriptive`
  (describes the fitting data; not predictive evidence).
- `mode = "validated"` runs nested cross-validation with fold-fitted
  preprocessing and returns held-out predictions, pooled and per-repeat
  performance, and held-out permutation importance labelled
  `cross_validated`. Binary outcomes report accuracy, balanced error rate,
  macro-F1, MCC, Brier, log-loss, and AUROC; gaussian outcomes report MAE,
  RMSE, and R². With `final_refit = TRUE` the result can be finalized into
  a portable model (see `stabl_finalize()`).

### `stabl_interpret()`

Builds the canonical interpretation table (coefficients, stability, effect
directions, correlation groups) for an exploratory analysis.

### Preprocessing

- `stabl_preprocess_fit()` / `stabl_preprocess_apply()` fit and apply
  fold-safe preprocessing so nothing learned on held-out data leaks into
  training. `stabl_preprocess_fit()` also accepts an unprepared
  `recipes::recipe` as `transforms`, bridging to a `recipes` workflow.

## Environment Stability

- `stabl_environment_stability()` runs the ICP-inspired transportability
  workflow across environments (batches, cohorts, sites); min-aggregation
  is the calibrated read within the documented WS-5 scope (small feature
  panels).
- `get_min_scores()` extracts the per-environment minimum scores.

## Reliability and Final Models

### `stabl_diagnose(fit, ...)`

Read-only reliability report for a fitted `"stabl_fit"`: machine-readable
`warning_codes` (FDP+ curve saturation/boundary/tie/low-decoy probes),
the `cutoff` block disclosing the selection rule, comparator, and
denominator conventions, the full `threshold_path`, and a
benchmark-provenanced `evidence_profile` for the fit's learner x decoy x
outcome x dimensionality regime. An empty `warning_codes` is not a
validity guarantee.

### `stabl_null_check(x, y, design, ...)`

Outcome-null validation: refits the same fit specification under an
explicit exchangeability design (`"independent"`, `"within_block"`,
`"whole_cluster"`, or `"simulated"`) and reports how often anything is
selected (`p_any_selection`). Invalid design/structure combinations are
refused fail-closed; replicate seeds are derived replayably from `seed`.

### `stabl_replication_profile(x, y, study_labels, ...)`

Cohort replication profile: pooled stability, leave-one-study-out training
replication, optional inverse-variance pooling of user-supplied per-study
effect estimates (Cochran's Q), and optional held-out downstream
performance. Preprocessing is recomputed inside each training partition;
no held-out cohort chooses thresholds.

### `stabl_finalize(development_result, x, y, preprocessing = NULL, ...)`

Packs a `"stabl_fit"` or a validated `"stabl_analysis"` (requires
`final_refit = TRUE`) into a portable `"stabl_predictor"` for binomial or
gaussian outcomes. Preprocessing is stored, not refit at predict time; no
training observations are serialized.

### `predict.stabl_predictor(object, newdata, type = NULL, extra_columns = ...)`

Predicts from a portable final model. `type` defaults to `"response"`;
binomial models also support `"class"` and `"probability"` (with the
stored event-class mapping). New data is validated internally against the
stored schema (`validate_newdata_schema()` is not exported): missing,
renamed, or duplicated columns always error; extra columns error unless
`extra_columns = "drop"`.

## Cooperative Fusion

Cooperative fusion is available through `stabl_multiomic_train_validate()` and
`stabl_multiomic_cv()` with `cooperative_fusion = TRUE`.

Key arguments:

- `rho`: numeric scalar or vector of non-negative cooperation strengths.
- `cooperation_selection`: only `"cv"` since 1.0.0; the `"validation"` mode
  tuned on the assessment set and was removed from the public API (F-CV-04).
- `cooperation_selector`: `"lambda.min"` or `"lambda.1se"`.
- `cooperation_type_measure`: tuning metric; defaults depend on `family`.
- `cooperation_nfolds`: finite integer-like scalar, at least 3.

Native cooperative fusion currently supports gaussian and binomial families.
Cox and Poisson cooperative fusion are intentionally rejected.

## Visualization

- `plot_stabl_path(object, title = "STABL Stability Path")` plots stability
  scores over the lambda path.
- `plot_fdr_graph(object, title = "FDP+ diagnostic", reference_fdp_plus = 0.05)`
  plots the FDP+ diagnostic over candidate thresholds; `reference_fdp_plus`
  draws a purely visual reference line (renamed from `fdr_target` in 1.0.0).
The ROC/PRC and per-feature distribution helpers (`plot_roc()`,
`plot_prc()`, `boxplot_features()`, `scatterplot_features()`) were
internalized in 1.0.0; use general plotting packages (the TCGA vignette
shows inline ggplot2 replacements). `save_stabl_results()` still uses the
internal implementations when writing report plots.

## Export

### `save_stabl_results()`

```r
save_stabl_results(
  object, path, x = NULL, y = NULL,
  figure_fmt = "pdf",
  new_hard_threshold = NULL,
  task_type = "binary",
  override = FALSE,
  what = c("all", "csv")
)
```

With the default `what = "all"`, writes score tables, selected-feature
tables, diagnostic plots, and feature-distribution plots for one
`stabl_fit` object.  `x` and `y` are required in this mode for feature
distribution plots.  `task_type` is one of `"binary"`, `"multiclass"`, or
`"regression"`.  Set `override = TRUE` to reuse an existing output
directory.

With `what = "csv"`, writes only the score tables below and requires
neither `x`/`y` nor ggplot2; an existing directory is reused without
error.

Score tables (both modes):

| File | When written |
|---|---|
| `STABL scores.csv` | Always. |
| `Max STABL scores.csv` | Always. |
| `STABL artificial scores.csv` | When artificial features were used. |
| `Max STABL artificial scores.csv` | When artificial features were used. |

Additional artifacts with `what = "all"`:

| File | When written |
|---|---|
| `FDR Graph.<fmt>` | When artificial features and FDP diagnostics are present. |
| `Stability Path.<fmt>` | Always. |
| `Selected Features/Selected features.csv` | Always. |
| `Selected Features/Feature distributions.<fmt>` | When at least one feature is selected. |

`export_stabl_to_csv()` was removed ahead of 1.0.0; use
`save_stabl_results(what = "csv")`.

## Selection Similarity Metrics

- `stabl_similarity()` computes scalar or pairwise similarities between
  selected-feature sets, with `method = "jaccard"` (default), `"adjusted"`,
  or `"pearson"`; `as_matrix = TRUE` returns the full pairwise matrix.
- `stabl_similarity_summary()` reduces pairwise values to a single
  median/mean statistic with a spread measure.

The eight 0.3.0 names (`jaccard_similarity()`,
`jaccard_similarity_matrix()`, `adjusted_similarity()`,
`adjusted_similarity_pairwise()`, `adjusted_similarity_summary()`,
`pearson_similarity()`, `pearson_similarity_pairwise()`, and
`pearson_similarity_summary()`) were removed ahead of 1.0.0.

The realized-performance benchmark helpers `fdr_similarity()`,
`tpr_similarity()`, and `fscore_similarity()` were internalized in 1.0.0
(they are benchmark metrics against known truth, not general
similarities).
