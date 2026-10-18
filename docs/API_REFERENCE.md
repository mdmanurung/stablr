# stablr API Reference

This file is a human-readable map of the exported `stablr` API. The canonical
function documentation is generated from roxygen comments in `R/` and lives in
`man/`; the documentation website groups the same functions through
`_pkgdown.yml`.

## Core STABL Engine

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

## Learner Adapters

- `make_glmnet_adapter()` creates a single-lambda glmnet learner adapter.
- `make_adaptive_lasso_adapter()` creates a single-lambda adaptive-lasso
  adapter.
- `make_sgl_adapter()` creates a single-lambda sparse-group-lasso adapter.

The core fitting path uses internal batch adapters for efficiency, but the
exported factories remain available for adapter-level workflows and tests.

## Artificial Features

- `make_artificial_features()` dispatches artificial-feature generation.
- `make_rp_features()` creates column-permuted decoys.
- `make_knockoff_features()` creates fixed-X knockoffs via the optional
  `knockoff` package.
- `make_knockoff_equi_features()` creates model-X equicorrelated Gaussian
  knockoffs via the optional `knockoff` package.
- `make_knockoff_mvr_features()` creates model-X MVR Gaussian knockoffs via
  the optional `knockoff` package.

Supported `artificial_type` values in fit functions are:

| Value | Meaning |
|---|---|
| `"random_permutation"` | Column-permuted decoys; no optional dependency. |
| `"knockoff"` | Fixed-X knockoffs. |
| `"knockoff_equi"` | Model-X equicorrelated Gaussian knockoffs. |
| `"knockoff_mvr"` | Model-X MVR Gaussian knockoffs. |
| `NULL` | No artificial features; FDP diagnostics are unavailable. |

## Input Validation and Bootstrapping

- `validate_sample_alignment()` validates single-matrix sample alignment.
- `validate_multiomic_inputs()` validates named multi-omic input lists.
- `classic_bootstrap_indices()` draws ordinary bootstrap/subsample indices.
- `group_bootstrap_indices()` draws group-aware bootstrap/subsample indices.

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
- `get_feature_names_out(object, new_hard_threshold = NULL)` returns selected
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

### `stacked_multi_omic()`

Performs random-search late-fusion weight selection for `task_type` values `"binary"`, `"regression"`, or `"multiclass"`.

The pinned cross-language solver evidence compares feature rankings and support,
not bit-identical coefficients. R and Python use different optimisation paths;
exact contracts are reserved for accumulation, FDP+, masks, metrics, and the
committed stacking-candidate matrix.

### `load_ool_data()`

Loads the bundled OOL example subset with aligned `cytof` and `proteomics`
matrices plus a named DOS outcome vector.

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
- `plot_roc(y_true, y_preds, title = "ROC Curve")` plots ROC curves for
  downstream binary predictors.
- `plot_prc(y_true, y_preds, show_iso = TRUE, title = "Precision-Recall Curve")`
  plots precision-recall curves.
- `boxplot_features(features, x, y, title = "Selected Features", ncol = 3L)`
  plots selected features by class.
- `scatterplot_features(features, x, y, title = "Selected Features", ncol =
  3L)` plots selected features against a continuous outcome.

ROC/PRC helpers require explicit 0/1 numeric or logical outcomes with both
classes present, plus finite predicted probabilities in `[0, 1]`. Character or
factor outcomes should be recoded by the caller so the positive class is
unambiguous. Tied prediction scores are aggregated at a single threshold before
ROC AUC and precision-recall average-precision computation, so metrics do not
depend on row order. Plot titles are validated as non-empty character scalars.
`ncol` and similar numeric controls are finite integer-like scalars.

## Export

### `export_stabl_to_csv(object, path)`

Writes STABL score matrices and max-score tables:

| File | When written |
|---|---|
| `STABL scores.csv` | Always. |
| `Max STABL scores.csv` | Always. |
| `STABL artificial scores.csv` | When artificial features were used. |
| `Max STABL artificial scores.csv` | When artificial features were used. |

### `save_stabl_results()`

```r
save_stabl_results(
  object, path, x, y,
  figure_fmt = "pdf",
  new_hard_threshold = NULL,
  task_type = "binary",
  override = FALSE
)
```

Writes score tables, selected-feature tables, diagnostic plots, and
feature-distribution plots for one `stabl_fit` object. `x` and `y` are required
for feature distribution plots. `task_type` is one of `"binary"`,
`"multiclass"`, or `"regression"`. Set `override = TRUE` to reuse an
existing output directory.

Main output artifacts:

| File | When written |
|---|---|
| `STABL scores.csv` | Always. |
| `Max STABL scores.csv` | Always. |
| `FDR Graph.<fmt>` | When artificial features and FDP diagnostics are present. |
| `Stability Path.<fmt>` | Always. |
| `Selected Features/Selected features.csv` | Always. |
| `Selected Features/Feature distributions.<fmt>` | When at least one feature is selected. |

## Selection Similarity Metrics

- `jaccard_similarity()` and `jaccard_matrix()`
- `adjusted_similarity()`, `adjusted_similarity_values()`, and
  `adjusted_similarity_measure()`
- `pearson_similarity()`, `pearson_similarity_values()`, and
  `pearson_similarity_measure()`
- `fdr_similarity()`, `tpr_similarity()`, and `fscore_similarity()`
