# Migrating from stablr 0.3.x to 1.0.0

stablr 1.0.0 is a **breaking** major release. It removes the deprecated
compatibility shims introduced at the 0.3.0 rename, deletes duplicate
beginner generics, internalizes learner-factory and utility-plotting
helpers that were never the extension contract used by `stabl_fit()`, and
retires the experimental oblique random-forest learner. Every removal has a
mechanical replacement; nothing about the statistical contracts of the
retained functions changes.

Existing release tags (<= 0.3.x) remain available for faithful historical
execution of archived scripts.

## 1. Removed deprecated aliases (call the canonical name)

| Removed (errors in 1.0.0) | Use instead |
|---|---|
| `classic_bootstrap_indices(...)` | `sample_bootstrap_indices(...)` |
| `get_feature_names_out(object, ...)` | `get_feature_names(object, ...)` |
| `jaccard_matrix(list_of_lists, ...)` | `jaccard_similarity_matrix(list_of_lists, ...)` |
| `adjusted_similarity_values(list_of_lists, n)` | `adjusted_similarity_pairwise(list_of_lists, n)` |
| `adjusted_similarity_measure(list_of_lists, n, stat)` | `adjusted_similarity_summary(list_of_lists, n, stat)` |
| `pearson_similarity_values(list_of_lists, d)` | `pearson_similarity_pairwise(list_of_lists, d)` |
| `pearson_similarity_measure(list_of_lists, d, stat)` | `pearson_similarity_summary(list_of_lists, d, stat)` |
| `stacked_multi_omic(...)` | `stabl_multiomic_stacked(...)` |

The Python-parity fixtures keep their historical file and key names
(`jaccard_matrix_remove_diag_rowmajor`, ...); only the R entry points moved.

## 2. Removed duplicate generics

| Removed | Use instead |
|---|---|
| `support(fit)` | `get_support(fit)` |
| `feature_names(fit)` | `get_feature_names(fit)` |

`stabl()` remains the beginner front door; only the duplicate S3 generic
layer was deleted.

## 3. Internalized learner factories (select learners by name)

The single-row adapter factories were never the dispatch contract used by
`stabl_fit()` (which uses internal batch adapters). They are no longer
exported:

`make_glmnet_adapter()`, `make_adaptive_lasso_adapter()`,
`make_sgl_adapter()`, `make_ranger_adapter()`, `make_ordinal_adapter()`,
`make_splsda_adapter()`.

**Migration:** pass the learner by name — `stabl_fit(x, y, lambda_grid =
..., base_learner = "rf")` (one of `"lasso"`, `"elastic_net"`,
`"adaptive_lasso"`, `"sparse_group_lasso"`, `"rf"`, `"ordinal"`,
`"splsda"`). The grid helpers `auto_lambda_grid()`,
`auto_rf_threshold_grid()`, `auto_ordinal_lambda_grid()`, and
`auto_splsda_keepx_grid()` remain exported. The factory implementations
remain in the package namespace for the test suite; a validated
`stabl_register_learner()`-style extension contract may return in a future
version.

## 4. Retired: oblique random forest

`base_learner = "oblique_rf"`, `make_aorsf_adapter()`, and
`auto_oblique_threshold_grid()` are removed. The oblique forest failed the
predeclared FDP+ calibration bar in the WS-3 benchmark (decoy importance is
inflated for oblique splits, so FDP+ is anti-conservative). `stabl_fit()`
now fails closed with an explanatory error. Use `base_learner = "rf"`. The
archived source and tests live in `attic/` for inspection; re-opening
requires predeclared positive new calibration evidence.

## 5. Internalized plotting and benchmark-metric utilities

| Internalized | Replacement |
|---|---|
| `plot_roc()`, `plot_prc()` | Any general plotting package (e.g. ggplot2); the vignettes show inline ROC/PRC construction |
| `boxplot_features()`, `scatterplot_features()` | Plain ggplot2 (`geom_boxplot()`/`geom_point()` + `facet_wrap()`) |
| `fdr_similarity()`, `tpr_similarity()`, `fscore_similarity()` | Internal benchmark metrics (realized proportions against known truth, not general similarities); still used by the shipped reproducibility script `inst/analysis/run_methodology_validation.R` via `stablr:::` |

The implementations remain in the namespace so published simulations stay
reproducible; they are simply no longer part of the public API.

## 6. Behavior removals (Phase 13 PR-A, restated)

- `stabl_fit(explore = ..., n_explore = ...)` removed; an empty selection
  stays empty. Use `get_importances()` for a descriptive ranking.
- `late_fusion_training = "python_legacy"` removed from production
  workflows (in-sample stacking biased evaluation); only leakage-safe
  `"oof"` stacking remains.
- `cooperation_selection = "validation"` removed; only `"cv"` (inner-fold
  tuning) is supported.
- `plot_fdr_graph(fdr_target = ...)` renamed to `reference_fdp_plus`
  (purely visual reference line; the FDP+ curve is a diagnostic, not a
  controlled FDR).

## 7. New in 1.0.0 (migration-relevant additions)

- `stabl_diagnose()`, `stabl_null_check()`, `stabl_replication_profile()` —
  reliability diagnostics, outcome-null validation, and cohort replication.
- Validated analysis supports `family = "gaussian"` end to end
  (MAE/RMSE/R² on held-out folds).
- `stabl_finalize()` + `predict()` produce a portable `stabl_predictor`
  with new-data schema validation (`validate_newdata_schema()`).

A machine-readable record of every removed symbol, its replacement, and its
reverse-dependency scan lives in `inst/release/removal-ledger.csv`.
