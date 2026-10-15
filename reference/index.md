# Package index

## Package Overview

- [`stablr`](https://mdmanurung.github.io/stablr/reference/stablr-package.md)
  [`stablr-package`](https://mdmanurung.github.io/stablr/reference/stablr-package.md)
  : stablr: Sparse and Reliable Biomarker Discovery in R

## Getting started

The one-stop-shop entry point for new users — a data.frame in, a
plain-language selection out.

- [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md) :
  STABL Stability Selection (Beginner Interface)
- [`print(`*`<stabl>`*`)`](https://mdmanurung.github.io/stablr/reference/print.stabl.md)
  : Print a Fitted STABL Object (Beginner Style)
- [`summary(`*`<stabl>`*`)`](https://mdmanurung.github.io/stablr/reference/summary.stabl.md)
  : Summarize a Fitted STABL Object (Beginner Style)
- [`plot(`*`<stabl>`*`)`](https://mdmanurung.github.io/stablr/reference/plot.stabl.md)
  : Plot Stability Paths for a Fitted STABL Object (Beginner Style)
- [`support()`](https://mdmanurung.github.io/stablr/reference/support.md)
  [`feature_names()`](https://mdmanurung.github.io/stablr/reference/support.md)
  : Support and Feature Names of a STABL Selection

## Core STABL Engine

Fit single-matrix STABL selectors and compute FDP+ thresholds (expert
interface).

- [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  [`print(`*`<stabl_fit>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  : Fit STABL (Stability-Penalized Feature Selection)
- [`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md)
  : Build a Data-Driven Lambda Grid
- [`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md)
  : Compute FDP+ (False Discovery Proportion Upper Bound)

## Knockoff Feature Construction

Factories that inject artificial features (Model-X knockoffs or random
permutations) for FDP+ calibration; plug into stabl_fit through the
knockoff/solver arguments.

- [`make_knockoff_features()`](https://mdmanurung.github.io/stablr/reference/make_knockoff_features.md)
  : Make Knockoff Artificial Features
- [`make_knockoff_equi_features()`](https://mdmanurung.github.io/stablr/reference/make_knockoff_equi_features.md)
  : Make Model-X Equicorrelated Knockoff Artificial Features
- [`make_knockoff_mvr_features()`](https://mdmanurung.github.io/stablr/reference/make_knockoff_mvr_features.md)
  : Make Model-X MVR Knockoff Artificial Features
- [`make_rp_features()`](https://mdmanurung.github.io/stablr/reference/make_rp_features.md)
  : Make Random-Permutation Artificial Features

## Learner Adapters

Public learner factories for glmnet-family and sparse-group backends.

- [`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md)
  : Build a glmnet Learner Adapter
- [`make_adaptive_lasso_adapter()`](https://mdmanurung.github.io/stablr/reference/make_adaptive_lasso_adapter.md)
  : Build an Adaptive Lasso Learner Adapter
- [`make_sgl_adapter()`](https://mdmanurung.github.io/stablr/reference/make_sgl_adapter.md)
  : Build a Sparse Group Lasso Learner Adapter

## Extended Learner Adapters

Extensions beyond the glmnet family. The random forest is calibrated
within a measured scope (WS-2); the oblique forest is fenced (WS-3); the
ordinal learner is experimental (WS-9).

- [`make_ranger_adapter()`](https://mdmanurung.github.io/stablr/reference/make_ranger_adapter.md)
  : Build a Random-Forest Learner Adapter (Calibrated, Scoped)
- [`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md)
  : Auto-Generate an Importance-Threshold Grid for Random-Forest STABL
  (Experimental)
- [`make_aorsf_adapter()`](https://mdmanurung.github.io/stablr/reference/make_aorsf_adapter.md)
  : Build an Oblique Random-Forest Learner Adapter (Fenced)
- [`auto_oblique_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_oblique_threshold_grid.md)
  : Auto-Generate an Importance-Threshold Grid for Oblique-Random-Forest
  STABL (Experimental)
- [`make_ordinal_adapter()`](https://mdmanurung.github.io/stablr/reference/make_ordinal_adapter.md)
  : Build an Ordinal-Regression Learner Adapter (Experimental)
- [`auto_ordinal_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_ordinal_lambda_grid.md)
  : Build a Data-Driven Lambda Grid for the Ordinal Learner
  (Experimental)

## Environment Stability

ICP-inspired transportability workflow across environments (batches,
cohorts, sites); promoted with a documented scope after the WS-5
transportability benchmark (min-aggregation is the calibrated read).

- [`stabl_environment_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_environment_stability.md)
  : Environment-Stability Workflow: Reliable and Transportable Feature
  Sets
- [`get_min_scores()`](https://mdmanurung.github.io/stablr/reference/get_min_scores.md)
  : Get Minimum-Across-Environment Stability Scores
- [`print(`*`<stabl_env_stability>`*`)`](https://mdmanurung.github.io/stablr/reference/print.stabl_env_stability.md)
  : Print an Environment-Stability Analysis

## Artificial Features

Decoy-feature generation for FDP+ calibration.

- [`make_artificial_features()`](https://mdmanurung.github.io/stablr/reference/make_artificial_features.md)
  : Dispatcher for Artificial Feature Generation

## Input Validation and Bootstrapping

- [`validate_sample_alignment()`](https://mdmanurung.github.io/stablr/reference/validate_sample_alignment.md)
  : Validate Sample Alignment Across Inputs
- [`validate_multiomic_inputs()`](https://mdmanurung.github.io/stablr/reference/validate_multiomic_inputs.md)
  : Validate Multi-Omic Input Contract
- [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md)
  : Classic Bootstrap Sampler Indices
- [`group_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/group_bootstrap_indices.md)
  : Group-Aware Bootstrap Sampler Indices

## Accessors

- [`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
  : Get the Feature Selection Mask from a Fitted STABL Object
- [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md)
  : Get the Names of Selected Features from a Fitted STABL Object
- [`get_fdrs()`](https://mdmanurung.github.io/stablr/reference/get_fdrs.md)
  : Get FDP+ Diagnostic Estimates Across the Swept Thresholds
- [`transform_stabl()`](https://mdmanurung.github.io/stablr/reference/transform_stabl.md)
  : Transform New Data to the Selected STABL Feature Set
- [`get_stabl_scores()`](https://mdmanurung.github.io/stablr/reference/get_stabl_scores.md)
  : Get the Full Stability Score Matrix from a Fitted STABL Object
- [`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md)
  : Get Per-Feature Importance Scores (Maximum Stability Score)
- [`get_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/get_lambda_grid.md)
  : Get the Fitted Lambda Grid
- [`get_cooperative_features()`](https://mdmanurung.github.io/stablr/reference/get_cooperative_features.md)
  : Get Cooperative-Fusion Selected Features
- [`get_cooperative_diagnostics()`](https://mdmanurung.github.io/stablr/reference/get_cooperative_diagnostics.md)
  : Get Cooperative-Fusion Tuning Diagnostics

## Multi-Omic Workflows

- [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  [`print(`*`<stabl_multiomic_fit>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  : Multi-Omic STABL Train/Validation Workflow
- [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  [`print(`*`<stabl_multiomic_cv>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  : Multi-Omic STABL Cross-Validation Workflow
- [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
  [`print(`*`<stabl_multiomic_nested_cv>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
  : Multi-Omic STABL Nested Cross-Validation
- [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)
  : Stacked Generalization Over Per-Omic Predictions
- [`stabl_multiomic_study_diagnostics()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_study_diagnostics.md)
  : Per-Study Batch Diagnostics for Multi-Omic STABL Inputs
- [`stabl_multiomic_study_split()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_study_split.md)
  : Study-Disjoint Train/Validation Split
- [`stabl_multiomic_scheme_comparison()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_scheme_comparison.md)
  : Compare Sample-I.I.D. and Study-Disjoint Nested CV
- [`stabl_multiomic_selection_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_selection_stability.md)
  : Selection Stability Across Study-Disjoint Folds
- [`load_ool_data()`](https://mdmanurung.github.io/stablr/reference/load_ool_data.md)
  : Load the Onset of Labor example dataset

## Analysis Workflow

Explicit exploratory/validated analysis with canonical evidence tables
and fold-safe preprocessing.

- [`stabl_analyze()`](https://mdmanurung.github.io/stablr/reference/stabl_analyze.md)
  [`print(`*`<stabl_analysis>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_analyze.md)
  [`summary(`*`<stabl_analysis>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_analyze.md)
  : STABL analysis with an explicit evidence mode
- [`stabl_interpret()`](https://mdmanurung.github.io/stablr/reference/stabl_interpret.md)
  : Interpret a fitted STABL selection
- [`stabl_preprocess_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit.md)
  : Fit a fold-safe preprocessing transform on training data
- [`stabl_preprocess_apply()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_apply.md)
  : Apply a fitted preprocessing transform
- [`stabl_preprocess_fit_recipes()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit_recipes.md)
  : Fit a preprocessing transform from a recipes recipe

## Visualization

- [`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md)
  : Plot the STABL Stability Path
- [`plot_fdr_graph()`](https://mdmanurung.github.io/stablr/reference/plot_fdr_graph.md)
  : Plot the FDP+ FDR Estimate Curve
- [`plot_roc()`](https://mdmanurung.github.io/stablr/reference/plot_roc.md)
  : Plot a ROC Curve
- [`plot_prc()`](https://mdmanurung.github.io/stablr/reference/plot_prc.md)
  : Plot a Precision-Recall Curve
- [`boxplot_features()`](https://mdmanurung.github.io/stablr/reference/boxplot_features.md)
  : Boxplots of Selected Features Grouped by Outcome
- [`scatterplot_features()`](https://mdmanurung.github.io/stablr/reference/scatterplot_features.md)
  : Scatterplots of Selected Features Against a Continuous Outcome

## Plot a fitted stabl_fit

S3 plot method for fitted STABL objects (stability path and FDP+ graph).

- [`plot(`*`<stabl_fit>`*`)`](https://mdmanurung.github.io/stablr/reference/plot.stabl_fit.md)
  : Plot a Fitted stabl_fit Object

## Summarise a fitted stabl_fit

S3 summary and print methods for fitted STABL objects.

- [`summary(`*`<stabl_fit>`*`)`](https://mdmanurung.github.io/stablr/reference/summary.stabl_fit.md)
  [`print(`*`<summary.stabl_fit>`*`)`](https://mdmanurung.github.io/stablr/reference/summary.stabl_fit.md)
  : Summarize a Fitted stabl_fit Object

## Export

- [`export_stabl_to_csv()`](https://mdmanurung.github.io/stablr/reference/export_stabl_to_csv.md)
  : Export STABL Stability Scores to CSV
- [`save_stabl_results()`](https://mdmanurung.github.io/stablr/reference/save_stabl_results.md)
  : Save All STABL Results to Disk

## Selection Similarity Metrics

- [`jaccard_similarity()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity.md)
  : Jaccard Similarity Between Two Feature Sets
- [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)
  : Pairwise Jaccard Matrix from a List of Feature Sets
- [`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md)
  : Adjusted Similarity Between Two Feature Sets
- [`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md)
  : Upper-Triangle Adjusted Similarity Values
- [`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
  : Summary Statistic of Adjusted Similarity Values
- [`pearson_similarity()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity.md)
  : Pearson-Corrected Similarity Between Two Feature Sets
- [`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)
  : Upper-Triangle Pearson Similarity Values
- [`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)
  : Summary Statistic of Pearson Similarity Values
- [`fdr_similarity()`](https://mdmanurung.github.io/stablr/reference/fdr_similarity.md)
  : FDR Between Two Feature Sets
- [`tpr_similarity()`](https://mdmanurung.github.io/stablr/reference/tpr_similarity.md)
  : TPR Between Two Feature Sets
- [`fscore_similarity()`](https://mdmanurung.github.io/stablr/reference/fscore_similarity.md)
  : F-Score Between Two Feature Sets

## Deprecated names

Old names kept as deprecation shims; they warn once and delegate to the
replacements above.

- [`jaccard_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_matrix.md)
  :

  Deprecated: use
  [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)

- [`adjusted_similarity_values()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_values.md)
  :

  Deprecated: use
  [`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md)

- [`adjusted_similarity_measure()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_measure.md)
  :

  Deprecated: use
  [`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)

- [`pearson_similarity_values()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_values.md)
  :

  Deprecated: use
  [`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md)

- [`pearson_similarity_measure()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_measure.md)
  :

  Deprecated: use
  [`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)

- [`get_feature_names_out()`](https://mdmanurung.github.io/stablr/reference/get_feature_names_out.md)
  :

  Deprecated: use
  [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md)

- [`classic_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/classic_bootstrap_indices.md)
  :

  Deprecated: use
  [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md)

- [`stacked_multi_omic()`](https://mdmanurung.github.io/stablr/reference/stacked_multi_omic.md)
  :

  Deprecated: use
  [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)
