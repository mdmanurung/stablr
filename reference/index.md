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

## Learner Grid Helpers

Data-driven threshold/lambda grid constructors for the non-glmnet base
learners. The random forest is calibrated within a measured scope
(WS-2); the ordinal learner is experimental (WS-9); sPLS-DA is
experimental and binomial-only. Learner factories are internal; select
learners via the `base_learner` argument of
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- [`auto_rf_threshold_grid()`](https://mdmanurung.github.io/stablr/reference/auto_rf_threshold_grid.md)
  : Auto-Generate an Importance-Threshold Grid for Random-Forest STABL
  (Experimental)
- [`auto_ordinal_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_ordinal_lambda_grid.md)
  : Build a Data-Driven Lambda Grid for the Ordinal Learner
  (Experimental)
- [`auto_splsda_keepx_grid()`](https://mdmanurung.github.io/stablr/reference/auto_splsda_keepx_grid.md)
  : Build an Automatic keepX Grid for the sPLS-DA Learner

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

## Reliability and Final Models

Post-selection reliability diagnostics, outcome-null validation, cohort
replication, and portable final models (Phase 13; see the
stablr-reliability article for a guided tour).

- [`stabl_diagnose()`](https://mdmanurung.github.io/stablr/reference/stabl_diagnose.md)
  : Diagnose selection reliability of a fitted STABL object

- [`print(`*`<stabl_diagnostics>`*`)`](https://mdmanurung.github.io/stablr/reference/print.stabl_diagnostics.md)
  : Print STABL Selection Diagnostics

- [`plot(`*`<stabl_diagnostics>`*`)`](https://mdmanurung.github.io/stablr/reference/plot.stabl_diagnostics.md)
  : Plot STABL selection diagnostics

- [`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md)
  : Outcome-Null Calibration Check for STABL

- [`print(`*`<stabl_null_validation>`*`)`](https://mdmanurung.github.io/stablr/reference/print.stabl_null_validation.md)
  : Print an Outcome-Null Validation Summary

- [`plot(`*`<stabl_null_validation>`*`)`](https://mdmanurung.github.io/stablr/reference/plot.stabl_null_validation.md)
  : Plot Outcome-Null Validation Results

- [`stabl_replication_profile()`](https://mdmanurung.github.io/stablr/reference/stabl_replication_profile.md)
  : Cohort Replication and Heterogeneity Profile

- [`print(`*`<stabl_replication_profile>`*`)`](https://mdmanurung.github.io/stablr/reference/print.stabl_replication_profile.md)
  : Print a Cohort Replication Profile

- [`stabl_finalize()`](https://mdmanurung.github.io/stablr/reference/stabl_finalize.md)
  [`predict(`*`<stabl_predictor>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_finalize.md)
  [`print(`*`<stabl_predictor>`*`)`](https://mdmanurung.github.io/stablr/reference/stabl_finalize.md)
  : Finalize a STABL selection into a portable predictive model

- [`validate_newdata_schema()`](https://mdmanurung.github.io/stablr/reference/validate_newdata_schema.md)
  :

  Validate new data against a `"stabl_predictor"` schema

## Visualization

- [`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md)
  : Plot the STABL Stability Path
- [`plot_fdr_graph()`](https://mdmanurung.github.io/stablr/reference/plot_fdr_graph.md)
  : Plot the FDP+ Diagnostic Curve

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
