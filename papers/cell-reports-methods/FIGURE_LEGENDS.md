# Figure and table legends

## Figure 1. `stablr` workflow for stability-based biomarker discovery

Schematic of the package workflow. A feature matrix and outcome are supplied to
`stabl_fit()` with a regularization grid, artificial-feature strategy, bootstrap
count, and random seed. Sparse models are fit across bootstrap samples, real and
artificial feature recurrence scores are compared through FDP+ calibration, and
selected features are returned through R-native accessors and export helpers.

Source file: `papers/application-note/artifacts/figure1_workflow_schematic.pdf`

## Figure 2. OOL proteomics parity validation

Onset of Labor proteomics validation at publication-scale settings. The release
analysis used 150 samples, 1,317 protein features, gaussian lasso models,
knockoff artificial features, 500 bootstraps, and 10 lambda values. The panel
combines the FDP+ curve and stability path used to inspect selected features.

Source file: `papers/application-note/artifacts/figure2_ool_validation.pdf`

## Table 1. Selection-level parity against the Python tutorial reference

Summary of OOL parity metrics. `stablr` selected 9 features and recovered all 7
Python tutorial selected features. The selected-set Jaccard similarity was
0.778, and Spearman correlation of maximum feature scores was 0.808.

Source file: `papers/application-note/artifacts/table1_ool_parity_metrics.csv`

## Optional supplementary table. Selected OOL features

Selected feature table from the publication-scale OOL parity run, including
feature names and release-generated scores.

Source file: `papers/application-note/artifacts/ool_publication_selected_features.csv`
