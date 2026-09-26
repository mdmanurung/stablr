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

## Table 1. OOL full-profile run against the Python tutorial reference (negative result)

Summary of a sealed full-profile OOL run. `stablr` selected 9 features and
recovered all 7 Python tutorial selected features, giving a selected-set Jaccard
similarity of 0.778 and a maximum-feature-score Spearman correlation of 0.785.
Both fall short of the frozen acceptance thresholds of 0.80 and 0.90, and the
FDP-threshold delta of 0.08 exceeds the 0.05 threshold. The run requested a
fixed-X knockoff construction that is infeasible when features outnumber
samples, so it fell back to random-permutation artificial features and its FDP+
was calibrated against a different null distribution than the reference. The
table is reported as a negative result; it does not support a parity claim.

Source file: `papers/application-note/artifacts/stablr-0.1.1-c59780093899/ool_parity_results.csv`

## Optional supplementary table. Selected OOL features

Selected feature table from the publication-scale OOL parity run, including
feature names and release-generated scores.

Source file: `papers/application-note/artifacts/ool_publication_selected_features.csv`
