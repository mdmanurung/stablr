# stablr: Sparse and Reliable Biomarker Discovery in R

`stablr` is an R-native implementation of STABL for sparse, stable
biomarker selection in high-dimensional clinical and omic datasets. The
package compiles C++ for its cooperative-learning backend and provides a
core bootstrap stability-selection engine, FDP+ threshold calibration
with artificial features, glmnet-family learner adapters, and multi-omic
workflows for per-omic, early-fusion, late-fusion, and
cooperative-fusion analyses.

## Main workflows

- [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  fits the core single-matrix STABL selector.

- [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  runs train/validation multi-omic workflows with optional early, late,
  and cooperative fusion.

- [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  runs outer cross-validation for named multi-omic inputs, preserving
  optional group structure.

## Learners and outcomes

The core selector supports lasso, elastic net, adaptive lasso, and
optional sparse group lasso backends. Supported `glmnet` families
include Gaussian, binomial, multinomial, and Cox where the selected
backend supports them.

## Outputs

Fitted objects expose stable S3 accessors for support masks, selected
feature names, stability scores, importances, cooperative-fusion
features, and cooperative diagnostics. Plotting and export helpers
provide stability paths, FDP+ curves, ROC/PRC plots, selected-feature
visualizations, and CSV output.

## See also

[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md),
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md),
[`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md),
[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md)

## Author

**Maintainer**: Mikaël Manurung <mikhael.manurung@gmail.com> \[copyright
holder\]

Authors:

- Mikaël Manurung <mikhael.manurung@gmail.com> \[copyright holder\]

Other contributors:

- Jean Hédou (Original STABL method (Hédou et al., Nat Biotechnol 2024))
  \[contributor\]

- Ivana Marić (Original STABL method (Hédou et al., Nat Biotechnol
  2024)) \[contributor\]

- Gregory Bellan (Original STABL method (Hédou et al., Nat Biotechnol
  2024)) \[contributor\]

- David Gfeller (Original STABL method (Hédou et al., Nat Biotechnol
  2024)) \[contributor\]

- Daisy Yi Ding (Cooperative learning / multiview method (Ding et al.,
  PNAS 2022)) \[contributor\]

- Shuangning Li (Cooperative learning / multiview method (Ding et al.,
  PNAS 2022)) \[contributor\]

- Balasubramanian Narasimhan (Cooperative learning / multiview method
  (Ding et al., PNAS 2022)) \[contributor\]

- Robert J. Tibshirani (Cooperative learning / multiview method (Ding et
  al., PNAS 2022)) \[contributor\]

- Trevor Hastie \[contributor\]

- Kenneth Tay \[contributor\]

- James Yang \[contributor\]

- Jonathan Taylor \[contributor\]
