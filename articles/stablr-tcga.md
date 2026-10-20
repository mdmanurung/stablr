# Research protocol: TCGA multi-omic biomarker discovery

## Contents

- 1 Scope and status of this protocol
- 2 Prerequisites
- 3 Loading and preparing the data
- 4 Per-omic STABL
- 5 Integrated multi-omic pipeline
- 6 Validation performance
- 7 Exporting results
- 8 Comparison with mixOmics (sgccaDA and DIABLO)
- 9 Session information

## 1 Scope and status of this protocol

> **Non-executable research protocol.** No TCGA result in this page is
> package or release evidence. Code chunks are set to `eval = FALSE`
> because the data are external and the fits are long-running. Execute
> the protocol only after pinning the dataset, candidate tarball,
> settings, and output manifest.

Breast cancer subtype is a molecular label, but it is also a test of
integration. The same tumour is measured by mRNA and miRNA assays. Some
signals should be visible in one view, some in the other, and some only
when the views are allowed to compete or cooperate.

This vignette translates the **mixOmics N-Integration Chapter 6** case
study into the `stablr` idiom. The original mixOmics workflow uses
sparse multi-block Discriminant Analysis (sgccaDA) to classify breast
cancer subtypes. The question addressed here is narrower: which mRNA and
miRNA features robustly distinguish Basal-like tumours from all other
subtypes.

The proposed analysis uses the `breast.TCGA` train/test split from
`mixOmics`, with compact settings shown only to specify the workflow. No
result is reported here. A publication analysis requires a separately
predeclared and hashed run.

**Dataset:** `breast.TCGA` from the `mixOmics` package.

- Training: mRNA and miRNA measurements on the same tumour samples.
- Test: held-out mRNA and miRNA measurements with matching feature sets.
- Outcome: PAM50 subtype binarised to **Basal vs. non-Basal**.
- Protein data is excluded to match the two-omic scope used here.

## 2 Prerequisites

The `mixOmics` package is required. The following chunk stops with an
installation hint when it is unavailable.

``` r

if (!requireNamespace("mixOmics", quietly = TRUE)) {
  stop(
    "This vignette requires mixOmics (Bioconductor). ",
    "Install with: BiocManager::install('mixOmics')"
  )
}
library(stablr)
library(mixOmics)
```

## 3 Loading and preparing the data

The train and test matrices are already aligned within the `mixOmics`
object. The mRNA and miRNA matrices are extracted, and the PAM50 subtype
is turned into a binary Basal vs non-Basal outcome.

The four matrices are extracted as follows.

``` r

data(breast.TCGA)

X_train_mrna  <- breast.TCGA$data.train$mrna
X_train_mirna <- breast.TCGA$data.train$mirna
X_test_mrna   <- breast.TCGA$data.test$mrna
X_test_mirna  <- breast.TCGA$data.test$mirna

cat("Training mRNA:  ", nrow(X_train_mrna),  "x", ncol(X_train_mrna),  "\n")
cat("Training miRNA: ", nrow(X_train_mirna), "x", ncol(X_train_mirna), "\n")
cat("Test mRNA:      ", nrow(X_test_mrna),   "x", ncol(X_test_mrna),   "\n")
cat("Test miRNA:     ", nrow(X_test_mirna),  "x", ncol(X_test_mirna),  "\n")
```

The output reports the number of samples and features in each of the
four matrices.

The outcome is then binarised.

``` r

subtype_train <- breast.TCGA$data.train$subtype
subtype_test  <- breast.TCGA$data.test$subtype

y_train <- as.integer(subtype_train == "Basal")
y_test  <- as.integer(subtype_test  == "Basal")

names(y_train) <- rownames(X_train_mrna)
names(y_test)  <- rownames(X_test_mrna)

cat("Training: ", sum(y_train), "Basal /", sum(y_train == 0), "non-Basal\n")
cat("Test:     ", sum(y_test),  "Basal /", sum(y_test == 0),  "non-Basal\n")
```

The output reports the number of Basal and non-Basal tumours in the
training and test splits.

## 4 Per-omic STABL

The first pass treats each assay as its own experiment. This establishes
whether mRNA and miRNA each contain stable Basal signal before the
integrated model is allowed to mix them. The vignette uses 30 bootstraps
and 12 lambda values to keep runtime bounded; 500-1000 bootstraps should
be used for final reporting.

The mRNA fit is run first.

``` r

lambda_mrna <- auto_lambda_grid(
  X_train_mrna, y_train,
  family   = "binomial",
  n_lambda = 12
)

fit_mrna <- stabl_fit(
  x               = X_train_mrna,
  y               = y_train,
  lambda_grid     = lambda_mrna,
  family          = "binomial",
  n_bootstraps    = 30L,
  artificial_type = "random_permutation",
  random_state    = 42L
)
fit_mrna
```

The printed summary reports the selected features and the stability
diagnostics of the mRNA fit.

The miRNA fit is run with the same settings.

``` r

lambda_mirna <- auto_lambda_grid(
  X_train_mirna, y_train,
  family   = "binomial",
  n_lambda = 12
)

fit_mirna <- stabl_fit(
  x               = X_train_mirna,
  y               = y_train,
  lambda_grid     = lambda_mirna,
  family          = "binomial",
  n_bootstraps    = 30L,
  artificial_type = "random_permutation",
  random_state    = 42L
)
fit_mirna
```

The printed summary reports the corresponding results for miRNA.

The selected feature names are printed for both assays.

``` r

cat("mRNA selected feature names:\n")
print(get_feature_names(fit_mrna))

cat("\nmiRNA selected feature names:\n")
print(get_feature_names(fit_mirna))
```

The two lists show which features survive stability selection in each
assay. Features selected in both assays are candidates for cross-omic
markers, whereas features specific to one assay reflect view-specific
signal.

The stability paths are plotted for both assays.

``` r

print(plot_stabl_path(fit_mrna,  title = "mRNA - stability path"))
print(plot_stabl_path(fit_mirna, title = "miRNA - stability path"))
```

Each plot shows the stability of every feature across the lambda grid
and should be inspected for features that remain stable over a wide
range of lambda values.

## 5 Integrated multi-omic pipeline

The integrated workflow repeats the per-omic fits inside a common
pipeline, adds early fusion on the concatenated matrix, and learns a
late-fusion validation predictor. In other words, it asks three related
questions: what is stable in each view, what is stable when all features
compete, and how much each view helps prediction on held-out tumours.

The integrated fit is run with early and late fusion enabled.

``` r

lambda_list <- list(
  mrna  = lambda_mrna,
  mirna = lambda_mirna
)

multi_fit <- stabl_multiomic_train_validate(
  x_train_list    = list(mrna = X_train_mrna,  mirna = X_train_mirna),
  x_valid_list    = list(mrna = X_test_mrna,   mirna = X_test_mirna),
  y_train         = y_train,
  y_valid         = y_test,
  lambda_grid     = lambda_list,
  family          = "binomial",
  n_bootstraps    = 30L,
  artificial_type = "random_permutation",
  random_state    = 42L,
  early_fusion    = TRUE,
  late_fusion     = TRUE,
  n_iter_lf       = 500L
)
multi_fit
```

The printed object summarises the per-omic fits and the fusion branches.

The selected features and the late-fusion weights are printed as
follows.

``` r

cat("mRNA (integrated):\n")
print(get_feature_names(multi_fit$fits$mrna))

cat("\nmiRNA (integrated):\n")
print(get_feature_names(multi_fit$fits$mirna))

cat("\nEarly fusion selected feature names:\n")
print(get_feature_names(multi_fit$early_fusion$fit))

cat("\nLate fusion omic weights:\n")
print(multi_fit$late_fusion$weights)
```

The output lists the features selected in each view under the integrated
pipeline, the features selected by early fusion, and the weights that
the late fusion assigns to the two omics.

The stability paths of the integrated per-omic fits are plotted as
follows.

``` r

print(plot_stabl_path(multi_fit$fits$mrna,  title = "Integrated - mRNA"))
print(plot_stabl_path(multi_fit$fits$mirna, title = "Integrated - miRNA"))
```

These paths can be compared with the per-omic paths in Section 4 to
determine whether integration changes the stable feature sets.

## 6 Validation performance

Feature stability is only one axis of evidence. The held-out test set
asks whether the selected feature names and fusion weights carry
predictive information outside the training data.

The late-fusion predictions are evaluated on the test set with ROC and
PRC curves, a confusion matrix, and derived error rates.

``` r

lf_preds <- multi_fit$late_fusion$valid_predictions

# ROC/PRC curves computed inline (the package's plot_roc/plot_prc helpers
# were internalized in 1.0.0; any general plotting package serves here).
ord <- order(lf_preds, decreasing = TRUE)
y_ord <- y_test[ord]
n_pos <- sum(y_ord == 1); n_neg <- sum(y_ord == 0)
tpr <- cumsum(y_ord == 1) / n_pos
fpr <- cumsum(y_ord == 0) / n_neg
roc_df <- data.frame(fpr = c(0, fpr), tpr = c(0, tpr))
prc_df <- data.frame(
  recall    = c(0, tpr),
  precision = c(1, cumsum(y_ord == 1) / seq_along(y_ord))
)

ggplot2::ggplot(roc_df, ggplot2::aes(fpr, tpr)) +
  ggplot2::geom_line() +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::labs(title = "Late fusion - ROC curve (test set)",
                x = "False positive rate", y = "True positive rate")
ggplot2::ggplot(prc_df, ggplot2::aes(recall, precision)) +
  ggplot2::geom_line() +
  ggplot2::labs(title = "Late fusion - PRC curve (test set)",
                x = "Recall", y = "Precision")

pred_class <- as.integer(lf_preds >= 0.5)
conf_mat <- table(
  Predicted = factor(pred_class, levels = c(0, 1)),
  Actual    = factor(y_test, levels = c(0, 1))
)
print(conf_mat)

sensitivity <- conf_mat["1", "1"] / sum(conf_mat[, "1"])
specificity <- conf_mat["0", "0"] / sum(conf_mat[, "0"])
ber <- 1 - (sensitivity + specificity) / 2

cat(sprintf("\nSensitivity (Basal):     %.3f\n", sensitivity))
cat(sprintf("Specificity (non-Basal): %.3f\n", specificity))
cat(sprintf("Balanced Error Rate:     %.3f\n", ber))
```

The ROC and PRC curves summarise the discrimination of Basal against
non-Basal tumours. The printed sensitivity, specificity, and balanced
error rate quantify performance at the 0.5 prediction threshold.

The distributions of the selected mRNA features are plotted by outcome
class.

``` r

mrna_feats <- get_feature_names(multi_fit$fits$mrna)
if (length(mrna_feats) > 0) {
  # Per-feature class distributions (boxplot_features() was internalized
  # in 1.0.0; plain ggplot2 serves the same purpose).
  plot_df <- data.frame(
    value   = as.vector(X_train_mrna[, mrna_feats]),
    feature = rep(mrna_feats, each = nrow(X_train_mrna)),
    class   = factor(rep(y_train, times = length(mrna_feats)))
  )
  print(
    ggplot2::ggplot(plot_df, ggplot2::aes(class, value)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~ feature, scales = "free_y",
                          ncol = min(3L, length(mrna_feats))) +
      ggplot2::labs(title = "Top mRNA features - training set",
                    x = "Class", y = "Expression")
  )
}
```

The boxplots show the separation of Basal and non-Basal tumours for each
selected mRNA feature in the training set.

## 7 Exporting results

For a vignette, exports go to
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) so that rendering
does not leave generated files in the package source tree. In a real
analysis, these CSVs are the bridge to notebooks, review tables, and
downstream biological interpretation.

The stability-score tables of both integrated per-omic fits are exported
as follows.

``` r

tcga_out_dir <- file.path(tempdir(), "stablr_tcga_results")
export_stabl_to_csv(multi_fit$fits$mrna,  path = file.path(tcga_out_dir, "mrna"))
export_stabl_to_csv(multi_fit$fits$mirna, path = file.path(tcga_out_dir, "mirna"))
list.files(tcga_out_dir, recursive = TRUE)
```

The file listing confirms that the CSV files were written under the
temporary output directory.

## 8 Comparison with mixOmics (sgccaDA and DIABLO)

The mixOmics ecosystem and `stablr` both address multi-omic
classification, but they optimise **different objectives**. Treating
them as interchangeable classifiers misses the point: one searches for
sparse latent factors that separate classes; the other searches for
**named features whose selection frequency survives resampling and
passes an FDP-aware threshold**.

### 8.1 Scientific question

| Question | mixOmics (sgccaDA / DIABLO) | stablr |
|----|----|----|
| What is optimised? | Sparse loadings on latent components that maximise block-wise covariance with the outcome | Resampling stability of individual feature coefficients |
| What is the deliverable? | Components, block loadings, [`selectVar()`](https://rdrr.io/pkg/mixOmics/man/selectVar.html) lists tied to components | Explicit feature names, stability scores, FDP+ threshold |
| Integration mechanism | Joint multiblock PLS / sGCCA with a design matrix coupling blocks | Per-omic STABL, early/late fusion, or cooperative fusion (`rho`) |
| FDR / stability calibration | Sparsity via `keepX` tuning; no bootstrap stability or knockoff FDP+ | Knockoffs or random artificial features; FDP+ threshold sweep |
| Prediction | Native (`predict`, `centroids.dist`) | Requires a downstream refit on selected feature names (e.g. multinomial glmnet) |

### 8.2 Objective strengths and limitations

**mixOmics DIABLO is preferred when:**

- The goal is **integrative discrimination** with interpretable latent
  structure across correlated blocks.
- A mature tuning workflow (`tune.block.splsda`, design-matrix priors)
  and rich visualisation (`plotIndiv`, `plotVar`, circos plots) are
  desired.
- Classification (or block-supervised reduction) is the primary endpoint
  and component loadings are an acceptable feature representation.

**stablr is preferred when:**

- The goal is a **reproducible sparse biomarker list** with an FDP+
  diagnostic under the artificial-feature assumptions. The fitted
  threshold minimises the diagnostic curve; it is not a universal
  false-discovery guarantee.
- Gaussian, binomial, multinomial, or Cox STABL are needed on the same
  API, with grouped bootstrap for repeated measures.
- An sPLS-DA-style classification learner is wanted inside the
  stability-selection machinery (`base_learner = "splsda"`; see
  `base_learner = "splsda"`).
- Integration strategies (per-omic, early, late, cooperative) are to be
  compared under one stability-selection framework.
- The glmnet / resampling idiom is preferred over PLS components.

**Neither replaces the other.** On TCGA, DIABLO may win on outer-fold
BER when blocks share strong latent structure; stablr may surface
features with higher cross-fold recurrence that are not top DIABLO
loadings. Quantitative head-to-head results live in
[`vignette("stablr-tcga-nestedcv")`](https://mdmanurung.github.io/stablr/articles/stablr-tcga-nestedcv.md);
the full benchmark protocol is documented in
`inst/analysis/BENCHMARK_PLAN.md`.

| Aspect | sgccaDA (Ch. 6) | DIABLO | stablr |
|----|----|----|----|
| Method | Sparse CCA + LDA | `block.splsda` + design matrix | Bootstrap-stability glmnet |
| Feature selection | L1 on components | `keepX` per block | FDP+ stability threshold |
| Integration | Joint latent space | Correlated latent components | Early, late, or cooperative fusion |
| Primary readout | BER + loadings | BER + `selectVar` | Recurrence + stability scores |
| Survival / Cox | Other mixOmics modules | Classification focus | Native `family = "cox"` |

Both sgccaDA and DIABLO target compact mRNA and miRNA signatures for
subtype classification. STABL’s distinctive contribution is FDP-aware
stability selection without assuming that omics must share a single
latent space.

## 9 Session information

The session information is recorded for reproducibility.

``` r

sessionInfo()
```
