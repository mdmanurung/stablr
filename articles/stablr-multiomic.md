# Multi-Omic Biomarker Discovery with stablr

> **Workflow tutorial, not publication evidence.** This executable
> vignette uses 150 training samples, 21 aligned validation samples, and
> the first 100 features from each of the CyTOF and Proteomics assays.
> Compact lambda grids, small bootstrap counts, and random-permutation
> decoys keep package builds bounded and catch API drift without
> optional knockoff dependencies.

## Overview

A multi-omic dataset is a set of partial views. One assay sees
immune-cell state, another sees circulating proteins, and the phenotype
is the same person seen through both. The statistical question is not
only “which features are stable?” but also “does a feature remain useful
when another view is allowed to speak?”

This vignette uses a real Onset of Labor (OOL) dataset bundled with
`stablr`. It contains two omic layers, CyTOF and Proteomics, measured in
pregnant women. The outcome is days-to-onset-of-labor (DOS), a
continuous regression target. The settings are intentionally bounded for
package builds while preserving the shape of a real train/validation
analysis.

We will:

1.  Load the bundled example data with
    [`load_ool_data()`](https://mdmanurung.github.io/stablr/reference/load_ool_data.md).
2.  Run per-omic STABL — a useful baseline before integration.
3.  Run the integrated multi-omic workflow with **early fusion** and
    **late fusion**.
4.  Compare selected feature names across strategies.
5.  Export results.

## Learning objectives

By the end, you should be able to:

1.  Understand the data structure returned by
    [`load_ool_data()`](https://mdmanurung.github.io/stablr/reference/load_ool_data.md).
2.  Run STABL independently on each omic and interpret the outputs.
3.  Explain the difference between early fusion and late fusion.
4.  Run
    [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
    and navigate its result list.
5.  Compare which features appear across integration strategies.
6.  Export results for downstream use.

``` r

library(stablr)
```

------------------------------------------------------------------------

## 1. Load the data

[`load_ool_data()`](https://mdmanurung.github.io/stablr/reference/load_ool_data.md)
reads the bundled OOL subset directly from the package. No file paths
are needed. It returns a named list with two elements:

- **`$x_list`**: a named list of feature matrices, one per omic layer.
  Each matrix has samples in rows and features in columns, with matching
  row names.
- **`$y`**: a named numeric vector of outcomes (DOS in days), aligned to
  the rows of every matrix in `$x_list`.

**Relation to the Python tutorial notebook:** this is a bounded package
example, not a full reproduction of
`Notebook examples/Tutorial Notebook.ipynb`. The notebook uses the full
OOL-CV training files (`CyTOF`, `Proteomics`, and `Metabolomics`) and
the COVID-19 proteomics classification dataset. It also applies a
preprocessing pipeline before STABL: zero-variance filtering,
high-missingness feature filtering, median imputation, and
standardisation. Here,
[`load_ool_data()`](https://mdmanurung.github.io/stablr/reference/load_ool_data.md)
ships 150 training samples and 21 aligned validation samples for `CyTOF`
and `Proteomics` only, with the first 100 features from each assay,
already aligned and without missing values. These are quick-build
settings, not the predeclared full OOL parity analysis. Use
`docs/PYTHON_TO_R_MAPPING.md` when the goal is to mirror the tutorial
notebook’s datasets and preprocessing choices.

``` r

ool_train <- load_ool_data(split = "train")
ool_valid <- load_ool_data(split = "valid")

# Training set dimensions
lapply(ool_train$x_list, dim)
#> $cytof
#> [1] 150 100
#> 
#> $proteomics
#> [1] 150 100
length(ool_train$y)
#> [1] 150

# Validation set dimensions
lapply(ool_valid$x_list, dim)
#> $cytof
#> [1]  21 100
#> 
#> $proteomics
#> [1]  21 100
```

The `lapply(..., dim)` output confirms the basic geometry: the omics
share rows but differ in columns. Row names must match across views.
`stablr` checks this alignment and stops early if it detects a mismatch.

``` r

# Outcome distribution (DOS = days before onset of labor)
summary(ool_train$y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -113.00  -70.75  -35.50  -43.39  -15.00    0.00
```

The [`summary()`](https://rdrr.io/r/base/summary.html) output shows the
spread of DOS. A wide range is expected; a strong skew or narrow range
changes how much prediction signal one should expect from molecular
measurements.

------------------------------------------------------------------------

## 2. Per-omic STABL

Before integrating views, it is worth listening to each one alone. A
per-omic fit is both a sanity check and a baseline:

1.  **Sanity check** — does each omic contain any signal at all? If
    STABL selects zero features from an omic, that layer may not
    contribute useful information to the integrated analysis.
2.  **Baseline** — the per-omic biomarker lists can be compared to the
    integrated results to see which features survive the joint analysis.

We use compact 6-point lambda grids and small bootstrap counts
throughout this vignette. Increase to 500-1000 bootstraps and broader
grids for publication-quality results.

When each fit prints, look first at the number of selected feature names
and the FDP estimate. The ideal result is not simply a long list. It is
a compact list whose stability is difficult for decoys to explain.

``` r

lambda_cytof <- auto_lambda_grid(
  ool_train$x_list$cytof, ool_train$y,
  family = "gaussian", n_lambda = 6
)
```

``` r

fit_cytof <- stabl_fit(
  x               = ool_train$x_list$cytof,
  y               = ool_train$y,
  lambda_grid     = lambda_cytof,
  family          = "gaussian",
  n_bootstraps    = 10L,
  artificial_type = "random_permutation",
  random_state    = 42L
)
fit_cytof
#> <stabl_fit>
#>   Features in:      100 
#>   Features selected: 0 
#>   Min FDP+:         1.044 
#>   FDP+ threshold:   1 
#>   Artificial:       random_permutation
```

``` r

proteomics_train_scaled <- scale(ool_train$x_list$proteomics)
lambda_prot <- auto_lambda_grid(
  proteomics_train_scaled, ool_train$y,
  family = "gaussian", n_lambda = 6
)

fit_prot <- stabl_fit(
  x               = proteomics_train_scaled,
  y               = ool_train$y,
  lambda_grid     = lambda_prot,
  family          = "gaussian",
  n_bootstraps    = 10L,
  artificial_type = "random_permutation",
  random_state    = 42L
)
fit_prot
#> <stabl_fit>
#>   Features in:      100 
#>   Features selected: 0 
#>   Min FDP+:         1.01 
#>   FDP+ threshold:   1 
#>   Artificial:       random_permutation
```

[`get_feature_names_out()`](https://mdmanurung.github.io/stablr/reference/get_feature_names_out.md)
returns the feature names that passed the stability threshold in each
omic. Use
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
when you need the named logical mask for subsetting. Non-overlapping
lists suggest view-specific information. Overlap, when the same molecule
is measured in both layers, is stronger evidence because the signal has
appeared twice through different measurement systems.

``` r

cat("CyTOF selected feature names:\n")
#> CyTOF selected feature names:
print(get_feature_names_out(fit_cytof))
#> Warning: `get_feature_names_out()` is deprecated; use `get_feature_names()`
#> instead.
#> character(0)

cat("\nProteomics selected feature names:\n")
#> 
#> Proteomics selected feature names:
print(get_feature_names_out(fit_prot))
#> Warning: `get_feature_names_out()` is deprecated; use `get_feature_names()`
#> instead.
#> character(0)
```

In the stability path, each line is one feature. Features that hold high
selection frequencies across a moderate penalty range are the most
convincing. Features that appear only when the penalty is very weak are
easier to dismiss as opportunistic.

``` r

plot_stabl_path(fit_cytof,     title = "CyTOF - stability path")
```

![Stability paths for CyTOF and
Proteomics](stablr-multiomic_files/figure-html/plot-per-omic-1.png)

``` r

plot_stabl_path(fit_prot,      title = "Proteomics - stability path")
```

![Stability paths for CyTOF and
Proteomics](stablr-multiomic_files/figure-html/plot-per-omic-2.png)

------------------------------------------------------------------------

## 3. Integrated multi-omic workflow

[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
puts the train/validation analysis into one call. It accepts matched
omic matrices and evaluates two complementary ways to let the views
interact:

| Strategy | How it works | Best when |
|----|----|----|
| **Early fusion** | Concatenates all omic matrices into one wide matrix, then runs a *single* STABL fit across all features jointly. | Omics are measured on the same samples; you want a unified biomarker list. |
| **Late fusion** | Fits one STABL + downstream model per omic separately, learns stacking weights from OOF predictions by default, and applies the weighted blend to full-refit validation predictions. | You want a held-out evaluation of an omic-weighted prediction rule. |

Key arguments:

| Argument | Purpose |
|----|----|
| `x_train_list` / `x_valid_list` | Named lists of feature matrices; names must match between train and valid. |
| `lambda_grid` | Named list of lambda grids, one per omic (names must match `x_train_list`). |
| `early_fusion` / `late_fusion` | Toggle each strategy on or off independently. |
| `late_fusion_training` | Use leakage-safe `"oof"` (default), or explicit `"python_legacy"` compatibility mode. |
| `late_fusion_nfolds` | Number of stacking folds; five folds plus the full refit cost roughly six complete fits. |
| `n_iter_lf` | Number of optimisation iterations for the late-fusion weight search. |

The printed `multi_fit` object summarises every active branch: per-omic
selected-feature counts, whether validation data are present,
early-fusion branch presence, and the late-fusion score. In OOF mode,
assessment samples and their groups never enter selector or
downstream-model fitting; use `$oof_predictions` and `$oof_score`.
Explicit `"python_legacy"` mode instead exposes `$in_sample_predictions`
and `$in_sample_score`. The older `$train_predictions` and `$score`
fields are exact 0.1.x compatibility aliases for the active mode. Binary
and regression late fusion use parity-preserving batched stacking;
multiclass stacking remains scalar by design to preserve
probability-normalisation and tie behavior.

``` r

# Build per-omic lambda grids
lambda_list <- list(
  cytof      = auto_lambda_grid(ool_train$x_list$cytof,
                                ool_train$y, family = "gaussian", n_lambda = 6),
  proteomics = auto_lambda_grid(ool_train$x_list$proteomics,
                                ool_train$y, family = "gaussian", n_lambda = 6)
)
```

When you plan to cross-validate with
[`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md),
prefer `lambda_grid = "auto"` or derive grids from training data only. A
grid computed on the full dataset leaks mild information into the
cross-validated estimates; `"auto"` sidesteps this by building each
fold’s grid from its training data alone.

``` r

multi_fit <- stabl_multiomic_train_validate(
  x_train_list    = ool_train$x_list,
  y_train         = ool_train$y,
  lambda_grid     = lambda_list,
  x_valid_list    = ool_valid$x_list,
  y_valid         = ool_valid$y,
  family          = "gaussian",
  n_bootstraps    = 10L,
  artificial_type = "random_permutation",
  early_fusion    = TRUE,
  late_fusion     = TRUE,
  n_iter_lf       = 100L,
  random_state    = 42L
)

multi_fit
#> <stabl_multiomic_fit>
#>   Omics:           2 (cytof, proteomics)
#>   Per-omic selected features:
#>     cytof: 0
#>     proteomics: 0
#>   Validation data:  yes 
#>   Early fusion:     yes (0 features selected) 
#>   Late fusion:      yes (score = 0.1901)
```

------------------------------------------------------------------------

## 4. Exploring results

`multi_fit` is a list because the workflow has several branches. The
main elements are:

- **`$fits`**: named list of `stabl_fit` objects, one per omic. These
  are the per-omic STABL results produced *within* the joint pipeline.
- **`$early_fusion`**: result of the early-fusion step (if
  `early_fusion = TRUE`), containing `$fit` (a `stabl_fit` on the
  concatenated matrix).
- **`$late_fusion`**: result of the late-fusion step (if
  `late_fusion = TRUE`), containing `$valid_predictions` (validation
  predictions) and `$weights` (training-stack model weights).

``` r

# Per-omic selected feature names from the joint pipeline
cat("CyTOF (from multi-omic pipeline):\n")
#> CyTOF (from multi-omic pipeline):
print(get_feature_names_out(multi_fit$fits$cytof))
#> Warning: `get_feature_names_out()` is deprecated; use `get_feature_names()`
#> instead.
#> character(0)

cat("\nProteomics (from multi-omic pipeline):\n")
#> 
#> Proteomics (from multi-omic pipeline):
print(get_feature_names_out(multi_fit$fits$proteomics))
#> Warning: `get_feature_names_out()` is deprecated; use `get_feature_names()`
#> instead.
#> character(0)
```

``` r

# Early fusion: joint selection across all omics simultaneously
if (!is.null(multi_fit$early_fusion)) {
  cat("Early fusion selected feature names:\n")
  print(get_feature_names_out(multi_fit$early_fusion$fit))
}
#> Early fusion selected feature names:
#> Warning: `get_feature_names_out()` is deprecated; use `get_feature_names()`
#> instead.
#> character(0)
```

Early fusion can disagree with the per-omic lists because the joint
model sees all features at once. A feature that is useful alone may
become redundant when another omic carries the same information. That is
a property of the joint question, not a failure of the method.

``` r

# Late fusion: weighted blend learned on training predictions and applied to validation samples
if (!is.null(multi_fit$late_fusion)) {
  cat("Late fusion validation predictions (first 6):\n")
  print(head(multi_fit$late_fusion$valid_predictions))
}
#> Late fusion validation predictions (first 6):
#>  004_31_A  004_34_B  015_30_A  015_34_B  015_36_C  047_26_A 
#> -43.39333 -43.39333 -43.39333 -43.39333 -43.39333 -43.39333
```

`$valid_predictions` is a numeric vector of predicted DOS values for the
held-out validation set, one prediction per sample. Compare these to
`ool_valid$y` to assess predictive performance. The `$weights` element
records the omic weights chosen by the training-prediction stacking
objective; a small weight means that view received little weight in this
fitted prediction rule, not that the assay is biologically
uninformative.

``` r

plot_stabl_path(multi_fit$fits$cytof,
                title = "CyTOF (multi-omic pipeline)")
```

![CyTOF stability path from multi-omic
fit](stablr-multiomic_files/figure-html/stability-paths-1.png)

``` r

plot_stabl_path(multi_fit$fits$proteomics,
                title = "Proteomics (multi-omic pipeline)")
```

![CyTOF stability path from multi-omic
fit](stablr-multiomic_files/figure-html/stability-paths-2.png)

------------------------------------------------------------------------

## 5. Export results

[`save_stabl_results()`](https://mdmanurung.github.io/stablr/reference/save_stabl_results.md)
writes the result of a single `stabl_fit` object to disk:

| File | Content |
|----|----|
| `STABL scores.csv` | Full stability-score matrix. |
| `Max STABL scores.csv` | Per-feature maximum stability scores. |
| `FDR Graph.<fmt>` | FDP diagnostic plot when artificial features are used. |
| `Stability Path.<fmt>` | Stability path plot. |
| `Selected Features/Selected features.csv` | Features above the threshold. |
| `Selected Features/Feature distributions.<fmt>` | Per-feature outcome distributions when features are selected. |

Set `override = TRUE` to overwrite an existing directory, which is
useful when an analysis is being rerun with more bootstraps or a revised
grid.

``` r

out_dir <- file.path(tempdir(), "stablr_ool_cytof")

save_stabl_results(
  object    = multi_fit$fits$cytof,
  path      = out_dir,
  x         = ool_train$x_list$cytof,
  y         = ool_train$y,
  task_type  = "regression",
  figure_fmt = "png",
  override   = TRUE
)

list.files(out_dir, recursive = TRUE)
#> [1] "FDR Graph.png"                          
#> [2] "Max STABL artificial scores.csv"        
#> [3] "Max STABL scores.csv"                   
#> [4] "Selected Features/Selected features.csv"
#> [5] "Stability Path.png"                     
#> [6] "STABL artificial scores.csv"            
#> [7] "STABL scores.csv"
```

------------------------------------------------------------------------

## Next steps

- Use
  [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
  for a fully cross-validated evaluation when no pre-defined validation
  split is available.
- Pass `groups` to respect subject-level grouping in bootstrap
  resampling and CV fold construction.
- See
  [`vignette("stablr-intro")`](https://mdmanurung.github.io/stablr/articles/stablr-intro.md)
  for single-omic usage and alternative learners.
- See
  [`vignette("stablr-cooperative")`](https://mdmanurung.github.io/stablr/articles/stablr-cooperative.md)
  for cross-view cooperative fusion (`rho` tuning) when blocks should
  agree without full concatenation.
- For integrative **classification** benchmarks against mixOmics DIABLO
  on TCGA, see
  [`vignette("stablr-tcga-nestedcv")`](https://mdmanurung.github.io/stablr/articles/stablr-tcga-nestedcv.md)
  and `inst/analysis/BENCHMARK_PLAN.md`.
- For **selection parity** with the original Python STABL
  implementation, see `docs/PYTHON_TO_R_MAPPING.md`.
