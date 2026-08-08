# stablr: Sparse and Reliable Biomarker Discovery in R

`stablr` is an R-native implementation of STABL for sparse, stable biomarker
selection in high-dimensional clinical and omic data. It ports the
parity-critical STABL semantics from the Python implementation while exposing
R-native S3 objects, `glmnet`-ecosystem learners, multi-omic workflows,
visualization helpers, and CSV/disk export utilities.

The package compiles C++ but has no Python or tidymodels runtime dependency.

## Feature Summary

- Core STABL selector: `stabl_fit()`
- FDP+ thresholding with random-permutation and knockoff artificial features
- Lasso, elastic net, adaptive lasso, and optional sparse group lasso learners
- Gaussian, binomial, multinomial, and Cox outcome support where the backend
  supports the family
- Classic and group-aware bootstrap sampling with reproducible seeds
- Multi-omic train/validation and outer-CV workflows
- Early fusion, late fusion, and optional cooperative fusion
- S3 accessors, plotting helpers, export helpers, and reproducibility metrics

## Installation

### System requirements

`stablr` compiles a C++17 native cooperative-learning backend at install time.
You need a C++17-capable compiler toolchain:

- **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Rtools42 or later)
- **macOS:** Xcode Command Line Tools (`xcode-select --install`)
- **Linux:** GCC ≥ 7 or Clang ≥ 5 (usually pre-installed; install via `apt install r-base-dev` or equivalent)

### Install from GitHub

```r
# install.packages("remotes")
remotes::install_github("mdmanurung/stablr")
```

### Install from a local clone

From the repository root:

```r
install.packages("devtools")
devtools::install_local(".")
```

### Optional packages

Optional functionality uses optional packages:

- `ggplot2`: plotting helpers and saved result figures
- `knockoff`: `artificial_type = "knockoff"`, `"knockoff_equi"`, or
  `"knockoff_mvr"`
- `future`, `furrr`: parallel bootstrap execution
- `sparsegl`: `base_learner = "sparse_group_lasso"`
- `mixOmics`: TCGA vignette dataset

## Quick Start

```r
library(stablr)

set.seed(42)
n <- 80
p <- 20
x <- matrix(
  rnorm(n * p),
  nrow = n,
  dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p)))
)
y <- setNames(1.2 * x[, 1] - 0.8 * x[, 2] + rnorm(n), rownames(x))

lambda_grid <- auto_lambda_grid(x, y, family = "gaussian", n_lambda = 10)

fit <- stabl_fit(
  x = x,
  y = y,
  lambda_grid = lambda_grid,
  family = "gaussian",
  n_bootstraps = 50L,
  artificial_type = "random_permutation",
  random_state = 42L
)

get_feature_names_out(fit)
head(sort(get_importances(fit), decreasing = TRUE))
```

Use larger bootstrap counts and broader lambda grids for final analyses. The
quick-start vignette intentionally uses small settings to keep package builds
fast.

## Main Workflows

### Single-Omic STABL

Use `stabl_fit()` for one feature matrix. Inputs are aligned by sample names,
not position, so `rownames(x)` must match the names or row names of `y`.

Important arguments:

- `lambda_grid`: a data frame with a `lambda` column, or `"auto"`
- `base_learner`: `"lasso"`, `"elastic_net"`, `"adaptive_lasso"`, or
  `"sparse_group_lasso"`
- `family`: `"gaussian"`, `"binomial"`, `"multinomial"`, or `"cox"`
- `artificial_type`: `"random_permutation"`, `"knockoff"`,
  `"knockoff_equi"`, `"knockoff_mvr"`, or `NULL`
- `groups`: optional named group vector for grouped bootstrap sampling
- `random_state`: top-level seed for reproducible artificial features,
  bootstrap indices, and learner calls

### Multi-Omic Workflows

Use `stabl_multiomic_train_validate()` for named omic lists with optional
validation data. Enable `early_fusion`, `late_fusion`, and
`cooperative_fusion` independently. Use `stabl_multiomic_cv()` when no fixed
validation split is available.

Cooperative-fusion results can be inspected with:

```r
get_cooperative_features(multi_fit)
get_cooperative_diagnostics(multi_fit)
```

Cooperative fusion is built into `stablr` (gaussian and binomial families in
native v1). Validation-based cooperative tuning is not supported for Cox models;
Cox cooperative fusion is deferred to a future release.

## API Reference

Core fitting:

- `stabl_fit()`, `auto_lambda_grid()`
- `make_glmnet_adapter()`, `make_adaptive_lasso_adapter()`, `make_sgl_adapter()`

Input validation and bootstrapping:

- `validate_sample_alignment()`, `validate_multiomic_inputs()`
- `classic_bootstrap_indices()`, `group_bootstrap_indices()`

Artificial features and FDP+:

- `make_artificial_features()` with `"random_permutation"`, `"knockoff"`,
  `"knockoff_equi"`, or `"knockoff_mvr"`; `compute_fdp_plus()`

Accessors:

- `get_support()`, `get_feature_names_out()`, `transform_stabl()`,
  `get_stabl_scores()`, `get_importances()`
- `get_cooperative_features()`, `get_cooperative_diagnostics()`

Multi-omic workflows:

- `stabl_multiomic_train_validate()`, `stabl_multiomic_cv()`
- `stacked_multi_omic()`, `load_ool_data()`

Visualization and export:

- `plot_stabl_path()`, `plot_fdr_graph()`, `plot_roc()`, `plot_prc()`
- `boxplot_features()`, `scatterplot_features()`
- `export_stabl_to_csv()`, `save_stabl_results()`

Selection reproducibility metrics:

- `jaccard_similarity()`, `jaccard_matrix()`
- `adjusted_similarity()`, `adjusted_similarity_values()`,
  `adjusted_similarity_measure()`
- `pearson_similarity()`, `pearson_similarity_values()`,
  `pearson_similarity_measure()`
- `fdr_similarity()`, `tpr_similarity()`, `fscore_similarity()`

## Interpretation limits

The fitted threshold is the first minimizer of the observed FDP+ diagnostic
curve. A plotted `fdr_target` is a visual reference, not a fitting target, and
neither is a universal false-discovery guarantee. Interpretation depends on the
artificial-feature assumptions. Knockoff modes fail closed above 3,000
predictors because a globally valid chunked construction has not been
established; random-permutation controls remain available at larger dimensions.

Late-fusion weights default to leakage-safe OOF training. Empty selections or
downstream fit failures use fold-training means/event priors/class priors and
record the reason. The historical in-sample algorithm is available only through
`late_fusion_training = "python_legacy"`. Cross-language solver validation uses
ranking/support contracts; it does not claim bit-identical coefficients.

## Vignettes

Canonical source vignettes live in `vignettes/`:

- `stablr-intro.Rmd`: quick simulated-data start
- `stablr-multiomic.Rmd`: real OOL multi-omic train/validation workflow
- `stablr-tcga.Rmd`: TCGA Breast Cancer multi-omic workflow
- `stablr-tcga-nestedcv.Rmd`: TCGA nested CV and DIABLO comparison workflow
- `stablr-cooperative.Rmd`: cooperative fusion with the built-in multiview engine
- `stablr-advanced.Rmd`: Cox, multinomial, knockoffs, grouped bootstrap, metrics, export, outer CV

The Python-to-R API mapping is maintained as developer documentation in
`docs/PYTHON_TO_R_MAPPING.md`, not as a source vignette.

The `doc/` directory is generated by `R CMD build` or
`devtools::build_vignettes()` and is ignored in this source tree; edit
`vignettes/` as the canonical source.

## Documentation Website

Build the pkgdown site from the repository root:

```bash
conda run -n R4_51 Rscript -e "pkgdown::build_site('.', install = FALSE)"
```

The site reference index is configured in `_pkgdown.yml`.

## Application Note reproduction

Publication-scale OOL parity metrics (Bioinformatics Advances track):

```bash
# From stablr-experiments/ root (requires sample-data/Onset of Labor)
conda run -n R4_51 Rscript analysis/generate_publication_parity_table.R

# Smoke check on bundled extdata only
conda run -n R4_51 Rscript analysis/generate_publication_parity_table.R --allow-bundled
```

## Citation

If you use STABL, cite:

Hédou, J., Marić, I., Bellan, G. et al. Discovery of sparse, reliable omic
biomarkers with Stabl. Nat Biotechnol (2024).
https://doi.org/10.1038/s41587-023-02033-x
