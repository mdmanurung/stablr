# stablr: Sparse and Reliable Biomarker Discovery in R

`stablr` is an R-native implementation of STABL for sparse, stable
biomarker selection in high-dimensional clinical and omic data. It ports
the parity-critical STABL semantics from the Python implementation while
exposing R-native S3 objects, `glmnet`-ecosystem learners, multi-omic
workflows, visualization helpers, and CSV/disk export utilities.

The package compiles C++ but has no Python or tidymodels runtime
dependency.

## Feature Summary

- Core STABL selector:
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
- FDP+ thresholding with random-permutation and knockoff artificial
  features
- Lasso, elastic net, adaptive lasso, and optional sparse group lasso
  learners
- Gaussian, binomial, multinomial, and Cox outcome support where the
  backend supports the family
- Classic and group-aware bootstrap sampling with reproducible seeds
- Multi-omic train/validation and outer-CV workflows
- Early fusion, late fusion, and optional cooperative fusion
- S3 accessors, plotting helpers, export helpers, and reproducibility
  metrics

## Installation

### System requirements

`stablr` compiles a C++17 native cooperative-learning backend at install
time. You need a C++17-capable compiler toolchain:

- **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  (Rtools42 or later)
- **macOS:** Xcode Command Line Tools (`xcode-select --install`)
- **Linux:** GCC ≥ 7 or Clang ≥ 5 (usually pre-installed; install via
  `apt install r-base-dev` or equivalent)

### Install from GitHub

``` r

# install.packages("remotes")
remotes::install_github("mdmanurung/stablr")
```

### Install from a local clone

From the repository root:

``` r

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

## Quick start

**New to stablr? Use
[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md).**
Give it a table of features and an outcome; it runs STABL with
statistically safe defaults (model-X knockoff artificial features, FDP+
calibration, strict fallback, automatic lambda grid) and prints a
plain-language answer:

``` r

library(stablr)

set.seed(42)
df <- data.frame(m1 = rnorm(80), m2 = rnorm(80), m3 = rnorm(80))
df$response <- 1.2 * df$m1 - 0.8 * df$m2 + rnorm(80)

fit <- stabl(df, y = "response", n_bootstraps = 1000L, random_state = 42L)
fit                    # plain-language answer
feature_names(fit)     # the selected features
plot(fit)              # stability paths + FDP+ curve
```

**Need full control? Use
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).**
The expert interface exposes all arguments — lambda grids, bootstrap
design, thresholds, learner adapters, sample grouping — and every
accessor and workflow applies to both `stabl` and `stabl_fit` objects.
See
[`?stabl_fit`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
and the vignettes.

## Five-minute example (expert interface)

``` r

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

get_feature_names(fit)
head(sort(get_importances(fit), decreasing = TRUE))
```

This is an executable orientation example, not release or publication
evidence. Use larger bootstrap counts and broader lambda grids for a
predeclared final analysis, and inspect the FDP+ diagnostic rather than
interpreting the selected set alone.

## Main Workflows

### Single-Omic STABL

Use
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
for one feature matrix. Inputs are aligned by sample names, not
position, so `rownames(x)` must match the names or row names of `y`.

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

Use
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
for named omic lists with optional validation data. Enable
`early_fusion`, `late_fusion`, and `cooperative_fusion` independently.
Use
[`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
when no fixed validation split is available.

Cooperative-fusion results can be inspected with:

``` r

get_cooperative_features(multi_fit)
get_cooperative_diagnostics(multi_fit)
```

Cooperative fusion is built into `stablr` (gaussian and binomial
families in native v1). Validation-based cooperative tuning is not
supported for Cox models; Cox cooperative fusion is deferred to a future
release.

## Extensions beyond Python STABL

The core STABL method remains parity-tested against the Python STABL
implementation (Hédou et al. 2024, Nat Biotechnol). Everything below is
an R-native extension and is labelled as such.

- **Cooperative fusion** (Ding et al. 2022, PNAS): multi-view learning
  that encourages omic blocks measured on shared samples to agree,
  useful when views are complementary rather than redundant; shipped on
  the vendored multiview engine (gaussian and binomial families).
- **Experimental random-forest learner** (stablr extension, validation
  pending): a `ranger`-based base learner candidate for non-linear
  signals; results are exploratory until its FDP+ calibration study
  completes.
- **Environment-stability workflow**: fit STABL per environment and
  aggregate selection frequencies by their minimum across environments,
  for feature sets that are both reliable and transportable.

## API Reference

Beginner entry point:

- [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md) —
  one-stop-shop wrapper with safe defaults;
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`support()`](https://mdmanurung.github.io/stablr/reference/support.md),
  [`feature_names()`](https://mdmanurung.github.io/stablr/reference/support.md)
  methods

Core fitting:

- [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
  [`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md)
- [`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
  [`make_adaptive_lasso_adapter()`](https://mdmanurung.github.io/stablr/reference/make_adaptive_lasso_adapter.md),
  [`make_sgl_adapter()`](https://mdmanurung.github.io/stablr/reference/make_sgl_adapter.md)

Input validation and bootstrapping:

- [`validate_sample_alignment()`](https://mdmanurung.github.io/stablr/reference/validate_sample_alignment.md),
  [`validate_multiomic_inputs()`](https://mdmanurung.github.io/stablr/reference/validate_multiomic_inputs.md)
- [`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md),
  [`group_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/group_bootstrap_indices.md)

Artificial features and FDP+:

- [`make_artificial_features()`](https://mdmanurung.github.io/stablr/reference/make_artificial_features.md)
  with `"random_permutation"`, `"knockoff"`, `"knockoff_equi"`, or
  `"knockoff_mvr"`;
  [`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md)

Accessors:

- [`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md),
  [`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md),
  [`transform_stabl()`](https://mdmanurung.github.io/stablr/reference/transform_stabl.md),
  [`get_stabl_scores()`](https://mdmanurung.github.io/stablr/reference/get_stabl_scores.md),
  [`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md)
- [`get_cooperative_features()`](https://mdmanurung.github.io/stablr/reference/get_cooperative_features.md),
  [`get_cooperative_diagnostics()`](https://mdmanurung.github.io/stablr/reference/get_cooperative_diagnostics.md)

Multi-omic workflows:

- [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md),
  [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
- [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md),
  [`load_ool_data()`](https://mdmanurung.github.io/stablr/reference/load_ool_data.md)

Visualization and export:

- [`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md),
  [`plot_fdr_graph()`](https://mdmanurung.github.io/stablr/reference/plot_fdr_graph.md),
  [`plot_roc()`](https://mdmanurung.github.io/stablr/reference/plot_roc.md),
  [`plot_prc()`](https://mdmanurung.github.io/stablr/reference/plot_prc.md)
- [`boxplot_features()`](https://mdmanurung.github.io/stablr/reference/boxplot_features.md),
  [`scatterplot_features()`](https://mdmanurung.github.io/stablr/reference/scatterplot_features.md)
- [`export_stabl_to_csv()`](https://mdmanurung.github.io/stablr/reference/export_stabl_to_csv.md),
  [`save_stabl_results()`](https://mdmanurung.github.io/stablr/reference/save_stabl_results.md)

Selection reproducibility metrics:

- [`jaccard_similarity()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity.md),
  [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)
- [`adjusted_similarity()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity.md),
  [`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md),
  [`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md)
- [`pearson_similarity()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity.md),
  [`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md),
  [`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md)
- [`fdr_similarity()`](https://mdmanurung.github.io/stablr/reference/fdr_similarity.md),
  [`tpr_similarity()`](https://mdmanurung.github.io/stablr/reference/tpr_similarity.md),
  [`fscore_similarity()`](https://mdmanurung.github.io/stablr/reference/fscore_similarity.md)

Deprecated names (still work; warn once and delegate):
[`get_feature_names_out()`](https://mdmanurung.github.io/stablr/reference/get_feature_names_out.md),
[`stacked_multi_omic()`](https://mdmanurung.github.io/stablr/reference/stacked_multi_omic.md),
[`jaccard_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_matrix.md),
[`classic_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/classic_bootstrap_indices.md),
[`adjusted_similarity_values()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_values.md),
[`adjusted_similarity_measure()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_measure.md),
[`pearson_similarity_values()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_values.md),
[`pearson_similarity_measure()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_measure.md)
— see `NEWS.md` for the full rename table.

## Interpretation limits

The fitted threshold is the first minimizer of the observed FDP+
diagnostic curve. A plotted `fdr_target` is a visual reference, not a
fitting target, and neither is a universal false-discovery guarantee.
Interpretation depends on the artificial-feature assumptions. Knockoff
modes fail closed above 3,000 predictors because a globally valid
chunked construction has not been established; random-permutation
controls remain available at larger dimensions.

Late-fusion weights default to leakage-safe OOF training. Empty
selections or downstream fit failures use fold-training means/event
priors/class priors and record the reason. The historical in-sample
algorithm is available only through
`late_fusion_training = "python_legacy"`. Read
`oof_predictions`/`oof_score` in OOF mode and
`in_sample_predictions`/`in_sample_score` in legacy mode;
`train_predictions` and `score` are exact 0.1.x compatibility aliases.
Cross-language solver validation uses ranking/support contracts; it does
not claim bit-identical coefficients.

## Vignettes

Canonical source vignettes live in `vignettes/`:

- `stablr-intro.Rmd`: quick simulated-data start
- `stablr-multiomic.Rmd`: executable OOL workflow tutorial using 150
  training samples, 21 aligned validation samples, and 100 features per
  assay
- `stablr-tcga.Rmd`: non-executable TCGA workflow research protocol
- `stablr-tcga-nestedcv.Rmd`: non-executable TCGA/DIABLO benchmark
  protocol; no external results are bundled
- `stablr-cooperative.Rmd`: cooperative fusion with the built-in
  multiview engine
- `stablr-advanced.Rmd`: Cox, multinomial, knockoffs, grouped bootstrap,
  metrics, export, outer CV

The machine-readable audit of all 42 exports and 19 registered S3
methods is maintained in `inst/release/public-interface-audit.csv`.

The Python-to-R API mapping is maintained as developer documentation in
`docs/PYTHON_TO_R_MAPPING.md`, not as a source vignette.

Contributor documentation: the release-evidence glossary and validation
vocabulary live in
[docs/CONTEXT.md](https://mdmanurung.github.io/stablr/docs/CONTEXT.md).

The `doc/` directory is generated by `R CMD build` or
`devtools::build_vignettes()` and is ignored in this source tree; edit
`vignettes/` as the canonical source.

## Documentation Website

Build the pkgdown site from the repository root:

``` bash
conda run -n R4_51 Rscript -e "pkgdown::build_site('.', install = FALSE)"
```

The site reference index is configured in `_pkgdown.yml`.

## Application workflow status

The bundled OOL vignette is a workflow tutorial. It is not
application-note evidence. Release claims require a separate run from an
exact, hashed candidate source tarball and the predeclared full OOL
inputs and gates; mutable-checkout execution does not qualify.

``` bash
# From stablr-experiments/ after building the candidate tarball. The output
# directory must not exist.
RUNTIME=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/lib/R/bin/Rscript
"$RUNTIME" --vanilla analysis/generate_publication_parity_table.R \
  --candidate-tarball /absolute/path/to/stablr_0.1.1.tar.gz \
  --runtime "$RUNTIME" \
  --outdir /absolute/new/path/ool-candidate-evidence \
  --profile full
```

The frozen full-candidate protocol currently produces a negative result:
support Jaccard 0.778 (required 0.80), score Spearman 0.785 (required
0.90), and FDP-threshold delta 0.08 (maximum 0.05). The lambda and
training-only preprocessing gates pass. The requested fixed knockoff
construction is infeasible for the 150-by-1,317 matrix (`n > p` is
required), so the typed, recorded fallback is `random_permutation` and
the actual-generator gate also fails. The sealed negative artifact is
retained; `stablr` 0.1.1 is therefore not ready for release review.

## Citation

If you use STABL, cite:

Hédou, J., Marić, I., Bellan, G. et al. Discovery of sparse, reliable
omic biomarkers with Stabl. Nat Biotechnol (2024).
<https://doi.org/10.1038/s41587-023-02033-x>
