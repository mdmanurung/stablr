# stablr: Sparse and Reliable Biomarker Discovery in R

### R-native STABL: sparse, stable biomarker selection with false-discovery calibration

The `stablr` package implements STABL (Stability-based Thresholding and
Aggregation for Biomarker discovery via LASSO) for high-dimensional
clinical and omic data. It resamples the data, fits a sparse learner
across a penalty grid, records how often each feature is selected, and
calibrates those selection frequencies against artificial decoy
features. The package ports the parity-critical STABL semantics from the
Python implementation while exposing R-native S3 objects,
`glmnet`-ecosystem learners, multi-omic workflows, visualization
helpers, and CSV/disk export utilities. It compiles C++ but has no
Python or tidymodels runtime dependency.

Major functionality of `stablr`:

- [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
                      One-stop-shop selector: a data.frame in, a
  plain-language selection out
- [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
                     Expert interface: full control of lambda grid,
  bootstrap design, thresholds, and learners
- [`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md)
                FDP+ thresholding with random-permutation and knockoff
  artificial features
- [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  Multi-omic train/validation with early, late, and cooperative fusion
- [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
      Outer CV with sample-i.i.d., grouped, or study-disjoint folds
- [`stabl_multiomic_scheme_comparison()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_scheme_comparison.md)
  Compare CV schemes arm-by-arm on held-out error
- [`stabl_environment_stability()`](https://mdmanurung.github.io/stablr/reference/stabl_environment_stability.md)
  Transportable selections across environments (experimental)
- [`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)
     Selection-reproducibility metrics

## Resources

- Hédou, J., Marić, I., Bellan, G. et al. Discovery of sparse, reliable
  omic biomarkers with Stabl. *Nat Biotechnol* (2024).
  [doi:10.1038/s41587-023-02033-x](https://doi.org/10.1038/s41587-023-02033-x)
- Ding, Y. D., Li, S., Narasimhan, B., Tibshirani, R. J. Cooperative
  learning for multi-view analysis. *PNAS* (2022).
  [doi:10.1073/pnas.2202113119](https://doi.org/10.1073/pnas.2202113119)
- [Documentation website](https://mdmanurung.github.io/stablr/) and the
  Vignettes menu above for executable tutorials
- [Release governance and review
  status](https://github.com/mdmanurung/stablr/blob/main/docs/RELEASE_STATUS.md)

## Motivation

A high-dimensional biomarker screen begins with an uncomfortable fact:
there are many more tempting variables than there are samples. A sparse
model can still choose a short list, but a single fit does not tell
whether the list is sturdy or merely lucky, and a selection that
reproduces only inside the sampling scheme that produced it has little
value for a new cohort or a new lab.

STABL turns the one-shot question into a repeated experiment. Bootstrap
resampling paired with artificial decoy features converts selection
frequency into a calibrated decision: a feature is retained only when it
out-competes decoys that carry no signal, and the fitted threshold is
the first minimizer of the observed FDP+ diagnostic curve. `stablr`
extends the same discipline across omic views (early, late, and
cooperative fusion) and across environments and studies (grouped and
study-disjoint resampling), so that reliability and transportability are
measured rather than assumed.

## stablr workflow

![stablr workflow](reference/figures/workflow.png)

## Technical intro

`stablr` is R-native without a Python runtime. The core STABL method is
parity-tested against the Python STABL implementation (Hédou et
al. 2024); cross-language solver validation uses ranking/support
contracts and does not claim bit-identical coefficients. The
cooperative-learning backend is compiled C++17 (vendored multiview
engine), and learners come from the `glmnet` ecosystem: lasso, elastic
net, adaptive lasso, and optional sparse group lasso, with gaussian,
binomial, multinomial, and Cox outcomes where the backend supports the
family.

The interface is two-tier.
[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md) runs
STABL with statistically safe defaults — model-X knockoff artificial
features, FDP+ calibration, strict fallback, automatic lambda grid — and
prints a plain-language answer:

``` r

library(stablr)
fit <- stabl(df, y = "response", n_bootstraps = 1000L, random_state = 42L)
fit                    # plain-language answer
feature_names(fit)     # the selected features
plot(fit)              # stability paths + FDP+ curve
```

[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
exposes every argument — lambda grids, bootstrap design, thresholds,
learner adapters, sample grouping — and every accessor and workflow
applies to both objects. See
[`?stabl_fit`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
and the vignettes.

Everything beyond core STABL is labelled as an extension: cooperative
fusion (Ding et al. 2022) ships on the vendored multiview engine for
gaussian and binomial families (Cox deferred); the random-forest learner
is experimental until its FDP+ calibration study completes; the
environment-stability workflow aggregates per-environment selections by
their minimum frequency.

Interpretation limits are part of the contract. The fitted threshold is
the first minimizer of the observed FDP+ curve; a plotted `fdr_target`
is a visual reference, not a fitting target, and neither is a universal
false-discovery guarantee. Knockoff modes fail closed above 3,000
predictors because a globally valid chunked construction has not been
established; random-permutation controls remain available at larger
dimensions.

## Install

`stablr` compiles a C++17 native cooperative-learning backend at install
time. You need a C++17-capable compiler toolchain:

- **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  (Rtools42 or later)
- **macOS:** Xcode Command Line Tools (`xcode-select --install`)
- **Linux:** GCC ≥ 7 or Clang ≥ 5 (usually pre-installed; install via
  `apt install r-base-dev` or equivalent)

### Install

``` r

# Install from GitHub
remotes::install_github("mdmanurung/stablr")

# Or from a local clone, from the repository root:
devtools::install_local(".")
```

### Dependencies

Optional functionality uses optional packages, in case these aren’t
installed automatically:

``` r

install.packages(c("ggplot2", "knockoff", "future", "furrr", "sparsegl"))
```

- `ggplot2`: plotting helpers and saved result figures
- `knockoff`: `artificial_type = "knockoff"`, `"knockoff_equi"`, or
  `"knockoff_mvr"`
- `future`, `furrr`: parallel bootstrap execution
- `sparsegl`: `base_learner = "sparse_group_lasso"`
- `mixOmics`: TCGA vignette dataset
