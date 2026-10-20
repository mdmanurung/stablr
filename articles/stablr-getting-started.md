# Getting started with stabl

## Contents

- 1 Introduction
- 2 A worked example
- 3 Defaults chosen by
  [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
- 4 Leakage-safe preprocessing
- 5 The expert interface
- 6 Migrating from pre-1.0.0 versions

## 1 Introduction

[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md) is
the high-level entry point to **stablr**. Given a table of features and
an outcome, it runs STABL with statistically safe defaults and returns
an object that states, in plain language, which features were selected
and how reliable the selection is estimated to be.

The full contract is as follows.

``` r

fit <- stabl(my_data, y = "response")
fit
```

Everything else in the package is available when more control is
required; see Section 5, The expert interface.

## 2 A worked example

The example uses a small simulated proteomic-style dataset: 60 samples
and 8 features, of which two truly drive the outcome.

``` r

set.seed(1L)
df <- data.frame(
  m1 = rnorm(60), m2 = rnorm(60), m3 = rnorm(60), m4 = rnorm(60),
  m5 = rnorm(60), n1 = rnorm(60), n2 = rnorm(60), n3 = rnorm(60)
)
df$response <- 0.9 * df$m1 + 0.9 * df$m2 + rnorm(60, sd = 0.4)
```

STABL is fitted with a reduced number of bootstraps so that the vignette
knits quickly; the default (`n_bootstraps = 1000`) should be used in
practice.

``` r

library(stablr)
fit <- stabl(df, y = "response", n_bootstraps = 200L, random_state = 1L)
```

Printing the fit gives the plain-language answer.

``` r

fit
#> <stabl stability selection>
#>   outcome family: gaussian | base learner: lasso | bootstraps: 200
#>   3 feature(s) selected at FDP+ estimate 0.333:
#>   m1, m2, n3
#>   use summary() for details, plot() for stability paths,
#>   get_feature_names() / get_support() for the selection,
#>   get_fdrs() for the FDP+ curve.
```

The printed output reports the selected features together with the
estimated reliability of the selection.

The selected features are also available as a character vector.

``` r

get_feature_names(fit)
#> [1] "m1" "m2" "n3"
```

[`summary()`](https://rdrr.io/r/base/summary.html) adds the stability
threshold and the FDP+ estimate.

``` r

summary(fit)
#> <stabl stability selection summary>
#>   features considered : 8
#>   bootstraps          : 200
#>   stability threshold : 0.990
#>   FDP+ estimate       : 0.333
#>   selected (3)       : m1, m2, n3
```

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
draws the stability paths and the FDP+ curve.

``` r

plot(fit)
```

The figure shows the stability paths of the features across the lambda
grid together with the FDP+ curve used for calibration.

## 3 Defaults chosen by `stabl()`

| Decision | Default chosen | Why |
|----|----|----|
| Artificial features | model-X knockoffs (`knockoff_equi`) | conservative decoys for tree learners; measured calibration is learner- and regime-dependent — see the pairing table in [`vignette("stablr-intro")`](https://mdmanurung.github.io/stablr/articles/stablr-intro.md) before trusting FDP+ |
| Calibration | FDP+ on (no fixed `hard_threshold`) | the threshold is chosen to control the estimated false discovery proportion rather than fixed by hand |
| Knockoff failure | `knockoff_fallback = "stop"` | a silent fallback would silently change what the FDP+ estimates mean |
| Lambda grid | automatic | a sensible grid is derived from the data |

Any expert argument — `artificial_proportion`, `sample_fraction`,
`sample_groups`, a custom adapter — is passed through `...`.

``` r

fit <- stabl(df, y = "response", base_learner = "elastic_net",
             artificial_proportion = 0.5)
```

## 4 Leakage-safe preprocessing

The worked example above fits on the full dataset. That is correct for
exploratory selection — the selection describes the data you have. The
moment you want held-out predictions or a train/test evaluation, one
rule applies to every *learned* transform (centering, scaling,
imputation, variance filtering):

**Fit the transform on the training samples only; apply the stored
training parameters to the test samples.**

Two shortcuts both corrupt the analysis. Fitting the transform on the
pooled train + test data leaks test-distribution information into
training — the test means and variances shift the fitted parameters, and
the held-out estimate is optimistically biased. Fitting a separate
transform on the test data is not leakage, but it places the test
samples in a different feature space, so a model trained on the training
transform mispredicts them.

[`stabl_preprocess_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit.md)
and
[`stabl_preprocess_apply()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_apply.md)
implement the rule as a fit/apply pair.

``` r

x <- as.matrix(df[, 1:8])
idx <- seq_len(40)  # first 40 samples train, remaining 20 test

pp <- stabl_preprocess_fit(x[idx, ], transforms = c("center", "scale"))
x_train <- stabl_preprocess_apply(pp, x[idx, ])
x_test  <- stabl_preprocess_apply(pp, x[-idx, ])
```

`pp` stores the training means and standard deviations, and `x_test` is
scaled with those values — not with its own. You can see the difference:
the test matrix is centred with the training means, so its own column
means are not exactly zero.

``` r

round(colMeans(x_test), 3)
#>     m1     m2     m3     m4     m5     n1     n2     n3 
#>  0.053  0.003  0.830 -0.396  0.087 -0.142  0.182  0.588
```

The same rule is applied automatically inside every fold of
`stabl_analyze(mode = "validated")` and inside the multi-omic
cross-validation workflows, so a validated analysis cannot leak through
preprocessing. For the multi-view case — including why the validation
matrix must be transformed with the training parameters before
late-fusion prediction — see
[`vignette("stablr-multiomic")`](https://mdmanurung.github.io/stablr/articles/stablr-multiomic.md).

## 5 The expert interface

[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md) is a
thin wrapper around
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
the expert interface with full control over all arguments (lambda grids,
bootstrap design, thresholds, learner adapters, grouping). Every
[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
result is also a valid `stabl_fit` object, so the entire expert toolkit
applies unchanged.

``` r

# accessors (identical on stabl and stabl_fit objects)
head(get_stabl_scores(fit))
#>     [,1]  [,2]  [,3]  [,4]  [,5] [,6] [,7]  [,8] [,9] [,10] [,11] [,12] [,13]
#> m1 0.510 0.875 0.975 0.995 1.000 1.00 1.00 1.000 1.00 1.000 1.000 1.000 1.000
#> m2 0.465 0.900 0.995 1.000 1.000 1.00 1.00 1.000 1.00 1.000 1.000 1.000 1.000
#> m3 0.000 0.000 0.000 0.000 0.000 0.01 0.07 0.155 0.26 0.400 0.540 0.650 0.755
#> m4 0.000 0.000 0.000 0.005 0.010 0.02 0.11 0.200 0.32 0.445 0.535 0.620 0.715
#> m5 0.000 0.020 0.060 0.155 0.390 0.58 0.68 0.755 0.80 0.830 0.835 0.845 0.850
#> n1 0.000 0.000 0.010 0.035 0.115 0.24 0.35 0.405 0.51 0.580 0.670 0.715 0.755
#>    [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21]
#> m1 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
#> m2 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
#> m3 0.815 0.865 0.900 0.915 0.945 0.950 0.970 0.975
#> m4 0.790 0.840 0.885 0.920 0.945 0.960 0.970 0.985
#> m5 0.850 0.895 0.895 0.900 0.930 0.950 0.960 0.980
#> n1 0.825 0.860 0.870 0.910 0.935 0.955 0.965 0.980
get_fdrs(fit)
#>         0      0.01      0.02      0.03      0.04      0.05      0.06      0.07 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.08      0.09       0.1      0.11      0.12      0.13      0.14      0.15 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.16      0.17      0.18      0.19       0.2      0.21      0.22      0.23 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.24      0.25      0.26      0.27      0.28      0.29       0.3      0.31 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.32      0.33      0.34      0.35      0.36      0.37      0.38      0.39 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>       0.4      0.41      0.42      0.43      0.44      0.45      0.46      0.47 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.48      0.49       0.5      0.51      0.52      0.53      0.54      0.55 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.56      0.57      0.58      0.59       0.6      0.61      0.62      0.63 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.64      0.65      0.66      0.67      0.68      0.69       0.7      0.71 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.72      0.73      0.74      0.75      0.76      0.77      0.78      0.79 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>       0.8      0.81      0.82      0.83      0.84      0.85      0.86      0.87 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.88      0.89       0.9      0.91      0.92      0.93      0.94      0.95 
#> 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 1.1250000 
#>      0.96      0.97      0.98      0.99 
#> 1.0000000 0.8750000 0.4000000 0.3333333

# the multiomic workflows, metrics, and plotting helpers all accept
# stabl_fit objects — including the one you just built with stabl()
class(fit)
#> [1] "stabl"     "stabl_fit"
```

The printed accessors show the stability scores, the FDP+ summary, and
the class of the fitted object, which confirms that a
[`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
result can be treated as a `stabl_fit` object. Read
[`?stabl_fit`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
for the full expert interface, and the other vignettes for multi-omic
workflows, nested cross-validation, and environment stability. For what
“Python parity” does and does not establish, see
[`vignette("stablr-python-parity")`](https://mdmanurung.github.io/stablr/articles/stablr-python-parity.md);
for a guide to the stability-path and FDP+ diagnostic plots, see
[`vignette("stablr-reading-plots")`](https://mdmanurung.github.io/stablr/articles/stablr-reading-plots.md).

## 6 Migrating from pre-1.0.0 versions

The 1.0.0 contract removes the pre-0.3.0 compatibility names outright
(`get_feature_names_out()`, `stacked_multi_omic()`, `jaccard_matrix()`,
`adjusted_similarity_values()`, `adjusted_similarity_measure()`,
`pearson_similarity_values()`, `pearson_similarity_measure()`,
`classic_bootstrap_indices()`); calling them now raises an error. Use
the canonical names —
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md),
[`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md),
[`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md),
[`adjusted_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_pairwise.md),
[`adjusted_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/adjusted_similarity_summary.md),
[`pearson_similarity_pairwise()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_pairwise.md),
[`pearson_similarity_summary()`](https://mdmanurung.github.io/stablr/reference/pearson_similarity_summary.md),
[`sample_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/sample_bootstrap_indices.md)
— and see `NEWS.md` and `docs/MIGRATION_0.3_TO_1.0.md` for the full
mapping.

The 1.0.0 contract also removes a small set of legacy behaviors (Phase
13 PR-A): the `explore`/`n_explore` arguments of
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
(an empty selection now stays empty instead of being relabelled as a
top-N “explore” set), the in-sample
`late_fusion_training = "python_legacy"` mode, and
`cooperation_selection = "validation"` (both biased performance
estimates). The beginner generics `support()` and `feature_names()` are
removed in favour of the canonical accessors
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md)
and
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md),
and the experimental oblique random-forest learner
(`base_learner = "oblique_rf"`) is retired (use `"rf"`); see the
migration guide for the complete list.
