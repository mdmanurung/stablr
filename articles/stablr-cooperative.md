# Cooperative Fusion for Multi-Omic Biomarker Discovery

## Overview

Early fusion asks all features to compete in one wide matrix. Late
fusion lets each omic make its own prediction and combines the answers
afterward. Cooperative fusion asks a subtler question: can the views be
fitted together without erasing their identities?

In `stablr`, cooperative fusion is built in natively (vendored from the
[multiview](https://cran.r-project.org/package=multiview)
cooperative-learning engine). Each view keeps its own coefficient path,
but the paths are encouraged to agree through a cooperation strength
`rho`. The omics are neither fully pooled nor fully isolated. Native v1
supports **gaussian** and **binomial** families.

**When should you use cooperative fusion?**

| Strategy | How it works | Best when |
|----|----|----|
| **Per-omic** | One STABL fit per omic, independently. | Baseline; or when views are likely orthogonal. |
| **Early fusion** | Concatenate all omics, one joint STABL fit. | Omics are tightly correlated; you want a single feature list. |
| **Late fusion** | Per-omic fits, then learn a weighted blend of predictions. | Each omic may predict independently; you want omic weights. |
| **Cooperative fusion** | Joint multiview fit with cooperation strength `rho`; each view has its own coefficient path. | Omics share some signal structure; you want cross-view feature selection to be coherent. |

The cooperation strength `rho` controls how much the views are
encouraged to share signal. At `rho = 0`, cooperative fusion reduces to
independent fits. As `rho` grows, the views are penalised toward a
common prediction. `stablr` tunes `rho` over a user-supplied grid.

## Learning objectives

By the end, you should be able to:

1.  Explain what `rho` controls in cooperative fusion and how to choose
    a grid.
2.  Run
    [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
    with `cooperative_fusion = TRUE` and navigate the
    `$cooperative_fusion` result element.
3.  Read the [`print()`](https://rdrr.io/r/base/print.html) output for
    cooperative fusion details.
4.  Understand which `cooperation_selection` and `cooperation_selector`
    combinations are supported and why.
5.  Compare cooperative fusion feature selections to early and late
    fusion on the same dataset.
6.  Use
    [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
    to evaluate cooperative fusion across outer folds.

------------------------------------------------------------------------

## 1. Load the data

We use the same Onset of Labor (OOL) dataset as
[`vignette("stablr-multiomic")`](https://mdmanurung.github.io/stablr/articles/stablr-multiomic.md):
CyTOF immune profiling and Proteomics measured on the same women, with
days-to-onset-of-labor (DOS) as the continuous outcome. The examples use
compact grids so the vignette remains buildable; enlarge the bootstrap
count and `rho` grid for a final analysis.

``` r

library(stablr)
```

``` r

ool_train <- load_ool_data(split = "train")
ool_valid <- load_ool_data(split = "valid")

# Shape of each omic layer
lapply(ool_train$x_list, dim)
#> $cytof
#> [1] 150 100
#> 
#> $proteomics
#> [1] 150 100
```

------------------------------------------------------------------------

## 2. Build lambda grids

Cooperative fusion still needs a lambda grid for each omic. We build one
grid per view with
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md)
and pass the grids as a named list. Twelve lambda values are enough to
demonstrate the workflow without making each cooperative fit too
expensive.

``` r

lambda_list <- list(
  cytof      = auto_lambda_grid(
    ool_train$x_list$cytof, ool_train$y,
    family = "gaussian", n_lambda = 12
  ),
  proteomics = auto_lambda_grid(
    ool_train$x_list$proteomics, ool_train$y,
    family = "gaussian", n_lambda = 12
  )
)

# Number of lambda values per omic
vapply(lambda_list, nrow, integer(1L))
#>      cytof proteomics 
#>         12         12
```

------------------------------------------------------------------------

## 3. Running cooperative fusion

Activate cooperative fusion with `cooperative_fusion = TRUE` and a `rho`
grid. `stablr` fits one multiview model per `rho` value using shared
inner CV folds, then picks the best `rho` by the chosen
`cooperation_type_measure`. This is **sequential** tuning over `rho`
(one full CV path per value), not a single joint CV over `(rho, lambda)`
pairs.

Key arguments:

| Argument | Default | Purpose |
|----|----|----|
| `rho` | `0` | Scalar or vector of cooperation strengths to try. |
| `cooperation_selection` | `"cv"` | `"cv"` (inner CV) or `"validation"` (held-out set). |
| `cooperation_selector` | `"lambda.min"` | Which lambda to pick: `"lambda.min"` or `"lambda.1se"` (CV only). |
| `cooperation_type_measure` | `"default"` | Tuning metric; defaults to `"mse"` for Gaussian. |
| `cooperation_nfolds` | `5L` | Number of inner CV folds when `cooperation_selection = "cv"`. |

``` r

cf_fit <- stabl_multiomic_train_validate(
  x_train_list          = ool_train$x_list,
  y_train               = ool_train$y,
  lambda_grid           = lambda_list,
  x_valid_list          = ool_valid$x_list,
  y_valid               = ool_valid$y,
  family                = "gaussian",
  n_bootstraps          = 30L,
  artificial_type       = "random_permutation",
  random_state          = 42L,
  cooperative_fusion    = TRUE,
  rho                   = c(0, 0.25, 0.5),
  cooperation_selection = "cv",
  cooperation_selector  = "lambda.min",
  cooperation_nfolds    = 3L
)
```

[`print()`](https://rdrr.io/r/base/print.html) shows every active
branch. When cooperative fusion is present, a `Cooperative fusion:`
block appears at the bottom:

``` r

cf_fit
#> <stabl_multiomic_fit>
#>   Omics:           2 (cytof, proteomics)
#>   Per-omic selected features:
#>     cytof: 0
#>     proteomics: 2
#>   Validation data:  yes 
#>   Early fusion:     no 
#>   Late fusion:      no 
#>   Cooperative fusion:
#>     selection:     cv
#>     rho (chosen):  0
#>     selector:      lambda.min
#>     type.measure:  mse
#>     score:         812.1337
#>     selected feats: 26 (across views)
```

**Reading the cooperative fusion block:**

- **`selection`** — how `rho` was tuned: `"cv"` here.
- **`rho (chosen)`** — the `rho` value selected by cross-validation. A
  value of `0` means the specified metric and folds preferred
  independent fits within the supplied grid. A larger value means that
  the same tuning surface preferred some coupling between views; it is
  not by itself evidence of biological cooperation.
- **`selector`** — which lambda rule was applied (`lambda.min`).
- **`type.measure`** — the CV metric used for tuning (e.g., `mse`).
- **`score`** — the CV score at the chosen `rho` and `lambda`.
- **`selected feats`** — total number of features selected across all
  views.

------------------------------------------------------------------------

## 4. Navigating `$cooperative_fusion`

The `cooperative_fusion` element contains the selected model and the
tuning surface that led to it.

``` r

cf <- cf_fit$cooperative_fusion

# Chosen hyperparameters
cat("Chosen rho:     ", cf$rho, "\n")
#> Chosen rho:      0
cat("Chosen lambda:  ", round(cf$selected_lambda, 4), "\n")
#> Chosen lambda:   2.4293
cat("Selector:       ", cf$selector, "\n")
#> Selector:        lambda.min
cat("CV score (MSE): ", round(cf$score, 4), "\n")
#> CV score (MSE):  812.1337
```

### Selected features per view

`$selected_features` is a named list with one character vector per omic:

``` r

cat("CyTOF features selected by cooperative fusion:\n")
#> CyTOF features selected by cooperative fusion:
print(cf$selected_features$cytof)
#> [1] "Granulocytes_S6_unstim" "Bcells_S6_unstim"       "MDSCs_CREB_unstim"     
#> [4] "DCs_CREB_unstim"        "pDCs_STAT3_unstim"      "pDCs_ERK_unstim"       
#> [7] "intMCs_CREB_unstim"

cat("\nProteomics features selected by cooperative fusion:\n")
#> 
#> Proteomics features selected by cooperative fusion:
print(cf$selected_features$proteomics)
#>  [1] "Mcl.1"               "HEMK2"               "CD59"               
#>  [4] "S100A6"              "SECTM1"              "CSH"                
#>  [7] "NEGR1"               "FSTL1"               "Osteopontin"        
#> [10] "SMOC1"               "protein.Z.inhibitor" "ISLR2"              
#> [13] "HXK2"                "SEM5A"               "ApoM"               
#> [16] "AMGO2"               "C1QR1"               "H2B2E"              
#> [19] "Myostatin"
```

### Predictions

`$train_predictions` and `$valid_predictions` hold the multiview
predictions from the chosen model:

``` r

cat("First 6 training predictions:\n")
#> First 6 training predictions:
print(head(cf$train_predictions))
#>  001_26_A  001_33_B  001_35_C  003_25_A  003_30_B  003_32_C 
#> -70.90491 -11.58418 -10.49450 -52.14157 -27.99975 -23.70652

cat("\nFirst 6 validation predictions:\n")
#> 
#> First 6 validation predictions:
print(head(cf$valid_predictions))
#>  004_31_A  004_34_B  015_30_A  015_34_B  015_36_C  047_26_A 
#> -56.56061 -46.02405 -44.23919 -11.32462 -33.57264 -48.94823
```

### Tuning diagnostics

`$diagnostics` records the CV score at every `rho` in the grid, so the
tuning choice can be inspected rather than treated as a black box:

``` r

# One row per rho value
print(cf$diagnostics[, c("rho", "lambda", "metric_value", "selected")])
#>    rho   lambda metric_value selected
#> 1 0.00 2.429267     812.1337     TRUE
#> 2 0.25 1.754128     816.5125    FALSE
#> 3 0.50 2.016820     842.6442    FALSE
```

The `selected` column marks the winning `rho`. If the surface is flat,
with similar `metric_value` across `rho`, the omics are not gaining much
from cooperation and `rho = 0` is a defensible choice.

------------------------------------------------------------------------

## 5. Policy constraints

Some combinations are rejected on purpose. These guards prevent analyses
that look reasonable in code but do not have a valid tuning
interpretation.

### Constraint 1 - `lambda.1se` requires CV mode

`lambda.1se` applies the “one-standard-error” rule and is only
well-defined when cross-validation is used to estimate the error
surface. Requesting it with `cooperation_selection = "validation"`
raises an error:

``` r

stabl_multiomic_train_validate(
  x_train_list          = ool_train$x_list,
  y_train               = ool_train$y,
  lambda_grid           = lambda_list,
  x_valid_list          = ool_valid$x_list,
  y_valid               = ool_valid$y,
  family                = "gaussian",
  n_bootstraps          = 3L,
  artificial_type       = NULL,
  cooperative_fusion    = TRUE,
  rho                   = 0,
  cooperation_selection = "validation",
  cooperation_selector  = "lambda.1se"
)
#> Error:
#> ! `cooperation_selector = 'lambda.1se'` is only available when `cooperation_selection = 'cv'`.
```

### Constraint 2 - Cox and Poisson families are not supported in native v1

Native cooperative fusion in `stablr` currently supports **gaussian**
and **binomial** families only. Poisson and Cox cooperative paths are
deferred to a future release:

``` r

library(survival)

set.seed(42)
n <- 20L
ids <- paste0("s", seq_len(n))
x_a <- matrix(rnorm(n * 4), nrow = n,
               dimnames = list(ids, paste0("a", seq_len(4))))
x_b <- matrix(rnorm(n * 4), nrow = n,
               dimnames = list(ids, paste0("b", seq_len(4))))
y_surv <- Surv(time = rexp(n), event = rbinom(n, 1L, 0.7))
rownames(y_surv) <- ids

stabl_multiomic_train_validate(
  x_train_list          = list(omic_a = x_a, omic_b = x_b),
  y_train               = y_surv,
  lambda_grid           = data.frame(lambda = c(0.2, 0.1)),
  family                = "cox",
  n_bootstraps          = 2L,
  artificial_type       = NULL,
  cooperative_fusion    = TRUE,
  rho                   = 0
)
#> Error:
#> ! `cooperative_fusion = TRUE` currently supports family = 'gaussian' or 'binomial'.
```

------------------------------------------------------------------------

## 6. Built-in cooperative engine

Cooperative fusion is built into `stablr`; no separate `multiview`
install is required. The engine is vendored from CRAN `multiview` v1.0
(GPL-2; see `inst/COPYING.cooperative`) and is enabled with
`cooperative_fusion = TRUE`.

Setting `cooperative_fusion = FALSE`, the default, runs only the
per-omic STABL paths and never invokes the cooperative engine.

------------------------------------------------------------------------

## 7. Comparing cooperative, early, and late fusion

The most useful comparison is often side by side. Run all three
strategies in one call, then compare the selected features:

``` r

all_fit <- stabl_multiomic_train_validate(
  x_train_list          = ool_train$x_list,
  y_train               = ool_train$y,
  lambda_grid           = lambda_list,
  x_valid_list          = ool_valid$x_list,
  y_valid               = ool_valid$y,
  family                = "gaussian",
  n_bootstraps          = 30L,
  artificial_type       = "random_permutation",
  random_state          = 42L,
  early_fusion          = TRUE,
  late_fusion           = TRUE,
  n_iter_lf             = 250L,
  cooperative_fusion    = TRUE,
  rho                   = c(0, 0.25, 0.5),
  cooperation_selection = "cv",
  cooperation_selector  = "lambda.min",
  cooperation_nfolds    = 3L
)
```

``` r

# print() now shows all four branches
all_fit
#> <stabl_multiomic_fit>
#>   Omics:           2 (cytof, proteomics)
#>   Per-omic selected features:
#>     cytof: 0
#>     proteomics: 2
#>   Validation data:  yes 
#>   Early fusion:     yes (2 features selected) 
#>   Late fusion:      yes (score = 0.1867) 
#>   Cooperative fusion:
#>     selection:     cv
#>     rho (chosen):  0
#>     selector:      lambda.min
#>     type.measure:  mse
#>     score:         812.1337
#>     selected feats: 26 (across views)
```

Compare feature sets across strategies:

``` r

ef_features  <- all_fit$early_fusion$selected_features
cf_features  <- unique(unlist(all_fit$cooperative_fusion$selected_features,
                              use.names = FALSE))
per_omic_features <- unique(unlist(all_fit$selected_features,
                                   use.names = FALSE))

cat("Early fusion selected:      ", length(ef_features), "features\n")
#> Early fusion selected:       2 features
cat("Cooperative fusion selected:", length(cf_features), "features\n")
#> Cooperative fusion selected: 26 features
cat("Per-omic union selected:    ", length(per_omic_features), "features\n")
#> Per-omic union selected:     2 features

cat("\nIn cooperative but not early fusion:\n")
#> 
#> In cooperative but not early fusion:
print(setdiff(cf_features, ef_features))
#>  [1] "Granulocytes_S6_unstim" "Bcells_S6_unstim"       "MDSCs_CREB_unstim"     
#>  [4] "DCs_CREB_unstim"        "pDCs_STAT3_unstim"      "pDCs_ERK_unstim"       
#>  [7] "intMCs_CREB_unstim"     "Mcl.1"                  "CD59"                  
#> [10] "S100A6"                 "SECTM1"                 "CSH"                   
#> [13] "NEGR1"                  "FSTL1"                  "Osteopontin"           
#> [16] "SMOC1"                  "protein.Z.inhibitor"    "ISLR2"                 
#> [19] "HXK2"                   "SEM5A"                  "ApoM"                  
#> [22] "AMGO2"                  "C1QR1"                  "Myostatin"

cat("\nIn early fusion but not cooperative:\n")
#> 
#> In early fusion but not cooperative:
print(setdiff(ef_features, cf_features))
#> character(0)
```

**Interpreting differences:** early fusion sees all features
simultaneously in one penalised model; cooperative fusion fits per-view
paths that are encouraged to agree. Features unique to cooperative
fusion may carry cross-view signal that is diluted in the concatenated
early-fusion space.

### Relation to mixOmics DIABLO

Cooperative fusion is the closest `stablr` analogue to multiblock
integration, but it is **not** equivalent to DIABLO:

|  | Cooperative fusion (`rho`) | DIABLO |
|----|----|----|
| Optimisation | Penalised regression paths per view, coupled by cooperation penalty | Sparse PLS components with block design matrix |
| Output | Per-view coefficient paths and selected features | Latent components + block loadings |
| Tuning | `rho`, `lambda`, `cooperation_nfolds` | `ncomp`, `keepX`, design matrix entries |
| Families (v1) | gaussian, binomial | Classification (sPLS-DA) |
| Stability / FDP | STABL bootstrap + knockoff layer | Component sparsity only |

Use cooperative fusion when you want **stable sparse features per omic**
with cross-view agreement. Use DIABLO when you want **joint latent
factors** for classification and mixOmics visualisation. Head-to-head
TCGA benchmarks that include early fusion (and future cooperative arms)
are described in `inst/analysis/BENCHMARK_PLAN.md`.

------------------------------------------------------------------------

## 8. Cross-validation with cooperative fusion

Use
[`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
when no predefined validation split is available. Cooperative fusion
parameters propagate to every outer fold. This is a realistic analysis
pattern, but it multiplies cost by the number of outer folds, so the
code is shown but not evaluated during vignette builds:

``` r

# Long-running vignette branch; exercised by test-multiomic-workflows.R.
cv_fit <- stabl_multiomic_cv(
  x_list                = ool_train$x_list,
  y                     = ool_train$y,
  lambda_grid           = lambda_list,
  v                     = 3L,
  family                = "gaussian",
  n_bootstraps          = 20L,
  artificial_type       = "random_permutation",
  random_state          = 42L,
  cooperative_fusion    = TRUE,
  rho                   = c(0, 0.3),
  cooperation_selection = "cv",
  cooperation_selector  = "lambda.min",
  cooperation_nfolds    = 3L
)

cv_fit
```

The `$diagnostics` data frame includes cooperative-specific columns
(`cooperative_rho`, `cooperative_lambda`, `cooperative_score`, etc.),
making it possible to inspect how `rho` was selected in each outer fold:

``` r

# Output contract exercised by test-multiomic-workflows.R.
print(cv_fit$diagnostics[, c("fold", "omic", "n_selected",
                              "cooperative_rho", "cooperative_score")])
```

------------------------------------------------------------------------

## 9. Choosing `rho` and design guidance

The useful `rho` grid is usually modest. Treat it as a scientific
sensitivity analysis: did the views need to cooperate, and if so, how
strongly?

| Decision | Guidance |
|----|----|
| **Grid density** | Start with `c(0, 0.1, 0.3, 0.5)`. Add finer points only if the CV surface has a clear non-zero minimum. |
| **Grid upper bound** | Values \> 0.5 can over-penalise disagreement between views; rarely beneficial unless views are near-identical. |
| **`lambda.min` vs `lambda.1se`** | `lambda.1se` gives a sparser model at the cost of a higher CV error; prefer `lambda.min` when prediction is the goal. |
| **`cooperation_nfolds`** | 5 (default) is suitable for n \>= 50. Use 3 for smaller datasets. Must be \>= 3. |
| **`cooperation_selection = "validation"`** | Use when you have a true held-out set and want to tune `rho` and `lambda` together on it. |
| **Native v1 families** | Cooperative fusion supports `family = "gaussian"` and `"binomial"` only. |
| **Binomial tuning metrics** | Try `cooperation_type_measure = "auc"` or `"class"` when the outcome is binary. |
| **rho = 0 wins** | If the chosen `rho` is always 0, the chosen metric prefers independent fits under this grid and fold design. Treat that as a prompt for sensitivity analysis, not proof that the omics share no useful structure. |

------------------------------------------------------------------------

## Next steps

- See
  [`vignette("stablr-multiomic")`](https://mdmanurung.github.io/stablr/articles/stablr-multiomic.md)
  for per-omic, early, and late fusion on the same OOL dataset.
- See
  [`vignette("stablr-intro")`](https://mdmanurung.github.io/stablr/articles/stablr-intro.md)
  for single-omic STABL and learner options.
- Pass `groups` to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  to prevent leakage in longitudinal designs; grouped bootstrapping and
  outer-fold construction both respect this argument.
- For binomial cooperative fusion, set `family = "binomial"` and
  consider `cooperation_type_measure = "auc"` or `"class"`.
- See
  [`vignette("stablr-advanced")`](https://mdmanurung.github.io/stablr/articles/stablr-advanced.md)
  for Cox survival STABL, knockoffs, grouped bootstrap, reproducibility
  metrics, export helpers, and outer CV.
- For mixOmics DIABLO comparison on TCGA, see
  [`vignette("stablr-tcga-nestedcv")`](https://mdmanurung.github.io/stablr/articles/stablr-tcga-nestedcv.md).
