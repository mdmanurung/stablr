# Cooperative Fusion for Multi-Omic Biomarker Discovery

## Contents

- 1 Overview
- 2 Learning objectives
- 3 Loading the data
- 4 Building the lambda grids
- 5 Running cooperative fusion
- 6 Navigating the cooperative fusion element
- 7 Policy constraints
- 8 Built-in cooperative engine
- 9 Comparing cooperative, early, and late fusion
- 10 Cross-validation with cooperative fusion
- 11 Choosing rho and design guidance
- 12 Next steps

## 1 Overview

Early fusion asks all features to compete in one wide matrix, and late
fusion lets each omic make its own prediction before the answers are
combined. Cooperative fusion addresses the remaining case: the views are
fitted together without erasing their identities.

In `stablr`, cooperative fusion is built in natively (vendored from the
[multiview](https://cran.r-project.org/package=multiview)
cooperative-learning engine). Each view keeps its own coefficient path,
but the paths are encouraged to agree through a cooperation strength
`rho`. The omics are neither fully pooled nor fully isolated. Native v1
supports **gaussian** and **binomial** families.

Table 1 summarises when each strategy is appropriate.

| Strategy | How it works | Best when |
|----|----|----|
| **Per-omic** | One STABL fit per omic, independently. | Baseline; or when views are likely orthogonal. |
| **Early fusion** | Concatenate all omics, one joint STABL fit. | Omics are tightly correlated; a single feature list is desired. |
| **Late fusion** | Per-omic fits, then learn a weighted blend of predictions. | Each omic may predict independently; omic weights are desired. |
| **Cooperative fusion** | Joint multiview fit with cooperation strength `rho`; each view has its own coefficient path. | Omics share some signal structure; cross-view feature selection should be coherent. |

The cooperation strength `rho` controls how much the views are
encouraged to share signal. At `rho = 0`, cooperative fusion reduces to
independent fits. As `rho` grows, the views are penalised toward a
common prediction. `stablr` tunes `rho` over a user-supplied grid.

## 2 Learning objectives

The vignette covers the following topics.

1.  What `rho` controls in cooperative fusion and how to choose a grid.
2.  Running
    [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
    with `cooperative_fusion = TRUE` and navigating the
    `$cooperative_fusion` result element.
3.  Reading the [`print()`](https://rdrr.io/r/base/print.html) output
    for cooperative fusion details.
4.  Which `cooperation_selection` and `cooperation_selector`
    combinations are supported and why.
5.  Comparing cooperative fusion feature selections to early and late
    fusion on the same dataset.
6.  Using
    [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
    to evaluate cooperative fusion across outer folds.

## 3 Loading the data

The same Onset of Labor (OOL) dataset as
[`vignette("stablr-multiomic")`](https://mdmanurung.github.io/stablr/articles/stablr-multiomic.md)
is used: CyTOF immune profiling and Proteomics measured on the same
women, with days-to-onset-of-labor (DOS) as the continuous outcome. The
examples use compact grids so that the vignette remains buildable; the
bootstrap count and the `rho` grid should be enlarged for a final
analysis.

The package is loaded as follows.

``` r

library(stablr)
```

The training and validation splits are then read, and the dimensions of
each omic layer are inspected.

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

The output reports the dimensions of the CyTOF and proteomics matrices.
The two views contain the same samples but different numbers of
features.

## 4 Building the lambda grids

Cooperative fusion still needs a lambda grid for each omic. One grid per
view is built with
[`auto_lambda_grid()`](https://mdmanurung.github.io/stablr/reference/auto_lambda_grid.md),
and the grids are passed as a named list. Twelve lambda values are
enough to demonstrate the workflow without making each cooperative fit
too expensive.

The grids are constructed as follows.

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

Both views return 12 lambda values.

## 5 Running cooperative fusion

Cooperative fusion is activated with `cooperative_fusion = TRUE` and a
`rho` grid. `stablr` fits one multiview model per `rho` value using
shared inner CV folds, then picks the best `rho` by the chosen
`cooperation_type_measure`. This is **sequential** tuning over `rho`
(one full CV path per value), not a single joint CV over `(rho, lambda)`
pairs.

The key arguments are listed in the table below.

| Argument | Default | Purpose |
|----|----|----|
| `rho` | `0` | Scalar or vector of cooperation strengths to try. |
| `cooperation_selection` | `"cv"` | `"cv"` (inner CV) or `"validation"` (held-out set). |
| `cooperation_selector` | `"lambda.min"` | Which lambda to pick: `"lambda.min"` or `"lambda.1se"` (CV only). |
| `cooperation_type_measure` | `"default"` | Tuning metric; defaults to `"mse"` for Gaussian. |
| `cooperation_nfolds` | `5L` | Number of inner CV folds when `cooperation_selection = "cv"`. |

The cooperative fit is run with three `rho` values and inner
cross-validation.

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

The [`print()`](https://rdrr.io/r/base/print.html) method shows every
active branch. When cooperative fusion is present, a
`Cooperative fusion:` block appears at the bottom.

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
#>     rho (chosen):  0.5
#>     selector:      lambda.min
#>     type.measure:  mse
#>     score:         669.5853
#>     selected feats: 90 (across views)
```

The block reports the following fields.

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

## 6 Navigating the cooperative fusion element

The `cooperative_fusion` element contains the selected model and the
tuning surface that led to it.

The chosen hyperparameters are printed as follows.

``` r

cf <- cf_fit$cooperative_fusion

# Chosen hyperparameters
cat("Chosen rho:     ", cf$rho, "\n")
#> Chosen rho:      0.5
cat("Chosen lambda:  ", round(cf$selected_lambda, 4), "\n")
#> Chosen lambda:   0.1795
cat("Selector:       ", cf$selector, "\n")
#> Selector:        lambda.min
cat("CV score (MSE): ", round(cf$score, 4), "\n")
#> CV score (MSE):  669.5853
```

The output reports the chosen `rho`, the chosen lambda, the selector
rule, and the CV score (MSE) of the selected model.

### 6.1 Selected features per view

`$selected_features` is a named list with one character vector per omic.

The selected features are printed for both views.

``` r

cat("CyTOF features selected by cooperative fusion:\n")
#> CyTOF features selected by cooperative fusion:
print(cf$selected_features$cytof)
#>  [1] "Granulocytes_CREB_unstim"     "Granulocytes_STAT5_unstim"   
#>  [3] "Granulocytes_p38_unstim"      "Granulocytes_STAT3_unstim"   
#>  [5] "Granulocytes_S6_unstim"       "Granulocytes_MAPKAPK2_unstim"
#>  [7] "Bcells_S6_unstim"             "Bcells_IkB_unstim"           
#>  [9] "Bcells_MAPKAPK2_unstim"       "Bcells_NFkB_unstim"          
#> [11] "cMCs_MAPKAPK2_unstim"         "CCR2poscMCs_IkB_unstim"      
#> [13] "MDSCs_CREB_unstim"            "MDSCs_p38_unstim"            
#> [15] "MDSCs_S6_unstim"              "DCs_STAT5_unstim"            
#> [17] "DCs_IkB_unstim"               "DCs_MAPKAPK2_unstim"         
#> [19] "mDCs_STAT3_unstim"            "mDCs_NFkB_unstim"            
#> [21] "pDCs_STAT5_unstim"            "pDCs_STAT3_unstim"           
#> [23] "pDCs_NFkB_unstim"             "pDCs_ERK_unstim"             
#> [25] "intMCs_CREB_unstim"           "intMCs_S6_unstim"            
#> [27] "intMCs_IkB_unstim"            "intMCs_MAPKAPK2_unstim"      
#> [29] "intMCs_NFkB_unstim"           "intMCs_ERK_unstim"           
#> [31] "CCR2posintMCs_CREB_unstim"

cat("\nProteomics features selected by cooperative fusion:\n")
#> 
#> Proteomics features selected by cooperative fusion:
print(cf$selected_features$proteomics)
#>  [1] "NSE"                        "PIAS4"                     
#>  [3] "IRF1"                       "c.Jun"                     
#>  [5] "Mcl.1"                      "c.Myc"                     
#>  [7] "SMAD3"                      "PDGFRA"                    
#>  [9] "IL.12"                      "STAT1"                     
#> [11] "Osteocalcin"                "IL.5"                      
#> [13] "GPDA"                       "HEMK2"                     
#> [15] "YKL.40"                     "Apo.L1"                    
#> [17] "CD59"                       "FABPL"                     
#> [19] "GDF.11"                     "BTC"                       
#> [21] "HIF.1a"                     "S100A6"                    
#> [23] "RSPO3"                      "PSP"                       
#> [25] "VEGF.D"                     "CSH"                       
#> [27] "EFNB1"                      "SNP25"                     
#> [29] "NEGR1"                      "BCL6"                      
#> [31] "FSTL1"                      "Osteopontin"               
#> [33] "CD177"                      "CHKB"                      
#> [35] "SMOC1"                      "protein.Z.inhibitor"       
#> [37] "FLRT2"                      "ISLR2"                     
#> [39] "DSC2"                       "LDLR"                      
#> [41] "HXK2"                       "HXK1"                      
#> [43] "SEM5A"                      "PIANP"                     
#> [45] "S100A4"                     "RNF43"                     
#> [47] "ZNRF3"                      "ApoM"                      
#> [49] "IFNA7"                      "AMGO2"                     
#> [51] "C1QR1"                      "H2B2E"                     
#> [53] "H2A3"                       "H31"                       
#> [55] "IL.1F8"                     "IL.1F6"                    
#> [57] "Ephrin.A3"                  "X14.3.3.protein.beta.alpha"
#> [59] "Myostatin"
```

The two lists give the CyTOF and proteomics features selected by the
cooperative model. Features that appear in both lists carry cross-view
support.

### 6.2 Predictions

`$train_predictions` and `$valid_predictions` hold the multiview
predictions from the chosen model.

The first six predictions from each split are printed.

``` r

cat("First 6 training predictions:\n")
#> First 6 training predictions:
print(head(cf$train_predictions))
#>  001_26_A  001_33_B  001_35_C  003_25_A  003_30_B  003_32_C 
#> -87.57170 -36.82602 -30.89577 -53.20318 -29.72827 -27.61738

cat("\nFirst 6 validation predictions:\n")
#> 
#> First 6 validation predictions:
print(head(cf$valid_predictions))
#>   004_31_A   004_34_B   015_30_A   015_34_B   015_36_C   047_26_A 
#>  30.525422  32.567613 -23.970698 658.466946  -3.683433  -5.073059
```

The two vectors hold the predicted DOS values for the training samples
and for the held-out validation samples.

### 6.3 Tuning diagnostics

`$diagnostics` records the CV score at every `rho` in the grid, so the
tuning choice can be inspected rather than treated as a black box.

One row per `rho` value is printed.

``` r

# One row per rho value
print(cf$diagnostics[, c("rho", "lambda", "metric_value", "selected")])
#>    rho    lambda metric_value selected
#> 1 0.00 3.8680805     784.1156    FALSE
#> 2 0.25 0.1713800     700.7326    FALSE
#> 3 0.50 0.1795404     669.5853     TRUE
```

The `selected` column marks the winning `rho`. If the surface is flat,
with similar `metric_value` across `rho`, the omics are not gaining much
from cooperation and `rho = 0` is a defensible choice.

The same tuning surface is easier to inspect as a plot: a flat curve
means the choice of `rho` is not doing much work on these data.

``` r

ggplot2::ggplot(cf$diagnostics, ggplot2::aes(rho, metric_value)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(ggplot2::aes(shape = selected), size = 3) +
  ggplot2::scale_x_continuous(breaks = cf$diagnostics$rho) +
  ggplot2::labs(title = "Cooperative tuning surface",
                y = paste("CV", cf$type_measure), shape = "selected") +
  ggplot2::theme_minimal()
```

![Cooperative tuning surface: CV metric versus
rho](stablr-cooperative_files/figure-html/rho-surface-1.png)

## 7 Policy constraints

Some combinations are rejected on purpose. These guards prevent analyses
that look reasonable in code but do not have a valid tuning
interpretation.

### 7.1 `lambda.1se` requires CV mode

`lambda.1se` applies the “one-standard-error” rule and is only
well-defined when cross-validation is used to estimate the error
surface. Requesting it with `cooperation_selection = "validation"`
raises an error.

The following call raises the corresponding error.

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

The error confirms that `lambda.1se` cannot be combined with
validation-based selection.

### 7.2 Cox and Poisson families are not supported in native v1

Native cooperative fusion in `stablr` currently supports **gaussian**
and **binomial** families only. Poisson and Cox cooperative paths are
deferred to a future release.

A Cox fit raises the corresponding error.

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

The call fails because `family = "cox"` is not available for cooperative
fusion in native v1.

## 8 Built-in cooperative engine

Cooperative fusion is built into `stablr`; no separate `multiview`
install is required. The engine is vendored from CRAN `multiview` v1.0
(GPL-2; see `inst/COPYING.cooperative`) and is enabled with
`cooperative_fusion = TRUE`.

Setting `cooperative_fusion = FALSE`, the default, runs only the
per-omic STABL paths and never invokes the cooperative engine.

## 9 Comparing cooperative, early, and late fusion

The most useful comparison is side by side. All three strategies are run
in one call, and the selected features are then compared.

The combined fit is specified as follows.

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

The [`print()`](https://rdrr.io/r/base/print.html) method now shows all
four branches.

``` r

# print() now shows all four branches
all_fit
#> <stabl_multiomic_fit>
#>   Omics:           2 (cytof, proteomics)
#>   Per-omic selected features:
#>     cytof: 0
#>     proteomics: 2
#>   Validation data:  yes 
#>   Early fusion:     yes (7 features selected) 
#>   Late fusion:      yes (score = 0.2565) 
#>   Cooperative fusion:
#>     selection:     cv
#>     rho (chosen):  0.5
#>     selector:      lambda.min
#>     type.measure:  mse
#>     score:         669.5853
#>     selected feats: 90 (across views)
```

The output contains the per-omic fits, the early-fusion branch, the
late-fusion branch, and the `Cooperative fusion:` block.

The feature sets of the three strategies are compared as follows.

``` r

ef_features  <- all_fit$early_fusion$selected_features
cf_features  <- unique(unlist(all_fit$cooperative_fusion$selected_features,
                              use.names = FALSE))
per_omic_features <- unique(unlist(all_fit$selected_features,
                                   use.names = FALSE))

cat("Early fusion selected:      ", length(ef_features), "features\n")
#> Early fusion selected:       7 features
cat("Cooperative fusion selected:", length(cf_features), "features\n")
#> Cooperative fusion selected: 90 features
cat("Per-omic union selected:    ", length(per_omic_features), "features\n")
#> Per-omic union selected:     2 features

cat("\nIn cooperative but not early fusion:\n")
#> 
#> In cooperative but not early fusion:
print(setdiff(cf_features, ef_features))
#>  [1] "Granulocytes_CREB_unstim"     "Granulocytes_STAT5_unstim"   
#>  [3] "Granulocytes_p38_unstim"      "Granulocytes_STAT3_unstim"   
#>  [5] "Granulocytes_MAPKAPK2_unstim" "Bcells_S6_unstim"            
#>  [7] "Bcells_IkB_unstim"            "Bcells_MAPKAPK2_unstim"      
#>  [9] "Bcells_NFkB_unstim"           "cMCs_MAPKAPK2_unstim"        
#> [11] "CCR2poscMCs_IkB_unstim"       "MDSCs_CREB_unstim"           
#> [13] "MDSCs_p38_unstim"             "MDSCs_S6_unstim"             
#> [15] "DCs_STAT5_unstim"             "DCs_IkB_unstim"              
#> [17] "DCs_MAPKAPK2_unstim"          "mDCs_STAT3_unstim"           
#> [19] "mDCs_NFkB_unstim"             "pDCs_STAT5_unstim"           
#> [21] "pDCs_STAT3_unstim"            "pDCs_NFkB_unstim"            
#> [23] "pDCs_ERK_unstim"              "intMCs_CREB_unstim"          
#> [25] "intMCs_S6_unstim"             "intMCs_IkB_unstim"           
#> [27] "intMCs_MAPKAPK2_unstim"       "intMCs_NFkB_unstim"          
#> [29] "intMCs_ERK_unstim"            "CCR2posintMCs_CREB_unstim"   
#> [31] "NSE"                          "PIAS4"                       
#> [33] "IRF1"                         "c.Jun"                       
#> [35] "Mcl.1"                        "c.Myc"                       
#> [37] "SMAD3"                        "PDGFRA"                      
#> [39] "IL.12"                        "STAT1"                       
#> [41] "Osteocalcin"                  "IL.5"                        
#> [43] "GPDA"                         "YKL.40"                      
#> [45] "Apo.L1"                       "FABPL"                       
#> [47] "GDF.11"                       "BTC"                         
#> [49] "HIF.1a"                       "S100A6"                      
#> [51] "RSPO3"                        "PSP"                         
#> [53] "VEGF.D"                       "CSH"                         
#> [55] "EFNB1"                        "SNP25"                       
#> [57] "BCL6"                         "FSTL1"                       
#> [59] "Osteopontin"                  "CD177"                       
#> [61] "CHKB"                         "SMOC1"                       
#> [63] "protein.Z.inhibitor"          "FLRT2"                       
#> [65] "ISLR2"                        "DSC2"                        
#> [67] "LDLR"                         "HXK2"                        
#> [69] "HXK1"                         "SEM5A"                       
#> [71] "PIANP"                        "S100A4"                      
#> [73] "RNF43"                        "ZNRF3"                       
#> [75] "IFNA7"                        "AMGO2"                       
#> [77] "C1QR1"                        "H2A3"                        
#> [79] "H31"                          "IL.1F8"                      
#> [81] "IL.1F6"                       "Ephrin.A3"                   
#> [83] "X14.3.3.protein.beta.alpha"   "Myostatin"

cat("\nIn early fusion but not cooperative:\n")
#> 
#> In early fusion but not cooperative:
print(setdiff(ef_features, cf_features))
#> [1] "Granulocytes_IkB_unstim"
```

The output reports the number of features selected by each strategy,
followed by the features unique to cooperative fusion and the features
unique to early fusion. Early fusion sees all features simultaneously in
one penalised model, whereas cooperative fusion fits per-view paths that
are encouraged to agree. Features unique to cooperative fusion may carry
cross-view signal that is diluted in the concatenated early-fusion
space.

### 9.1 Relation to mixOmics DIABLO

Cooperative fusion is the closest `stablr` analogue to multiblock
integration, but it is **not** equivalent to DIABLO:

|  | Cooperative fusion (`rho`) | DIABLO |
|----|----|----|
| Optimisation | Penalised regression paths per view, coupled by cooperation penalty | Sparse PLS components with block design matrix |
| Output | Per-view coefficient paths and selected features | Latent components + block loadings |
| Tuning | `rho`, `lambda`, `cooperation_nfolds` | `ncomp`, `keepX`, design matrix entries |
| Families (v1) | gaussian, binomial | Classification (sPLS-DA) |
| Stability / FDP | STABL bootstrap + knockoff layer | Component sparsity only |

Cooperative fusion is appropriate when **stable sparse features per
omic** with cross-view agreement are the goal. DIABLO is appropriate
when **joint latent factors** for classification and mixOmics
visualisation are the goal. For a like-for-like classification
comparison, `stablr` also ships an sPLS-DA base learner
(`base_learner = "splsda"`, backed by mixOmics; see
[`?make_splsda_adapter`](https://mdmanurung.github.io/stablr/reference/make_splsda_adapter.md)),
which runs the same stability-selection machinery on top of the sPLS-DA
learner. Head-to-head TCGA benchmarks that include early fusion (and
future cooperative arms) are described in
`inst/analysis/BENCHMARK_PLAN.md`.

## 10 Cross-validation with cooperative fusion

[`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
is used when no predefined validation split is available. Cooperative
fusion parameters propagate to every outer fold. This is a realistic
analysis pattern, but it multiplies cost by the number of outer folds,
so the code is shown but not evaluated during vignette builds.

The outer-CV fit is specified as follows.

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
making it possible to inspect how `rho` was selected in each outer fold.

The cooperative columns are printed as follows.

``` r

# Output contract exercised by test-multiomic-workflows.R.
print(cv_fit$diagnostics[, c("fold", "omic", "n_selected",
                              "cooperative_rho", "cooperative_score")])
```

## 11 Choosing `rho` and design guidance

The useful `rho` grid is usually modest. The grid should be treated as a
scientific sensitivity analysis: whether the views needed to cooperate,
and if so, how strongly.

| Decision | Guidance |
|----|----|
| **Grid density** | Start with `c(0, 0.1, 0.3, 0.5)`. Add finer points only if the CV surface has a clear non-zero minimum. |
| **Grid upper bound** | Values \> 0.5 can over-penalise disagreement between views; rarely beneficial unless views are near-identical. |
| **`lambda.min` vs `lambda.1se`** | `lambda.1se` gives a sparser model at the cost of a higher CV error; prefer `lambda.min` when prediction is the goal. |
| **`cooperation_nfolds`** | 5 (default) is suitable for n \>= 50. Use 3 for smaller datasets. Must be \>= 3. |
| **`cooperation_selection = "validation"`** | Use when a true held-out set is available and `rho` and `lambda` are to be tuned together on it. |
| **Native v1 families** | Cooperative fusion supports `family = "gaussian"` and `"binomial"` only. |
| **Binomial tuning metrics** | Try `cooperation_type_measure = "auc"` or `"class"` when the outcome is binary. |
| **rho = 0 wins** | If the chosen `rho` is always 0, the chosen metric prefers independent fits under this grid and fold design. Treat that as a prompt for sensitivity analysis, not proof that the omics share no useful structure. |

## 12 Next steps

The following resources extend the workflows described above.

- See
  [`vignette("stablr-multiomic")`](https://mdmanurung.github.io/stablr/articles/stablr-multiomic.md)
  for per-omic, early, and late fusion on the same OOL dataset.
- See
  [`vignette("stablr-intro")`](https://mdmanurung.github.io/stablr/articles/stablr-intro.md)
  for single-omic STABL and learner options.
- Pass `groups_train` to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  (or `groups` to
  [`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md))
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
