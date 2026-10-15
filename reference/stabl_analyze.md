# STABL analysis with an explicit evidence mode

Runs STABL feature selection in one of two explicitly declared modes:
`"exploratory"` (full-data selection with descriptive, apparent
interpretation) or `"validated"` (binary classification with the whole
pipeline — preprocessing, selection, tuning and prediction — isolated
inside nested cross-validation training partitions).

## Usage

``` r
stabl_analyze(
  x,
  y,
  family = c("binomial", "gaussian"),
  mode,
  preprocessing = NULL,
  resampling = c("nested_cv"),
  event_class = NULL,
  outer_v = 5L,
  outer_repeats = 1L,
  inner_v = 5L,
  stratified = TRUE,
  study_labels = NULL,
  lambda_grid = "auto",
  n_bootstraps = 100L,
  artificial_type = "random_permutation",
  hard_threshold = NULL,
  metric = c("ber", "accuracy"),
  cor_method = c("pearson", "spearman"),
  cor_threshold = 0.7,
  n_repeats = 20L,
  loss = NULL,
  random_state = NULL,
  final_refit = FALSE,
  ...
)

# S3 method for class 'stabl_analysis'
print(x, ...)

# S3 method for class 'stabl_analysis'
summary(object, ...)
```

## Arguments

- x:

  Numeric feature matrix (samples × features) with row names (sample
  IDs) and column names (feature names).

- y:

  Outcome aligned to the rows of `x`. Named vector or factor whose names
  cover exactly the row names of `x`. For `family = "binomial"`, exactly
  two classes; for `family = "gaussian"`, numeric.

- family:

  `"binomial"` or `"gaussian"`. Validated mode currently supports
  `"binomial"` only.

- mode:

  **Required.** `"exploratory"` or `"validated"`. There is no default:
  the evidence mode must be chosen explicitly.

- preprocessing:

  Optional fold-safe preprocessing specification: `NULL` (identity), a
  character vector of transforms for
  [`stabl_preprocess_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit.md)
  (`"zv"`, `"impute_median"`, `"center"`, `"scale"`), or an unprepped
  [`recipes::recipe()`](https://recipes.tidymodels.org/reference/recipe.html).
  In validated mode the transform is fitted anew on every inner- and
  outer-training split; in exploratory mode it is fitted once on the
  full data (an apparent, fitting-data transform). Fitted
  `stabl_preprocessing` objects are rejected.

- resampling:

  Resampling scheme for validated mode; currently only `"nested_cv"`
  (via
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)).
  Ignored in exploratory mode.

- event_class:

  For binary outcomes, the level treated as the event for probabilities,
  odds ratios and probability metrics; defaults to the second level of
  the outcome factor.

- outer_v, outer_repeats, inner_v, stratified, study_labels,
  lambda_grid, n_bootstraps, artificial_type, hard_threshold, metric:

  Validated-mode resampling and selection settings, passed to
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md).

- cor_method, cor_threshold:

  Correlation method and \|correlation\| threshold for redundancy
  groups. In validated mode, groups are derived on each fold's
  transformed training data only.

- n_repeats:

  Number of permutation repeats for importance (per fold in validated
  mode).

- loss:

  Loss for permutation importance: `"squared_error"` (gaussian),
  `"brier"` or `"log_loss"` (binomial). Positive importance means
  `permuted_loss - baseline_loss > 0`; negative values are retained.

- random_state:

  Seed. Fold assignment, selection, and per-fold importance seeds are
  all derived from it deterministically.

- final_refit:

  Logical; validated mode only. When `TRUE`, an explicitly labelled
  full-development refit (selection + interpretation on all data) is
  added for future inference. It is a different model from the
  cross-validation models and is not part of performance estimation.

- ...:

  Additional arguments passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  (and, in validated mode, to
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)).

- object:

  A `"stabl_analysis"` object (for `summary`).

## Value

An object of class `"stabl_analysis"` with a `mode` field. Exploratory
mode returns the full-data `fit` (a `"stabl_fit"`), its `interpretation`
(a `"stabl_interpretation"`, evidence type
`"post_selection_descriptive"`), the fitted `preprocessing` transform,
and `provenance`. Validated mode returns the underlying
`"stabl_multiomic_nested_cv"` object in `nested`; fold-level held-out
`predictions` (one row per sample per repeat, with `event_probability`
for the named `event_class`); `performance` (pooled and per-repeat
accuracy, BER, macro-F1, MCC, Brier score, clipped log loss and AUROC);
fold-level `selected_features`; `heldout_importance` (per-feature and
per-group permutation importance computed with each fold's own predictor
on that fold's held-out rows, with correlation groups derived from
training data only); and optionally `final_refit`. Validated evidence is
labelled `"cross_validated"` with evaluation scope
`"held_out_assessment"`.

## Functions

- `print(stabl_analysis)`: Print a concise summary of a
  `"stabl_analysis"` object; invisibly returns `x`.

- `summary(stabl_analysis)`: Return the key evidence table invisibly:
  the interpretation evidence for exploratory mode, or the held-out
  performance and importance summary for validated mode.

## Examples

``` r
set.seed(2L)
ids <- paste0("s", 1:80)
x <- matrix(rnorm(80 * 8), 80, 8,
            dimnames = list(ids, paste0("f", 1:8)))
y <- setNames(factor(rep(c("A", "B"), each = 40L)), ids)
x[, 1] <- x[, 1] + 1.5 * as.numeric(y == "B")
x[, 2] <- x[, 2] - 1.2 * as.numeric(y == "B")

exploratory <- suppressWarnings(stabl_analyze(
  x, y, family = "binomial", mode = "exploratory",
  lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
  n_bootstraps = 10L, random_state = 2L, n_repeats = 3L
))
summary(exploratory)
#> <stabl_analysis: exploratory >
#>   family:      binomial
#>   event class: B
#>   preprocessing: identity
#>   selected:     3 of 8 features
#>   evidence:    post_selection_descriptive (fitting_data; apparent)
#> 
#> Selected features (post-selection descriptive):
#>   feature stability coefficient effect_direction correlation_group
#> 1      f1       1.0   1.0698420         positive                 1
#> 2      f2       1.0  -0.8424217         negative                 2
#> 7      f7       0.7   0.5137596         positive                 3
#>   permutation_importance
#> 1             0.10213965
#> 2             0.05213274
#> 7             0.01112795

validated <- suppressWarnings(stabl_analyze(
  x, y, family = "binomial", mode = "validated",
  preprocessing = c("center", "scale"),
  lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
  outer_v = 3L, inner_v = 3L, n_bootstraps = 10L,
  random_state = 2L, n_repeats = 3L
))
validated$performance$pooled
#>    n n_events accuracy balanced_error_rate  macro_f1       mcc     brier
#> 1 80       40      0.7                 0.3 0.6992481 0.4020151 0.1853578
#>   log_loss    auroc
#> 1 0.549004 0.811875
```
