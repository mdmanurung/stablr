# Multi-Omic STABL Nested Cross-Validation

Runs a repeated outer cross-validation loop for generalisation
assessment. Within each outer training split, candidate STABL workflows
are compared using an inner cross-validation loop. The selected
candidate is then refit on the full outer-training split and evaluated
once on the held-out outer fold.

## Usage

``` r
stabl_multiomic_nested_cv(
  x_list,
  y,
  candidates = NULL,
  lambda_grid = "auto",
  outer_v = 5L,
  outer_repeats = 1L,
  inner_v = 5L,
  stratified = TRUE,
  strata = NULL,
  strata_bins = 5L,
  metric = c("ber", "accuracy"),
  family = "multinomial",
  n_bootstraps = 100L,
  artificial_type = "random_permutation",
  hard_threshold = NULL,
  random_state = NULL,
  n_lambda = 30L,
  l1_ratio = NULL,
  workers = 1L,
  cv_workers = 1L,
  study_labels = NULL,
  preprocessing = NULL,
  event_class = NULL,
  return_probabilities = FALSE,
  keep_fold_predictors = FALSE,
  ...
)

# S3 method for class 'stabl_multiomic_nested_cv'
print(x, ...)
```

## Arguments

- x_list:

  Named list of omic matrices or data frames with identical sample row
  names.

- y:

  Named factor or character vector aligned to `x_list` rows.

- candidates:

  Optional list of candidate definitions. Each candidate is a list with
  `name` and `blocks`; `blocks` names one or more omics from `x_list`.
  When `NULL`, one candidate is created for each omic plus one
  early-fusion candidate using all omics.

- lambda_grid:

  Either `"auto"`, a single lambda-grid data frame, or a named list of
  lambda grids. Named entries can match omic names or candidate names.

- outer_v:

  Number of outer folds per repeat.

- outer_repeats:

  Number of repeated outer CV rounds.

- inner_v:

  Number of inner folds used for candidate selection.

- stratified:

  Logical. When `TRUE` (default), outer and inner folds preserve class
  proportions as closely as possible.

- strata:

  Optional named vector used to stratify folds instead of `y`.
  Categorical values are used directly. Numeric values are binned into
  quantile groups before fold assignment. Cannot be combined with
  `study_labels`.

- strata_bins:

  Number of quantile bins for numeric `strata`.

- metric:

  Candidate-selection metric, `"ber"` or `"accuracy"`.

- family:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).
  Defaults to `"multinomial"`.

- n_bootstraps:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- artificial_type:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- hard_threshold:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- random_state:

  Integer or `NULL`; top-level seed. Default: `NULL`.

- n_lambda:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  when `lambda_grid = "auto"`.

- l1_ratio:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  when `lambda_grid = "auto"`. Use this with
  `base_learner = "elastic_net"` to generate alpha-aware train-fold
  grids.

- workers:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  for bootstrap-level parallelism.

- cv_workers:

  Number of outer folds to evaluate in parallel. Uses
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix-like systems and falls back to sequential execution on Windows.

- study_labels:

  Optional named vector assigning each sample to a study (cohort, site,
  or batch); names must cover exactly the sample IDs. When `NULL`
  (default), folds are assigned sample-wise exactly as in prior releases
  (`fold_scheme = "sample_iid"`). When supplied, the CV unit becomes the
  study (`fold_scheme = "study_disjoint"`): outer folds assign whole
  studies (with `outer_v` equal to the number of studies this is
  leave-one-study-out), inner candidate-selection folds are
  study-disjoint within each outer training split, and
  `stratified = TRUE` is reinterpreted at the study level through a
  deterministic largest-first greedy pass that minimises
  class-proportion imbalance across folds (sample-level stratification
  is impossible under study disjointness). Study disjointness concerns
  the CV folds only: within-training bootstrap resampling is untouched,
  so `bootstrap_strata` (sampling design) and SGL `feature_groups`
  (feature grouping) remain orthogonal to `study_labels` (CV unit).

- preprocessing:

  Optional fold-safe preprocessing specification: `NULL` (default;
  identity, legacy behaviour), a character vector of transforms for
  [`stabl_preprocess_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit.md)
  (`"zv"`, `"impute_median"`, `"center"`, `"scale"`), or an unprepped
  [`recipes::recipe()`](https://recipes.tidymodels.org/reference/recipe.html)
  (requires the `recipes` package). The transform is fitted anew on
  every inner- and outer-training split and applied unchanged to the
  corresponding held-out rows; a fitted `stabl_preprocessing` object is
  rejected because it would leak information across folds. Recipe steps
  see the block-prefixed candidate columns (`block::feature`).

- event_class:

  For binary outcomes, the level treated as the event when
  `return_probabilities = TRUE`; defaults to the second level of `y`.

- return_probabilities:

  Logical; when `TRUE` and `y` is binary, each outer held-out prediction
  gains an `event_probability` column giving P(`event_class`) from the
  fold-trained selected-feature predictor, and the object records the
  `event_class`. Hard-class columns are unchanged.

- keep_fold_predictors:

  Logical; when `TRUE`, each outer fold's fitted downstream predictor
  (with a `predict_event_prob()` closure for binary outcomes) and its
  fold-fitted preprocessing transform are retained in `fold_results`,
  enabling held-out permutation importance with the same fold-trained
  predictor. Defaults to `FALSE` (legacy behaviour).

- ...:

  Additional arguments passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- x:

  A numeric matrix or `data.frame` (samples \\\times\\ features) with
  row names used as sample IDs.

## Value

An object of class `"stabl_multiomic_nested_cv"` containing fold
definitions, inner candidate diagnostics, outer held-out predictions,
selected features, and aggregate performance. The object also records
the sampling unit as `fold_scheme` (`"sample_iid"` or
`"study_disjoint"`) and the validated `study_labels`; under
`fold_scheme = "study_disjoint"` each outer/inner fold additionally
carries `study_ids` (per-fold study membership) and `class_counts`
(per-fold held-out class counts). When `study_labels` is supplied,
`performance_by_study` holds per-study outer held-out accuracy and
balanced error rate (pooled across folds and repeats; `NA` balanced
error rate for studies whose pooled held-out predictions contain a
single outcome class).

## Methods (by generic)

- `print(stabl_multiomic_nested_cv)`: Print a concise summary of a
  `stabl_multiomic_nested_cv` object; invisibly returns `x`.

## Note

Cooperative fusion (`cooperative_fusion = TRUE`) is not supported in
nested CV. Use
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
or
[`stabl_multiomic_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_cv.md)
for cooperative workflows. The held-out metrics estimate predictive
performance for the declared candidate workflow; they do not establish
dataset-wide FDP control. Late-fusion stacking is not a candidate type
in this helper.

Batch harmonization (e.g. ComBat) fitted on the full dataset before
cross-validation leaks study-level information across folds and inflates
performance estimates. Fit harmonization within training folds, or
interpret all reported metrics as conditional on the harmonization step;
study-disjoint folds do not repair leakage from pooled preprocessing.

## Examples

``` r
set.seed(3L)
ids <- paste0("s", 1:18)
y <- setNames(factor(rep(c("A", "B", "C"), each = 6L)), ids)
signal <- stats::model.matrix(~ y - 1)
x <- matrix(
  rnorm(18 * 6, sd = 0.3), 18, 6,
  dimnames = list(ids, paste0("g", 1:6))
)
x[, 1:3] <- x[, 1:3] + signal
nested <- suppressWarnings(stabl_multiomic_nested_cv(
  x_list = list(rna = x), y = y,
  candidates = list(list(name = "rna", blocks = "rna")),
  lambda_grid = data.frame(lambda = c(0.2, 0.1)),
  outer_v = 3L, inner_v = 3L, n_bootstraps = 3L,
  artificial_type = NULL, hard_threshold = 1e-9,
  sample_fraction = 1, random_state = 3L
))
nested$performance$balanced_error_rate
#> [1] 0.2222222

# Study-disjoint folds: each block of four samples is treated as one
# study; outer CV is leave-one-study-out.
set.seed(4L)
ids2 <- paste0("s", 1:20)
y2 <- setNames(factor(rep(c("A", "B"), times = 10L)), ids2)
x2 <- matrix(
  rnorm(20 * 6, sd = 0.3), 20, 6,
  dimnames = list(ids2, paste0("g", 1:6))
)
x2[, 1:2] <- x2[, 1:2] + stats::model.matrix(~ y2 - 1)
study_labels <- setNames(rep(paste0("study", 1:5), each = 4L), ids2)
nested2 <- suppressWarnings(stabl_multiomic_nested_cv(
  x_list = list(rna = x2), y = y2,
  candidates = list(list(name = "rna", blocks = "rna")),
  lambda_grid = data.frame(lambda = c(0.2, 0.1)),
  outer_v = 5L, inner_v = 2L, n_bootstraps = 3L,
  artificial_type = NULL, hard_threshold = 1e-9,
  sample_fraction = 1, random_state = 4L,
  study_labels = study_labels
))
nested2$fold_scheme
#> [1] "study_disjoint"
nested2$outer_folds[[1L]]$study_ids
#> [1] "study5"
```
