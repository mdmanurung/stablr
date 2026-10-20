# Multi-Omic STABL Train/Validation Workflow

Fits
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
independently on each omic block from a named list, then returns
per-omic fitted objects and selected-feature matrices for downstream
composition. Optional early-fusion, late-fusion, and cooperative-fusion
branches are additive to the per-omic STABL results.

## Usage

``` r
stabl_multiomic_train_validate(
  x_train_list,
  y_train,
  lambda_grid,
  x_valid_list = NULL,
  y_valid = NULL,
  groups_train = NULL,
  base_learner = "lasso",
  family = "gaussian",
  n_bootstraps = 100L,
  artificial_type = "random_permutation",
  hard_threshold = NULL,
  stratify_bootstrap = FALSE,
  bootstrap_strata_train = NULL,
  l1_ratio = NULL,
  random_state = NULL,
  early_fusion = FALSE,
  late_fusion = FALSE,
  late_fusion_training = "oof",
  late_fusion_nfolds = 5L,
  n_iter_lf = 10000L,
  cooperative_fusion = FALSE,
  rho = NULL,
  cooperation_selection = "cv",
  cooperation_selector = c("lambda.min", "lambda.1se"),
  cooperation_type_measure = "default",
  cooperation_nfolds = 5L,
  ...
)

# S3 method for class 'stabl_multiomic_fit'
print(x, ...)
```

## Arguments

- x_train_list:

  Named list of training omic tables (`data.frame` or numeric matrix),
  each with row names as sample IDs.

- y_train:

  Named outcome vector for training samples.

- lambda_grid:

  Either a shared lambda grid (`data.frame` or `"auto"`) used for all
  omics, or a named list mapping each omic name to its own lambda grid.
  When `late_fusion = TRUE`, the same full-training grid is used inside
  the OOF stacking procedure, so a grid computed on the full data leaks
  mild information into the OOF stacking weights; prefer `"auto"` or a
  grid computed from training data only.

- x_valid_list:

  Optional named list of validation omic tables. When supplied, names
  must match `x_train_list` names.

- y_valid:

  Optional named outcome vector for validation samples. When supplied
  together with `x_valid_list`, sample alignment is validated.

- groups_train:

  Optional named grouping vector for training samples. When supplied,
  grouped bootstrap sampling is used in each per-omic fit.

- base_learner:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- family:

  Character; `glmnet` response family (`"gaussian"`, `"binomial"`,
  `"multinomial"`, or `"cox"`), or `"ordinal"` for
  `base_learner = "ordinal"`. `base_learner = "splsda"` requires
  `"binomial"`. Default: `"gaussian"`.

- n_bootstraps:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- artificial_type:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- hard_threshold:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- stratify_bootstrap:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).
  When `TRUE`, per-omic and early-fusion STABL fits draw
  class-stratified bootstrap subsamples.

- bootstrap_strata_train:

  Optional categorical bootstrap stratification design for training
  samples, forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- l1_ratio:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  when `lambda_grid = "auto"`. Use this with
  `base_learner = "elastic_net"` to generate alpha-aware auto grids.

- random_state:

  Integer or `NULL`; top-level seed. Default: `NULL`.

- early_fusion:

  Logical. When `TRUE`, a single
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  is run on the column-bound concatenation of all omic matrices in
  addition to the per-omic fits. Results are returned in the
  `early_fusion` field.

- late_fusion:

  Logical. When `TRUE`, a downstream predictor is fitted per omic on its
  selected features, and
  [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)
  combines the per-omic predictions. If no features are selected for an
  omic, late fusion falls back to class priors for multinomial tasks or
  the train-set mean for other tasks. The `task_type` (binary,
  regression, or multiclass) is inferred from `family`. Results are
  returned in the `late_fusion` field.

- late_fusion_training:

  How stacking weights are trained. Only `"oof"` (the default) is
  supported: STABL selection and the downstream model are refit inside
  every stacking fold. The historical in-sample `"python_legacy"`
  algorithm was removed from the public API in 1.0.0 because it trains
  the stacking weights on the same data used to fit the base models
  (optimistic); it survives only as the unexported
  `.late_fusion_multiomic_fit()` for Python-parity characterization
  tests.

- late_fusion_nfolds:

  Number of stacking folds used in OOF mode.

- n_iter_lf:

  Number of random weight draws passed to
  [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)
  during late fusion. Ignored when `late_fusion = FALSE`. Naming note:
  [`stabl_multiomic_stacked()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_stacked.md)
  calls this same quantity `n_iter`; the multiomic workflows expose it
  as `n_iter_lf` so that the `...` channel (forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md))
  stays free of name collisions. `n_iter` itself is not a
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  argument and is therefore rejected if passed through `...`.

- cooperative_fusion:

  Logical. When `TRUE`, fit a built-in cooperative learning branch
  (vendored multiview engine) in addition to the existing per-omic STABL
  fits. Native v1 supports `family = "gaussian"` and `"binomial"` only.

- rho:

  Numeric scalar or vector of non-negative cooperation strengths. When
  `NULL`, defaults to `0`.

- cooperation_selection:

  Character scalar. Only `"cv"` is supported: tuning over `rho` uses
  shared inner fold assignments. The `"validation"` option (tuning `rho`
  and `lambda` on the supplied validation set) was removed from the
  public API in 1.0.0 because it tunes on the assessment data and biases
  the performance estimate (F-CV-04); the branch survives only in
  unexported form for characterization tests.

- cooperation_selector:

  Character scalar. Selection rule for the cooperative `lambda`. Both
  `"lambda.min"` and `"lambda.1se"` are available with
  `cooperation_selection = "cv"`.

- cooperation_type_measure:

  Character scalar controlling the cooperative tuning metric. Supported
  values follow the active `family` and the multiview CV API.

- cooperation_nfolds:

  Number of inner folds used when `cooperation_selection = "cv"`.

- ...:

  Additional arguments forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- x:

  A numeric matrix or `data.frame` (samples \\\times\\ features) with
  row names used as sample IDs.

## Value

A named list with class `"stabl_multiomic_fit"` containing:

- `fits`:

  Named list of per-omic `stabl_fit` objects. When `late_fusion = TRUE`,
  the top-level `fits`, `selected_features`, `selected_train`, and
  `selected_valid` come from the late-fusion full refit (its
  `full_per_omic`, which is stripped from the returned `late_fusion`
  list), not from a separate per-omic pass, so they can differ from the
  `late_fusion = FALSE` result on identical input and seed.

- `selected_features`:

  Named list of selected feature names.

- `selected_train`:

  Named list of training matrices restricted to selected features
  (possibly 0-column).

- `selected_valid`:

  Named list of validation matrices restricted to selected features, or
  `NULL` when no validation input is provided.

- `early_fusion`:

  `NULL` when `early_fusion = FALSE`. Otherwise a list with `fit`,
  `selected_features`, `selected_train`, and `selected_valid` for the
  concatenated single-STABL run.

- `late_fusion`:

  `NULL` when `late_fusion = FALSE`. Otherwise a list with `weights`,
  `valid_predictions`, and OOF training results (`oof_predictions`,
  `oof_score`). (The removed `"python_legacy"` mode additionally
  provided `in_sample_predictions` and `in_sample_score`.)
  `oof_predictions` is a *base-model* out-of-fold matrix; the stacked
  `oof_score` is computed with stacking weights tuned on that same OOF
  matrix (maximum over `n_iter_lf` random weight draws), so it is an
  optimistic, in-sample-for-the-weight-step statistic and NOT a
  full-stack OOF prediction. Use outer cross-validation
  ([`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md))
  for an honest performance estimate. `train_predictions` and `score`
  remain exact 0.1.x compatibility aliases for the active mode-specific
  fields. `task_type` is always present; multinomial tasks also include
  `levels`, `log_loss`, and classification metrics.

- `cooperative_fusion`:

  Present only when `cooperative_fusion = TRUE`. A list containing the
  selected multiview fit, chosen `rho` and `lambda`, selected features
  per view, train/validation predictions, and tuning diagnostics.

## Functions

- `print(stabl_multiomic_fit)`: Print a concise summary of a
  `stabl_multiomic_fit` object; invisibly returns `x`.

## Examples

``` r
set.seed(1L)
ids <- paste0("s", 1:12)
x_list <- list(
  rna = matrix(rnorm(12 * 4), 12, 4,
               dimnames = list(ids, paste0("g", 1:4))),
  protein = matrix(rnorm(12 * 3), 12, 3,
                   dimnames = list(ids, paste0("p", 1:3)))
)
y <- setNames(rnorm(12), ids)
fit <- stabl_multiomic_train_validate(
  x_list, y, lambda_grid = data.frame(lambda = c(0.2, 0.1)),
  n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
  random_state = 1L
)
names(fit$fits)
#> [1] "rna"     "protein"
```
