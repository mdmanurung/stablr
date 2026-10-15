# Multi-Omic STABL Cross-Validation Workflow

Builds deterministic fold splits over a named multi-omic input, fits
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
on each training fold, and returns fold-wise selection diagnostics
together with selected train/validation matrices for downstream
inspection.

## Usage

``` r
stabl_multiomic_cv(
  x_list,
  y,
  lambda_grid,
  v = 5L,
  groups = NULL,
  base_learner = "lasso",
  family = "gaussian",
  n_bootstraps = 100L,
  artificial_type = "random_permutation",
  hard_threshold = NULL,
  stratify_bootstrap = FALSE,
  bootstrap_strata = NULL,
  l1_ratio = NULL,
  random_state = NULL,
  early_fusion = FALSE,
  late_fusion = FALSE,
  late_fusion_training = c("oof", "python_legacy"),
  late_fusion_nfolds = 5L,
  n_iter_lf = 10000L,
  cooperative_fusion = FALSE,
  rho = NULL,
  cooperation_selection = c("cv", "validation"),
  cooperation_selector = c("lambda.min", "lambda.1se"),
  cooperation_type_measure = "default",
  cooperation_nfolds = 5L,
  ...
)

# S3 method for class 'stabl_multiomic_cv'
print(x, ...)
```

## Arguments

- x_list:

  Named list of omic tables (`data.frame` or numeric matrix), each with
  row names as sample IDs.

- y:

  Named outcome vector.

- lambda_grid:

  Either a shared lambda grid (`data.frame` or `"auto"`) used for all
  omics, or a named list mapping each omic name to its own lambda grid.
  Prefer `"auto"` (each fold's grid is derived from its training data)
  or a grid computed from training data only: a grid computed on the
  full data leaks mild information into the cross-validated estimates.

- v:

  Number of folds. Must be at least 2.

- groups:

  Optional named grouping vector. When supplied, all samples in the same
  group are assigned to the same assessment fold.

- base_learner:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- family:

  Character; `glmnet` response family (`"gaussian"`, `"binomial"`,
  `"multinomial"`, or `"cox"`), or `"ordinal"` for
  `base_learner = "ordinal"`. Default: `"gaussian"`.

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

- bootstrap_strata:

  Optional categorical bootstrap stratification design forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  for each training fold.

- l1_ratio:

  Passed to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  when `lambda_grid = "auto"`.

- random_state:

  Optional integer seed used for deterministic fold assignment and
  forwarded to each per-fold
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  call.

- early_fusion:

  Logical. Forwarded to each per-fold
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  call.

- late_fusion:

  Logical. Forwarded to each per-fold
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
  call.

- late_fusion_training:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- late_fusion_nfolds:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- n_iter_lf:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- cooperative_fusion:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- rho:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- cooperation_selection:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- cooperation_selector:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- cooperation_type_measure:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- cooperation_nfolds:

  Forwarded to
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- ...:

  Additional arguments forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- x:

  A numeric matrix or `data.frame` (samples \\\times\\ features) with
  row names used as sample IDs.

## Value

A list with class `"stabl_multiomic_cv"` containing:

- `folds`:

  List of fold descriptors with `fold`, `train_ids`, and `valid_ids`.

- `fold_results`:

  Named list of per-fold `stabl_multiomic_fit` results.

- `diagnostics`:

  Data frame with one row per fold/omic and columns `fold`, `omic`,
  `n_selected`, `threshold`, and `max_score`.

- `fold_scheme`:

  `"group_disjoint"` when `groups` is supplied (all samples of a group
  share one assessment fold — pass study labels here for study-disjoint
  CV), otherwise `"sample_iid"`.

## Details

Optional early-fusion, late-fusion, and cooperative-fusion branches are
forwarded to each fold-specific
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
call.

## Methods (by generic)

- `print(stabl_multiomic_cv)`: Print a concise summary of a
  `stabl_multiomic_cv` object; invisibly returns `x`.

## Note

Batch harmonization (e.g. ComBat) fitted on the full dataset before
cross-validation leaks study-level information across folds and inflates
performance estimates. Fit harmonization within training folds, or
interpret all reported metrics as conditional on the harmonization step;
group-disjoint folds do not repair leakage from pooled preprocessing.

## Examples

``` r
set.seed(2L)
ids <- paste0("s", 1:12)
x_list <- list(
  rna = matrix(rnorm(12 * 4), 12, 4,
               dimnames = list(ids, paste0("g", 1:4))),
  protein = matrix(rnorm(12 * 3), 12, 3,
                   dimnames = list(ids, paste0("p", 1:3)))
)
y <- setNames(rnorm(12), ids)
cv <- stabl_multiomic_cv(
  x_list, y, lambda_grid = data.frame(lambda = 0.1), v = 2L,
  n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.5,
  random_state = 2L
)
cv$diagnostics
#>          fold    omic n_selected threshold max_score
#> Fold1.1 Fold1     rna          1       0.5         1
#> Fold1.2 Fold1 protein          1       0.5         1
#> Fold2.1 Fold2     rna          1       0.5         1
#> Fold2.2 Fold2 protein          1       0.5         1
```
