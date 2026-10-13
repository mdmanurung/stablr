# Stacked Generalization Over Per-Omic Predictions

Finds optimal omic weights by random search, matching the Python
`stabl_multiomic_stacked` algorithm from
`stabl/stacked_generalization.py`. Missing values in `predictions` are
handled per-row: rows with all-NA predictions receive `NA` in the
stacked output.

## Usage

``` r
stabl_multiomic_stacked(
  predictions,
  y,
  task_type = c("binary", "regression", "multiclass"),
  n_iter = 10000L,
  random_state = NULL
)
```

## Arguments

- predictions:

  For `task_type = "binary"` or `"regression"`, a `data.frame` or
  numeric matrix with one column per omic and one row per sample. For
  `task_type = "multiclass"`, a named list of class-probability
  matrices/data frames, one per omic, with samples in rows and classes
  in columns.

- y:

  Numeric outcome vector. When `y` is named, prediction rows must have
  unique sample IDs with the same set of names and `y` is aligned to
  their order. Positional alignment is used only when both sides are
  unnamed. Mixed, duplicated, or mismatched identity states are errors.
  For `task_type = "binary"` values must be `0`/`1`; for
  `task_type = "multiclass"`, values are coerced to a factor with levels
  matching the probability columns, and every label must be present in
  those columns.

- task_type:

  `"binary"` (maximise AUC), `"regression"` (maximise R²), or
  `"multiclass"` (minimise multiclass log loss).

- n_iter:

  Number of random weight draws to try. Mirrors `n_iter` in the Python
  implementation (default `10000`). The multiomic workflow functions
  expose this same quantity as `n_iter_lf`; see
  [`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md).

- random_state:

  Optional integer seed for reproducibility. Saves and restores the RNG
  state on exit so caller state is unaffected.

## Value

A named list with three elements:

- `predictions`:

  `data.frame` with one column per omic plus a
  `"Stacked Gen. Predictions"` column of the best weighted average.

- `weights`:

  `data.frame` with one row per omic and a column `Associated_weight`
  containing the optimal weight.

- `score`:

  Best score achieved (AUC, R², or negative log loss).

## Details

For binary and regression tasks, candidate weighted predictions are
evaluated in chunks while preserving the scalar algorithm's RNG order,
strict `score > best_score` incumbent update, and first-maximum tie
policy. The multiclass path remains scalar because exact
probability-normalization and tie parity are prioritized for this
release.

## Examples

``` r
set.seed(1L)
# Simulate two-omic predictions for 20 binary samples
y    <- rep(c(0L, 1L), 10)
pred <- data.frame(
  omic_a = y + rnorm(20, sd = 0.5),
  omic_b = y + rnorm(20, sd = 1.5)
)
res <- stabl_multiomic_stacked(pred, y, task_type = "binary",
                               n_iter = 200L, random_state = 42L)
res$weights
#>        Associated_weight
#> omic_a          9.782264
#> omic_b          1.174874
res$score
#> [1] 0.93
```
