# Interpret a fitted STABL selection

Builds a canonical feature-evidence table for a fitted
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
or [`stabl()`](https://mdmanurung.github.io/stablr/reference/stabl.md)
object **without re-running STABL**. For each feature the table records
the stability score, the selection flag, a conditional signed effect
from a separately fitted downstream model, a reproducible
correlation/redundancy group, and a seeded permutation importance.

## Usage

``` r
stabl_interpret(
  object,
  x,
  y,
  family = c("gaussian", "binomial"),
  assessment = FALSE,
  cor_method = c("pearson", "spearman"),
  cor_threshold = 0.7,
  n_repeats = 20L,
  loss = NULL,
  random_state = NULL,
  event_class = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `"stabl_fit"` or `"stabl"` object.

- x:

  Feature matrix (samples × features) with row and column names. Must
  contain every selected feature; extra columns are ignored.

- y:

  Outcome aligned to the rows of `x` (named vector or factor).

- family:

  `"gaussian"` or `"binomial"`; the downstream effect model.

- assessment:

  Logical; set `TRUE` only when `x`/`y` are genuinely held out from the
  fitting procedure. Changes the `evaluation_scope` label, nothing else.

- cor_method:

  `"pearson"` or `"spearman"`.

- cor_threshold:

  Numeric in `[0, 1]`; \|correlation\| at which two features join one
  redundancy group.

- n_repeats:

  Number of permutation repeats for importance.

- loss:

  Loss for permutation importance: `"squared_error"` (gaussian),
  `"brier"` or `"log_loss"` (binomial).

- random_state:

  Seed for permutation importance.

- event_class:

  For binomial, the level treated as the event; defaults to the second
  level of the outcome factor.

- ...:

  Reserved for future use.

## Value

An object of class `"stabl_interpretation"`: a list with `evidence` (the
canonical data.frame), `effects_model` (the downstream fit),
`correlation` (method, threshold, groups), `importance` (per-repeat
matrix and settings), `event_class`, `family`, and `provenance` (call,
seed, scopes). Missing or unestimable entries are explicit `NA` with a
reason column.

## Details

**Evidence labels.** All effects and importances are computed on the
data supplied as `x`/`y`. When that data is the fitting data, every row
is labelled `post_selection_descriptive` (effects) and `apparent`
(importance): the same samples informed both selection and the
downstream fit, so the quantities are optimistic and carry no
inferential meaning. Supplying genuinely held-out `x`/`y` relabels the
importance scope, but the caller is responsible for that data being
untouched by fitting — this function cannot verify it.

**What each column is.** `stability` is the max-over-lambda selection
frequency from
[`get_importances()`](https://mdmanurung.github.io/stablr/reference/get_importances.md);
it is not a coefficient magnitude or a prediction importance.
`coefficient` / `standardized_effect` / `effect_direction` come from a
downstream [`stats::glm()`](https://rdrr.io/r/stats/glm.html) fit on the
selected features only and are conditional, descriptive associations.
`correlation_group` is a reproducible high-\|correlation\| cluster; it
is increase in loss when the feature is permuted (positive = helpful);
negative values are retained.

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
y <- setNames(0.8 * x[, "f1"] + rnorm(60), rownames(x))
fit <- stabl_fit(x, y,
                 lambda_grid = data.frame(lambda = c(0.3, 0.1, 0.05)),
                 n_bootstraps = 10L, artificial_type = "random_permutation",
                 random_state = 1L)
insight <- stabl_interpret(fit, x, y, family = "gaussian", n_repeats = 3L)
summary(insight)
#>   feature selected stability selection_threshold coefficient
#> 1      f1     TRUE       1.0                 0.8   0.6703844
#> 2      f2     TRUE       1.0                 0.8  -0.2830241
#> 3      f6     TRUE       0.9                 0.8   0.2003799
#>   standardized_effect effect_direction  effect_measure correlation_group
#> 1           0.5732918         positive mean_difference                 1
#> 2          -0.2587687         negative mean_difference                 2
#> 3           0.1986708         positive mean_difference                 3
#>   permutation_importance importance_sd importance_n_repeats
#> 1             0.64144613    0.09960395                    3
#> 2             0.09023568    0.02097216                    3
#> 3             0.03652539    0.04748862                    3
#>                evidence_type evaluation_scope n_observations         model_id
#> 1 post_selection_descriptive     fitting_data             60 downstream_refit
#> 2 post_selection_descriptive     fitting_data             60 downstream_refit
#> 3 post_selection_descriptive     fitting_data             60 downstream_refit
#>   na_reason
#> 1      <NA>
#> 2      <NA>
#> 3      <NA>
```
