# Finalize a STABL selection into a portable predictive model

Takes a completed development result — a
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
object, or a
[`stabl_analyze()`](https://mdmanurung.github.io/stablr/reference/stabl_analyze.md)
result — and refits a downstream predictive engine on the selected panel
using the full development data. The result is a `"stabl_predictor"`
that can be saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and applied to new
data with predict() without carrying any training observations.

## Usage

``` r
stabl_finalize(
  development_result,
  x,
  y,
  preprocessing = NULL,
  event_class = NULL,
  random_state = NULL
)

# S3 method for class 'stabl_predictor'
predict(object, newdata, type = NULL, extra_columns = c("error", "drop"), ...)

# S3 method for class 'stabl_predictor'
print(x, ...)
```

## Arguments

- development_result:

  A `"stabl_fit"` from
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
  or a `"stabl_analysis"` from
  [`stabl_analyze()`](https://mdmanurung.github.io/stablr/reference/stabl_analyze.md).
  For a validated analysis, the object must have been produced with
  `final_refit = TRUE` (the full-development refit supplies the selected
  panel); the cross-validated performance summary is carried into the
  predictor's `development` field, explicitly labelled as distinct from
  the full-data refit stored here. For an exploratory analysis the
  full-data fit is used and the development summary is labelled apparent
  / post-selection.

- x:

  Development feature matrix (samples x features) with row and column
  names; the data on which the downstream engine and the preprocessing
  transform are fitted. Must contain every selected feature.

- y:

  Development outcome aligned to the rows of `x` (named, or length
  `nrow(x)`). Factor with exactly two levels for `family = "binomial"`,
  numeric for `family = "gaussian"`. The family is taken from the fitted
  object, not re-declared.

- preprocessing:

  Optional preprocessing specification: `NULL` (identity), a character
  vector of transforms for
  [`stabl_preprocess_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_fit.md)
  (`"zv"`, `"impute_median"`, `"center"`, `"scale"`), or an unprepped
  [`recipes::recipe()`](https://recipes.tidymodels.org/reference/recipe.html).
  Fitted here on the full development data — the correct choice for a
  final model, since no held-out claim is derived from this fit. The
  stored transform is applied unchanged at predict time.

- event_class:

  For binary outcomes, the level treated as the event for probability
  predictions; defaults to the second level of the outcome factor.

- random_state:

  Seed for the downstream engine's internal cross-validation (glmnet
  path). `NULL` leaves the ambient stream; set it for a
  bitwise-reproducible final model.

- object:

  A `"stabl_predictor"` from `stabl_finalize()`.

- newdata:

  Development-schema new data (see
  [`validate_newdata_schema()`](https://mdmanurung.github.io/stablr/reference/validate_newdata_schema.md)).

- type:

  Prediction type: `"response"` (numeric; the only type for gaussian),
  `"class"` (hard class; binomial default) or `"probability"` (event
  probability for the stored `event_class`; binomial only).

- extra_columns:

  Policy for unexpected columns; see
  [`validate_newdata_schema()`](https://mdmanurung.github.io/stablr/reference/validate_newdata_schema.md).

- ...:

  Ignored.

## Value

A `"stabl_predictor"` object with the packed `model`, `engine` tag,
`family`, outcome `levels` / `event_class`, `selected_features`,
`required_columns`, the fitted `preprocessing` transform, the `learner`
specification, `feature_provenance` (selection scores), a `development`
summary, the stablr version and the call. No training observations are
stored.

## Details

The downstream engine mirrors the nested-CV machinery: an empty panel or
degenerate outcome falls back to a deterministic intercept-only model
(training event prevalence / majority class for classification, training
mean for regression); a one-feature panel uses a glm (stored as bare
coefficients); two or more features use a seeded
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
fit. Only binomial and gaussian fits are supported; other families fail
clearly.

When `preprocessing` is `NULL`, new data must supply exactly the
selected panel. When a transform is fitted, new data must supply all
development columns, because the stored transform may depend on any of
them; the panel is selected after transformation.

## Functions

- `predict(stabl_predictor)`: Predict from a portable
  `"stabl_predictor"`. The stored preprocessing transform is applied
  unchanged (imputation and normalization parameters are never
  recomputed from new data) and columns are matched by name, so input
  column order does not matter.

- `print(stabl_predictor)`: Print a concise summary of a
  `"stabl_predictor"`; invisibly returns `x`.

## Examples

``` r
set.seed(2L)
ids <- paste0("s", 1:60)
x <- matrix(rnorm(60 * 6), 60, 6,
            dimnames = list(ids, paste0("f", 1:6)))
y <- setNames(2 * x[, 1] - x[, 2] + rnorm(60), ids)

fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.3, 0.2, 0.1)),
                 family = "gaussian", n_bootstraps = 20L, random_state = 2L)
final <- stabl_finalize(fit, x, y, preprocessing = c("center", "scale"),
                        random_state = 2L)
final
#> <stabl_predictor: gaussian >
#>   engine:      cvglmnet
#>   panel:       2 of 6 features
#>   selected:    f1, f2
#>   requires:    6 input column(s)
#>   preprocessing: center, scale
#>   development: stabl_fit (full-data refit; not cross-validated performance)
#>   stablr:      0.3.0
predict(final, x[1:3, ], type = "response")
#> [1]  0.1957291 -1.9545288  4.2163516
```
