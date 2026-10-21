# Fit a fold-safe preprocessing transform on training data

Learns preprocessing parameters from a training matrix only. The result
is a `"stabl_preprocessing"` object that
[`stabl_preprocess_apply()`](https://mdmanurung.github.io/stablr/reference/stabl_preprocess_apply.md)
applies to training or assessment data without refitting. When
`transforms` is an unprepped
[`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html),
the recipe is prepped on the training data only (wrapping
[`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html)).

## Usage

``` r
stabl_preprocess_fit(
  x_train,
  transforms = c("center", "scale"),
  zv_tol = 0,
  ...
)
```

## Arguments

- x_train:

  Numeric training matrix (samples × features) with column names.

- transforms:

  Character vector naming the transforms to fit, or an unprepped
  [`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html)
  object (requires the `recipes` package, in `Suggests`). Supported
  character values: `"zv"` (zero-variance filter), `"impute_median"`
  (median imputation), `"center"`, `"scale"`. Transforms are always
  applied in the fixed order `"zv"` -\> `"impute_median"` -\> `"center"`
  -\> `"scale"`, regardless of the order in `transforms`. `NULL` or
  `character(0)` fits an identity transform that leaves input unchanged.

- zv_tol:

  Non-negative numeric; a feature whose training variance is `<= zv_tol`
  is removed by `"zv"`. Ignored when `transforms` is a recipe.

- ...:

  Passed to
  [`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html)
  when `transforms` is a recipe; reserved for future use otherwise.

## Value

A `"stabl_preprocessing"` object recording the fitted parameters, the
retained feature names, and the feature mapping.

## Examples

``` r
x <- matrix(rnorm(40 * 5), 40, 5,
            dimnames = list(paste0("s", 1:40), paste0("f", 1:5)))
pp <- stabl_preprocess_fit(x, transforms = c("center", "scale"))
xt <- stabl_preprocess_apply(pp, x)

if (requireNamespace("recipes", quietly = TRUE)) {
  rec <- recipes::recipe(~ ., data = as.data.frame(x))
  rec <- recipes::step_center(rec, recipes::all_predictors())
  pp_rec <- stabl_preprocess_fit(x, transforms = rec)
  xt_rec <- stabl_preprocess_apply(pp_rec, x)
}
```
