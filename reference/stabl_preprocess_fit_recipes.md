# Fit a preprocessing transform from a recipes recipe

Wraps
[`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html)
on training data only. Requires the `recipes` package (in `Suggests`,
not `Imports`).

## Usage

``` r
stabl_preprocess_fit_recipes(recipe, x_train, ...)
```

## Arguments

- recipe:

  A
  [`recipes::recipe()`](https://recipes.tidymodels.org/reference/recipe.html)
  object.

- x_train:

  Training data frame or matrix.

- ...:

  Passed to
  [`recipes::prep()`](https://recipes.tidymodels.org/reference/prep.html).

## Value

A `"stabl_preprocessing_recipes"` object wrapping the prepped recipe.

## Examples

``` r
if (requireNamespace("recipes", quietly = TRUE)) {
  x <- matrix(rnorm(40 * 4), 40, 4,
              dimnames = list(paste0("s", 1:40), paste0("f", 1:4)))
  rec <- recipes::recipe(~ ., data = as.data.frame(x))
  rec <- recipes::step_center(rec, recipes::all_predictors())
  pp <- stabl_preprocess_fit_recipes(rec, x)
  xt <- stabl_preprocess_apply(pp, x)
}
```
