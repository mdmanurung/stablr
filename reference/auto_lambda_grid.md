# Build a Data-Driven Lambda Grid

Runs `glmnet` on the full training data to obtain a calibrated penalty
sequence, returning the result as a pre-expanded `data.frame` suitable
for use as the `lambda_grid` argument of
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Usage

``` r
auto_lambda_grid(x, y, family = "gaussian", n_lambda = 30L, l1_ratio = NULL)
```

## Arguments

- x:

  Numeric matrix of predictors (samples by features).

- y:

  Outcome vector. For `"gaussian"` provide a numeric vector; for
  `"binomial"`/`"multinomial"` provide a factor or 0/1 integer vector;
  for `"cox"` provide a
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object.

- family:

  Character; `glmnet` response family (`"gaussian"`, `"binomial"`,
  `"multinomial"`, or `"cox"`). Default `"gaussian"`.

- n_lambda:

  Positive integer; desired number of lambda values per alpha level.
  Actual length may be slightly shorter if `glmnet` detects saturation
  early. Default `30L`.

- l1_ratio:

  Numeric scalar, numeric vector, or `NULL`. When `NULL` (default), a
  pure lasso path (alpha = 1) is used. When a vector is supplied (e.g.
  `c(0.5, 0.75, 1.0)`), one path is fitted per value and the results are
  combined; an `alpha` column is added to the output.

## Value

A `data.frame` with at least a `lambda` column. An `alpha` column is
included whenever `l1_ratio` is not `NULL`.

## Details

Choosing lambda values manually is error-prone: too large a lambda
selects nothing; too small a lambda selects everything. This function
delegates the sequence computation to `glmnet`'s own warm-start path
algorithm, which guarantees the grid spans from near-zero sparsity to
near-full sparsity for the given data, family, and alpha.

For elastic-net models, supplying a vector of `l1_ratio` values causes
the function to fit a separate path per alpha, row-bind the resulting
grids, and add an `alpha` column so that
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md)
can dispatch correctly. This mirrors `auto_mode_lambda_grid()` in the
Python reference implementation.

## See also

[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
which calls this automatically when `lambda_grid = "auto"`.

## Examples

``` r
set.seed(1L)
x <- matrix(rnorm(40 * 6), 40, 6,
             dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
y <- setNames(rnorm(40), rownames(x))
grid <- auto_lambda_grid(x, y, family = "gaussian", n_lambda = 10L)
head(grid)
#>        lambda
#> 1 0.256317038
#> 2 0.092115567
#> 3 0.033104618
#> 4 0.011897183
#> 5 0.004275626
#> 6 0.001536580
```
