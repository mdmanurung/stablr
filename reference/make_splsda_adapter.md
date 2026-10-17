# Build a sPLS-DA Learner Adapter

Returns a closure that fits a
[`mixOmics::splsda`](https://rdrr.io/pkg/mixOmics/man/splsda.html) model
on a single bootstrap subsample and returns a logical selection mask
over the features. The closure mirrors the adapter contract of
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md):
the STABL bootstrap loop supplies the subsample and one lambda-grid row,
and the adapter answers with the selected-feature indicator.

## Usage

``` r
make_splsda_adapter(
  family = "binomial",
  ncomp = 2L,
  bootstrap_threshold = .SPLSDA_LOAD_TOL
)
```

## Arguments

- family:

  Character; must be `"binomial"`. sPLS-DA is a classification learner;
  there is no gaussian, cox, or ordinal variant.

- ncomp:

  Positive integer; number of latent components to retain in every fit.
  Default `2L`.

- bootstrap_threshold:

  Positive numeric; loading-magnitude cutoff above which a feature is
  counted as selected. Default `1e-10` (any non-zero sparse loading
  counts).

## Value

A function with signature
`function(x, y, lambda_val) -> logical vector of length ncol(x)` for use
inside
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

## Details

The regularisation parameter swept by the STABL lambda grid is sPLS-DA's
`keepX` (the number of features kept per component). Grid rows are read
from a `keepX` column when present and fall back to the `lambda` column
otherwise, so
[`auto_splsda_keepx_grid()`](https://mdmanurung.github.io/stablr/reference/auto_splsda_keepx_grid.md)
output and hand-built `data.frame(lambda = ...)` grids both work.

A feature is counted as selected in a bootstrap when the magnitude of
its loading on any retained component exceeds `bootstrap_threshold`.
Sparse PLS loadings are exactly zero for non-retained features, so this
is an exact selection indicator, not a soft threshold.

## See also

[`auto_splsda_keepx_grid()`](https://mdmanurung.github.io/stablr/reference/auto_splsda_keepx_grid.md),
[`make_glmnet_adapter()`](https://mdmanurung.github.io/stablr/reference/make_glmnet_adapter.md),
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)

## Examples

``` r
if (requireNamespace("mixOmics", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(120L), 30L, 20L)
  y <- factor(rep(c("0", "1"), each = 15L))
  adapter <- make_splsda_adapter()
  mask <- adapter(x, y, data.frame(lambda = 5))
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
#> selected: 10 of 20 features
```
