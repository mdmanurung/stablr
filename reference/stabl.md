# STABL Stability Selection (Beginner Interface)

One-stop-shop stability selection for biomarker discovery. Give it a
table of features and an outcome; it runs STABL with statistically safe
defaults and returns an object that tells you, in plain language, which
features were selected and how reliable the selection is.

## Usage

``` r
stabl(
  x,
  y,
  family = "gaussian",
  base_learner = "lasso",
  n_bootstraps = 1000L,
  random_state = NULL,
  ...
)
```

## Arguments

- x:

  A `data.frame` or numeric matrix of predictors (samples \\\times\\
  features). A `data.frame` may contain the outcome as a column (see
  `y`). Row names become sample IDs when present; otherwise IDs are
  generated as `s1, s2, ...`.

- y:

  The outcome. One of: a vector (numeric, character, or factor) aligned
  to `x`'s rows (by names when both are named, otherwise by position);
  or, when `x` is a `data.frame`, a length-1 character naming the
  outcome column, which is then removed from the predictors.

- family:

  Character; `"gaussian"` (default), `"binomial"`, `"multinomial"`, or
  `"cox"`. Use `"binomial"` for two-class outcomes, `"multinomial"` for
  more.

- base_learner:

  Character; `"lasso"` (default), `"elastic_net"`, `"adaptive_lasso"`,
  `"sparse_group_lasso"`, or `"rf"`.

- n_bootstraps:

  Integer; bootstrap iterations per lambda. Default `1000L` (the
  published STABL default); lower it only for quick exploration.

- random_state:

  Integer or `NULL`; seed for reproducibility.

- ...:

  Additional arguments forwarded to
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  (the expert interface), e.g. `artificial_proportion`,
  `sample_fraction`, `sample_groups`, `knockoff_fallback`.

## Value

An object of class `c("stabl", "stabl_fit")`. Print it for a
plain-language answer; use
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`get_support()`](https://mdmanurung.github.io/stablr/reference/get_support.md),
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md),
and
[`get_fdrs()`](https://mdmanurung.github.io/stablr/reference/get_fdrs.md)
for details. All `stabl_fit` accessors and workflows apply unchanged.

## Details

This is the recommended entry point for new users. It wraps the expert
interface
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
with opinionated defaults:

- `artificial_type = "knockoff_equi"` — model-X knockoff artificial
  features, the pairing validated for FDP+ calibration (the package
  default `"random_permutation"` is retained on the expert interface for
  Python-parity workflows but is anti-conservative for multivariate
  learners).

- FDP+ calibration on (no `hard_threshold`), strict knockoff fallback
  (`knockoff_fallback = "stop"`).

- Automatic lambda grid (always: `stabl()` fixes `lambda_grid` to
  `"auto"`; supply a manual grid via
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  if needed).

Everything else — other learners, bootstrap design, thresholds, groups —
is one `...` away: unnamed arguments in `...` are forwarded to
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md),
which documents the full expert interface.

## See also

[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
for the expert interface;
[`get_feature_names()`](https://mdmanurung.github.io/stablr/reference/get_feature_names.md)
for the selected features;
[`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md)
for the calibration machinery.

## Examples

``` r
set.seed(1L)
df <- data.frame(
  age  = rnorm(40),
  m1   = rnorm(40), m2 = rnorm(40), m3 = rnorm(40),
  m4   = rnorm(40), m5 = rnorm(40)
)
df$y <- 0.8 * df$m1 + 0.8 * df$m2 + rnorm(40, sd = 0.5)
# \donttest{
fit <- stabl(df, y = "y", n_bootstraps = 50L, random_state = 1L)
fit                      # plain-language print
#> <stabl stability selection>
#>   outcome family: gaussian | base learner: lasso | bootstraps: 50
#>   5 feature(s) selected at FDP+ estimate 1.000:
#>   m1, m2, m3, m4, m5
#>   use summary() for details, plot() for stability paths,
#>   get_feature_names() / get_support() for the selection,
#>   get_fdrs() for the FDP+ curve.
get_feature_names(fit)   # selected features
#> [1] "m1" "m2" "m3" "m4" "m5"
# }
```
