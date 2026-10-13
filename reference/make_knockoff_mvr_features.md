# Make Model-X MVR Knockoff Artificial Features

Generates **model-X MVR (minimum-variance-reconstructability)** knockoff
features. The S-matrix is solved via `solve_mvr()` (an R-native
coordinate- descent port of `knockpy.mrc._solve_mvr_ungrouped`), then
the knockoff sample is drawn with
`knockoff::create.gaussian(..., diag_s = S)`. This is a novel feature
exclusive to `stablr` — the Python STABL library does not implement MVR
knockoffs. Inputs above 3 000 predictors fail closed because global
exchangeability for a chunked construction has not been established. For
supported dimensions, MVR-solver failure falls back to equi;
`create.gaussian` failure falls back to random permutation.

## Usage

``` r
make_knockoff_mvr_features(
  x,
  n_injected,
  random_state = NULL,
  knockoff_fallback = c("stop", "allow")
)
```

## Arguments

- x:

  Numeric matrix of predictors (samples \\\times\\ features).

- n_injected:

  Integer; number of knockoff columns to select.

- random_state:

  Optional integer seed; passed to `solve_mvr()` for the
  coordinate-shuffle RNG.

## Value

Named list with elements `x_augmented`, `noise_col_indices`, and
`artificial_provenance`; see
[`make_rp_features()`](https://mdmanurung.github.io/stablr/reference/make_rp_features.md)
for details. The retained `artificial_provenance$mvr_chunking`
compatibility field always records that chunking was not applied.
