# Make Model-X Equicorrelated Knockoff Artificial Features

Generates **model-X equicorrelated** knockoff features via
`knockoff::create.gaussian(..., method = "equi")`. This matches the
`GaussianSampler(X, method='equicorrelated')` call used by the Python
STABL library, making it the parity-correct knockoff type for
cross-language comparisons. Inputs above 3 000 predictors fail closed
because a globally valid chunked construction has not been established.
Supported-dimensional constructor failures fall back to
random-permutation features.

## Usage

``` r
make_knockoff_equi_features(
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

  Optional integer seed.

## Value

Named list with elements `x_augmented`, `noise_col_indices`, and
`artificial_provenance`; see
[`make_rp_features()`](https://mdmanurung.github.io/stablr/reference/make_rp_features.md)
for details.
