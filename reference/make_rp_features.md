# Make Random-Permutation Artificial Features

Randomly selects `n_injected` columns from `x`, copies them, and
shuffles each copy independently. Mirrors the `"random_permutation"`
branch of `Stabl._make_artificial_features()` in the Python STABL
library.

## Usage

``` r
make_rp_features(x, n_injected)
```

## Arguments

- x:

  Numeric matrix of predictors (samples \\\times\\ features).

- n_injected:

  Integer; number of artificial columns to generate.

## Value

Named list:

- x_augmented:

  Original matrix with artificial columns appended.

- noise_col_indices:

  Integer vector (1-based) of original column indices selected as
  sources for the artificial block.

- artificial_provenance:

  List describing the requested generator, `actual_type` used by the
  selected artificial columns, all `actual_types`, and the complete
  ordered `fallback_history`.
