# Make Knockoff Artificial Features

Generates **fixed-X** knockoff features via
[`knockoff::create.fixed()`](https://rdrr.io/pkg/knockoff/man/create.fixed.html),
with up to 3 000 predictors. Larger inputs fail closed because the
former chunked construction did not establish a valid global knockoff
design. Falls back to random-permutation features when the
supported-dimensional knockoff constructor fails (e.g., rank-deficient
input).

## Usage

``` r
make_knockoff_features(
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
for details. `artificial_provenance$standardization` records the fixed-X
scale contract (`knockoff_rescaled_to_input_scale = TRUE` when the
fixed-X construction ran) and carries `standardized = FALSE` with NULL
`means`/`sds` when the random-permutation fallback was taken.

## Details

`create.fixed()` column-normalizes its input internally (scale only, no
centering) and returns knockoffs on that internal unit-norm scale. The
knockoffs are therefore rescaled by the input column norms so that the
shipped pair `(x, xk)` satisfies the fixed-X Gram condition on the
input's own scale: `t(xk) xk = t(x) x` and `t(x) xk = t(x) x - diag(s)`
with `s_j >= 0`. The normalization actually applied is recorded in
`artificial_provenance$standardization`.

## Caveat on bootstrap reuse

Fixed-X knockoffs are generated once on the full design, so their
exchangeability guarantee does not survive the row subsampling of the
bootstrap loop: a knockoff row is only exchangeable with the matching
original row of that same realized design. When correlation structure
matters (in particular with `sample_fraction < 1`), prefer the model-X
modes `"knockoff_equi"` or `"knockoff_mvr"`, whose knockoff rows are
i.i.d. draws from a fitted Gaussian model and remain valid under row
subsampling.
