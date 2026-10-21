# Dispatcher for Artificial Feature Generation

Selects and calls the appropriate artificial-feature generator based on
`type`, returning the augmented predictor matrix together with the
column indices of the injected noise block.

## Usage

``` r
make_artificial_features(
  x,
  n_injected,
  type,
  random_state = NULL,
  knockoff_fallback = c("stop", "allow")
)
```

## Arguments

- x:

  Numeric matrix of predictors (samples \\\times\\ features). Must have
  more columns than `n_injected` for random permutation; for knockoffs,
  constructor failure (e.g., rank-deficient input) is governed by
  `knockoff_fallback`: by default (`"stop"`) it is an error, and with
  `knockoff_fallback = "allow"` a fallback to random permutation is
  attempted. Fixed-X knockoffs that require row augmentation also fail
  (or fall back under `"allow"`) because augmented knockoff rows cannot
  be aligned with the original outcome. All knockoff modes fail closed
  above 3 000 predictors; random permutation has no such dimensional
  limit.

- n_injected:

  Positive integer; number of artificial columns to append. Typically
  `round(ncol(x) * artificial_proportion)` as computed in
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- type:

  Character string; one of `"random_permutation"`, `"knockoff"`,
  `"knockoff_equi"`, or `"knockoff_mvr"`.

- random_state:

  Optional integer; passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before any random
  operations for reproducibility. `NULL` leaves the RNG unchanged.
  Seeding happens exactly once in this dispatcher; downstream generators
  inherit the seeded RNG state and do not re-seed (audit M-5).

- knockoff_fallback:

  Character; the fallback policy when a knockoff constructor fails on a
  supported-dimensional input (audit A1, strict fallback policy).
  `"stop"` (the default) raises an error: a
  knockoff-to-random-permutation fallback silently changes the FDP+
  calibration target, because permutation features do not support FDP+
  calibration. `"allow"` restores the pre-0.2.2 behavior (warn and fall
  back to random-permutation features for the affected chunk; the
  `fallback_history` provenance records the event). The benign
  `knockoff_mvr` -\> `knockoff_equi` fallback is unaffected: both are
  model-X constructions with valid FDP+ calibration, so it remains a
  warning under either policy.

## Value

Named list with three elements:

- `x_augmented`:

  Numeric matrix of size (nrow(x)) \\\times\\ (ncol(x) + n_injected)
  with the artificial columns appended after the original features.

- `noise_col_indices`:

  Integer vector of length `n_injected` containing the 1-based indices
  into the **original** `x` columns (not into the artificial block) that
  identify which source features were used to build each artificial
  column. Used by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
  to look up sparse-group-lasso group memberships for the artificial
  block via `.append_noise_groups`.

- `artificial_provenance`:

  Additive metadata describing the requested artificial-feature type,
  `actual_type` used for scientific classification, all generator modes
  used by the selected artificial columns, and the ordered structured
  `fallback_history`. Per-chunk counts and reasons remain available for
  0.1.x compatibility.

## Details

Injecting artificial features supports STABL's FDP+ diagnostic
calibration: by mixing known-noise columns into the predictor matrix
alongside real features, STABL can empirically estimate how often a
variable of pure noise is selected at a given stability threshold. This
observed noise-selection rate drives the FDP+ estimate computed in
[`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md).
The current rule chooses its minimum and does not imply universal
error-rate control.

Four noise strategies are supported:

- `"random_permutation"`:

  Copies `n_injected` randomly chosen real columns and shuffles each
  copy independently, breaking all signal while preserving marginal
  distributions. Fast and broadly applicable.

- `"knockoff"`:

  Generates **fixed-X** knockoffs via
  [`knockoff::create.fixed()`](https://rdrr.io/pkg/knockoff/man/create.fixed.html),
  which preserve the second-moment (Gram) structure of the original
  design on the input's own scale: `t(xk) xk = t(x) x` and
  `t(x) xk = t(x) x - diag(s)` with `s >= 0` (0.2.2 scale contract).
  Kept for backward compatibility; note that fixed-X knockoffs are
  generated once on the full design, so their exchangeability does not
  survive the row subsampling of the bootstrap loop — prefer
  `"knockoff_equi"` or `"knockoff_mvr"` when correlation structure
  matters.

- `"knockoff_equi"`:

  Generates **model-X equicorrelated** knockoffs via
  `knockoff::create.gaussian(..., method = "equi")`. This matches the
  `GaussianSampler(method='equicorrelated')` call in Python STABL and is
  the parity-correct knockoff type for cross-language comparisons.

- `"knockoff_mvr"`:

  Generates **model-X MVR** (minimum-variance- reconstructability)
  knockoffs. The S-matrix is solved by `solve_mvr()` (an R-native port
  of `knockpy.mrc`); sampling uses
  `knockoff::create.gaussian(..., diag_s = S)`. This is a novel feature
  exclusive to `stablr`.

## See also

[`compute_fdp_plus()`](https://mdmanurung.github.io/stablr/reference/compute_fdp_plus.md)
which consumes the artificial-feature scores,
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
which calls this function automatically.

## Examples

``` r
set.seed(1L)
x <- matrix(
  rnorm(20 * 4), 20, 4,
  dimnames = list(paste0("s", 1:20), paste0("f", 1:4))
)
generated <- stablr:::make_artificial_features(
  x, n_injected = 2L, type = "random_permutation", random_state = 1L
)
dim(generated$x_augmented)
#> [1] 20  6
generated$artificial_provenance$actual_type
#> [1] "random_permutation"
try(stablr:::make_artificial_features(x, 0L, "random_permutation"))
#> Error : `n_injected` must be >= 1.
```
