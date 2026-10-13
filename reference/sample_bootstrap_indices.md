# Classic Bootstrap Sampler Indices

Draws a random set of row indices from a training set for a single STABL
bootstrap iteration. This is the standard (non-grouped) sampler used
when all samples are independent.

## Usage

``` r
sample_bootstrap_indices(
  y,
  n_subsamples,
  replace = TRUE,
  class_weights = NULL,
  stratify = FALSE,
  strata = NULL,
  random_state = NULL,
  seed = NULL
)
```

## Arguments

- y:

  Outcome vector whose length equals the number of training samples.
  Values are only examined to check class diversity and to apply
  `class_weights`; the type is arbitrary (numeric, factor, character).

- n_subsamples:

  Positive integer; number of row indices to draw. Must not exceed
  `length(y)` when `replace = FALSE`.

- replace:

  Logical; sample with replacement? Default `TRUE`.

- class_weights:

  Optional named numeric vector keyed by class labels (as produced by
  `as.character(y)`). When `NULL` all samples are equally likely.

- stratify:

  Logical; if `TRUE`, draw approximately the same class proportions as
  `y` by sampling within each outcome class. This is useful for small or
  imbalanced classification tasks where unstratified subsampling can
  produce very sparse minority-class bootstraps. Default `FALSE`
  preserves the original STABL sampling behavior.

- strata:

  Optional categorical stratification design. Provide a vector for one
  stratification factor, or a `data.frame`/matrix/list for a joint
  design across multiple factors. Bootstrap samples are drawn within the
  interaction of all supplied columns. When `NULL`, no stratification is
  used unless `stratify = TRUE`, in which case `y` is used as the
  strata.

- random_state:

  Optional integer; seed used for this draw. The function does **not**
  clobber the caller's global RNG state: the previous `.Random.seed` is
  restored on exit (via
  [`withr::local_seed()`](https://withr.r-lib.org/reference/with_seed.html)).

- seed:

  Deprecated; use `random_state` instead. If supplied, a deprecation
  warning is emitted and the value is used as `random_state` when the
  latter is `NULL`. It will be removed in a future release.

## Value

Integer vector of length `n_subsamples` containing 1-based row indices
into the training set.

## Details

STABL accumulates how often each feature is selected across many
bootstrap subsamples; these subsamples must be drawn consistently to
avoid bias. This function centralises that sampling so that
class-imbalance handling and the "degenerate bootstrap" guard (see
below) are applied uniformly.

**Class weighting:** When `class_weights` is supplied each sample's
inclusion probability is proportional to the weight of its class. This
is useful for imbalanced binary or multi-class outcomes where unweighted
subsampling would frequently produce single-class bootstraps.

**Degenerate bootstrap guard:** When the drawn subsample contains only
one unique class but the full outcome contains at least two, the
function retries (without a seed) until a class-diverse subsample is
found. This prevents downstream model failures in classifiers that
require at least two classes.

## See also

[`group_bootstrap_indices()`](https://mdmanurung.github.io/stablr/reference/group_bootstrap_indices.md)
for repeated-measures/grouped data,
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
which calls this sampler automatically.

## Examples

``` r
set.seed(42L)
y <- c(rep(0, 15), rep(1, 15))  # balanced binary outcome
idx <- sample_bootstrap_indices(y, n_subsamples = 20L, random_state = 1L)
table(y[idx])  # class distribution in the bootstrap
#> 
#>  0  1 
#> 11  9 

# Stratified sampling preserves class proportions
idx_str <- sample_bootstrap_indices(y, n_subsamples = 20L,
                                     stratify = TRUE, random_state = 1L)
table(y[idx_str])
#> 
#>  0  1 
#> 10 10 
```
