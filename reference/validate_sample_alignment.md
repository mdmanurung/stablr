# Validate Sample Alignment Across Inputs

Checks that the predictor matrix, outcome vector, and optional group
vector all refer to exactly the same set of samples (by name), and that
they can be safely aligned before modelling.

## Usage

``` r
validate_sample_alignment(x, y, groups = NULL)
```

## Arguments

- x:

  A `data.frame` or numeric matrix with non-empty, non-`NA` row names
  used as sample IDs.

- y:

  A named vector (or matrix-like object such as
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) with
  row names) whose names identify the samples. The set of names in `y`
  must be identical to the set of row names in `x`; order does not need
  to match.

- groups:

  Optional named vector where names are sample IDs and values are group
  memberships. When supplied, its name set must also match
  `rownames(x)`. Pass `NULL` to skip group validation.

## Value

Invisibly returns `TRUE` when all checks pass. Raises an informative
error as soon as the first violation is found.

## Details

STABL enforces strict name-based alignment (rather than positional
alignment) to mirror the pandas index-alignment semantics of the Python
reference implementation and to prevent silent sample-order bugs that
could corrupt stability scores or introduce outcome leakage.

This function is called automatically by
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
and
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md);
you only need to call it directly when building a custom pre-processing
step that receives `x` and `y` separately.

## See also

[`validate_multiomic_inputs()`](https://mdmanurung.github.io/stablr/reference/validate_multiomic_inputs.md)
for multi-omic list inputs.

## Examples

``` r
x <- matrix(
  1:12, 4, 3,
  dimnames = list(paste0("s", 1:4), paste0("f", 1:3))
)
y <- setNames(c(0, 1, 0, 1), c("s3", "s1", "s4", "s2"))
validate_sample_alignment(x, y)
try(validate_sample_alignment(x, setNames(y, paste0("v", 1:4))))
#> Error : Sample mismatch between `x` row names and `y` names.
```
