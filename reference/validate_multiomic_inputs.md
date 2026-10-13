# Validate Multi-Omic Input Contract

Enforces the canonical `stablr` input contract for multi-omic analyses:
a named list of omic tables with identical sample IDs and row order
across all views, plus strict alignment with the outcome and optional
group vectors.

## Usage

``` r
validate_multiomic_inputs(x_list, y, groups = NULL)
```

## Arguments

- x_list:

  Named list of `data.frame` or numeric matrix omic tables. Each table
  must have non-empty row names identifying samples. All tables must
  have the same set of row names **in the same order**.

- y:

  Named outcome vector aligned to the samples in `x_list`. The name set
  must match the row names of every omic table.

- groups:

  Optional named group vector. When supplied its names must match the
  sample IDs in `x_list`. Pass `NULL` to skip group validation.

## Value

Invisibly returns `TRUE` when all checks pass. Raises an informative
error as soon as the first violation is found.

## Details

Multi-omic analyses require that all omic matrices have exactly the same
rows in exactly the same order so that row-wise operations (bootstrap
sampling, cooperative learning) are coherent. This function catches
mismatches early and provides informative error messages naming the
offending omic layer, which is much easier to diagnose than silent
misalignment detected later in the modelling pipeline.

## See also

[`validate_sample_alignment()`](https://mdmanurung.github.io/stablr/reference/validate_sample_alignment.md)
for single-omic inputs,
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
which calls this automatically.

## Examples

``` r
ids <- paste0("s", 1:4)
x_list <- list(
  rna = matrix(1:12, 4, 3, dimnames = list(ids, paste0("g", 1:3))),
  protein = matrix(1:8, 4, 2, dimnames = list(ids, paste0("p", 1:2)))
)
y <- setNames(c(0, 1, 0, 1), ids)
validate_multiomic_inputs(x_list, y)
try(validate_multiomic_inputs(x_list[c("rna", "rna")], y))
#> Error : omic names must not contain duplicate values.
```
