# Print an Outcome-Null Validation Summary

Print an Outcome-Null Validation Summary

## Usage

``` r
# S3 method for class 'stabl_null_validation'
print(x, ...)
```

## Arguments

- x:

  A `stabl_null_validation` object from
  [`stabl_null_check()`](https://mdmanurung.github.io/stablr/reference/stabl_null_check.md).

- ...:

  Ignored; present for S3 compatibility.

## Value

The input `x`, invisibly.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
nv <- stabl_null_check(
  x, y,
  fit_spec = list(n_bootstraps = 10,
                  artificial_type = "random_permutation"),
  n_reps = 3, seed = 7
)
print(nv)
#> STABL outcome-null validation
#>   design: independent | truth: induced_global_null_permutation
#>   replicates: 3 successful / 0 failed (of 3 attempted)
#>   P(any selection on null): 0.000 | median panel size: 0
#>   reported FDP+ minima: median 1.12 (range 1.12-1.12); realized FDP on this null is 1 whenever the panel is non-empty
#>   no feature was ever selected on the null replicates
#>   note: permutation induces a global null; this quantifies null selection behavior,
#>   not FDP validity for partly non-null observational data.
```
