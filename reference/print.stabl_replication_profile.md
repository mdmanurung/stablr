# Print a Cohort Replication Profile

Print a Cohort Replication Profile

## Usage

``` r
# S3 method for class 'stabl_replication_profile'
print(x, ...)
```

## Arguments

- x:

  A `stabl_replication_profile` object.

- ...:

  Ignored; present for S3 compatibility.

## Value

The input `x`, invisibly.

## Examples

``` r
set.seed(1)
ids <- paste0("s", 1:60)
x <- matrix(rnorm(60 * 8), 60, 8,
            dimnames = list(ids, paste0("f", 1:8)))
y <- setNames(x[, 1] * 3 + rnorm(60, sd = 0.8), ids)
studies <- setNames(rep(paste0("study", 1:3), each = 20L), ids)
rp <- stabl_replication_profile(
  x, y, studies,
  fit_spec = list(n_bootstraps = 50,
                  artificial_type = "random_permutation"),
  seed = 5
)
print(rp)
#> STABL cohort replication profile
#>   studies: 3 | features: 8 | pooled support: 7
#>   training replication: median 0.67 (range 0.33-1.00) over features
#>   pooled panel LOSO retention:
#>     f1: 1.00
#>     f2: 0.67 (dropped on removal of: study1)
#>     f3: 1.00
#>     f5: 0.67 (dropped on removal of: study3)
#>     f6: 0.33 (dropped on removal of: study1,study3)
#>     f7: 1.00
#>     f8: 0.33 (dropped on removal of: study1,study3)
#>   note: four distinct estimands (pooled stability, training-cohort
#>   replication, direction replication, held-out performance);
#>   selection agreement is not causal agreement. Cohort-level
#>   observations only -- no binomial CIs over dependent folds.
```
