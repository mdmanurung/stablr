# Compare Sample-I.I.D. and Study-Disjoint Nested CV

Runs
[`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
twice on the same data — once with sample-i.i.d. folds and once with
study-disjoint folds — and summarises the gap between the two
generalization estimates. A large gap signals that the sample-i.i.d.
estimate is optimistic because held-out samples share study-level
structure (e.g. outcome prevalence or feature shifts) with the training
data.

## Usage

``` r
stabl_multiomic_scheme_comparison(x_list, y, study_labels, ...)
```

## Arguments

- x_list:

  Named list of omic matrices or data frames with identical sample row
  names.

- y:

  Named outcome vector (factor for classification, numeric for
  regression).

- study_labels:

  Named vector assigning each sample to a study; names must cover
  exactly the sample IDs. Required.

- ...:

  Additional arguments passed to
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
  (`candidates`, `lambda_grid`, `metric`, `family`, `n_bootstraps`,
  `artificial_type`, `hard_threshold`, `stratified`, `random_state`,
  `n_lambda`, `l1_ratio`, `outer_v`, `inner_v`, `cv_workers`,
  `workers`).

## Value

A list with `sample_iid` and `study_disjoint` (the two
`stabl_multiomic_nested_cv` objects), `comparison` (one row per scheme
with `balanced_error_rate`, `accuracy`, and `macro_f1`), and `ber_gap`
(study-disjoint minus sample-i.i.d. balanced error rate; positive values
mean the sample-i.i.d. estimate was optimistic).

## Details

Both arms receive identical arguments (including `random_state`), so the
only difference is the CV sampling unit. Under `stratified = TRUE` the
stratification is interpreted at the sample level for the sample-i.i.d.
arm and at the study level for the study-disjoint arm, matching the
semantics of each scheme.

## See also

[`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md),
[`stabl_multiomic_study_diagnostics()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_study_diagnostics.md)

## Examples

``` r
set.seed(6L)
ids <- paste0("s", 1:48)
y <- setNames(factor(rep(c("A", "B"), times = 24L)), ids)
x <- matrix(rnorm(48 * 4, sd = 0.5), 48, 4,
            dimnames = list(ids, paste0("g", 1:4)))
x[, 1:2] <- x[, 1:2] + 1.2 * (y == "A")
study_labels <- setNames(rep(paste0("study", 1:6), each = 8L), ids)
# study6 is outcome-skewed: sample-i.i.d. folds leak its composition
y[study_labels == "study6"] <-
  factor(rep(c("A", "A", "A", "A", "A", "A", "B", "B")), levels = levels(y))
cmp <- stabl_multiomic_scheme_comparison(
  x_list = list(rna = x), y = y, study_labels = study_labels,
  candidates = list(list(name = "rna", blocks = "rna")),
  lambda_grid = data.frame(lambda = c(0.2, 0.1)),
  outer_v = 2L, inner_v = 2L, n_bootstraps = 10L,
  artificial_type = NULL, hard_threshold = 1e-9,
  sample_fraction = 1, random_state = 6L
)
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
cmp$comparison
#>           scheme balanced_error_rate  accuracy  macro_f1
#> 1     sample_iid           0.1486014 0.8541667 0.8525669
#> 2 study_disjoint           0.2167832 0.7916667 0.7857143
cmp$ber_gap
#> [1] 0.06818182
```
