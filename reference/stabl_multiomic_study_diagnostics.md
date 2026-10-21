# Per-Study Batch Diagnostics for Multi-Omic STABL Inputs

Summarises each study's outcome composition and the location shift of
its features relative to the pooled rest of the data. The diagnostics
quantify how much study-level structure a sample-i.i.d. cross-validation
scheme would leak across folds, and feed the study-disjoint benchmark
(`benchmarks/b1_study_disjoint_cv.R`).

## Usage

``` r
stabl_multiomic_study_diagnostics(x_list, y, study_labels)
```

## Arguments

- x_list:

  Named list of omic matrices or data frames with identical sample row
  names.

- y:

  Named outcome vector aligned to `x_list` rows. Factors or characters
  are summarised as per-level counts and prevalences; numeric outcomes
  as per-study mean and SD.

- study_labels:

  Named vector assigning each sample to a study (cohort, site, or
  batch); names must cover exactly the sample IDs. Each study must
  contain at least 2 samples.

## Value

A list with two data frames: `study_summary` (one row per study with
sample counts and outcome composition) and `feature_shift` (one row per
study and omic block with `n_features`, `smd_mean_abs`, `smd_p95_abs`,
`smd_max_abs`, and `n_zero_sd`).

## Details

Feature shifts are reported as absolute standardized mean differences
(SMD) per feature: `(mean_study - mean_rest) / pooled_sd` with
`pooled_sd = sqrt((sd_study^2 + sd_rest^2) / 2)`. Features with zero
pooled SD are excluded from the summaries and counted in `n_zero_sd`.

## See also

[`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
for study-disjoint cross-validation and
[`stabl_multiomic_scheme_comparison()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_scheme_comparison.md)
for comparing study-disjoint validation schemes (both take
`study_labels` directly).

## Examples

``` r
set.seed(5L)
ids <- paste0("s", 1:12)
y <- setNames(factor(rep(c("A", "B"), times = 6L)), ids)
x <- matrix(rnorm(12 * 4), 12, 4, dimnames = list(ids, paste0("g", 1:4)))
study_labels <- setNames(rep(paste0("study", 1:3), each = 4L), ids)
# study1 carries a location shift on the first two features
x[study_labels == "study1", 1:2] <- x[study_labels == "study1", 1:2] + 2
diag <- stabl_multiomic_study_diagnostics(
  x_list = list(rna = x), y = y, study_labels = study_labels
)
diag$study_summary
#>    study n_samples n_A prop_A n_B prop_B
#> 1 study1         4   2    0.5   2    0.5
#> 2 study2         4   2    0.5   2    0.5
#> 3 study3         4   2    0.5   2    0.5
diag$feature_shift
#>    study block n_features smd_mean_abs smd_p95_abs smd_max_abs n_zero_sd
#> 1 study1   rna          4    1.1676817    1.665309    1.702477         0
#> 2 study2   rna          4    1.4616423    2.208855    2.344658         0
#> 3 study3   rna          4    0.6915646    1.270697    1.369367         0
#>   n_missing
#> 1         0
#> 2         0
#> 3         0
```
