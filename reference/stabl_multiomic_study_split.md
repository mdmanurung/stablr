# Study-Disjoint Train/Validation Split

Assigns whole studies to a training and a validation split so that
[`stabl_multiomic_train_validate()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_train_validate.md)
evaluates a cross-study deployment setting rather than a sample-i.i.d.
one. Studies are assigned largest-first by a deterministic greedy pass
targeting `prop_train` of the samples while keeping outcome class
proportions as balanced as possible across the two sides; a seeded
shuffle fixes the order of size ties.

## Usage

``` r
stabl_multiomic_study_split(
  x_list,
  y,
  study_labels,
  prop_train = 0.75,
  stratified = TRUE,
  random_state = NULL
)
```

## Arguments

- x_list:

  Named list of omic matrices or data frames with identical sample row
  names.

- y:

  Named outcome vector aligned to `x_list` rows.

- study_labels:

  Named vector assigning each sample to a study; names must cover
  exactly the sample IDs. Each study must contain at least 2 samples,
  and at least 2 studies are required.

- prop_train:

  Target proportion of samples in the training split (strictly between 0
  and 1).

- stratified:

  Logical. When `TRUE` (default), class-proportion imbalance is the
  secondary assignment criterion after the size target.

- random_state:

  Optional integer seed for deterministic splits.

## Value

A list with `train_ids`, `valid_ids` (sample IDs in the original row
order), `study_assignment` (named vector mapping each study to `"train"`
or `"valid"`), `class_counts` (per-side held-out outcome counts), and
`fold_scheme = "study_disjoint"`.

## Details

`prop_train` is a target, not a guarantee: study granularity governs the
achievable split, and extreme values can fail the fail-closed checks
below.

## Examples

``` r
set.seed(6L)
ids <- paste0("s", 1:12)
y <- setNames(factor(rep(c("A", "B"), times = 6L)), ids)
x <- matrix(rnorm(12 * 4), 12, 4, dimnames = list(ids, paste0("g", 1:4)))
study_labels <- setNames(rep(paste0("study", 1:3), each = 4L), ids)
split <- stabl_multiomic_study_split(
  x_list = list(rna = x), y = y, study_labels = study_labels,
  prop_train = 2 / 3, random_state = 6L
)
split$study_assignment
#>  study1  study2  study3 
#> "valid" "train" "train" 
```
