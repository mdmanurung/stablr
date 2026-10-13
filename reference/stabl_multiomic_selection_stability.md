# Selection Stability Across Study-Disjoint Folds

Quantifies how consistent the selected feature sets are across the outer
folds of a study-disjoint nested CV. Each fold's selection signature is
the union of features selected for the winning candidate across repeats;
pairwise Jaccard similarities between fold signatures summarise whether
the support is stable when whole studies are held out.

## Usage

``` r
stabl_multiomic_selection_stability(object)
```

## Arguments

- object:

  A `stabl_multiomic_nested_cv` object returned by
  [`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
  with `fold_scheme = "study_disjoint"`.

## Value

A list with `similarity_matrix` (fold x fold Jaccard matrix, rows and
columns in outer-fold order), `mean_pairwise_similarity`,
`min_pairwise_similarity`, `fold_selection_sizes` (named integer vector
of per-fold signature sizes), `n_repeats`, and `fold_scheme`. Render the
matrix with e.g. `stats::heatmap(result$similarity_matrix)`.

## See also

[`jaccard_similarity_matrix()`](https://mdmanurung.github.io/stablr/reference/jaccard_similarity_matrix.md)
for the underlying metric and
[`stabl_multiomic_nested_cv()`](https://mdmanurung.github.io/stablr/reference/stabl_multiomic_nested_cv.md)
for the study-disjoint scheme.

## Examples

``` r
set.seed(7L)
ids <- paste0("s", 1:48)
y <- setNames(factor(rep(c("A", "B"), times = 24L)), ids)
x <- matrix(rnorm(48 * 5, sd = 0.5), 48, 5,
            dimnames = list(ids, paste0("g", 1:5)))
x[, 1:2] <- x[, 1:2] + 1.2 * (y == "A")
study_labels <- setNames(rep(paste0("study", 1:6), each = 8L), ids)
nested <- suppressWarnings(stabl_multiomic_nested_cv(
  x_list = list(rna = x), y = y, study_labels = study_labels,
  candidates = list(list(name = "rna", blocks = "rna")),
  lambda_grid = data.frame(lambda = c(0.2, 0.1)),
  outer_v = 2L, inner_v = 2L, n_bootstraps = 10L,
  artificial_type = NULL, hard_threshold = 1e-9,
  sample_fraction = 1, random_state = 7L
))
stab <- stabl_multiomic_selection_stability(nested)
stab$similarity_matrix
#>               Repeat1_Fold1 Repeat1_Fold2
#> Repeat1_Fold1             1             1
#> Repeat1_Fold2             1             1
stab$mean_pairwise_similarity
#> [1] 1
```
