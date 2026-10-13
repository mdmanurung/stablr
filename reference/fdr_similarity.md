# FDR Between Two Feature Sets

Computes the False Discovery Rate when `list1` is treated as the
predicted selection and `list2` as the known ground-truth. This metric
is only meaningful in simulation benchmarks where the true signal
features are known in advance.

## Usage

``` r
fdr_similarity(list1, list2)
```

## Arguments

- list1:

  Predicted feature set (character or integer vector).

- list2:

  True feature set (ground truth).

## Value

Numeric scalar in \\\[0, 1\]\\. Returns 0 when the predicted set is
empty (no false discoveries possible).

## See also

[`tpr_similarity()`](https://mdmanurung.github.io/stablr/reference/tpr_similarity.md),
[`fscore_similarity()`](https://mdmanurung.github.io/stablr/reference/fscore_similarity.md)

## Examples

``` r
true_signal <- c("f1", "f2", "f3")
predicted   <- c("f1", "f2", "f4", "f5")  # f4, f5 are false discoveries
fdr_similarity(predicted, true_signal)      # 2/4 = 0.5
#> [1] 0.5
fdr_similarity(character(0), true_signal)   # 0 (empty prediction)
#> [1] 0
```
