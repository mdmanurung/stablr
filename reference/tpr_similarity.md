# TPR Between Two Feature Sets

Computes the True Positive Rate (sensitivity / recall) when `list1` is
the predicted selection and `list2` is the known ground-truth. Use
together with
[`fdr_similarity()`](https://mdmanurung.github.io/stablr/reference/fdr_similarity.md)
and
[`fscore_similarity()`](https://mdmanurung.github.io/stablr/reference/fscore_similarity.md)
to characterise the precision-recall trade-off in simulation benchmarks.

## Usage

``` r
tpr_similarity(list1, list2)
```

## Arguments

- list1:

  Predicted feature set (character or integer vector).

- list2:

  True feature set (ground truth).

## Value

Numeric scalar in \\\[0, 1\]\\. Returns 0 when the true set is empty.

## See also

[`fdr_similarity()`](https://mdmanurung.github.io/stablr/reference/fdr_similarity.md),
[`fscore_similarity()`](https://mdmanurung.github.io/stablr/reference/fscore_similarity.md)

## Examples

``` r
true_signal <- c("f1", "f2", "f3")
predicted   <- c("f1", "f2", "f4")  # misses f3
stablr:::tpr_similarity(predicted, true_signal)       # 2/3
#> [1] 0.6666667
stablr:::tpr_similarity(character(0), true_signal)    # 0
#> [1] 0
stablr:::tpr_similarity(c("f1","f2","f3"), true_signal)  # 1
#> [1] 1
```
