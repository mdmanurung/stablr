# Plot a Precision-Recall Curve

Generates a ggplot2 Precision-Recall Curve (PRC) from binary outcome
labels and predicted probabilities, optionally with iso-F1 contour
lines.

## Usage

``` r
plot_prc(y_true, y_preds, show_iso = TRUE, title = "Precision-Recall Curve")
```

## Arguments

- y_true:

  Numeric 0/1 or logical vector of binary outcomes (1/`TRUE` for the
  positive class, 0/`FALSE` for the negative class). Must contain both
  classes and be the same length as `y_preds`. Character or factor
  labels should be recoded explicitly so the positive class is
  unambiguous.

- y_preds:

  Numeric vector of predicted probabilities for the positive class.
  Values must be finite and in \\\[0, 1\]\\.

- show_iso:

  Logical; if `TRUE` (default), four iso-F1 contour lines at F1 = 0.2,
  0.4, 0.6, 0.8 are drawn as reference guides.

- title:

  Character scalar; plot title. Default `"Precision-Recall Curve"`.

## Value

A `ggplot` object. Average precision is shown as a caption.

## Details

Precision-recall curves are preferred over ROC curves when the positive
class is rare (class imbalance), because they focus on the model's
performance on positive predictions without being diluted by the large
number of true negatives. Use this plot after applying STABL feature
selection and fitting a downstream binary classifier. Tied prediction
scores are aggregated at a single threshold before stepwise average
precision computation, making the result invariant to row order and
avoiding optimistic linear interpolation across precision drops. Average
precision is shown in the caption.

## See also

[`plot_roc()`](https://mdmanurung.github.io/stablr/reference/plot_roc.md)
for the ROC curve alternative.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  y_true  <- sample(0:1, 50, replace = TRUE)
  y_preds <- pmin(pmax(y_true * 0.6 + rnorm(50, 0, 0.3), 0), 1)
  stablr:::plot_prc(y_true, y_preds)
}
```
