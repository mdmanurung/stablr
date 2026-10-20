# Plot a ROC Curve

Generates a ggplot2 ROC (Receiver Operating Characteristic) curve from
binary outcome labels and predicted probabilities, together with the
chance-level diagonal for reference.

## Usage

``` r
plot_roc(y_true, y_preds, title = "ROC Curve")
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

- title:

  Character scalar; plot title. Default `"ROC Curve"`.

## Value

A `ggplot` object. The AUC is shown as a caption.

## Details

ROC curves visualise the trade-off between sensitivity (TPR) and
specificity (1 - FPR) across all possible classification thresholds. Use
this plot to assess overall discriminative ability after applying STABL
feature selection and fitting a downstream classifier on the selected
features. Tied prediction scores are aggregated at a single threshold
before AUC computation, making the result invariant to row order. The
Area Under the Curve (AUC) is shown in the caption.

## See also

[`plot_prc()`](https://mdmanurung.github.io/stablr/reference/plot_prc.md)
for precision-recall curves (preferred when classes are severely
imbalanced).

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  y_true  <- sample(0:1, 50, replace = TRUE)
  y_preds <- pmin(pmax(y_true * 0.6 + rnorm(50, 0, 0.3), 0), 1)
  stablr:::plot_roc(y_true, y_preds)
}
```
