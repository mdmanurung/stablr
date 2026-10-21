# Save STABL Results to Disk

Writes STABL result artefacts into a single output directory. With the
default `what = "all"` this is the complete set — stability scores, the
list of selected features, the FDR diagnostic graph, the stability path
plot, and per-feature distribution plots — and is the recommended way to
persist and share results from a single
[`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md)
call. With `what = "csv"` only the tabular CSV exports are written.

## Usage

``` r
save_stabl_results(
  object,
  path,
  x = NULL,
  y = NULL,
  figure_fmt = "pdf",
  new_hard_threshold = NULL,
  task_type = "binary",
  override = FALSE,
  what = c("all", "csv")
)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object returned by
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- path:

  Character string; path to the output directory. Created recursively if
  it does not already exist. When `what = "all"`, raises an error if the
  directory exists and `override = FALSE`; when `what = "csv"`, an
  existing directory is reused without error.

- x:

  Numeric matrix of predictors used to fit `object`. Required when
  `what = "all"` (for per-feature distribution plots); must contain at
  least the selected feature columns. Ignored when `what = "csv"`.

- y:

  Outcome vector, factor, or
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object
  used to fit `object`. Required when `what = "all"`; ignored when
  `what = "csv"`.

- figure_fmt:

  Character; graphics device extension used when saving plots. Common
  choices: `"pdf"` (default, publication quality), `"png"`, `"svg"`.

- new_hard_threshold:

  Numeric in `(0, 1]` or `NULL`. When supplied, overrides the stored
  threshold for support extraction and stability-path plotting.

- task_type:

  Character; one of `"binary"` (default), `"multiclass"`, or
  `"regression"`. Determines which per-feature plot is generated:
  grouped boxplots for classification, LOESS scatter for regression.

- override:

  Logical; if `TRUE`, an existing directory with the same path is
  accepted and its contents may be overwritten. Default `FALSE`. Only
  used when `what = "all"`.

- what:

  Character; `"all"` (default) for the complete result set including
  plots, or `"csv"` for the CSV tables only. The `"csv"` mode requires
  neither `x`/`y` nor a ggplot2 installation.

## Value

Invisibly returns `path` as a normalised absolute string.

## Details

With `what = "all"` the function orchestrates several lower-level
helpers: the CSV export engine for tabular data,
[`plot_fdr_graph()`](https://mdmanurung.github.io/stablr/reference/plot_fdr_graph.md)
and
[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md)
for diagnostic plots, and
[`boxplot_features()`](https://mdmanurung.github.io/stablr/reference/boxplot_features.md)
or
[`scatterplot_features()`](https://mdmanurung.github.io/stablr/reference/scatterplot_features.md)
for per-feature distribution visualisations. All graphics are saved
using
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
in the format specified by `figure_fmt`.

CSV files written (both modes):

- `STABL scores.csv` — full p x L stability score matrix.

- `Max STABL scores.csv` — per-feature maximum stability score, sorted
  descending.

- `STABL artificial scores.csv` — artificial feature scores (only when
  artificial-feature scores are present, i.e.
  `object$stabl_scores_artificial_` is not `NULL`).

- `Max STABL artificial scores.csv` — artificial max scores, sorted
  descending (only when artificial features were used).

## See also

[`plot_stabl_path()`](https://mdmanurung.github.io/stablr/reference/plot_stabl_path.md),
[`plot_fdr_graph()`](https://mdmanurung.github.io/stablr/reference/plot_fdr_graph.md),
[`boxplot_features()`](https://mdmanurung.github.io/stablr/reference/boxplot_features.md)

## Examples

``` r
# \donttest{
set.seed(1L)
x <- matrix(rnorm(40 * 6), 40, 6,
             dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
y <- setNames(rnorm(40), rownames(x))
fit <- stabl_fit(x, y,
                 lambda_grid  = data.frame(lambda = c(0.3, 0.1, 0.05)),
                 n_bootstraps = 6L, hard_threshold = 0.3, random_state = 1L)

# CSV tables only (no ggplot2 required)
out_csv <- file.path(tempdir(), "stabl_csv_export")
save_stabl_results(fit, out_csv, what = "csv")
list.files(out_csv)
#> [1] "Max STABL artificial scores.csv" "Max STABL scores.csv"           
#> [3] "STABL artificial scores.csv"     "STABL scores.csv"               

# Complete result set
if (requireNamespace("ggplot2", quietly = TRUE)) {
  out_dir <- file.path(tempdir(), "stabl_results")
  save_stabl_results(fit, out_dir, x = x, y = y,
                     task_type = "regression", figure_fmt = "png")
  list.files(out_dir)
}
#> [1] "FDR Graph.png"                   "Max STABL artificial scores.csv"
#> [3] "Max STABL scores.csv"            "STABL artificial scores.csv"    
#> [5] "STABL scores.csv"                "Selected Features"              
#> [7] "Stability Path.png"             
# }
```
