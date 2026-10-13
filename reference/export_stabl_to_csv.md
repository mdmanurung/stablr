# Export STABL Stability Scores to CSV

Writes the full stability-score matrix, the per-feature maximum
stability scores, and (when artificial features were used) the
corresponding artificial-feature scores to CSV files in the directory
specified by `path`.

## Usage

``` r
export_stabl_to_csv(object, path)
```

## Arguments

- object:

  A fitted `"stabl_fit"` object from
  [`stabl_fit()`](https://mdmanurung.github.io/stablr/reference/stabl_fit.md).

- path:

  Character string; path to an existing (or to-be-created) directory.

## Value

Invisibly returns `path` (as a normalized string).

## Details

Files written:

- `STABL scores.csv` — full p x L stability score matrix.

- `Max STABL scores.csv` — per-feature maximum stability score, sorted
  descending.

- `STABL artificial scores.csv` — artificial feature scores (only when
  `object$artificial_type` is not `NULL`).

- `Max STABL artificial scores.csv` — artificial max scores, sorted
  descending (only when artificial features were used).

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
out_dir <- file.path(tempdir(), "stabl_csv_export")
export_stabl_to_csv(fit, out_dir)
list.files(out_dir)
#> [1] "Max STABL artificial scores.csv" "Max STABL scores.csv"           
#> [3] "STABL artificial scores.csv"     "STABL scores.csv"               
# }
```
