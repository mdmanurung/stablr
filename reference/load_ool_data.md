# Load the Onset of Labor example dataset

Reads the subsetted Onset of Labor (OOL) cytokine and proteomics data
that is bundled with the stablr package. The files contain 150 training
samples and up to 27 validation samples, with the first 100 feature
columns retained from each omic layer. Sample IDs are intersected across
all three sources (CyTOF, Proteomics, outcome) so the returned matrices
and vectors are fully aligned with no missing values.

## Usage

``` r
load_ool_data(split = c("train", "valid"))
```

## Arguments

- split:

  One of `"train"` (default) or `"valid"` to select the training or
  validation split.

## Value

A named list with three elements:

- `x_list`:

  A named list of numeric matrices, one per omic (`"cytof"`,
  `"proteomics"`). Rows are samples, columns are features.

- `y`:

  A named numeric vector of gestational-age-at-delivery offsets (days
  before onset of labor, DOS), aligned to the rows of `x_list`.

- `ids`:

  A character vector of sample IDs (same order as rows).

## Examples

``` r
ool <- load_ool_data()
dim(ool$x_list$cytof)      # 150 x 100
#> [1] 150 100
dim(ool$x_list$proteomics) # 150 x 100
#> [1] 150 100
length(ool$y)              # 150
#> [1] 150

ool_val <- load_ool_data(split = "valid")
dim(ool_val$x_list$cytof)  # up to 21 x 100 (intersection of omics)
#> [1]  21 100
```
