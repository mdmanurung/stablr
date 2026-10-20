#' Load the Onset of Labor example dataset
#'
#' Reads the subsetted Onset of Labor (OOL) cytokine and proteomics data that
#' is bundled with the stablr package. The files contain 150 training samples
#' and 21 validation samples (after ID intersection), with the first 100
#' feature columns retained from each omic layer.  Sample IDs are intersected
#' across all three sources (CyTOF, Proteomics, outcome) so the returned
#' matrices and vectors are fully aligned with no missing values.
#'
#' @param split One of `"train"` (default) or `"valid"` to select the training
#'   or validation split.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{`x_list`}{A named list of numeric matrices, one per omic
#'       (`"cytof"`, `"proteomics"`). Rows are samples, columns are features.}
#'     \item{`y`}{A named numeric vector of gestational-age-at-delivery offsets
#'       (days before onset of labor, DOS), aligned to the rows of `x_list`.}
#'     \item{`ids`}{A character vector of sample IDs (same order as rows).}
#'   }
#'
#' @examples
#' ool <- load_ool_data()
#' dim(ool$x_list$cytof)      # 150 x 100
#' dim(ool$x_list$proteomics) # 150 x 100
#' length(ool$y)              # 150
#'
#' ool_val <- load_ool_data(split = "valid")
#' dim(ool_val$x_list$cytof)  # up to 21 x 100 (intersection of omics)
#'
#' @export
load_ool_data <- function(split = c("train", "valid")) {
  split <- match.arg(split)

  pkg <- "stablr"
  if (split == "train") {
    cytof_file <- system.file("extdata", "ool_cytof_train.csv.gz",      package = pkg, mustWork = TRUE)
    prot_file  <- system.file("extdata", "ool_proteomics_train.csv.gz", package = pkg, mustWork = TRUE)
    dos_file   <- system.file("extdata", "ool_dos_train.csv",           package = pkg, mustWork = TRUE)
  } else {
    cytof_file <- system.file("extdata", "ool_cytof_valid.csv.gz",      package = pkg, mustWork = TRUE)
    prot_file  <- system.file("extdata", "ool_proteomics_valid.csv.gz", package = pkg, mustWork = TRUE)
    dos_file   <- system.file("extdata", "ool_dos_valid.csv",           package = pkg, mustWork = TRUE)
  }

  cytof_df <- read.csv(cytof_file, check.names = FALSE)
  prot_df  <- read.csv(prot_file,  check.names = FALSE)
  dos_df   <- read.csv(dos_file,   check.names = FALSE)

  # First column of each file is the sample ID
  cytof_ids <- cytof_df[[1]]
  prot_ids  <- prot_df[[1]]
  dos_ids   <- dos_df[[1]]

  # Fail closed on duplicate IDs (match() below would silently take the
  # first) and non-numeric feature columns (review A16).
  if (anyDuplicated(cytof_ids) || anyDuplicated(prot_ids) ||
      anyDuplicated(dos_ids)) {
    stop("Duplicate sample IDs found in at least one OOL data file.",
         call. = FALSE)
  }
  if (!all(vapply(cytof_df[, -1L, drop = FALSE], is.numeric, logical(1L))) ||
      !all(vapply(prot_df[, -1L, drop = FALSE], is.numeric, logical(1L)))) {
    stop("Non-numeric feature column found in an OOL data file.",
         call. = FALSE)
  }

  # Align to common sample IDs (intersection)
  common_ids <- Reduce(intersect, list(cytof_ids, prot_ids, dos_ids))
  if (length(common_ids) == 0L) {
    stop("No common sample IDs found across CyTOF, Proteomics, and DOS files.")
  }

  .align <- function(df, ids) {
    row_idx <- match(common_ids, ids)
    m <- as.matrix(df[row_idx, -1L, drop = FALSE])
    rownames(m) <- common_ids
    storage.mode(m) <- "double"
    m
  }

  cytof_mat <- .align(cytof_df, cytof_ids)
  prot_mat  <- .align(prot_df,  prot_ids)

  dos_idx <- match(common_ids, dos_ids)
  y_vec   <- setNames(as.numeric(dos_df[[2]][dos_idx]), common_ids)

  list(
    x_list = list(cytof = cytof_mat, proteomics = prot_mat),
    y      = y_vec,
    ids    = common_ids
  )
}
