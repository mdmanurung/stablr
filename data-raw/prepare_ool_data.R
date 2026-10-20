#!/usr/bin/env Rscript

# Rebuild the bounded package fixtures from the six Dryad study exports.
# Usage:
# Rscript data-raw/prepare_ool_data.R \
#   CyTOF.csv Proteomics.csv DOS.csv \
#   CyTOF_validation.csv Proteomics_validation.csv DOS_validation.csv OUT_DIR

.ool_source_sha256 <- c(
  cytof_train = "c17d4d2ad44a2a8b1da887a58b0cf190c1f74c23755f0f80c6231945f97f47f6",
  proteomics_train = "5037b8e2ff3c28c4d562bb32e6d04977c1532f73d659bdccaefa18271785dfba",
  dos_train = "11889faf2d1c3846025f57fa89e3a92f036aeff7162596a4d69f7eccf72695ec",
  cytof_valid = "912053964c201c215ce557ab9bc43488a74da3f9c1b2b331f5ec5b85de40ccfd",
  proteomics_valid = "55fb10bd1d2093d8bc593e7d4509fb51389f8882b842db1e8bf8139ff6584c71",
  dos_valid = "f1bf4ed0dc1f577fb59e44b4a9014bcc4b7516677705c524790cadfc226c8ab6"
)

.sha256_file <- function(path) {
  # NOTE: openssl::sha256(file(path)) silently hashes only the first chunk
  # of the connection for files larger than a few hundred bytes (verified
  # 2026-10-20 against sha256sum and Python hashlib).  Prefer the system
  # tool; the openssl fallback reads the bytes explicitly with readBin().
  tool <- Sys.which("sha256sum")
  if (nzchar(tool)) {
    return(strsplit(system2(tool, path, stdout = TRUE), "[[:space:]]+")[[1L]][[1L]])
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    bytes <- readBin(path, "raw", n = file.size(path))
    return(as.character(openssl::sha256(bytes)))
  }
  stop("Need openssl or sha256sum to verify OOL inputs.", call. = FALSE)
}

prepare_ool_data <- function(inputs, out_dir, verify_hashes = TRUE,
                             max_feature_columns = 100L) {
  required <- names(.ool_source_sha256)
  if (is.null(names(inputs)) || !setequal(names(inputs), required)) {
    stop("`inputs` must name the six train/validation source tables.", call. = FALSE)
  }
  inputs <- inputs[required]
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) stop("Missing OOL source file(s): ", paste(missing, collapse = ", "))
  if (isTRUE(verify_hashes)) {
    observed <- vapply(inputs, .sha256_file, character(1L))
    if (!identical(unname(observed), unname(.ool_source_sha256))) {
      bad <- names(observed)[observed != .ool_source_sha256]
      stop("OOL source checksum mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
    }
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write_features <- function(src, dst) {
    x <- utils::read.csv(src, check.names = FALSE)
    if (!identical(names(x)[[1L]], "ID")) stop(src, " must start with ID.")
    keep <- seq_len(min(ncol(x), max_feature_columns + 1L))
    con <- gzfile(dst, "wt")
    on.exit(close(con), add = TRUE)
    utils::write.csv(x[, keep, drop = FALSE], con, row.names = FALSE)
  }
  write_outcome <- function(src, dst) {
    x <- utils::read.csv(src, check.names = FALSE)
    if (!identical(names(x), c("ID", "DOS"))) {
      stop(src, " must contain exactly ID and DOS.", call. = FALSE)
    }
    utils::write.csv(x, dst, row.names = FALSE)
  }

  write_features(inputs[["cytof_train"]], file.path(out_dir, "ool_cytof_train.csv.gz"))
  write_features(inputs[["cytof_valid"]], file.path(out_dir, "ool_cytof_valid.csv.gz"))
  write_features(inputs[["proteomics_train"]], file.path(out_dir, "ool_proteomics_train.csv.gz"))
  write_features(inputs[["proteomics_valid"]], file.path(out_dir, "ool_proteomics_valid.csv.gz"))
  write_outcome(inputs[["dos_train"]], file.path(out_dir, "ool_dos_train.csv"))
  write_outcome(inputs[["dos_valid"]], file.path(out_dir, "ool_dos_valid.csv"))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 7L) {
    stop("Expected six source CSV files followed by OUT_DIR.", call. = FALSE)
  }
  prepare_ool_data(
    setNames(args[seq_len(6L)], names(.ool_source_sha256)), args[[7L]]
  )
}
