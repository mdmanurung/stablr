#!/usr/bin/env Rscript

# Build the PBMC 3k vignette fixtures from the cEBMF generation outputs.
# Upstream step: data-raw/prepare_pbmc3k_fixtures.py downloads the 10x
# Genomics PBMC 3k filtered matrix, selects HVGs, computes Louvain cluster
# labels, and fits cEBMF (cebmf_torch, CPU, seeded).  This script assembles
# the R fixtures and their SHA-256 manifest.
# Usage:
#   Rscript data-raw/prepare_pbmc3k_fixtures.R [GEN_DIR] [OUT_DIR]
# Defaults: GEN_DIR = /workspace/pbmc3k/out, OUT_DIR = inst/extdata/pbmc3k

args <- commandArgs(trailingOnly = TRUE)
gen_dir <- if (length(args) >= 1L) args[[1L]] else "/workspace/pbmc3k/out"
out_dir <- if (length(args) >= 2L) args[[2L]] else "inst/extdata/pbmc3k"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(Matrix))

# ── counts fixture (cells x HVGs, RAW counts) ───────────────────────────────
# Raw counts: the metacells vignette normalizes (CP10k + log1p) and aggregates
# live; raw integers also compress far better than normalized doubles.
counts <- readMM(file.path(gen_dir, "counts_hvg_raw.mtx"))
barcodes <- scan(file.path(gen_dir, "barcodes.txt"), what = character(), quiet = TRUE)
genes <- read.csv(file.path(gen_dir, "genes_hvg.txt"), header = FALSE,
                  col.names = c("ensembl", "symbol", "symbol_unique"))
counts <- as(counts, "CsparseMatrix")
dimnames(counts) <- list(barcodes, genes$symbol_unique)
stopifnot(nrow(counts) == length(barcodes), ncol(counts) == nrow(genes))
saveRDS(counts, file.path(out_dir, "pbmc3k_counts_hvg.rds"))

# ── cluster labels ───────────────────────────────────────────────────────────
clusters <- read.csv(file.path(gen_dir, "clusters.csv"))
stopifnot(identical(clusters$barcode, barcodes))
cluster_f <- factor(paste0("c", clusters$cluster))
names(cluster_f) <- clusters$barcode
saveRDS(cluster_f, file.path(out_dir, "pbmc3k_clusters.rds"))

# ── cEBMF fixture (factor scores + gene loadings) ───────────────────────────
L <- as.matrix(read.csv(file.path(gen_dir, "cebmf_L.csv")))
rownames(L) <- barcodes
F_mat <- read.csv(file.path(gen_dir, "cebmf_F.csv"))
stopifnot(identical(F_mat$symbol, genes$symbol_unique))
cebmf <- list(
  L = L,                       # cells x K cEBMF factor scores
  F = as.matrix(F_mat[, -1L]), # genes x K cEBMF loadings
  K = ncol(L),
  K_converged_note = "cebmf_torch cEBMF(maxit = 50), CPU, seed 20261020"
)
rownames(cebmf$F) <- F_mat$symbol
saveRDS(cebmf, file.path(out_dir, "pbmc3k_cebmf.rds"))

# ── provenance: metadata copy + SHA-256 manifest ────────────────────────────
file.copy(file.path(gen_dir, "run_metadata.json"),
          file.path(out_dir, "run_metadata.json"), overwrite = TRUE)

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
  stop("Need openssl or sha256sum.", call. = FALSE)
}

fixtures <- file.path(out_dir, c("pbmc3k_counts_hvg.rds", "pbmc3k_clusters.rds",
                                 "pbmc3k_cebmf.rds", "run_metadata.json"))
manifest <- vapply(fixtures, .sha256_file, character(1L))
writeLines(paste(manifest, basename(fixtures), sep = "  "),
           file.path(out_dir, "SHA256SUMS"))
cat("fixtures written:\n")
print(file.info(fixtures)[, "size", drop = FALSE])
cat(paste(manifest, basename(fixtures), sep = "  "), sep = "\n")
