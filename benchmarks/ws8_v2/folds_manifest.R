# WS-8 v2 — Fold definition table (protocol sections 6/10 conformance).
#
# Writes results/folds.csv: the donor-level assignment of the 8 fixed outer
# LOCO folds (one held-out cohort per fold, alphabetical; the training set is
# the complement).  Deterministic given the cohort vector; no fitting here.
# Per-method per-fold manifests (ID hashes, drop sets, panel sizes, seed) are
# written by run.R as fold_manifest.csv in each run directory.
#
# Usage: Rscript benchmarks/ws8_v2/folds_manifest.R   # run from repo root

source("benchmarks/ws8_v2/loco.R")

in_dir <- Sys.getenv("WS8_V2_IN", file.path("benchmarks/ws8_v2/results/m1"))
out_dir <- Sys.getenv("WS8_V2_OUT", file.path("benchmarks/ws8_v2/results"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
meta <- obj$meta
folds <- make_loco_folds(meta$donor, meta$dataset)
write.csv(folds, file.path(out_dir, "folds.csv"), row.names = FALSE)
cat(sprintf("folds.csv: %d donors x %d folds; test sizes: %s\n",
            nrow(folds), length(unique(folds$fold)),
            paste(sprintf("%s=%d", names(table(folds$fold)),
                          as.integer(table(folds$fold))), collapse = ", ")))
cat("Done. Output in ", out_dir, "\n")
