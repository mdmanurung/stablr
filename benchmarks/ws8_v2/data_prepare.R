# WS-8 v2 M1 — Rebuild donor-level tables from the raw portal file.
#
# Reconstructs the donor x cell-type count matrix and donor metadata directly
# from `all_cells_metadata.csv.gz`, enforcing the predeclared `Bad quality`
# exclusion (the v1 derivation script was never archived; this is the
# reproducible replacement).  Uses data.table for fast aggregation of the
# ~13.1M-row file.  Writes a data manifest with hashes, row counts, and
# exclusion logs.
#
# Usage: Rscript benchmarks/ws8_v2/data_prepare.R   # run from the repo root
#   env WS8_DATA_DIR  directory holding all_cells_metadata.csv.gz (default /workspace)
#   env WS8_V2_OUT    output dir (default benchmarks/ws8_v2/results/m1)
#
# Application case study — NOT calibration evidence.

suppressPackageStartupMessages(library(data.table))

data_dir <- Sys.getenv("WS8_DATA_DIR", "/workspace")
out_dir  <- Sys.getenv("WS8_V2_OUT", file.path("benchmarks/ws8_v2/results/m1"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- file.path(data_dir, "all_cells_metadata.csv.gz")
stopifnot(file.exists(raw))
md5 <- unname(tools::md5sum(raw))
cat(sprintf("input: %s\nmd5:   %s\n", raw, md5))
expected_md5 <- "50651968855bf6ca3efcdf199e2a9fd5"
if (md5 != expected_md5) warning("MD5 mismatch: expected ", expected_md5, ", got ", md5)

# ---- Read only the needed columns --------------------------------------------
keep_cols <- c("donor_id", "age", "sex", "ethnicity", "CMV", "Dataset", "annotation")
dt <- fread(raw, select = keep_cols, showProgress = FALSE)
cat(sprintf("cell rows read: %d\n", nrow(dt)))
cat("columns:", paste(names(dt), collapse = ", "), "\n")

# ---- Donor x cell-type counts (all annotation levels, incl. Bad quality) -----
counts <- dt[, .(n = .N), by = .(donor_id, annotation)]
setnames(counts, c("donor", "celltype", "n"))

# ---- One metadata row per donor ----------------------------------------------
meta <- unique(dt[, .(donor_id, age, sex, ethnicity, CMV, Dataset)],
               by = "donor_id")
# a donor should map to exactly one metadata combination; flag conflicts
n_meta_per_donor <- dt[, .(n_combo = uniqueN(paste(age, sex, ethnicity, CMV, Dataset))),
                       by = donor_id]
conflicts <- n_meta_per_donor[n_combo > 1L]
if (nrow(conflicts) > 0L) {
  warning(nrow(conflicts), " donors have inconsistent metadata rows; keeping first.")
}
setnames(meta, c("donor", "age", "sex", "ethnicity", "cmv", "dataset"))
meta[, age := as.numeric(age)]
cat(sprintf("donors: %d; donor x celltype pairs: %d; metadata conflicts: %d\n",
            nrow(meta), nrow(counts), nrow(conflicts)))
rm(dt); gc(verbose = FALSE)

# ---- Annotation levels and the Bad-quality exclusion --------------------------
ann_levels <- sort(unique(counts$celltype))
cat(sprintf("distinct annotation levels: %d\n", length(ann_levels)))
bad <- grep("Bad quality", ann_levels, value = TRUE)
cat("Bad-quality level(s):", paste(bad, collapse = ", "), "\n")
ct_levels <- setdiff(ann_levels, bad)
cat(sprintf("cell-type features after exclusion: %d\n", length(ct_levels)))

# ---- Build donor x cell-type matrix (Bad quality excluded from num + denom) ---
donor_levels <- sort(intersect(meta$donor, unique(counts$donor)))
ct_mat <- matrix(0L, nrow = length(donor_levels), ncol = length(ct_levels),
                 dimnames = list(donor_levels, ct_levels))
keep_ct <- counts$celltype %in% ct_levels
idx <- cbind(match(counts$donor[keep_ct], donor_levels),
             match(counts$celltype[keep_ct], ct_levels))
ct_mat[idx] <- counts$n[keep_ct]
total <- rowSums(ct_mat)  # denominator = total typed (non-Bad) cells per donor

# ---- Missingness log ----------------------------------------------------------
is_missing <- function(v) is.na(v) | v %in% c("", "No_info", "NA", "na", "Unknown")
miss <- data.frame(
  variable = c("age", "sex", "ethnicity", "cmv"),
  n_missing = c(sum(is.na(meta$age)), sum(is_missing(meta$sex)),
                sum(is_missing(meta$ethnicity)), sum(is_missing(meta$cmv))),
  stringsAsFactors = FALSE
)
miss_by_cohort <- meta[, .(
  n = .N,
  age_na = sum(is.na(age)),
  sex_na = sum(is_missing(sex)),
  eth_na = sum(is_missing(ethnicity)),
  cmv_na = sum(is_missing(cmv))
), by = dataset]
setorder(miss_by_cohort, dataset)

# ---- Persist ------------------------------------------------------------------
saveRDS(list(counts = as.data.frame(counts), meta = as.data.frame(meta),
             ct_mat = ct_mat, total = total, ct_levels = ct_levels,
             bad_levels = bad),
        file.path(out_dir, "donor_tables.rds"))
write.csv(miss, file.path(out_dir, "missingness_overall.csv"), row.names = FALSE)
write.csv(as.data.frame(miss_by_cohort),
          file.path(out_dir, "missingness_by_cohort.csv"), row.names = FALSE)

manifest <- data.frame(
  field = c("input_file", "input_md5", "cell_rows", "n_donors",
            "n_annotation_levels", "bad_quality_levels", "n_celltype_features",
            "n_donor_celltype_pairs", "metadata_conflicts", "prepared_at"),
  value = c(basename(raw), md5, nrow(counts) == 0L, nrow(meta),
            length(ann_levels), paste(bad, collapse = ";"), length(ct_levels),
            nrow(counts), nrow(conflicts),
            format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  stringsAsFactors = FALSE
)
manifest$value[manifest$field == "cell_rows"] <- format(
  sum(counts$n), scientific = FALSE)
write.csv(manifest, file.path(out_dir, "data_manifest.csv"), row.names = FALSE)

cat("\n== M1 data_prepare summary ==\n")
cat("total typed cells (excl Bad quality):", format(sum(ct_mat), big.mark = ","), "\n")
cat("donors:", length(donor_levels), "\n")
cat("\nMissingness by cohort:\n")
print(as.data.frame(miss_by_cohort), row.names = FALSE)
cat("\nCohort sizes:\n")
print(sort(table(meta$dataset), decreasing = TRUE))
cat("\nDone. Outputs in ", out_dir, "\n")
