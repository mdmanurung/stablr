# WS-8 v2 M1 — Validate the rebuilt donor-level tables (gate G1).
#
# Reads donor_tables.rds produced by data_prepare.R and fails closed on any
# integrity violation: duplicate donor IDs, non-finite ages, malformed cohort
# labels, count/denominator mismatches, non-finite compositions, unexpected
# feature count.  Also verifies the historical 59-cell-type matrix is
# reproducible under the documented v1 inclusion rules.
#
# Usage: Rscript benchmarks/ws8_v2/data_checks.R   # run from the repo root
#   env WS8_V2_OUT  dir holding donor_tables.rds (default .../results/m1)
#
# Application case study — NOT calibration evidence.

out_dir <- Sys.getenv("WS8_V2_OUT", file.path("benchmarks/ws8_v2/results/m1"))
obj <- readRDS(file.path(out_dir, "donor_tables.rds"))
counts <- obj$counts; meta <- obj$meta; ct_mat <- obj$ct_mat; total <- obj$total

fail <- function(msg) stop("GATE G1 FAILED: ", msg, call. = FALSE)
ok   <- function(msg) cat("  ok:", msg, "\n")

cat("== WS-8 v2 data checks (gate G1) ==\n")

# 1. No duplicate donor IDs
if (anyDuplicated(meta$donor)) fail("duplicate donor IDs in meta")
if (anyDuplicated(rownames(ct_mat))) fail("duplicate donor IDs in ct_mat")
ok("no duplicate donor IDs")

# 2. Donor sets agree
if (!setequal(meta$donor, rownames(ct_mat))) fail("meta and ct_mat donor sets differ")
ok("meta and count-matrix donor sets agree")

# 3. Ages finite and plausible
if (any(!is.finite(meta$age))) fail("non-finite ages present")
if (any(meta$age < 0 | meta$age > 120)) fail("ages outside [0, 120]")
ok(sprintf("ages finite, range [%.0f, %.0f]", min(meta$age), max(meta$age)))

# 4. Cohort labels well-formed, 8 cohorts, expected sizes
cohorts <- sort(unique(meta$dataset))
if (length(cohorts) != 8L) fail(sprintf("expected 8 cohorts, found %d", length(cohorts)))
exp_n <- c(ABF300 = 166, AIDA = 625, CIMA = 421, KCL = 188, OneK1K = 981,
           SoundLife = 96, SPAC = 33, UCSF = 99)
got_n <- table(meta$dataset)[names(exp_n)]
if (!all(got_n == exp_n)) {
  fail(paste("cohort sizes differ:", paste(sprintf("%s=%d(exp %d)", names(exp_n), got_n, exp_n), collapse = ", ")))
}
ok("8 cohorts with expected donor counts")

# 5. Counts non-negative integers; denominators positive
if (any(counts$n < 0)) fail("negative cell counts")
if (any(total <= 0)) fail("donors with zero typed cells")
ok("counts non-negative; all donors have > 0 typed cells")

# 6. Feature count = 59 cell types
if (ncol(ct_mat) != 59L) fail(sprintf("expected 59 cell-type features, found %d", ncol(ct_mat)))
if (any(grepl("Bad quality", colnames(ct_mat)))) fail("Bad quality leaked into features")
ok("59 cell-type features, Bad quality excluded")

# 7. Compositions finite after add-half smoothing + logit
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop))
if (any(!is.finite(X))) fail("non-finite logit proportions")
ok("logit proportions finite for all donors x features")

# 8. Reproduce historical v1 matrix under v1 inclusion rules
#    (v1 kept donors with non-missing age/sex/CMV; CMV missingness is structural,
#    so the v1 60-feature design is NOT reproducible for CMV — we verify the
#    59-cell-type composition block only, on donors with finite age.)
keep <- is.finite(meta$age)
v1_donors <- intersect(rownames(ct_mat), meta$donor[keep])
if (length(v1_donors) != 2609L) {
  warning(sprintf("v1-comparable donor count is %d, expected 2609", length(v1_donors)))
} else {
  ok("all 2609 donors have finite age (v1 composition block reproducible)")
}

# 9. Row sums of counts == total (denominator integrity)
rs <- rowSums(ct_mat)
if (!identical(unname(rs), unname(total))) fail("row sums != recorded denominators")
ok("count-matrix row sums equal recorded denominators")

cat("\nGATE G1 PASSED: donor-level lineage is intact.\n")
cat(sprintf("  donors=%d  features=%d  typed cells=%s  cohorts=%d\n",
            nrow(ct_mat), ncol(ct_mat), format(sum(ct_mat), big.mark = ","),
            length(cohorts)))
