# WS-8 v2 M6 — Biological interpretation (secondary).
#
# Reads the committed v2 results and the M1 donor tables.  Reports: (1) the
# eight a priori reference populations' held-out directional evidence; (2) the
# post-hoc candidate list; (3) the GZMK+/GZMB+ CD8 Tem ratio (secondary,
# predefined); (4) per-cohort heterogeneity of the core features.  No causal or
# clinical claims.  The eight reference populations are an a priori concordance
# set, NOT the paper's complete clock (Table S1 unavailable).
#
# Usage: Rscript benchmarks/ws8_v2/biology.R   # run from the repo root
# Application case study — NOT calibration evidence.

base <- "benchmarks/ws8_v2/results"
out_dir <- file.path(base, "biology")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(base, "m1", "donor_tables.rds"))
ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- meta[match(rownames(X), meta$donor), ]
y <- as.numeric(meta$age); cohort <- meta$dataset
cohorts <- sort(unique(cohort))

# Per-cohort Spearman rho for every feature
rho_mat <- sapply(cohorts, function(cv) {
  idx <- cohort == cv
  apply(X[idx, ], 2, function(f) suppressWarnings(cor(f, y[idx], method = "spearman")))
})
rownames(rho_mat) <- colnames(X)
write.csv(round(rho_mat, 4), file.path(out_dir, "per_cohort_rho.csv"))

# ---- 1. Reference populations -------------------------------------------------
ref <- c("CT:CD8+ T Naive", "CT:CD8+ RTE", "CT:CD4+ RTE", "CT:CD4+ Temra",
         "CT:MAIT", "CT:CD8+ Temra", "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")
ref <- intersect(ref, rownames(rho_mat))
ref_tab <- data.frame(
  feature = ref,
  mean_rho = round(rowMeans(rho_mat[ref, ]), 3),
  min_rho  = round(apply(rho_mat[ref, ], 1, min), 3),
  max_rho  = round(apply(rho_mat[ref, ], 1, max), 3),
  n_pos    = rowSums(rho_mat[ref, ] > 0),
  n_neg    = rowSums(rho_mat[ref, ] < 0),
  sign_consistent = apply(rho_mat[ref, ], 1, function(r) all(r > 0) || all(r < 0)),
  stringsAsFactors = FALSE)
write.csv(ref_tab, file.path(out_dir, "reference_populations.csv"), row.names = FALSE)
cat("== Reference populations (per-cohort Spearman with age) ==\n")
print(ref_tab, row.names = FALSE)

# ---- 2. Post-hoc candidates ---------------------------------------------------
cand <- c("CT:CD4+ T HLA-DR+ memory", "CT:CD8+ IFN+ cells", "CT:CD8+ HLA-DR+",
          "CT:CD8+ Tcm CCR4+")
cand <- intersect(cand, rownames(rho_mat))
cand_tab <- data.frame(
  feature = cand,
  mean_rho = round(rowMeans(rho_mat[cand, ]), 3),
  sign_consistent = apply(rho_mat[cand, ], 1, function(r) all(r > 0) || all(r < 0)),
  stringsAsFactors = FALSE)
write.csv(cand_tab, file.path(out_dir, "candidate_features.csv"), row.names = FALSE)
cat("\n== Post-hoc candidates (hypothesis-generating only) ==\n")
print(cand_tab, row.names = FALSE)

# ---- 3. GZMK+/GZMB+ CD8 Tem ratio (secondary, predefined) ---------------------
if (all(c("CD8+ Tem GZMK+", "CD8+ Tem GZMB+") %in% colnames(ct_mat))) {
  gzk <- ct_mat[, "CD8+ Tem GZMK+"]; gzb <- ct_mat[, "CD8+ Tem GZMB+"]
  ratio <- log1p((gzk + 0.5) / (gzb + 0.5))
  ratio_rho <- sapply(cohorts, function(cv) {
    idx <- cohort == cv
    suppressWarnings(cor(ratio[idx], y[idx], method = "spearman"))
  })
  ratio_tab <- data.frame(cohort = cohorts, ratio_rho = round(ratio_rho, 3),
                          stringsAsFactors = FALSE)
  write.csv(ratio_tab, file.path(out_dir, "gzmk_gzmb_ratio.csv"), row.names = FALSE)
  cat("\n== GZMK+/GZMB+ CD8 Tem ratio vs age (per-cohort Spearman) ==\n")
  print(ratio_tab, row.names = FALSE)
}

# ---- 4. Heterogeneity of the core ---------------------------------------------
# I^2-like summary per reference feature: between-cohort variance of rho.
het <- data.frame(
  feature = ref,
  var_between = round(apply(rho_mat[ref, ], 1, var), 4),
  sd_between  = round(apply(rho_mat[ref, ], 1, sd), 3),
  stringsAsFactors = FALSE)
write.csv(het, file.path(out_dir, "heterogeneity.csv"), row.names = FALSE)
cat("\n== Between-cohort heterogeneity (reference features) ==\n")
print(het, row.names = FALSE)

cat("\nDone. Outputs in ", out_dir, "\n")
