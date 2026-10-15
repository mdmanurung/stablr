# WS-8 v2 M5 / S0 — Within-cohort age-permutation null probe.
#
# Shuffles age WITHIN each cohort (preserving cell features, cohort structure,
# correlations, and each cohort's age distribution), refits the full M4
# environment-aware STABL selector on the permuted data, and records the native
# panel size and score distribution.  Repeated over independent seeds to build
# an empirical null for the number of selected features and the max min-score.
#
# This is a diagnostic null, NOT a calibrated FDR guarantee.  Real-data and
# null curves are kept separate.
#
# Usage: Rscript inst/benchmarks/ws8_v2/null_probe.R
#   env WS8_V2_IN    dir holding donor_tables.rds (default .../results/m1)
#   env WS8_V2_OUT   output dir (default .../results/null)
#   env WS8_N_PERM   number of permutations (default 10; smoke before scaling)
#   env WS8_PERM_START  first permutation index (default 1) for chunking
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("inst/benchmarks/ws8_v2/selectors.R")

in_dir  <- Sys.getenv("WS8_V2_IN", file.path("inst/benchmarks/ws8_v2/results/m1"))
out_dir <- Sys.getenv("WS8_V2_OUT", file.path("inst/benchmarks/ws8_v2/results/null"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
n_perm    <- as.integer(Sys.getenv("WS8_N_PERM", "10"))
perm_start <- as.integer(Sys.getenv("WS8_PERM_START", "1"))

# ---- Design matrix ------------------------------------------------------------
obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- meta[match(rownames(X), meta$donor), ]
y <- setNames(as.numeric(meta$age), meta$donor)
cohort <- setNames(as.character(meta$dataset), meta$donor)
cat(sprintf("null probe: %d donors x %d features, %d cohorts, %d permutations (start %d)\n",
            nrow(X), ncol(X), length(unique(cohort)), n_perm, perm_start))

# ---- Permutation loop ---------------------------------------------------------
rows <- list()
for (pi in seq_len(n_perm)) {
  perm_id <- perm_start + pi - 1L
  set.seed(20261013 + perm_id)  # independent stream per permutation
  t0 <- Sys.time()
  # shuffle age WITHIN each cohort
  y_perm <- y
  for (cv in unique(cohort)) {
    idx <- which(cohort == cv)
    y_perm[idx] <- sample(y[idx])
  }
  # refit the full M4 selector on permuted data (full atlas, all 8 cohorts)
  es <- stabl_environment_stability(
    x = X, y = y_perm, environments = cohort, lambda_grid = "auto",
    base_learner = "rf", artificial_type = "knockoff_equi",
    n_bootstraps = 200L, sample_fraction = 1, replace = TRUE,
    random_state = 20261013 + perm_id)
  mins <- get_min_scores(es)
  rows[[pi]] <- data.frame(
    perm = perm_id,
    n_selected = sum(mins >= 0.5),
    max_min_score = max(mins),
    mean_min_score = mean(mins),
    n_ge_0.3 = sum(mins >= 0.3),
    stringsAsFactors = FALSE)
  # per-feature null scores for distribution
  write.csv(data.frame(perm = perm_id, feature = names(mins),
                       min_score = as.numeric(mins), stringsAsFactors = FALSE),
            file.path(out_dir, sprintf("null_scores_perm%03d.csv", perm_id)),
            row.names = FALSE)
  cat(sprintf("perm %3d: n_selected=%2d  max_min=%.3f  (%.1f min)\n",
              perm_id, sum(mins >= 0.5), max(mins),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

null_df <- do.call(rbind, rows)
write.csv(null_df, file.path(out_dir, "null_results.csv"), row.names = FALSE)

cat("\n== S0 null summary (M4 native panel, min-agg >= 0.5) ==\n")
cat(sprintf("permutations: %d\n", nrow(null_df)))
cat(sprintf("n_selected: median %.1f, range [%d, %d]\n",
            median(null_df$n_selected), min(null_df$n_selected), max(null_df$n_selected)))
cat(sprintf("max min-score: median %.3f, range [%.3f, %.3f]\n",
            median(null_df$max_min_score), min(null_df$max_min_score), max(null_df$max_min_score)))
cat("\nReal-data reference: 24 features selected, max min-score 1.000.\n")
cat("Done. Outputs in ", out_dir, "\n")
