# WS-8 v3 — Within-cohort age-permutation null probe for the NEW STABL-machinery
# arms (M7 mint_stabl, M8 contrast_stabl), mirroring benchmarks/ws8_v2/null_probe.R.
#
# Shuffles age WITHIN each cohort (preserving features, cohort structure, and
# each cohort's age distribution), refits the selector on the FULL atlas, and
# records native panel size and score distribution.  10 permutations, locked
# n_boot = 200 (same as real fits — no reduced-budget deviation).
#
# Also fits the REAL full-atlas M7/M8 once (perm_id = 0): doubles as the
# full-atlas panel needed for the paper panel-recovery co-primary and as the
# real-data comparator for the null.
#
# Diagnostic null, NOT a calibrated FDR guarantee.
#
# Usage: Rscript null_probe_v3.R
#   env WS8_NULL_METHODS  comma subset of M7,M8 (default both)
#   env WS8_N_PERM        number of permutations (default 10)
#   env WS8_V3_NULL_OUT   output dir (default benchmarks/ws8_v3/results/null_v3)
#
# Application case study — NOT calibration evidence.

# Path config (repo-relative defaults; override via env vars).
v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
v2_dir <- Sys.getenv("WS8_V2_DIR", "benchmarks/ws8_v2")
in_dir <- Sys.getenv("WS8_V2_IN", file.path(v2_dir, "results", "m1"))
suppressPackageStartupMessages({
  library(stablr); library(glmnet); library(matrixStats); library(mixOmics)
})
engine_dir <- Sys.getenv("WS8_ENGINE", v2_dir)
source(file.path(engine_dir, "metrics.R"))
source(file.path(engine_dir, "loco.R"))
source(file.path(engine_dir, "selectors.R"))
source(file.path(v3_dir, "selectors_v3.R"))

out_dir <- Sys.getenv("WS8_V3_NULL_OUT", file.path(v3_dir, "results", "null_v3"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
methods <- strsplit(Sys.getenv("WS8_NULL_METHODS", "M7,M8"), ",")[[1]]
n_perm <- as.integer(Sys.getenv("WS8_N_PERM", "10"))

# ---- Design matrix (identical to run_v3.R) ------------------------------------
obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- meta[match(rownames(X), meta$donor), ]
y <- setNames(as.numeric(meta$age), meta$donor)
cohort <- setNames(as.character(meta$dataset), meta$donor)

# Same transform code path as the real fits; depends only on X, so the drop
# set is identical across permutations by construction.
recipe <- fit_fold_transform(X, cohort)
Xf <- apply_fold_transform(recipe, X)
cat(sprintf("null probe v3: %d x %d features, %d cohorts, methods=%s, perms=%d\n",
            nrow(Xf), ncol(Xf), length(unique(cohort)),
            paste(methods, collapse = ","), n_perm))

summ_rows <- list()
score_rows <- list()

fit_one <- function(m, yv, seed) {
  sel <- ws8v3_selectors[[m]]
  t0 <- Sys.time()
  res <- sel(Xf, yv, cohort, n_boot = 200L, seed = seed)
  rt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  list(res = res, runtime = rt)
}

for (m in methods) {
  # --- real full-atlas fit (perm 0) ---
  r0 <- fit_one(m, y, 20261013L)
  sc <- r0$res$scores
  summ_rows[[length(summ_rows) + 1L]] <- data.frame(
    method = m, perm = 0L, n_selected = length(r0$res$native),
    max_score = max(sc), mean_score = mean(sc), n_ge_0.3 = sum(sc >= 0.3),
    thr_star = unname(r0$res$extra$thr_star),
    min_fdr = unname(r0$res$extra$min_fdr), runtime_min = r0$runtime,
    stringsAsFactors = FALSE)
  score_rows[[length(score_rows) + 1L]] <- data.frame(
    method = m, perm = 0L, feature = names(sc), score = as.numeric(sc),
    native = names(sc) %in% r0$res$native, stringsAsFactors = FALSE)
  cat(sprintf("%s REAL: native=%d max=%.3f thr=%.2f (%.1f min)\n",
              m, length(r0$res$native), max(sc),
              unname(r0$res$extra$thr_star), r0$runtime))

  # --- permutations ---
  for (pi in seq_len(n_perm)) {
    set.seed(20261013L + pi)
    y_perm <- y
    for (cv in unique(cohort)) {
      idx <- which(cohort == cv)
      y_perm[idx] <- sample(y[idx])
    }
    rp <- fit_one(m, y_perm, 20261013L + pi)
    scp <- rp$res$scores
    summ_rows[[length(summ_rows) + 1L]] <- data.frame(
      method = m, perm = pi, n_selected = length(rp$res$native),
      max_score = max(scp), mean_score = mean(scp), n_ge_0.3 = sum(scp >= 0.3),
      thr_star = unname(rp$res$extra$thr_star),
      min_fdr = unname(rp$res$extra$min_fdr), runtime_min = rp$runtime,
      stringsAsFactors = FALSE)
    score_rows[[length(score_rows) + 1L]] <- data.frame(
      method = m, perm = pi, feature = names(scp), score = as.numeric(scp),
      native = names(scp) %in% rp$res$native, stringsAsFactors = FALSE)
    cat(sprintf("%s perm %2d: native=%2d max=%.3f thr=%.2f (%.1f min)\n",
                m, pi, length(rp$res$native), max(scp),
                unname(rp$res$extra$thr_star), rp$runtime))
  }
}

summ_df <- do.call(rbind, summ_rows)
score_df <- do.call(rbind, score_rows)
write.csv(summ_df, file.path(out_dir, "null_results_v3.csv"), row.names = FALSE)
write.csv(score_df, file.path(out_dir, "null_scores_v3.csv"), row.names = FALSE)

cat("\n== v3 null summary (native panel sizes; perm 0 = real) ==\n")
for (m in methods) {
  sub <- summ_df[summ_df$method == m, ]
  real_n <- sub$n_selected[sub$perm == 0]
  null_n <- sub$n_selected[sub$perm > 0]
  cat(sprintf("%s: real native=%d | null median %.1f range [%d, %d] | null max_score median %.3f\n",
              m, real_n, median(null_n), min(null_n), max(null_n),
              median(sub$max_score[sub$perm > 0])))
}
cat("Done. Outputs in ", out_dir, "\n")
