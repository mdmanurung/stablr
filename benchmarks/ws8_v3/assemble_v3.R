# WS-8 v3 — Assemble candidate benchmark: primary replication, contrasts,
# paper panel-recovery co-primary, reference populations, prediction, nulls.
#
# Inputs (all produced on the pinned stablr commit a8f3a58):
#   results/eq_m012, eq_m4, eq_m35        v2 arms rerun like-for-like
#   results/run_v3_m6, run_v3_m7, run_v3_m8, run_v3_m9   new arms
#   results/full_atlas/full_atlas_scores.csv             full-atlas panels
#   results/null_v3/                       M7/M8 nulls + real full-atlas fits
#   v2_results/run/                        historical v2 (older package commits)
#
# Outputs: <out>/{tables,figures} (default benchmarks/ws8_v3/results; env WS8_V3_FIGOUT)
#
# Application case study — NOT calibration evidence.

# Path config (repo-relative defaults; override via env vars).
v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
v2_dir <- Sys.getenv("WS8_V2_DIR", "benchmarks/ws8_v2")
in_dir <- Sys.getenv("WS8_V2_IN", file.path(v2_dir, "results", "m1"))
suppressPackageStartupMessages({
  library(ggplot2); library(data.table)
})
engine_dir <- Sys.getenv("WS8_ENGINE", v2_dir)
source(file.path(engine_dir, "metrics.R"))
source(file.path(engine_dir, "loco.R"))   # for training_signs

res_dir <- file.path(v3_dir, "results")
out <- Sys.getenv("WS8_V3_FIGOUT", file.path(v3_dir, "results"))
dir.create(file.path(out, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, "figures"), showWarnings = FALSE, recursive = TRUE)

METHODS <- c("M0","M1","M2","M3","M4","M5","M6","M7","M7c","M8","M9","M9s","M9b")
METHOD_LABELS <- c(
  M0 = "M0 RE-meta ranking", M1 = "M1 cohort vote", M2 = "M2 pooled lasso",
  M3 = "M3 pooled STABL", M4 = "M4 env-STABL (min)", M5 = "M5 env-STABL (vote)",
  M6 = "M6 het-filtered RE-meta", M7 = "M7 MINT-STABL",
  M7c = "M7c MINT tuned (no stability)", M8 = "M8 contrast-STABL",
  M9 = "M9 ComBat+STABL", M9s = "M9s ComBat std-only",
  M9b = "M9b ComBat(no-age)+STABL")
NEW_ARMS <- c("M6", "M7", "M7c", "M8", "M9", "M9s", "M9b")

# Paper's reduced clock panel (20 populations; authors' repo notebook cell 32,
# all verified to map onto the harmonized 59 features).  Concordance target,
# NOT ground truth.
P20 <- c("CT:CD4+ RTE", "CT:CD4+ Th2", "CT:CD4+ T HLA-DR+ memory",
         "CT:CD4+ Treg memory", "CT:CD8+ RTE", "CT:CD8+ T Naive",
         "CT:CD8+ Tcm CCR4+", "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+",
         "CT:CD8+ HLA-DR+", "CT:CD8+ Tmem KLRC2+", "CT:gd T naive",
         "CT:gd T Vd1 GZMK+", "CT:gd T Vd1 GZMB+", "CT:gd T Vd2 GZMK+",
         "CT:NK CD56bright", "CT:Myeloid Non-classical monocytes",
         "CT:CD8+ Temra", "CT:CD4+ Temra", "CT:CD4+ Terminal effector")

# 8 curated reference populations, sign-consistent across all 8 cohorts (v2).
REFPOPS <- c("CT:CD8+ T Naive" = -1L, "CT:CD8+ RTE" = -1L, "CT:CD4+ RTE" = -1L,
             "CT:CD4+ Temra" = 1L, "CT:MAIT" = -1L, "CT:CD8+ Temra" = 1L,
             "CT:CD8+ Tem GZMK+" = 1L, "CT:CD8+ Tem GZMB+" = 1L)
GZMK_GZMB <- c("CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")

# ---- 1. Combine LOCO replication results --------------------------------------
run_dirs <- c("eq_m012", "eq_m4", "eq_m35", "run_v3_m6", "run_v3_m7",
              "run_v3_m8", "run_v3_m9", "run_v3_m9b")
read_piece <- function(dir, file) {
  f <- file.path(res_dir, dir, file)
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}
rep_df <- do.call(rbind, lapply(run_dirs, read_piece, file = "heldout_replication.csv"))
feat_df <- do.call(rbind, lapply(run_dirs, read_piece, file = "feature_scores.csv"))
met_df <- do.call(rbind, lapply(run_dirs, read_piece, file = "prediction_metrics.csv"))
fman_df <- do.call(rbind, lapply(run_dirs, read_piece, file = "fold_manifest.csv"))
rep_df <- rep_df[order(rep_df$method, rep_df$fold, rep_df$K), ]
stopifnot(!any(duplicated(rep_df[, c("method", "fold", "K")])))
write.csv(rep_df, file.path(out, "tables", "heldout_replication_v3.csv"),
          row.names = FALSE)

Ks <- c(10, 5, 15, 20)
folds <- sort(unique(rep_df$fold))
macro <- sapply(Ks, function(K) {
  sub <- rep_df[rep_df$K == K, ]
  sapply(METHODS, function(m) mean(sub$R[sub$method == m]))
})
rownames(macro) <- METHODS; colnames(macro) <- paste0("K", Ks)
write.csv(round(macro, 4), file.path(out, "tables", "macro_R_by_K_v3.csv"))
cat("== PRIMARY (pinned package, like-for-like): macro signed-rho R(K) ==\n")
print(round(macro, 3))

k10 <- rep_df[rep_df$K == 10, ]
getR <- function(m) sapply(folds, function(f) k10$R[k10$method == m & k10$fold == f])
wide <- t(sapply(METHODS, getR)); colnames(wide) <- folds
wide <- cbind(macro_R10 = round(rowMeans(wide), 3), round(wide, 3))
write.csv(wide, file.path(out, "tables", "per_cohort_replication_K10_v3.csv"))

# ---- 2. Version-drift reference (v2 historical vs v3 pinned) -------------------
v2 <- read.csv(file.path(v2_dir, "results", "run", "heldout_replication.csv"))
v2 <- v2[v2$K == 10, ]
drift <- data.frame(
  method = paste0("M", 0:5),
  v2_macro_R10 = sapply(paste0("M", 0:5), function(m) mean(v2$R[v2$method == m])),
  v3_macro_R10 = sapply(paste0("M", 0:5), function(m) mean(k10$R[k10$method == m])))
drift$delta <- round(drift$v3_macro_R10 - drift$v2_macro_R10, 4)
write.csv(drift, file.path(out, "tables", "version_drift_reference.csv"),
          row.names = FALSE)
cat("\n== Package version drift (v2 commits -> pinned a8f3a58), macro R@10 ==\n")
print(drift)

# ---- 3. Paired contrasts at K=10 (v2 methodology: t CI + exact sign-flip) -----
sign_flip_p <- function(d) {
  n <- length(d); obs <- mean(d)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- apply(signs, 1, function(s) mean(s * d))
  mean(abs(null) >= abs(obs))
}
contrast_row <- function(a, b, confirmatory) {
  d <- getR(a) - getR(b)
  n <- length(d)
  tc <- stats::qt(0.975, n - 1) * stats::sd(d) / sqrt(n)
  data.frame(contrast = paste0(a, " - ", b),
             macro_diff = round(mean(d), 4),
             ci_low = round(mean(d) - tc, 4), ci_high = round(mean(d) + tc, 4),
             spread_low = round(unname(quantile(d, 0.025)), 4),
             spread_high = round(unname(quantile(d, 0.975)), 4),
             sign_flip_p = round(sign_flip_p(d), 4), n_cohorts = n,
             confirmatory = confirmatory, stringsAsFactors = FALSE)
}
cl <- list()
for (m in c("M6", "M7", "M8", "M9")) {           # predeclared candidate contrasts
  cl[[length(cl) + 1L]] <- contrast_row(m, "M0", TRUE)
  cl[[length(cl) + 1L]] <- contrast_row(m, "M4", TRUE)
}
cl[[length(cl) + 1L]] <- contrast_row("M4", "M0", TRUE)   # v2 confirmatory, rerun
cl[[length(cl) + 1L]] <- contrast_row("M4", "M3", TRUE)   # v2 confirmatory, rerun
cl[[length(cl) + 1L]] <- contrast_row("M7", "M7c", FALSE) # stability contribution
cl[[length(cl) + 1L]] <- contrast_row("M9", "M9s", FALSE) # transform sensitivity
cl[[length(cl) + 1L]] <- contrast_row("M9", "M3", FALSE)  # ComBat contribution
cl[[length(cl) + 1L]] <- contrast_row("M9b", "M9", FALSE) # age-protection sensitivity
cl[[length(cl) + 1L]] <- contrast_row("M7c", "M0", FALSE)
contrasts <- do.call(rbind, cl)
write.csv(contrasts, file.path(out, "tables", "method_comparisons_v3.csv"),
          row.names = FALSE)
cat("\n== Contrasts (paired cohort-level difference in R@10) ==\n")
print(contrasts, row.names = FALSE)

# ---- 4. Co-primary: recovery of the paper's P20 clock panel --------------------
# Per-fold: panel = top-K at K=10 (and K=20 sensitivity) from feature_scores.
rec_rows <- list()
for (K in c(10, 20)) {
  fk <- feat_df[feat_df$K == K, ]
  for (m in METHODS) {
    for (f in folds) {
      panel <- fk$feature[fk$method == m & fk$fold == f]
      hit <- intersect(panel, P20)
      rec_rows[[length(rec_rows) + 1L]] <- data.frame(
        method = m, fold = f, K = K, n_hit = length(hit),
        recovery = length(hit) / length(P20), precision = length(hit) / K,
        stringsAsFactors = FALSE)
    }
  }
}
rec_df <- do.call(rbind, rec_rows)
write.csv(rec_df, file.path(out, "tables", "panel_recovery_per_fold_v3.csv"),
          row.names = FALSE)
rec_macro <- aggregate(cbind(recovery, precision) ~ method + K, rec_df, mean)
rec_macro <- rec_macro[order(rec_macro$K, -rec_macro$recovery), ]
write.csv(rec_macro, file.path(out, "tables", "panel_recovery_macro_v3.csv"),
          row.names = FALSE)
cat("\n== CO-PRIMARY: macro P20 recovery (LOCO top-K panels) ==\n")
print(rec_macro[rec_macro$K == 10, ], row.names = FALSE)

# ---- 5. Full-atlas panels + reference populations + GZMK/GZMB ------------------
fa <- do.call(rbind, lapply(c("full_atlas", "full_atlas_fast", "full_atlas_m9b"), function(d) {
  f <- file.path(res_dir, d, "full_atlas_scores.csv")
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}))
null_sc <- read.csv(file.path(res_dir, "null_v3", "null_scores_v3.csv"))
fa_m78 <- null_sc[null_sc$perm == 0, c("method", "feature", "score", "native")]
fa_m78$rank <- ave(fa_m78$score, fa_m78$method,
                   FUN = function(s) rank(-s, ties.method = "min"))
fa_all <- rbind(fa[, c("method", "feature", "score", "rank", "native")], fa_m78)

# full-atlas training signs (for direction concordance)
obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
prop <- (obj$ct_mat + 0.5) / (obj$total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- obj$meta[match(rownames(X), obj$meta$donor), ]
y <- setNames(as.numeric(meta$age), meta$donor)
cohort <- setNames(as.character(meta$dataset), meta$donor)
recipe <- fit_fold_transform(X, cohort)
Xf <- apply_fold_transform(recipe, X)
fa_signs <- training_signs(Xf, y, cohort)

fa_rows <- list()
for (m in METHODS) {
  sub <- fa_all[fa_all$method == m, ]
  if (nrow(sub) == 0) next
  top20 <- sub$feature[sub$rank <= 20]
  nat <- sub$feature[sub$native]
  ref_in_top <- intersect(names(REFPOPS), top20)
  ref_in_nat <- intersect(names(REFPOPS), nat)
  sign_ok_top <- sum(fa_signs[ref_in_top] == REFPOPS[ref_in_top])
  fa_rows[[length(fa_rows) + 1L]] <- data.frame(
    method = m,
    p20_native = length(intersect(nat, P20)),
    p20_top20 = length(intersect(top20, P20)),
    p20_recovery_top20 = length(intersect(top20, P20)) / 20,
    refpop_top20 = length(ref_in_top), refpop_native = length(ref_in_nat),
    refpop_sign_concordant_top20 = sign_ok_top,
    gzmk_rank = sub$rank[sub$feature == "CT:CD8+ Tem GZMK+"],
    gzmb_rank = sub$rank[sub$feature == "CT:CD8+ Tem GZMB+"],
    gzmk_native = "CT:CD8+ Tem GZMK+" %in% nat,
    gzmb_native = "CT:CD8+ Tem GZMB+" %in% nat,
    n_native = length(nat), stringsAsFactors = FALSE)
}
fa_summary <- do.call(rbind, fa_rows)
write.csv(fa_summary, file.path(out, "tables", "full_atlas_summary_v3.csv"),
          row.names = FALSE)
write.csv(fa_all, file.path(out, "tables", "full_atlas_scores_v3.csv"),
          row.names = FALSE)
cat("\n== Full-atlas summary ==\n")
print(fa_summary, row.names = FALSE)

# ---- 6. Prediction summary ------------------------------------------------------
met_macro <- aggregate(cbind(mae, rmse, rank_cor) ~ method + panel, met_df,
                       FUN = function(x) mean(x, na.rm = TRUE),
                       na.action = na.pass)
met_macro <- met_macro[order(met_macro$panel, met_macro$mae), ]
write.csv(met_macro, file.path(out, "tables", "prediction_summary_v3.csv"),
          row.names = FALSE)

# ---- 7. Null summary (M7/M8) ----------------------------------------------------
null_f <- file.path(res_dir, "null_v3", "null_results_v3.csv")
if (file.exists(null_f)) {
  nd <- read.csv(null_f)
  null_sum <- do.call(rbind, lapply(split(nd, nd$method), function(sub) {
    real <- sub[sub$perm == 0, ]; nul <- sub[sub$perm > 0, ]
    data.frame(method = sub$method[1], real_native = real$n_selected,
               null_native_median = median(nul$n_selected),
               null_native_max = max(nul$n_selected),
               real_max_score = real$max_score,
               null_max_score_max = max(nul$max_score),
               emp_p_native = (1 + sum(nul$n_selected >= real$n_selected)) /
                 (1 + nrow(nul)),
               stringsAsFactors = FALSE)
  }))
  write.csv(null_sum, file.path(out, "tables", "null_summary_v3.csv"),
            row.names = FALSE)
  cat("\n== Null probe summary ==\n"); print(null_sum, row.names = FALSE)
}

# ---- 8. Runtime ---------------------------------------------------------------
rt <- aggregate(runtime_min ~ method, fman_df, sum)
write.csv(rt, file.path(out, "tables", "runtime_v3.csv"), row.names = FALSE)

cat("\nAssembly done. Tables in ", file.path(out, "tables"), "\n")
