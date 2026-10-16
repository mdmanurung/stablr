# WS-8 v2 S2a — assemble the SPAC-excluded (7-cohort) LOCO results.
#
# Merges the per-method heldout_replication.csv files from
# results/s2a_nospac/run_m012 (M0-M2), run_m3, run_m4, run_m5 into one table
# and computes the same paired contrasts as the primary assembly (t-CI on
# n-1 = 6 df; exact sign-flip p), plus a side-by-side with the full-atlas
# amended results (results/run/).  SPAC is excluded throughout: every
# remaining cohort has n > p, so no knockoff PD projection is involved.
#
# Usage: Rscript benchmarks/ws8_v2/s2a_assemble.R   # run from the repo root
# Application case study — NOT calibration evidence.

base <- "benchmarks/ws8_v2/results"
s2a <- file.path(base, "s2a_nospac")
out_dir <- s2a

read_rep <- function(dir) {
  f <- file.path(s2a, dir, "heldout_replication.csv")
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE)
}
rep_df <- do.call(rbind, lapply(c("run_m012", "run_m3", "run_m4", "run_m5"),
                                read_rep))
rep_df <- rep_df[order(rep_df$method, rep_df$fold, rep_df$K), ]
write.csv(rep_df, file.path(out_dir, "heldout_replication_s2a.csv"),
          row.names = FALSE)

methods <- c("M0", "M1", "M2", "M3", "M4", "M5")
Ks <- c(10, 5, 15, 20)
macro <- sapply(Ks, function(K) {
  sub <- rep_df[rep_df$K == K, ]
  sapply(methods, function(m) mean(sub$R[sub$method == m]))
})
rownames(macro) <- methods; colnames(macro) <- paste0("K", Ks)
cat("== S2a (SPAC-excluded, 7 folds): macro R(K) ==\n")
print(round(macro, 3))

sign_flip_p <- function(d) {
  n <- length(d); obs <- mean(d)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  mean(abs(apply(signs, 1, function(s) mean(s * d))) >= abs(obs))
}
k10 <- rep_df[rep_df$K == 10, ]
folds <- sort(unique(k10$fold))
getR <- function(m) sapply(folds, function(f) k10$R[k10$method == m & k10$fold == f])

contrasts <- do.call(rbind, lapply(c("M0", "M1", "M2", "M3", "M5"), function(m) {
  d <- getR("M4") - getR(m)
  n <- length(d)
  tc <- stats::qt(0.975, n - 1) * stats::sd(d) / sqrt(n)
  data.frame(contrast = paste0("M4 - ", m), macro_diff = round(mean(d), 4),
             ci_low = round(mean(d) - tc, 4), ci_high = round(mean(d) + tc, 4),
             sign_flip_p = round(sign_flip_p(d), 4), n_cohorts = n,
             stringsAsFactors = FALSE)
}))
write.csv(contrasts, file.path(out_dir, "method_comparisons_s2a.csv"),
          row.names = FALSE)
cat("\n== S2a contrasts (paired over 7 cohorts) ==\n")
print(contrasts, row.names = FALSE)

# ---- Side-by-side with the full-atlas amended results -------------------------
full <- read.csv(file.path(base, "run", "heldout_replication.csv"),
                 stringsAsFactors = FALSE)
full10 <- tapply(full$R[full$K == 10], full$method[full$K == 10], mean)
sbs <- data.frame(method = methods,
                  full_atlas_R10 = round(as.numeric(full10[methods]), 3),
                  nospac_R10 = round(macro[methods, "K10"], 3),
                  delta = round(macro[methods, "K10"] - as.numeric(full10[methods]), 3))
write.csv(sbs, file.path(out_dir, "s2a_vs_full_atlas.csv"), row.names = FALSE)
cat("\n== Full atlas (8 cohorts) vs SPAC-excluded (7 cohorts), macro R@10 ==\n")
print(sbs, row.names = FALSE)
cat("\nDone. Outputs in", out_dir, "\n")
