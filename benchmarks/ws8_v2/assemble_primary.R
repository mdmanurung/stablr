# WS-8 v2 M3 — Assemble the combined primary comparison across all methods.
#
# Merges the per-method heldout_replication.csv files (M0/M1/M2 from run_m012b,
# M3 from run_m3b, M4 from run_m4b, M5 from run_m5b — the Amendment-1 amended
# re-run) into one table and computes the predeclared confirmatory contrasts:
# M4 vs M0 and M4 vs M3 at K=10, with cohort-level uncertainty (8
# generalization units) and an exact paired sign-flip permutation p-value.
# All other contrasts are exploratory.
#
# Usage: Rscript benchmarks/ws8_v2/assemble_primary.R   # run from the repo root
# Application case study — NOT calibration evidence.

base <- "benchmarks/ws8_v2/results"
out_dir <- file.path(base, "run")   # combined output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_rep <- function(dir, methods) {
  f <- file.path(base, dir, "heldout_replication.csv")
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d[d$method %in% methods, ]
}

rep_df <- do.call(rbind, list(
  read_rep("run_m012b", c("M0", "M1", "M2")),
  read_rep("run_m3b", "M3"),
  read_rep("run_m4b", "M4"),
  read_rep("run_m5b", "M5")
))
rep_df <- rep_df[order(rep_df$method, rep_df$fold, rep_df$K), ]
write.csv(rep_df, file.path(out_dir, "heldout_replication.csv"), row.names = FALSE)

# ---- Macro R(K) per method -----------------------------------------------------
methods <- c("M0", "M1", "M2", "M3", "M4", "M5")
Ks <- c(10, 5, 15, 20)
macro <- sapply(Ks, function(K) {
  sub <- rep_df[rep_df$K == K, ]
  sapply(methods, function(m) mean(sub$R[sub$method == m]))
})
rownames(macro) <- methods; colnames(macro) <- paste0("K", Ks)
cat("== PRIMARY: macro signed-rho replication R(K) ==\n")
print(round(macro, 3))

# ---- Confirmatory contrasts at K=10: M4 vs M0, M4 vs M3 ------------------------
# Uncertainty basis (protocol section 7): the 8 cohorts are the generalization
# units.  ci_low/ci_high are a paired t interval for the MEAN cohort difference
# (t on n-1 = 7 df; assumes approximate symmetry of the cohort differences;
# n = 8 is imprecise — read with the exact sign-flip p).  spread_low/spread_high
# are the 2.5-97.5% quantiles of the 8 paired differences (a descriptive spread
# of cohort effects, NOT a confidence interval for the mean).
sign_flip_p <- function(d) {
  # exact paired sign-flip permutation test on the mean difference
  n <- length(d); obs <- mean(d)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- apply(signs, 1, function(s) mean(s * d))
  mean(abs(null) >= abs(obs))
}
k10 <- rep_df[rep_df$K == 10, ]
folds <- sort(unique(k10$fold))
getR <- function(m) sapply(folds, function(f) k10$R[k10$method == m & k10$fold == f])

contrasts <- do.call(rbind, lapply(c("M0", "M1", "M2", "M3", "M5"), function(m) {
  d <- getR("M4") - getR(m)
  n <- length(d)
  tc <- stats::qt(0.975, n - 1) * stats::sd(d) / sqrt(n)
  data.frame(
    contrast = paste0("M4 - ", m),
    macro_diff = round(mean(d), 4),
    ci_low = round(mean(d) - tc, 4),
    ci_high = round(mean(d) + tc, 4),
    spread_low = round(unname(quantile(d, 0.025)), 4),
    spread_high = round(unname(quantile(d, 0.975)), 4),
    sign_flip_p = round(sign_flip_p(d), 4),
    n_cohorts = n,
    confirmatory = m %in% c("M0", "M3"),
    stringsAsFactors = FALSE)
}))
write.csv(contrasts, file.path(out_dir, "method_comparisons.csv"), row.names = FALSE)
cat("\n== M4 vs comparators (paired cohort-level difference in R@10) ==\n")
print(contrasts, row.names = FALSE)

# ---- Per-cohort R@10 table (all methods) ---------------------------------------
wide <- t(sapply(methods, function(m) getR(m)))
colnames(wide) <- folds
wide <- cbind(macro_R10 = round(rowMeans(wide), 3), round(wide, 3))
write.csv(wide, file.path(out_dir, "per_cohort_replication_K10.csv"))
cat("\n== Per-cohort R@10 ==\n")
print(wide)
cat("\nDone. Combined outputs in ", out_dir, "\n")
