# WS-8 v2 M3 — Assemble the combined primary comparison across all methods.
#
# Merges the per-method heldout_replication.csv files (M0/M1/M2 from run_m012,
# M3 from the main-run log, M4 from run_m4, M5 from run_m5) into one table and
# computes the predeclared confirmatory contrasts: M4 vs M0 and M4 vs M3 at
# K=10, with cohort-level uncertainty (8 generalization units) and an exact
# paired sign-flip permutation p-value.  All other contrasts are exploratory.
#
# Usage: Rscript inst/benchmarks/ws8_v2/assemble_primary.R
# Application case study — NOT calibration evidence.

base <- "inst/benchmarks/ws8_v2/results"
out_dir <- file.path(base, "run")   # combined output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_rep <- function(dir, methods) {
  f <- file.path(base, dir, "heldout_replication.csv")
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d[d$method %in% methods, ]
}

rep_df <- do.call(rbind, list(
  read_rep("run_m012", c("M0", "M1", "M2")),
  read_rep("run_m4", "M4"),
  read_rep("run_m5", "M5")
))
# M3 comes from the main-run log (its CSV was not written before the run was
# stopped after completing all 8 folds).  Hard-code the verified per-fold values.
m3_vals <- c(ABF300 = 0.199, AIDA = 0.119, CIMA = 0.208, KCL = 0.157,
             OneK1K = 0.083, SoundLife = 0.230, SPAC = -0.059, UCSF = 0.199)
m3 <- expand.grid(method = "M3", fold = names(m3_vals), K = c(10, 5, 15, 20),
                  stringsAsFactors = FALSE)
m3 <- m3[m3$K == 10, ]
m3$R <- as.numeric(m3_vals[m3$fold])
m3$sign_agreement <- NA_real_; m3$n_panel <- 10L
m3 <- m3[, names(rep_df)]
rep_df <- rbind(rep_df, m3)
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
  data.frame(
    contrast = paste0("M4 - ", m),
    macro_diff = round(mean(d), 4),
    ci_low = round(unname(quantile(d, 0.025)), 4),
    ci_high = round(unname(quantile(d, 0.975)), 4),
    sign_flip_p = round(sign_flip_p(d), 4),
    n_cohorts = length(d),
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
