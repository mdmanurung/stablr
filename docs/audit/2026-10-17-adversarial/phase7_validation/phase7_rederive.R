# Phase 7 — independent re-derivation of the aging case-study saturation chain
# from the committed ws8_v2 artifacts (compendium repo), plus vignette-claim
# spot checks. Adversarial audit: numbers are recomputed from CSVs, never
# copied from the claims ledger.
#
# Artifacts: /workspace/stablr-reproducibility/results/ws8_v2/
# Vignette:  vignettes/stablr-aging-case-study.Rmd (stablr repo)
# Ledger:    benchmarks/ws8_v2/claims_ledger.md (V8, V12, V13; V1-V7, V11, V15)

base <- "/workspace/stablr-reproducibility/results/ws8_v2"

cat("==================================================================\n")
cat("1. V8 — S0 within-cohort age-permutation null (8 cohorts)\n")
cat("==================================================================\n")
nr <- read.csv(file.path(base, "null", "null_results.csv"))
cat("permutations committed:", nrow(nr), "\n")
cat("perms selecting 59/59 :", sum(nr$n_selected == 59L), "\n")
cat("max_min_score: median", round(median(nr$max_min_score), 3),
    " range [", round(min(nr$max_min_score), 3), ",",
    round(max(nr$max_min_score), 3), "]\n")
# Clopper-Pearson lower bound, x = n = 50:
# two-sided (alpha/2 lower tail): p^(n) = 0.025  -> 0.025^(1/50)
# one-sided (alpha lower tail):   p^(n) = 0.05   -> 0.05^(1/50)
cat("CP lower bound, two-sided 0.025^(1/50):", round(0.025^(1/50), 4),
    " (ledger V8 quotes 0.929)\n")
cat("CP lower bound, one-sided  0.05^(1/50):", round(0.05^(1/50), 4),
    " (ledger V12 quotes 0.942)\n")
rr <- read.csv(file.path(base, "null", "real_reference_scores.csv"))
cat("same-design real-data reference: ", sum(rr$native), "/ 59 native",
    " (ledger V8: 35/59)\n")
p1 <- read.csv(file.path(base, "null", "null_scores_perm001.csv"))
cat("perm001 per-feature min_score: median", round(median(p1$min_score), 3),
    " min", round(min(p1$min_score), 3), "\n")

cat("\n")
cat("==================================================================\n")
cat("2. V12 — S2a SPAC-excluded null (7 cohorts)\n")
cat("==================================================================\n")
files <- list.files(file.path(base, "s2a_nospac", "null"),
                    pattern = "^null_scores_perm.*\\.csv$", full.names = TRUE)
cat("permutation files:", length(files), "\n")
sat <- 0L
maxes <- numeric(length(files))
for (i in seq_along(files)) {
  d <- read.csv(files[i])
  maxes[i] <- max(d$min_score)
  sat <- sat + (sum(d$min_score >= 0.5) == 59L)
}
cat("perms selecting 59/59 :", sat, "\n")
cat("max min-score: median", round(median(maxes), 3),
    " range [", round(min(maxes), 3), ",", round(max(maxes), 3), "]\n")
ref2 <- read.csv(file.path(base, "s2a_nospac", "null_ref",
                           "real_reference_scores.csv"))
cat("real reference (no SPAC):", sum(ref2$native), "/ 59 native",
    " (ledger V12: 30/59)\n")

cat("\n")
cat("==================================================================\n")
cat("3. V13 — S1 synthetic saturation (4 truth regimes x 5 replicates)\n")
cat("==================================================================\n")
s1 <- read.csv(file.path(base, "s1", "s1_summary.csv"))
m45 <- s1[s1$method %in% c("M4", "M5"), ]
cat("M4/M5 n_native by scenario (must be 59 everywhere):\n")
print(tapply(m45$n_native, m45$scenario, unique))
cat("emp_fdp_native range:", round(min(m45$emp_fdp_native), 2), "-",
    round(max(m45$emp_fdp_native), 2), " (ledger V13: 0.83-1.00)\n")
m1 <- s1[s1$method == "M1", ]
cat("M1 n_native by scenario (calibration contrast):\n")
print(m1[, c("scenario", "n_native", "emp_fdp_native")])

cat("\n")
cat("==================================================================\n")
cat("4. Vignette section 10.3 — macro R@10 and confirmatory contrasts\n")
cat("==================================================================\n")
hr <- read.csv(file.path(base, "run", "heldout_replication.csv"))
k10 <- hr[hr$K == 10, ]
mac <- tapply(k10$R, k10$method, mean)
cat("macro R@10 (vignette table: M0 0.395, M1 0.386, M5 0.321,",
    "M4 0.286, M3 0.142, M2 0.118):\n")
print(round(mac[order(-mac)], 3))
mc <- read.csv(file.path(base, "run", "method_comparisons.csv"))
cat("M4 contrasts (vignette: M4-M3 +0.145 [+0.083,+0.206] p .0078;",
    "M4-M0 -0.108 [-0.174,-0.043] p .016):\n")
print(mc[mc$contrast %in% c("M4 - M3", "M4 - M0", "M4 - M1", "M4 - M2", "M4 - M5"),
         c("contrast", "macro_diff", "ci_low", "ci_high", "sign_flip_p",
           "confirmatory")])

cat("\n")
cat("==================================================================\n")
cat("5. Vignette section 10.5 — prediction panel summary\n")
cat("==================================================================\n")
pp <- read.csv(file.path(base, "run", "prediction_panel_summary.csv"))
cat("(vignette: topk 12.19/13.56, all59 12.35/14.35, median 20.11/22.74;",
    "slopes_po 0.32-0.68)\n")
print(pp[, c("panel", "macro_mae", "micro_mae", "slope_po_min", "slope_po_max")])

cat("\n")
cat("==================================================================\n")
cat("6. Vignette section 4 — v1 feature claims (ws8_feature_stats.csv)\n")
cat("==================================================================\n")
fs <- read.csv(file.path("/workspace/integrate/stablr/benchmarks",
                         "ws8_feature_stats.csv"))
cat("selected:", sum(fs$selected), "of", nrow(fs), "(vignette: 24 of 60)\n")
pre <- fs[fs$paper == TRUE, ]
cat("predeclared min_scores:",
    paste(sort(pre$min_score, decreasing = TRUE), collapse = ", "), "\n")
cat("(vignette Table 2: 1.000, 1.000, 0.995, 0.755, 0.740, 0.535, 0.505, 0.445)\n")
cat("predeclared selected:", sum(pre$selected), "of", nrow(pre),
    "(vignette: 7 of 8); sign_k values:", paste(unique(pre$sign_k), collapse = ","),
    "(vignette: all 8/8)\n")
sel <- fs[fs$selected == TRUE, ]
cat("consistent features with sign_p < 0.05:", sum(sel$sign_p < 0.05),
    "(vignette: 11 of 24)\n")

cat("\nALL RE-DERIVATIONS COMPLETE\n")
