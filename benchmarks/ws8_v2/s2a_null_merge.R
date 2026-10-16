# WS-8 v2 S2a — merge the SPAC-excluded null-probe chunks and summarize.
#
# Combines null_c1/ (perms 1-25) and null_c2/ (perms 26-50) into
# s2a_nospac/null/, alongside the same-design real-data reference fit
# (null_ref/).  Mirrors the full-atlas S0 analysis: selection count per
# permutation, Clopper-Pearson bound on the saturation probability, and the
# real-data comparator.  Answers the S2a falsification hook for V8: does the
# null saturation persist when SPAC (the rank-deficient, PD-projected cohort)
# is removed?
#
# Usage: Rscript benchmarks/ws8_v2/s2a_null_merge.R   # from repo root

base <- "benchmarks/ws8_v2/results/s2a_nospac"
out_dir <- file.path(base, "null")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Merge chunk outputs ------------------------------------------------------
chunks <- c("null_c1", "null_c2")
res <- do.call(rbind, lapply(chunks, function(cc)
  read.csv(file.path(base, cc, "null_results.csv"), stringsAsFactors = FALSE)))
res <- res[order(res$perm), ]
stopifnot(nrow(res) == 50L, identical(res$perm, 1:50))
write.csv(res, file.path(out_dir, "null_results.csv"), row.names = FALSE)

for (cc in chunks) {
  fs <- list.files(file.path(base, cc), pattern = "^null_scores_perm.*[.]csv$",
                   full.names = TRUE)
  invisible(file.copy(fs, file.path(out_dir, basename(fs))))
}

# Real-data reference (same design, SPAC-excluded)
ref <- read.csv(file.path(base, "null_ref", "real_reference_scores.csv"),
                stringsAsFactors = FALSE)
file.copy(file.path(base, "null_ref", "real_reference_scores.csv"),
          file.path(out_dir, "real_reference_scores.csv"), overwrite = TRUE)
if (file.exists(file.path(base, "null_ref", "real_reference_env_scores.csv")))
  file.copy(file.path(base, "null_ref", "real_reference_env_scores.csv"),
            file.path(out_dir, "real_reference_env_scores.csv"),
            overwrite = TRUE)

# ---- Summary ------------------------------------------------------------------
n_sat <- sum(res$n_selected == 59L)
cp <- stats::binom.test(n_sat, nrow(res), p = 0.5, alternative = "greater")
real_n <- sum(ref$native)
cat(sprintf("S2a null (SPAC-excluded): %d/%d permutations select 59/59 features\n",
            n_sat, nrow(res)))
cat(sprintf("max min-score: median %.3f, range [%.3f, %.3f]\n",
            median(res$max_min_score), min(res$max_min_score),
            max(res$max_min_score)))
cat(sprintf("Clopper-Pearson 95%% lower bound for saturation probability: %.3f\n",
            cp$conf.int[1]))
cat(sprintf("real-data reference (SPAC-excluded): %d/59 native\n", real_n))

summ <- data.frame(
  design = "SPAC-excluded (7 cohorts, all n > p; no PD projection)",
  n_perm = nrow(res), n_saturated = n_sat,
  median_max_min_score = round(median(res$max_min_score), 4),
  cp_lower_95 = round(cp$conf.int[1], 4),
  real_native_n = real_n, stringsAsFactors = FALSE)
write.csv(summ, file.path(out_dir, "null_summary_s2a.csv"), row.names = FALSE)
cat("\nDone. Outputs in", out_dir, "\n")
