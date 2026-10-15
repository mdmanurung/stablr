# WS-8 review follow-up 4: paired per-donor MAE differences between the
# all-feature and selected-panel LOCO models.
#
# For each donor, d = |error all-feature| - |error selected-panel|.  A
# positive d means the selected panel has the smaller error for that donor.
# We report the per-cohort mean of d with a 95% bootstrap interval (2,000
# donor-level resamples inside the cohort, matching the interval convention
# of ws8_figures.R), plus a pooled mean with a cohort-stratified bootstrap.
# Reads the objects saved by ws8_aging_showcase.R; no refitting.
#
# Note on model comparability: the all-feature model includes COV:sex; the
# 24-feature panel does not (sex was not selected).  The paired comparison
# therefore mixes a change of feature set with a change of covariate set.
#
# Usage: run from the repo root, e.g.
#   Rscript benchmarks/ws8_review_paired_mae.R
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)

obj_path <- Sys.getenv("WS8_OBJECTS", "/workspace/ws8_out/ws8_objects.rds")
out_dir  <- Sys.getenv("WS8_OUT_DIR", "benchmarks")
obj <- readRDS(obj_path)

la <- obj$loco_all   # donor-level: cohort, n, mae, slope, pred, obs
ls <- obj$loco_sel
stopifnot(identical(la$cohort, ls$cohort), identical(la$obs, ls$obs))

set.seed(20261013)
B <- 2000

cohorts <- unique(la$cohort)
rows <- lapply(cohorts, function(cv) {
  i  <- la$cohort == cv
  d  <- abs(la$pred[i] - la$obs[i]) - abs(ls$pred[i] - ls$obs[i])
  bs <- replicate(B, mean(d[sample.int(length(d), length(d), replace = TRUE)]))
  data.frame(cohort = cv, n = length(d),
             mean_ae_all = mean(abs(la$pred[i] - la$obs[i])),
             mean_ae_sel = mean(abs(ls$pred[i] - ls$obs[i])),
             mean_diff   = mean(d),
             ci_low  = unname(quantile(bs, 0.025)),
             ci_high = unname(quantile(bs, 0.975)),
             stringsAsFactors = FALSE)
})
paired <- do.call(rbind, rows)

# Pooled mean difference, cohort-stratified bootstrap.
d_all <- abs(la$pred - la$obs) - abs(ls$pred - ls$obs)
bs_pool <- replicate(B, {
  idx <- unlist(lapply(cohorts, function(cv) {
    i  <- which(la$cohort == cv)
    sample(i, length(i), replace = TRUE)
  }))
  mean(d_all[idx])
})
pooled <- data.frame(cohort = "pooled", n = length(d_all),
                     mean_ae_all = mean(abs(la$pred - la$obs)),
                     mean_ae_sel = mean(abs(ls$pred - ls$obs)),
                     mean_diff = mean(d_all),
                     ci_low = unname(quantile(bs_pool, 0.025)),
                     ci_high = unname(quantile(bs_pool, 0.975)),
                     stringsAsFactors = FALSE)

out <- rbind(pooled, paired)
write.csv(out, file.path(out_dir, "ws8_paired_mae_diff.csv"), row.names = FALSE)
cat("Paired per-donor MAE difference (all-feature minus selected; positive =\nselected panel better), 95% bootstrap CI, 2,000 resamples:\n")
print(out, row.names = FALSE)
