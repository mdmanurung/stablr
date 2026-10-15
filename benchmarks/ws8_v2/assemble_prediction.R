# WS-8 v2 M4 — Assemble the secondary predictive benchmark (protocol sec. 5).
#
# Reads the M4 run's prediction_metrics.csv and computes: macro AND micro
# (donor-weighted) MAE per panel; paired per-cohort MAE differences for the
# predeclared comparisons (top-K vs all-59, top-K vs median, native vs all-59)
# with paired t intervals (7 df; n = 8 cohorts, imprecise), descriptive spread
# quantiles, and exact sign-flip p-values; and a calibration summary (both
# slopes and bias ranges).  Negative differences favor the first-named panel.
#
# Usage: Rscript benchmarks/ws8_v2/assemble_prediction.R  # from repo root
#   env WS8_V2_PRED_DIR  dir holding prediction_metrics.csv
#                        (default .../results/run_m4b)
#   env WS8_V2_OUT       output dir (default .../results/run)
#
# Application case study — NOT calibration evidence.

base <- "benchmarks/ws8_v2/results"
pred_dir <- Sys.getenv("WS8_V2_PRED_DIR", file.path(base, "run_m4b"))
out_dir  <- Sys.getenv("WS8_V2_OUT", file.path(base, "run"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pm <- read.csv(file.path(pred_dir, "prediction_metrics.csv"),
               stringsAsFactors = FALSE)
panels <- c("topk", "native", "all59", "median")
folds <- sort(unique(pm$fold))

sign_flip_p <- function(d) {
  n <- length(d); obs <- mean(d)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- apply(signs, 1, function(s) mean(s * d))
  mean(abs(null) >= abs(obs))
}

# ---- Panel summary: macro + micro MAE, calibration ranges ---------------------
summ <- do.call(rbind, lapply(panels, function(pn) {
  d <- pm[pm$panel == pn, ]
  data.frame(
    panel = pn,
    macro_mae = round(mean(d$mae), 3),
    micro_mae = round(sum(d$mae * d$n) / sum(d$n), 3),
    macro_rmse = round(mean(d$rmse), 3),
    bias_min = round(min(d$bias), 2), bias_max = round(max(d$bias), 2),
    slope_po_min = round(min(d$slope_po), 3),
    slope_po_max = round(max(d$slope_po), 3),
    slope_op_min = round(min(d$slope_op), 3),
    slope_op_max = round(max(d$slope_op), 3),
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(out_dir, "prediction_panel_summary.csv"),
          row.names = FALSE)
cat("== M4 prediction panels: macro/micro MAE and calibration ranges ==\n")
print(summ, row.names = FALSE)

# ---- Paired per-cohort MAE contrasts ------------------------------------------
get_mae <- function(pn) sapply(folds, function(f) pm$mae[pm$panel == pn & pm$fold == f])
cmp <- do.call(rbind, lapply(list(c("topk", "all59"), c("topk", "median"),
                                  c("native", "all59")), function(pr) {
  d <- get_mae(pr[1]) - get_mae(pr[2])
  n <- length(d)
  tc <- stats::qt(0.975, n - 1) * stats::sd(d) / sqrt(n)
  data.frame(
    contrast = paste(pr[1], "-", pr[2]),
    macro_diff = round(mean(d), 4),
    ci_low = round(mean(d) - tc, 4), ci_high = round(mean(d) + tc, 4),
    spread_low = round(unname(quantile(d, 0.025)), 4),
    spread_high = round(unname(quantile(d, 0.975)), 4),
    sign_flip_p = round(sign_flip_p(d), 4),
    n_cohorts = n, stringsAsFactors = FALSE)
}))
write.csv(cmp, file.path(out_dir, "prediction_comparisons.csv"),
          row.names = FALSE)
cat("\n== Paired per-cohort MAE differences (negative favors first panel) ==\n")
print(cmp, row.names = FALSE)
cat("\nDone. Outputs in ", out_dir, "\n")
