# WS-8 v2 S1 — assemble synthetic-calibration results (followup_s1_s2.md).
#
# Collects the per-replicate run.R outputs under results/s1/ and evaluates
# each selector against the KNOWN truth (s1_truth.csv):
#   - macro R@10 (endpoint behaviour under known truth)
#   - top-10 recovery: |panel ∩ true| / min(10, n_true) per fold
#   - sign accuracy: train_sign == sign(true beta) on recovered features
#   - native panel size (selection-count saturation readout, s1c)
#   - empirical FDP: |native ∩ untrue| / max(|native|, 1) per fold
# Empirical FDP against the simulated truth is reported instead of the
# fitted FDP+ estimate (run.R does not persist fit-level FDP+); in
# simulations the truth-based quantity is the stronger readout.
#
# Usage: Rscript benchmarks/ws8_v2/s1_assemble.R   # from repo root
# Application case study — simulation calibration, NOT application evidence.

in_root  <- Sys.getenv("WS8_S1_RES", file.path("benchmarks/ws8_v2/results/s1"))
out_dir  <- in_root
truth    <- read.csv("benchmarks/ws8_v2/results/s1_truth.csv",
                     stringsAsFactors = FALSE)

dirs <- list.dirs(in_root, recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("_rep[0-9]+$", dirs)]
stopifnot(length(dirs) > 0)

rows <- list()
for (d in dirs) {
  tag <- basename(d)                       # e.g. s1a_hom_rep1
  scenario <- sub("_rep[0-9]+$", "", tag)
  rep <- as.integer(sub(".*_rep", "", tag))
  f_rep <- file.path(d, "heldout_replication.csv")
  f_scr <- file.path(d, "training_scores.csv")
  f_fea <- file.path(d, "feature_scores.csv")
  if (!file.exists(f_rep) || !file.exists(f_scr)) next
  rr <- read.csv(f_rep, stringsAsFactors = FALSE)
  sc <- read.csv(f_scr, stringsAsFactors = FALSE)
  fe <- if (file.exists(f_fea)) read.csv(f_fea, stringsAsFactors = FALSE) else NULL

  # Truth is per scenario x rep (base_beta redrawn each rep); filter both.
  # s1c_null has zero true rows: guard the paste0 prefix, which would
  # otherwise recycle character(0) to the length-1 string "CT:".
  t_idx <- truth$scenario == scenario & truth$rep == rep
  true_f <- truth$feature[t_idx]
  true_f <- if (length(true_f)) paste0("CT:", true_f) else character(0)
  true_sign <- setNames(sign(truth$base_beta[t_idx]), true_f)

  for (m in unique(rr$method)) {
    scr_m <- sc[sc$method == m, ]
    for (fd in sort(unique(rr$fold))) {
      r10 <- rr$R[rr$method == m & rr$K == 10 & rr$fold == fd]
      if (!length(r10)) next
      nat <- scr_m$feature[scr_m$fold == fd & scr_m$native]
      n_native <- length(nat)
      emp_fdp <- if (n_native) sum(!(nat %in% true_f)) / n_native else NA_real_
      rec <- NA_real_; sacc <- NA_real_
      if (!is.null(fe) && length(true_f)) {
        panel <- fe$feature[fe$method == m & fe$fold == fd & fe$K == 10]
        rec <- sum(panel %in% true_f) / min(10, length(true_f))
        got <- intersect(panel, true_f)
        if (length(got)) {
          ts <- fe$train_sign[fe$method == m & fe$fold == fd &
                                fe$K == 10][match(got, panel)]
          sacc <- mean(ts == true_sign[got])
        }
      }
      rows[[length(rows) + 1]] <- data.frame(
        scenario = scenario, rep = rep, method = m, fold = fd,
        R10 = r10, n_native = n_native, emp_fdp_native = emp_fdp,
        top10_recovery = rec, sign_accuracy = sacc, stringsAsFactors = FALSE)
    }
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "s1_by_fold.csv"), row.names = FALSE)

# ---- Summary: per scenario x method, mean over folds then reps ---------------
agg <- aggregate(cbind(R10, n_native, emp_fdp_native, top10_recovery,
                       sign_accuracy) ~ scenario + method, res, mean,
                 na.action = na.pass)
se  <- aggregate(R10 ~ scenario + method, res,
                 function(x) stats::sd(x) / sqrt(length(x)))
names(se)[3] <- "R10_se"
agg <- merge(agg, se, by = c("scenario", "method"))
agg <- agg[order(agg$scenario, -agg$R10), ]
write.csv(agg, file.path(out_dir, "s1_summary.csv"), row.names = FALSE)
agg[] <- lapply(agg, function(x) if (is.numeric(x)) round(x, 3) else x)
print(agg)
cat("\nDone. Outputs in", out_dir, "\n")
