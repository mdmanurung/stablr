# WS-8 v2 S2b — aggregation-strength sensitivity (followup_s1_s2.md).
#
# Re-aggregates the PERSISTED per-cohort M4 selection frequencies
# (results/run_m4b/training_env_detail.csv) at cross-cohort quantiles
# q in {0, 0.10, 0.25} instead of refitting: the STABL fits are
# seed-deterministic (V10), so quantile re-aggregation of the committed
# env scores is bit-identical to refitting M4q10/M4q25 via run.R.
# Only the aggregation rule changes; the locked tie-break is kept.
#
# Validation: q = 0 (min) must reproduce the committed M4 per-fold R@10
# exactly; the script stops with an error otherwise.
#
# Usage: Rscript benchmarks/ws8_v2/s2b_aggregation_quantile.R   # from repo root
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("benchmarks/ws8_v2/metrics.R")
source("benchmarks/ws8_v2/loco.R")
source("benchmarks/ws8_v2/selectors.R")

in_dir  <- Sys.getenv("WS8_V2_IN", file.path("benchmarks/ws8_v2/results/m1"))
det_f   <- Sys.getenv("WS8_S2B_DETAIL",
                      file.path("benchmarks/ws8_v2/results/run_m4b",
                                "training_env_detail.csv"))
out_dir <- Sys.getenv("WS8_V2_OUT",
                      file.path("benchmarks/ws8_v2/results/s2b_aggregation"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
Ks <- c(10L, 5L, 15L, 20L)
qs <- c(`M4q00` = 0, `M4q10` = 0.10, `M4q25` = 0.25)

# ---- Design matrix (identical construction to run.R) --------------------------
obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
prop <- (obj$ct_mat + 0.5) / (obj$total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- obj$meta[match(rownames(X), obj$meta$donor), ]
stopifnot(identical(meta$donor, rownames(X)))
y <- as.numeric(meta$age); cohort <- meta$dataset; donor <- meta$donor

# ---- Persisted M4 per-cohort selection frequencies, pivoted per fold ---------
det <- read.csv(det_f, stringsAsFactors = FALSE)
stopifnot(all(det$method == "M4"), all(det$value_type == "sel_freq"))
fold_ids <- sort(unique(cohort))

esc_by_fold <- lapply(fold_ids, function(h) {
  d <- det[det$fold == h, ]
  esc <- tapply(d$value, list(d$feature, d$cohort), identity)
  esc <- as.matrix(esc)
  esc[order(rownames(esc)), , drop = FALSE]
})
names(esc_by_fold) <- fold_ids

# ---- Evaluate each quantile rule on each fold --------------------------------
rows <- list()
for (h in fold_ids) {
  te <- cohort == h; tr <- !te
  recipe <- fit_fold_transform(X[tr, , drop = FALSE], cohort[tr])
  Xtr <- apply_fold_transform(recipe, X[tr, , drop = FALSE])
  Xte <- apply_fold_transform(recipe, X[te, , drop = FALSE])
  ytr <- setNames(y[tr], donor[tr]); yte <- setNames(y[te], donor[te])
  ctr <- setNames(as.character(cohort[tr]), donor[tr])
  signs <- training_signs(Xtr, ytr, ctr)
  rho_te <- sapply(colnames(Xte), function(f)
    suppressWarnings(stats::cor(Xte[, f], yte, method = "spearman")))
  names(rho_te) <- colnames(Xte)

  esc <- esc_by_fold[[h]]
  # env scores were fit on the transformed training design; align by name
  esc <- esc[rownames(esc) %in% colnames(Xtr), , drop = FALSE]
  for (lbl in names(qs)) {
    agg <- aggregate_env_scores(esc, qs[[lbl]])
    for (K in Ks) {
      panel <- head(agg$ranked, K)
      rr <- signed_rho_replication(panel, signs, rho_te)
      rows[[length(rows) + 1]] <- data.frame(
        method = lbl, fold = h, K = K, R = rr$R,
        sign_agreement = rr$sign_agreement, n_panel = length(panel),
        n_native = length(agg$native), stringsAsFactors = FALSE)
    }
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "s2b_quantile_replication.csv"),
          row.names = FALSE)

# ---- Validation: q = 0 must reproduce committed M4 R@10 bit-exactly -----------
ref <- read.csv(file.path("benchmarks/ws8_v2/results/run_m4b",
                          "heldout_replication.csv"))
ref <- ref[ref$method == "M4", ]
chk <- merge(res[res$method == "M4q00", c("fold", "K", "R")],
             ref[, c("fold", "K", "R")], by = c("fold", "K"),
             suffixes = c("_s2b", "_committed"))
max_diff <- max(abs(chk$R_s2b - chk$R_committed))
cat(sprintf("validation q=0 vs committed M4: max|diff| = %g over %d rows\n",
            max_diff, nrow(chk)))
if (max_diff > 1e-12) stop("S2b chain does not reproduce committed M4; aborting.")

# ---- Summary: macro R@10 per rule + paired contrasts vs min -------------------
k10 <- res[res$K == 10, ]
macro <- tapply(k10$R, k10$method, mean)
cat("\n== S2b macro R@10 by aggregation rule ==\n")
print(round(macro, 4))
nat <- tapply(k10$n_native, k10$method, mean)
cat("\n== mean native panel size (score >= 0.5) ==\n")
print(round(nat, 1))

sign_flip_p <- function(d) {
  n <- length(d); obs <- mean(d)
  signs_m <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  mean(abs(apply(signs_m, 1, function(s) mean(s * d))) >= abs(obs))
}
contr <- do.call(rbind, lapply(c("M4q10", "M4q25"), function(m) {
  d <- k10$R[k10$method == m][order(k10$fold[k10$method == m])] -
       k10$R[k10$method == "M4q00"][order(k10$fold[k10$method == "M4q00"])]
  n <- length(d); tc <- stats::qt(0.975, n - 1) * stats::sd(d) / sqrt(n)
  data.frame(contrast = paste0(m, " - M4(min)"), macro_diff = round(mean(d), 4),
             ci_low = round(mean(d) - tc, 4), ci_high = round(mean(d) + tc, 4),
             sign_flip_p = round(sign_flip_p(d), 4), n_cohorts = n)
}))
write.csv(contr, file.path(out_dir, "s2b_quantile_contrasts.csv"),
          row.names = FALSE)
cat("\n== paired contrasts vs min-aggregation (R@10) ==\n")
print(contr, row.names = FALSE)
cat("\nDone. Outputs in", out_dir, "\n")
