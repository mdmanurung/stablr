# WS-8 v2 M3/M4 — Driver: primary replication benchmark + secondary prediction.
#
# Builds the 59-feature composition design from the M1 donor tables, runs all
# six selectors (M0–M5) through the leakage-proof LOCO engine across the 8
# folds, computes the primary signed-rho replication endpoint (matched K=10,
# sensitivity {5,15,20}) and the secondary ridge prediction (top-K and native
# panels, plus median-age and all-59 baselines), and writes the required
# output tables.  No held-out label influences selection, signs, or tuning.
#
# Usage: Rscript benchmarks/ws8_v2/run.R   # run from the repository root
#   env WS8_V2_IN   dir holding donor_tables.rds (default .../results/m1)
#   env WS8_V2_OUT  output dir (default .../results/run)
#   env WS8_METHODS comma subset of M0..M5 (default all)
#
# Amendment 1: persists full per-fold training scores (training_scores.csv),
# per-cohort audit detail (training_env_detail.csv), and a per-fold manifest
# (fold_manifest.csv) per protocol sections 6/10; empty native panels fall
# back to the training-MEDIAN-age baseline per protocol section 5.
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)
suppressPackageStartupMessages(library(glmnet))
source("benchmarks/ws8_v2/metrics.R")
source("benchmarks/ws8_v2/loco.R")
source("benchmarks/ws8_v2/selectors.R")

in_dir  <- Sys.getenv("WS8_V2_IN", file.path("benchmarks/ws8_v2/results/m1"))
out_dir <- Sys.getenv("WS8_V2_OUT", file.path("benchmarks/ws8_v2/results/run"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
methods <- strsplit(Sys.getenv("WS8_METHODS", "M0,M1,M2,M3,M4,M5"), ",")[[1]]
K_primary <- 10L; K_sens <- c(5L, 15L, 20L)
# S1 simulations only: WS8_N_BOOT shrinks the STABL bootstrap count (default
# 200, the locked value).  Applies only to selectors that take n_boot.
n_boot <- as.integer(Sys.getenv("WS8_N_BOOT", "200"))

# ---- Design matrix (59 cell-type logit-proportions) ---------------------------
obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop))
colnames(X) <- paste0("CT:", colnames(X))
meta <- meta[match(rownames(X), meta$donor), ]
stopifnot(identical(meta$donor, rownames(X)))
y <- as.numeric(meta$age); cohort <- meta$dataset; donor <- meta$donor
covar <- data.frame(sex = meta$sex, ethnicity = meta$ethnicity,
                    row.names = meta$donor)
cat(sprintf("design: %d donors x %d features, %d cohorts\n",
            nrow(X), ncol(X), length(unique(cohort))))

# ---- Run all methods x folds --------------------------------------------------
fold_ids <- sort(unique(cohort))
rep_rows <- list(); pred_rows <- list(); feat_rows <- list(); met_rows <- list()
score_rows <- list(); detail_rows <- list(); fman_rows <- list()

for (m in methods) {
  sel_fn <- ws8v2_selectors[[m]]
  if (n_boot != 200L && "n_boot" %in% names(formals(sel_fn))) {
    base_fn <- sel_fn
    sel_fn <- function(X, y, cohort, covar = NULL)
      base_fn(X, y, cohort, covar = covar, n_boot = n_boot)
  }
  cat(sprintf("\n=== method %s ===\n", m))
  for (h in fold_ids) {
    t0 <- Sys.time()
    fr <- run_one_fold(h, X, y, cohort, donor, sel_fn,
                       K = max(K_primary, max(K_sens)), covar = covar)

    # --- Persist full training scores + audit detail (protocol 6/10) ---
    for (f in names(fr$scores)) {
      score_rows[[length(score_rows) + 1]] <- data.frame(
        method = m, fold = h, feature = f, score = unname(fr$scores[f]),
        native = f %in% fr$panel_native, stringsAsFactors = FALSE)
    }
    if (!is.null(fr$sel_full$detail)) {
      dm <- fr$sel_full$detail
      for (cn in colnames(dm)) {
        detail_rows[[length(detail_rows) + 1]] <- data.frame(
          method = m, fold = h, feature = rownames(dm), cohort = cn,
          value = dm[, cn], value_type = fr$sel_full$detail_type,
          stringsAsFactors = FALSE)
      }
    }
    hash_ids <- function(ids) {  # no extra deps: tools::md5sum on a tempfile
      tf <- tempfile(); on.exit(unlink(tf)); writeLines(sort(ids), tf)
      unname(tools::md5sum(tf))
    }
    fman_rows[[length(fman_rows) + 1]] <- data.frame(
      method = m, fold = h, n_train = length(fr$train_ids),
      n_test = length(fr$test_ids),
      train_ids_md5 = hash_ids(fr$train_ids),
      test_ids_md5  = hash_ids(fr$test_ids),
      n_dropped_features = length(fr$recipe$drop),
      drop_set = paste(fr$recipe$drop, collapse = ";"),
      n_native = length(fr$panel_native), n_topk = length(fr$panel_topk),
      seed = 20261013, stringsAsFactors = FALSE)

    # --- Primary: signed-rho replication at each K ---
    # held-out Spearman rho per feature (evaluation only)
    rho_te <- sapply(colnames(fr$Xte), function(f)
      suppressWarnings(stats::cor(fr$Xte[, f], fr$yte, method = "spearman")))
    names(rho_te) <- colnames(fr$Xte)
    R_at_K <- setNames(numeric(length(c(K_primary, K_sens))),
                       as.character(c(K_primary, K_sens)))
    for (K in c(K_primary, K_sens)) {
      panel <- head(fr$panel_topk, K)
      rr <- signed_rho_replication(panel, fr$signs, rho_te)
      R_at_K[as.character(K)] <- rr$R
      rep_rows[[length(rep_rows) + 1]] <- data.frame(
        method = m, fold = h, K = K, R = rr$R,
        sign_agreement = rr$sign_agreement, n_panel = length(panel),
        stringsAsFactors = FALSE)
      for (f in panel) {
        feat_rows[[length(feat_rows) + 1]] <- data.frame(
          method = m, fold = h, K = K, feature = f,
          train_sign = unname(fr$signs[f]), heldout_rho = unname(rho_te[f]),
          signed_rho = unname(rr$per_feature[f]),
          native = f %in% fr$panel_native, stringsAsFactors = FALSE)
      }
    }

    # --- Secondary: ridge prediction (top-K and native) + baselines ---
    fit_pred <- function(feats) {
      # < 2 features: cv.glmnet needs >= 2 columns; fall back to the
      # training-median baseline (same rule as the empty panel, sec. 5).
      # No committed result is affected: every real-atlas panel had >= 2.
      if (length(feats) < 2L) return(rep(median(fr$ytr), length(fr$yte)))
      cv <- train_ridge_cohort_cv(fr$Xtr[, feats, drop = FALSE], fr$ytr,
                                  fr$cohort_tr)
      as.numeric(predict(cv, newx = fr$Xte[, feats, drop = FALSE], s = "lambda.min"))
    }
    pred_topk   <- fit_pred(head(fr$panel_topk, K_primary))
    pred_native <- if (length(fr$panel_native) > 0L) fit_pred(fr$panel_native)
                   else rep(median(fr$ytr), length(fr$yte))  # zero-panel: median (sec. 5)
    pred_all    <- fit_pred(colnames(fr$Xtr))
    pred_median <- rep(median(fr$ytr), length(fr$yte))

    for (kind in c("topk", "native", "all59", "median")) {
      pr <- switch(kind, topk = pred_topk, native = pred_native,
                   all59 = pred_all, median = pred_median)
      ae <- abs(pr - fr$yte)
      met_rows[[length(met_rows) + 1]] <- data.frame(
        method = m, fold = h, panel = kind, n = length(fr$yte),
        mae = mean(ae), rmse = sqrt(mean((pr - fr$yte)^2)),
        bias = mean(pr - fr$yte),
        rank_cor = suppressWarnings(stats::cor(pr, fr$yte, method = "spearman")),
        slope_po = unname(coef(lm(pr ~ fr$yte))[2]),   # predicted on observed
        slope_op = unname(coef(lm(fr$yte ~ pr))[2]),   # observed on predicted
        stringsAsFactors = FALSE)
      pred_rows[[length(pred_rows) + 1]] <- data.frame(
        method = m, fold = h, panel = kind, donor = fr$test_ids,
        obs = fr$yte, pred = pr, stringsAsFactors = FALSE)
    }
    cat(sprintf("  %-9s R@10=%.3f  topK MAE=%.2f  native MAE=%.2f (%.1f min)\n",
                h, unname(R_at_K[as.character(K_primary)]),
                mean(abs(pred_topk - fr$yte)), mean(abs(pred_native - fr$yte)),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

rep_df  <- do.call(rbind, rep_rows)
pred_df <- do.call(rbind, pred_rows)
feat_df <- do.call(rbind, feat_rows)
met_df  <- do.call(rbind, met_rows)
score_df  <- do.call(rbind, score_rows)
fman_df   <- do.call(rbind, fman_rows)
detail_df <- if (length(detail_rows) > 0L) do.call(rbind, detail_rows) else NULL

write.csv(rep_df,  file.path(out_dir, "heldout_replication.csv"), row.names = FALSE)
write.csv(pred_df, file.path(out_dir, "predictions.csv"), row.names = FALSE)
write.csv(feat_df, file.path(out_dir, "feature_scores.csv"), row.names = FALSE)
write.csv(met_df,  file.path(out_dir, "prediction_metrics.csv"), row.names = FALSE)
write.csv(score_df, file.path(out_dir, "training_scores.csv"), row.names = FALSE)
write.csv(fman_df,  file.path(out_dir, "fold_manifest.csv"), row.names = FALSE)
if (!is.null(detail_df)) {
  write.csv(detail_df, file.path(out_dir, "training_env_detail.csv"),
            row.names = FALSE)
}

# ---- Primary summary: macro R(K) per method -----------------------------------
cat("\n== PRIMARY: macro signed-rho replication R(K) ==\n")
for (K in c(K_primary, K_sens)) {
  sub <- rep_df[rep_df$K == K, ]
  cat(sprintf("K=%2d: ", K))
  for (m in methods) {
    cat(sprintf("%s=%.3f  ", m, mean(sub$R[sub$method == m])))
  }
  cat("\n")
}

# ---- Method comparison: paired cohort-level differences at K=10 ---------------
k10 <- rep_df[rep_df$K == K_primary, ]
cmp <- do.call(rbind, lapply(setdiff(methods, "M4"), function(m) {
  d <- k10$R[k10$method == "M4"] - k10$R[k10$method == m][match(
    k10$fold[k10$method == "M4"], k10$fold[k10$method == m])]
  data.frame(contrast = paste0("M4-", m), macro_diff = mean(d),
             ci_low = unname(quantile(d, 0.025)), ci_high = unname(quantile(d, 0.975)),
             n_cohorts = length(d), stringsAsFactors = FALSE)
}))
write.csv(cmp, file.path(out_dir, "method_comparisons.csv"), row.names = FALSE)
cat("\n== M4 vs comparators (macro paired difference in R@10) ==\n")
print(cmp, row.names = FALSE)

# ---- Run manifest --------------------------------------------------------------
manifest <- list(
  commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA),
  r_version = R.version.string,
  methods = methods, K_primary = K_primary, K_sensitivity = K_sens,
  n_donors = nrow(X), n_features = ncol(X), seed = 20261013,
  run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "run_manifest.json"))

cat("\nDone. Outputs in ", out_dir, "\n")
