# WS-8 v3 — Driver: candidate arms (M6–M9) through the locked v2 LOCO engine.
#
# Adapted from benchmarks/ws8_v2/run.R.  Differences:
#   * methods come from selectors_v3.R (M6, M7, M7c, M8, M9, M9s); M4 can be
#     included for harness-equivalence verification.
#   * per-(method, fold) checkpoint RDS in CKPT_DIR; reruns skip completed cells.
#   * a selector may return `test_transform` (and `train_transformed`): applied
#     to Xte/Xtr ONLY for the ridge prediction step (M9).  The replication
#     endpoint always uses the unadjusted held-out matrix (within-cohort
#     Spearman is invariant to the per-feature monotone ComBat adjustment).
#   * selector `extra` named numerics are persisted long-format in extras.csv.
#
# Usage: Rscript run_v3.R
#   env WS8_METHODS  comma subset (default "M6,M7c,M8,M9,M9s"; M7 separate)
#   env WS8_V3_OUT   output dir (default benchmarks/ws8_v3/results/run_v3)
#   env CKPT_DIR     checkpoint dir (default benchmarks/ws8_v3/results/checkpoints)
#
# Application case study — NOT calibration evidence.

# Path config (repo-relative defaults; override via env vars).
v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
v2_dir <- Sys.getenv("WS8_V2_DIR", "benchmarks/ws8_v2")
in_dir <- Sys.getenv("WS8_V2_IN", file.path(v2_dir, "results", "m1"))
suppressPackageStartupMessages({
  library(stablr); library(glmnet); library(matrixStats)
})
engine_dir <- Sys.getenv("WS8_ENGINE", v2_dir)
source(file.path(engine_dir, "metrics.R"))
source(file.path(engine_dir, "loco.R"))
source(file.path(engine_dir, "selectors.R"))   # v2 selectors incl. M4 + helpers
source(file.path(v3_dir, "selectors_v3.R"))

out_dir <- Sys.getenv("WS8_V3_OUT", file.path(v3_dir, "results", "run_v3"))
ckpt_dir <- Sys.getenv("CKPT_DIR", file.path(v3_dir, "results", "checkpoints"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
methods <- strsplit(Sys.getenv("WS8_METHODS", "M6,M7c,M8,M9,M9s"), ",")[[1]]
K_primary <- 10L; K_sens <- c(5L, 15L, 20L)
n_boot <- as.integer(Sys.getenv("WS8_N_BOOT", "200"))

# ---- Design matrix (identical construction to v2 run.R) -----------------------
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

all_selectors <- c(ws8v2_selectors, ws8v3_selectors)
fold_ids <- sort(unique(cohort))

rep_rows <- list(); pred_rows <- list(); feat_rows <- list(); met_rows <- list()
score_rows <- list(); detail_rows <- list(); fman_rows <- list(); extra_rows <- list()

for (m in methods) {
  sel_fn <- all_selectors[[m]]
  stopifnot(!is.null(sel_fn))
  cat(sprintf("\n=== method %s ===\n", m))
  for (h in fold_ids) {
    ck <- file.path(ckpt_dir, sprintf("fold_%s_%s.rds", m, h))
    if (file.exists(ck)) {
      fr <- readRDS(ck)
      cat(sprintf("  %-9s [checkpoint]\n", h))
    } else {
      t0 <- Sys.time()
      fr <- run_one_fold(h, X, y, cohort, donor, sel_fn,
                         K = max(K_primary, max(K_sens)), covar = covar)
      fr$runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      saveRDS(fr, ck)
    }

    # --- Persist training scores / detail / extras (protocol 6/10 schema) ---
    for (f in names(fr$scores)) {
      score_rows[[length(score_rows) + 1L]] <- data.frame(
        method = m, fold = h, feature = f, score = unname(fr$scores[f]),
        native = f %in% fr$panel_native, stringsAsFactors = FALSE)
    }
    if (!is.null(fr$sel_full$detail)) {
      dm <- fr$sel_full$detail
      for (cn2 in colnames(dm)) {
        detail_rows[[length(detail_rows) + 1L]] <- data.frame(
          method = m, fold = h, feature = rownames(dm), cohort = cn2,
          value = dm[, cn2], value_type = fr$sel_full$detail_type,
          stringsAsFactors = FALSE)
      }
    }
    if (!is.null(fr$sel_full$extra)) {
      for (en in names(fr$sel_full$extra)) {
        ev <- fr$sel_full$extra[[en]]
        if (is.numeric(ev) && !is.null(names(ev))) {
          for (fn2 in names(ev)) {
            extra_rows[[length(extra_rows) + 1L]] <- data.frame(
              method = m, fold = h, field = en, feature = fn2,
              value = unname(ev[fn2]), stringsAsFactors = FALSE)
          }
        } else if (is.numeric(ev) && length(ev) == 1L) {
          extra_rows[[length(extra_rows) + 1L]] <- data.frame(
            method = m, fold = h, field = en, feature = "",
            value = unname(ev), stringsAsFactors = FALSE)
        } else if (is.numeric(ev)) {  # unnamed vector (e.g. fdp_curve)
          for (ii in seq_along(ev)) {
            extra_rows[[length(extra_rows) + 1L]] <- data.frame(
              method = m, fold = h, field = en, feature = as.character(ii - 1L),
              value = unname(ev[ii]), stringsAsFactors = FALSE)
          }
        } else if (is.logical(ev) && !is.null(names(ev))) {
          for (fn2 in names(ev)) {
            extra_rows[[length(extra_rows) + 1L]] <- data.frame(
              method = m, fold = h, field = en, feature = fn2,
              value = as.numeric(ev[fn2]), stringsAsFactors = FALSE)
          }
        }
      }
    }
    hash_ids <- function(ids) {
      tf <- tempfile(); on.exit(unlink(tf)); writeLines(sort(ids), tf)
      unname(tools::md5sum(tf))
    }
    fman_rows[[length(fman_rows) + 1L]] <- data.frame(
      method = m, fold = h, n_train = length(fr$train_ids),
      n_test = length(fr$test_ids),
      train_ids_md5 = hash_ids(fr$train_ids),
      test_ids_md5 = hash_ids(fr$test_ids),
      n_dropped_features = length(fr$recipe$drop),
      drop_set = paste(fr$recipe$drop, collapse = ";"),
      n_native = length(fr$panel_native), n_topk = length(fr$panel_topk),
      runtime_min = ifelse(is.null(fr$runtime_min), NA_real_, fr$runtime_min),
      seed = 20261013, stringsAsFactors = FALSE)

    # --- Primary: signed-rho replication (unadjusted held-out data) ---
    rho_te <- sapply(colnames(fr$Xte), function(f)
      suppressWarnings(stats::cor(fr$Xte[, f], fr$yte, method = "spearman")))
    names(rho_te) <- colnames(fr$Xte)
    for (K in c(K_primary, K_sens)) {
      panel <- head(fr$panel_topk, K)
      rr <- signed_rho_replication(panel, fr$signs, rho_te)
      rep_rows[[length(rep_rows) + 1L]] <- data.frame(
        method = m, fold = h, K = K, R = rr$R,
        sign_agreement = rr$sign_agreement, n_panel = length(panel),
        stringsAsFactors = FALSE)
      for (f in panel) {
        feat_rows[[length(feat_rows) + 1L]] <- data.frame(
          method = m, fold = h, K = K, feature = f,
          train_sign = unname(fr$signs[f]), heldout_rho = unname(rho_te[f]),
          signed_rho = unname(rr$per_feature[f]),
          native = f %in% fr$panel_native, stringsAsFactors = FALSE)
      }
    }

    # --- Secondary: ridge prediction (test_transform hook for M9) ---
    Xtr_p <- fr$Xtr; Xte_p <- fr$Xte
    if (!is.null(fr$sel_full$test_transform)) {
      Xtr_p <- fr$sel_full$train_transformed
      Xte_p <- fr$sel_full$test_transform(fr$Xte)
      stopifnot(identical(colnames(Xtr_p), colnames(Xte_p)))
    }
    fit_pred <- function(feats) {
      if (length(feats) < 2L) return(rep(median(fr$ytr), length(fr$yte)))
      cv <- train_ridge_cohort_cv(Xtr_p[, feats, drop = FALSE], fr$ytr,
                                  fr$cohort_tr)
      as.numeric(predict(cv, newx = Xte_p[, feats, drop = FALSE],
                         s = "lambda.min"))
    }
    pred_topk   <- fit_pred(head(fr$panel_topk, K_primary))
    pred_native <- if (length(fr$panel_native) > 0L) fit_pred(fr$panel_native)
                   else rep(median(fr$ytr), length(fr$yte))
    pred_all    <- fit_pred(colnames(Xtr_p))
    pred_median <- rep(median(fr$ytr), length(fr$yte))

    for (kind in c("topk", "native", "all59", "median")) {
      pr <- switch(kind, topk = pred_topk, native = pred_native,
                   all59 = pred_all, median = pred_median)
      ae <- abs(pr - fr$yte)
      met_rows[[length(met_rows) + 1L]] <- data.frame(
        method = m, fold = h, panel = kind, n = length(fr$yte),
        mae = mean(ae), rmse = sqrt(mean((pr - fr$yte)^2)),
        bias = mean(pr - fr$yte),
        rank_cor = suppressWarnings(stats::cor(pr, fr$yte, method = "spearman")),
        slope_po = unname(coef(lm(pr ~ fr$yte))[2]),
        slope_op = unname(coef(lm(fr$yte ~ pr))[2]),
        stringsAsFactors = FALSE)
      pred_rows[[length(pred_rows) + 1L]] <- data.frame(
        method = m, fold = h, panel = kind, donor = fr$test_ids,
        obs = fr$yte, pred = pr, stringsAsFactors = FALSE)
    }
    cat(sprintf("  %-9s R@10=%.3f  topK MAE=%.2f  native MAE=%.2f  nativeK=%d\n",
                h, rep_rows[[length(rep_rows) - 3L]]$R,
                mean(abs(pred_topk - fr$yte)), mean(abs(pred_native - fr$yte)),
                length(fr$panel_native)))
  }
}

rep_df <- do.call(rbind, rep_rows)
pred_df <- do.call(rbind, pred_rows)
feat_df <- do.call(rbind, feat_rows)
met_df <- do.call(rbind, met_rows)
score_df <- do.call(rbind, score_rows)
fman_df <- do.call(rbind, fman_rows)
detail_df <- if (length(detail_rows) > 0L) do.call(rbind, detail_rows) else NULL
extra_df <- if (length(extra_rows) > 0L) do.call(rbind, extra_rows) else NULL

write.csv(rep_df, file.path(out_dir, "heldout_replication.csv"), row.names = FALSE)
write.csv(pred_df, file.path(out_dir, "predictions.csv"), row.names = FALSE)
write.csv(feat_df, file.path(out_dir, "feature_scores.csv"), row.names = FALSE)
write.csv(met_df, file.path(out_dir, "prediction_metrics.csv"), row.names = FALSE)
write.csv(score_df, file.path(out_dir, "training_scores.csv"), row.names = FALSE)
write.csv(fman_df, file.path(out_dir, "fold_manifest.csv"), row.names = FALSE)
if (!is.null(detail_df))
  write.csv(detail_df, file.path(out_dir, "training_env_detail.csv"), row.names = FALSE)
if (!is.null(extra_df))
  write.csv(extra_df, file.path(out_dir, "extras.csv"), row.names = FALSE)

cat("\n== PRIMARY: macro signed-rho replication R(K) ==\n")
for (K in c(K_primary, K_sens)) {
  sub <- rep_df[rep_df$K == K, ]
  cat(sprintf("K=%2d: ", K))
  for (m in methods) cat(sprintf("%s=%.3f  ", m, mean(sub$R[sub$method == m])))
  cat("\n")
}

manifest <- list(
  code = "ws8_v3 run_v3.R", r_version = R.version.string,
  stablr_version = as.character(packageVersion("stablr")),
  methods = methods, K_primary = K_primary, K_sensitivity = K_sens,
  n_donors = nrow(X), n_features = ncol(X), seed = 20261013,
  run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "run_manifest.json"))
cat("\nDone. Outputs in ", out_dir, "\n")
