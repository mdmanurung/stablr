# WS-8 v3 addendum — M10: Invariant Causal Prediction (ICP) arm.
#
# Post-hoc addendum to the v3 candidate benchmark: direct evidence on ICP
# (Peters, Buhlmann & Meinshausen 2016) itself, which inspired the package's
# environment-stability workflow (M4/M5) but was not among the benchmarked
# arms.  Same locked LOCO harness, same primary endpoint (held-out signed-rho
# replication R@10), same seed policy.
#
# Selector contract (locked WS-8 v2):
#   fn(X_train, y_train, cohort_train, covar = NULL)
#     -> list(ranked, native, scores, detail = NULL, detail_type, ...)
#
# M10 design:
#   * InvariantCausalPrediction::ICP() with the 7 training cohorts as
#     environments (ExpInd), alpha = 0.05, package-default variable selection
#     ("boosting" for p > 8), stopIfEmpty = FALSE.
#   * Ranking: ICP p-value descending (large p = no evidence against the
#     variable's invariant role), tie-broken by |aggregate coefficient| desc,
#     then feature name (full determinism).
#   * native = accepted set (pvalue >= alpha).
#   * scores = named vector of ICP p-values (persisted by the harness).
#
# Application case study — NOT calibration evidence.

m10_icp_selector <- function(X, y, cohort, covar = NULL, alpha = 0.05,
                             seed = 20261013) {
  stopifnot(requireNamespace("InvariantCausalPrediction", quietly = TRUE))
  # Drop zero-variance columns (ICP/OLS cannot use them).
  v <- apply(X, 2, stats::var, na.rm = TRUE)
  keep <- names(v)[is.finite(v) & v > 1e-12]
  Xd <- X[, keep, drop = FALSE]  # matrix; ICP builds model.frame(Y ~ X - 1)
  storage.mode(Xd) <- "numeric"
  env <- setNames(as.integer(factor(cohort)), names(cohort))

  set.seed(seed)
  t0 <- Sys.time()
  fit <- InvariantCausalPrediction::ICP(
    X = Xd, Y = as.numeric(y), ExpInd = env, alpha = alpha,
    stopIfEmpty = FALSE, showAcceptedSets = FALSE, showCompletion = FALSE
  )
  runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  pv <- fit$pvalues
  nm <- if (!is.null(fit$colnames) && length(fit$colnames) == length(pv))
    as.character(fit$colnames) else keep
  stopifnot(length(nm) == length(pv))
  names(pv) <- nm
  pv[is.na(pv)] <- 0
  cf <- fit$maximinCoefficients  # p-vector aggregate coefficient estimate
  if (is.null(cf) || length(cf) != length(pv)) cf <- rep(0, length(pv))
  cf <- as.numeric(cf); cf[is.na(cf)] <- 0
  ord <- order(-pv, -abs(cf), nm)
  ranked <- nm[ord]
  native <- nm[pv >= alpha]

  list(
    ranked = ranked,
    native = native,
    scores = pv,
    detail = NULL, detail_type = NA_character_,
    extra = c(
      n_accepted = length(native),
      n_tested = length(pv),
      alpha = alpha,
      runtime_min = runtime_min
    )
  )
}

# ---- Standalone LOCO run for M10 (reuses the locked v2 engine) --------------
if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(matrixStats))
  engine_dir <- Sys.getenv("WS8_ENGINE", "/workspace/ws8v3/engine")
  source(file.path(engine_dir, "metrics.R"))
  source(file.path(engine_dir, "loco.R"))
  source(file.path(engine_dir, "selectors.R"))

  out_dir <- Sys.getenv("M10_OUT", "/workspace/ws8v3/results/run_v3_m10")
  ckpt_dir <- Sys.getenv("M10_CKPT", "/mnt/shared-workspace/ws8v3_ckpt_m10")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

  obj <- readRDS("/workspace/ws8v3/donor_tables.rds")
  ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
  prop <- (ct_mat + 0.5) / (total + 0.5)
  X <- log(prop / (1 - prop))
  colnames(X) <- paste0("CT:", colnames(X))
  meta <- meta[match(rownames(X), meta$donor), ]
  y <- as.numeric(meta$age); cohort <- meta$dataset; donor <- meta$donor

  K_primary <- 10L; K_sens <- c(5L, 15L, 20L)
  fold_ids <- sort(unique(cohort))
  rep_rows <- list(); feat_rows <- list(); fman_rows <- list(); extra_rows <- list()

  for (h in fold_ids) {
    ck <- file.path(ckpt_dir, sprintf("fold_%s_M10.rds", h))
    if (file.exists(ck)) {
      fr <- readRDS(ck); cat(sprintf("  %-9s [checkpoint]\n", h))
    } else {
      t0 <- Sys.time()
      fr <- run_one_fold(h, X, y, cohort, donor, m10_icp_selector,
                         K = max(K_primary, max(K_sens)))
      fr$runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      saveRDS(fr, ck)
    }
    rho_te <- sapply(colnames(fr$Xte), function(f)
      suppressWarnings(stats::cor(fr$Xte[, f], fr$yte, method = "spearman")))
    names(rho_te) <- colnames(fr$Xte)
    for (K in c(K_primary, K_sens)) {
      panel <- head(fr$panel_topk, K)
      rr <- signed_rho_replication(panel, fr$signs, rho_te)
      rep_rows[[length(rep_rows) + 1L]] <- data.frame(
        method = "M10", fold = h, K = K, R = rr$R,
        sign_agreement = rr$sign_agreement, n_panel = length(panel),
        stringsAsFactors = FALSE)
      for (f in panel) {
        feat_rows[[length(feat_rows) + 1L]] <- data.frame(
          method = "M10", fold = h, K = K, feature = f,
          train_sign = unname(fr$signs[f]), heldout_rho = unname(rho_te[f]),
          signed_rho = unname(rr$per_feature[f]),
          native = f %in% fr$panel_native, stringsAsFactors = FALSE)
      }
    }
    fman_rows[[length(fman_rows) + 1L]] <- data.frame(
      method = "M10", fold = h, n_train = length(fr$train_ids),
      n_test = length(fr$test_ids),
      n_dropped_features = length(fr$recipe$drop),
      n_native = length(fr$panel_native), n_topk = length(fr$panel_topk),
      runtime_min = fr$runtime_min, seed = 20261013,
      stringsAsFactors = FALSE)
    for (en in names(fr$sel_full$extra)) {
      extra_rows[[length(extra_rows) + 1L]] <- data.frame(
        method = "M10", fold = h, field = en, feature = "",
        value = unname(fr$sel_full$extra[[en]]), stringsAsFactors = FALSE)
    }
    cat(sprintf("  %-9s R@10 = %.3f  native = %d  (%.1f min)\n", h,
                rep_rows[[length(rep_rows)]]$R, length(fr$panel_native),
                fr$runtime_min))
  }

  rep_df <- do.call(rbind, rep_rows)
  feat_df <- do.call(rbind, feat_rows)
  fman_df <- do.call(rbind, fman_rows)
  extra_df <- do.call(rbind, extra_rows)
  write.csv(rep_df, file.path(out_dir, "heldout_replication_m10.csv"), row.names = FALSE)
  write.csv(feat_df, file.path(out_dir, "per_feature_m10.csv"), row.names = FALSE)
  write.csv(fman_df, file.path(out_dir, "fold_manifest_m10.csv"), row.names = FALSE)
  write.csv(extra_df, file.path(out_dir, "extras_m10.csv"), row.names = FALSE)

  macro <- function(K) mean(rep_df$R[rep_df$K == K])
  cat(sprintf("\n== M10 macro R@10 = %.3f (K5 %.3f, K15 %.3f, K20 %.3f) ==\n",
              macro(10), macro(5), macro(15), macro(20)))
  cat(sprintf("mean accepted-set size: %.1f\n", mean(fman_df$n_native)))
}
