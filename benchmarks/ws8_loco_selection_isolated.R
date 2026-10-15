# WS-8 PR5: selection-isolated LOCO for the aging case study.
#
# Protocol predeclared in ws8_showcase_predesign.md, Addendum 2 (2026-10-15),
# before this run.  For each held-out cohort: run stabl_environment_stability
# on the other seven cohorts (rf, knockoff_equi, 200 bootstraps,
# sample_fraction = 1, replace = TRUE, random_state = 20261013), apply the
# min-agg rule at >= 0.5 over the seven TRAINING cohorts, fit ridge
# (glmnet alpha = 0, nfolds = 5, lambda.min) on the within-fold panel, and
# predict the held-out cohort.  Fallback: a fold with an empty panel uses
# the all-feature model and is flagged.
#
# Usage: run from the repo root, e.g.
#   Rscript benchmarks/ws8_loco_selection_isolated.R
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)
suppressPackageStartupMessages(library(glmnet))

obj_path <- Sys.getenv("WS8_OBJECTS", "/workspace/ws8_out/ws8_objects.rds")
out_dir  <- Sys.getenv("WS8_OUT_DIR", "benchmarks")
obj <- readRDS(obj_path)
X <- obj$X; y_age <- obj$y_age; env <- obj$env
full_panel <- obj$sel_rf$feature[obj$sel_rf$selected_minagg]
cohorts <- sort(unique(env))

rows   <- list(); panels <- list(); fits <- list(); preds <- list()
for (h in cohorts) {
  t0 <- Sys.time()
  tr <- env != h; te <- env == h
  es <- stabl_environment_stability(
    x = X[tr, , drop = FALSE], y = y_age[tr], environments = env[tr],
    lambda_grid = "auto", base_learner = "rf",
    artificial_type = "knockoff_equi", n_bootstraps = 200L,
    sample_fraction = 1, replace = TRUE, random_state = 20261013
  )
  mins  <- get_min_scores(es)
  panel <- names(mins)[mins >= 0.5]
  fallback <- length(panel) == 0L
  if (fallback) panel <- colnames(X)
  cvfit <- cv.glmnet(X[tr, panel, drop = FALSE], y_age[tr], alpha = 0,
                     nfolds = 5)
  pr <- as.numeric(predict(cvfit, newx = X[te, panel, drop = FALSE],
                           s = "lambda.min"))
  rows[[h]] <- data.frame(
    cohort = h, n = sum(te), n_panel = length(panel), fallback = fallback,
    overlap_full_panel = length(intersect(panel, full_panel)),
    mae = mean(abs(pr - y_age[te])),
    slope = unname(coef(lm(pr ~ y_age[te]))[2]),
    stringsAsFactors = FALSE
  )
  panels[[h]]  <- data.frame(cohort = h, feature = panel,
                             stringsAsFactors = FALSE)
  fits[[h]]    <- es
  preds[[h]]   <- data.frame(cohort = h, obs = y_age[te], pred = pr,
                             stringsAsFactors = FALSE)
  cat(sprintf("%-9s n=%4d panel=%2d overlap=%2d mae=%.2f (%.1f min)\n",
              h, sum(te), length(panel),
              length(intersect(panel, full_panel)),
              mean(abs(pr - y_age[te])),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res    <- do.call(rbind, rows)
panels <- do.call(rbind, panels)
write.csv(res, file.path(out_dir, "ws8_loco_silo.csv"), row.names = FALSE)
write.csv(panels, file.path(out_dir, "ws8_loco_silo_features.csv"),
          row.names = FALSE)
saveRDS(list(fits = fits, preds = preds, full_panel = full_panel),
        file.path(out_dir, "ws8_loco_silo_objects.rds"))

cat("\nSelection-isolated LOCO (PR5):\n")
print(res, row.names = FALSE)
cat(sprintf("\npooled MAE: %.2f y (all-feature 12.14, post-selection panel 11.72)\n",
            mean(res$mae)))
