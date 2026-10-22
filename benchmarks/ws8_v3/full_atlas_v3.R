# WS-8 v3 — Full-atlas (all 8 cohorts) fits for the paper-comparison co-primary.
#
# Every method is fit once on the full 2,609-donor atlas (same fold-safe
# transform code path as the LOCO runs).  M5 is derived from the SAME
# stabl_environment_stability fit as M4 (identical args + seed => identical
# env_scores; deriving in-session avoids any cross-session determinism
# assumption).  M9s shares M9's selection (only the held-out transform mode
# differs, which is irrelevant for a full-atlas panel).  M7/M8 full-atlas fits
# live in null_probe_v3.R (perm 0 = real) and are merged at assembly.
#
# Usage: Rscript full_atlas_v3.R
#   env WS8_FA_METHODS  comma subset (default "M0,M1,M2,M3,M4,M5,M6,M7c,M9,M9s")
#   env WS8_FA_OUT      output dir (default benchmarks/ws8_v3/results/full_atlas)
#
# Application case study — NOT calibration evidence.

# Path config (repo-relative defaults; override via env vars).
v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
v2_dir <- Sys.getenv("WS8_V2_DIR", "benchmarks/ws8_v2")
in_dir <- Sys.getenv("WS8_V2_IN", file.path(v2_dir, "results", "m1"))
suppressPackageStartupMessages({
  library(stablr); library(glmnet); library(matrixStats); library(mixOmics)
})
engine_dir <- Sys.getenv("WS8_ENGINE", v2_dir)
source(file.path(engine_dir, "metrics.R"))
source(file.path(engine_dir, "loco.R"))
source(file.path(engine_dir, "selectors.R"))
source(file.path(v3_dir, "selectors_v3.R"))

out_dir <- Sys.getenv("WS8_FA_OUT", file.path(v3_dir, "results", "full_atlas"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
methods <- strsplit(Sys.getenv("WS8_FA_METHODS",
                               "M0,M1,M2,M3,M4,M5,M6,M7c,M9,M9s"), ",")[[1]]

obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- meta[match(rownames(X), meta$donor), ]
y <- setNames(as.numeric(meta$age), meta$donor)
cohort <- setNames(as.character(meta$dataset), meta$donor)
covar <- data.frame(sex = meta$sex, ethnicity = meta$ethnicity,
                    row.names = meta$donor)
recipe <- fit_fold_transform(X, cohort)
Xf <- apply_fold_transform(recipe, X)
cat(sprintf("full-atlas fits: %d x %d, %d cohorts; methods=%s\n",
            nrow(Xf), ncol(Xf), length(unique(cohort)),
            paste(methods, collapse = ",")))

score_rows <- list()
emit <- function(m, res, runtime) {
  sc <- res$scores
  score_rows[[length(score_rows) + 1L]] <<- data.frame(
    method = m, feature = names(sc), score = as.numeric(sc),
    rank = match(names(sc), res$ranked),
    native = names(sc) %in% res$native, runtime_min = runtime,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-4s native=%2d  (%.1f min)\n", m, length(res$native), runtime))
}

for (m in methods) {
  t0 <- Sys.time()
  if (m == "M5" && "M4" %in% methods) next          # derived from M4 below
  if (m == "M9s" && "M9" %in% methods) next         # shares M9 selection
  if (m == "M4" && "M5" %in% methods) {
    # one env-stability fit, both aggregations (same seed => same fits)
    es <- stabl_environment_stability(
      x = Xf, y = y, environments = cohort, lambda_grid = "auto",
      base_learner = "rf", artificial_type = "knockoff_equi",
      n_bootstraps = 200L, sample_fraction = 1, replace = TRUE,
      random_state = 20261013)
    rt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    mins <- get_min_scores(es)
    means <- setNames(rowMeans(es$env_scores), rownames(es$env_scores))
    ord <- order(-mins, -means[names(mins)], names(mins))
    emit("M4", list(ranked = names(mins)[ord],
                    native = names(mins)[mins >= 0.5],
                    scores = mins), rt)
    esc <- es$env_scores
    votes <- rowSums(esc >= 0.5); means5 <- rowMeans(esc)
    ord5 <- order(-votes, -means5, rownames(esc))
    emit("M5", list(ranked = rownames(esc)[ord5],
                    native = rownames(esc)[votes >= 5L],
                    scores = setNames(as.numeric(votes), rownames(esc))), 0)
    next
  }
  sel <- c(ws8v2_selectors, ws8v3_selectors)[[m]]
  stopifnot(!is.null(sel))
  res <- sel(Xf, y, cohort, covar = covar)
  rt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  emit(m, res, rt)
}

fa_df <- do.call(rbind, score_rows)
write.csv(fa_df, file.path(out_dir, "full_atlas_scores.csv"), row.names = FALSE)
cat("Done. Outputs in ", out_dir, "\n")
