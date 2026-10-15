# WS-8 v2 M5 / S0 — Real-data reference fit for the null probe.
#
# Runs the M4 environment-aware STABL selector on the FULL ATLAS (all 8
# cohorts, real age) through the identical code path as null_probe.R (same
# fold transform, same engine settings, same seed family).  This is the
# auditable real-data comparator for the S0 null: the null permutations and
# this reference share the 8-cohort full-atlas design, so the comparison of
# native panel size and score distributions is apples-to-apples (review
# finding 2/9).  The LOCO native panels (7-cohort training fits) remain
# available per fold in training_scores.csv.
#
# Usage: Rscript benchmarks/ws8_v2/real_reference.R   # run from repo root
#   env WS8_V2_IN    dir holding donor_tables.rds (default .../results/m1)
#   env WS8_V2_OUT   output dir (default .../results/null)
#
# Application case study — NOT calibration evidence.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("benchmarks/ws8_v2/metrics.R")
source("benchmarks/ws8_v2/loco.R")
source("benchmarks/ws8_v2/selectors.R")

in_dir  <- Sys.getenv("WS8_V2_IN", file.path("benchmarks/ws8_v2/results/m1"))
out_dir <- Sys.getenv("WS8_V2_OUT", file.path("benchmarks/ws8_v2/results/null"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Design matrix (identical construction to null_probe.R) -------------------
obj <- readRDS(file.path(in_dir, "donor_tables.rds"))
ct_mat <- obj$ct_mat; total <- obj$total; meta <- obj$meta
prop <- (ct_mat + 0.5) / (total + 0.5)
X <- log(prop / (1 - prop)); colnames(X) <- paste0("CT:", colnames(X))
meta <- meta[match(rownames(X), meta$donor), ]
y <- setNames(as.numeric(meta$age), meta$donor)
cohort <- setNames(as.character(meta$dataset), meta$donor)

recipe <- fit_fold_transform(X, cohort)
cat(sprintf("transform: %d features dropped (zero-variance within a cohort)\n",
            length(recipe$drop)))
Xf <- apply_fold_transform(recipe, X)

# ---- Full-atlas real M4 fit ---------------------------------------------------
t0 <- Sys.time()
es <- stabl_environment_stability(
  x = Xf, y = y, environments = cohort, lambda_grid = "auto",
  base_learner = "rf", artificial_type = "knockoff_equi",
  n_bootstraps = 200L, sample_fraction = 1, replace = TRUE,
  random_state = 20261013)
mins <- get_min_scores(es)
means <- setNames(rowMeans(es$env_scores), rownames(es$env_scores))

ref <- data.frame(
  feature = names(mins),
  min_score = as.numeric(mins),
  mean_score = as.numeric(means[names(mins)]),
  native = as.numeric(mins) >= 0.5,
  stringsAsFactors = FALSE)
write.csv(ref, file.path(out_dir, "real_reference_scores.csv"), row.names = FALSE)
write.csv(data.frame(feature = rownames(es$env_scores), es$env_scores,
                     row.names = NULL, check.names = FALSE),
          file.path(out_dir, "real_reference_env_scores.csv"), row.names = FALSE)

cat(sprintf("\n== Real-data full-atlas M4 reference (%.1f min) ==\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("native panel (min-agg >= 0.5): %d of %d features\n",
            sum(ref$native), nrow(ref)))
cat(sprintf("max min-score: %.3f\n", max(ref$min_score)))
cat("Done. Outputs in ", out_dir, "\n")
