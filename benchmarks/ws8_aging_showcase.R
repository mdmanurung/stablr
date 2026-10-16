# WS-8 — Application showcase: stabl on the Terekhova et al. (Immunity, 2026)
# human immune aging atlas.  Design predeclared in ws8_showcase_predesign.md
# (written before this run).  Application case study — NOT calibration evidence.

# Run against THIS checkout, not an installed copy (the ordinal routing and
# current RF adapter live here).
pkgload::load_all(".", compile = FALSE, quiet = TRUE)

set.seed(20261013)

data_dir <- Sys.getenv("WS8_DATA_DIR", "/workspace")
out_dir  <- Sys.getenv("WS8_OUT_DIR", "benchmarks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("== WS-8 step 1: load donor-level tables ==\n")
counts_raw <- read.delim(file.path(data_dir, "donor_celltype_counts.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)
key_parts <- strsplit(counts_raw[[1]], "|", fixed = TRUE)
counts <- data.frame(
  donor    = vapply(key_parts, `[`, character(1), 1),
  celltype = vapply(key_parts, `[`, character(1), 2),
  n        = as.integer(counts_raw[[2]]),
  stringsAsFactors = FALSE
)
meta_raw <- read.delim(file.path(data_dir, "donor_meta.tsv"), header = FALSE,
                       stringsAsFactors = FALSE)
stopifnot(ncol(meta_raw) >= 3)
# donor_meta.tsv rows: "META\t<donor>\t<Dataset|age|sex|ethnicity|CMV>"
meta_raw <- meta_raw[meta_raw[[1]] == "META", ]
parts <- strsplit(meta_raw[[3]], "|", fixed = TRUE)
meta <- data.frame(
  donor     = meta_raw[[2]],
  dataset   = vapply(parts, `[`, character(1), 1),
  age       = as.numeric(vapply(parts, `[`, character(1), 2)),
  sex       = vapply(parts, `[`, character(1), 3),
  ethnicity = vapply(parts, `[`, character(1), 4),
  cmv       = vapply(parts, `[`, character(1), 5),
  stringsAsFactors = FALSE
)
cat(sprintf("donors: %d; donor x celltype pairs: %d\n",
            nrow(meta), nrow(counts)))

# Cell-type matrix: donors x types, add-half smoothing + logit (base R).
ct_levels <- sort(unique(counts$celltype))
donor_levels <- intersect(meta$donor, unique(counts$donor))
ct_mat <- matrix(0L, nrow = length(donor_levels), ncol = length(ct_levels),
                 dimnames = list(donor_levels, ct_levels))
idx <- cbind(match(counts$donor, donor_levels),
             match(counts$celltype, ct_levels))
ct_mat[idx] <- counts$n
total <- rowSums(ct_mat)
prop  <- (ct_mat + 0.5) / (total + 0.5)
logit <- function(p) log(p / (1 - p))
X_ct  <- logit(prop)
colnames(X_ct) <- paste0("CT:", colnames(X_ct))
cat(sprintf("cell types: %d\n", ncol(X_ct)))

# Covariate features: sex, CMV, ancestry (levels with < 25 donors -> Other).
meta <- meta[match(rownames(X_ct), meta$donor), , drop = FALSE]
stopifnot(identical(meta$donor, rownames(X_ct)))

sex_f   <- ifelse(is.na(meta$sex), NA, as.integer(meta$sex == "Female"))
cmv_f   <- ifelse(is.na(meta$cmv), NA,
                  as.integer(meta$cmv %in% c("Positive", "positive", "TRUE")))
eth     <- meta$ethnicity
eth[is.na(eth)] <- "Unknown"
tab_eth <- table(eth)
small   <- names(tab_eth)[tab_eth < 25]
eth[eth %in% small] <- "Other"
eth_f   <- factor(eth)
cat(sprintf("ancestry levels: %s\n",
            paste(sprintf("%s(%d)", names(tab_eth), tab_eth),
                  collapse = ", ")))

keep <- !is.na(meta$age) & !is.na(sex_f) & !is.na(cmv_f)
cat(sprintf("dropping %d donors with missing age/sex/CMV\n", sum(!keep)))
X_ct   <- X_ct[keep, , drop = FALSE]
meta   <- meta[keep, , drop = FALSE]
sex_f  <- sex_f[keep]
cmv_f  <- cmv_f[keep]
eth_f  <- droplevels(eth_f[keep])
y_age  <- setNames(meta$age, rownames(X_ct))
env    <- setNames(meta$dataset, rownames(X_ct))

keep_levels <- setdiff(levels(eth_f), "Other")
eth_mat <- model.matrix(~ 0 + eth_f)[, paste0("eth_f", keep_levels),
                                    drop = FALSE]
colnames(eth_mat) <- paste0("ANC:", keep_levels)

X <- cbind(X_ct,
           data.frame("COV:sex" = sex_f, "COV:CMV" = cmv_f,
                      check.names = FALSE),
           eth_mat)
X <- as.matrix(X)
storage.mode(X) <- "numeric"
cat(sprintf("final design: %d donors x %d features; cohorts: %s\n",
            nrow(X), ncol(X),
            paste(sort(unique(env)), collapse = ", ")))

# Drop features that are constant within any cohort: a constant column has
# no within-cohort age signal, and its zero variance makes the equicorrelated
# knockoff construction numerically infeasible (singular covariance). Only
# covariate columns (COV:CMV, ancestry dummies) are affected; no cell-type
# feature is constant in any cohort.
zv_features <- unique(unlist(lapply(sort(unique(env)), function(e) {
  Xe <- X[env == e, , drop = FALSE]
  v <- apply(Xe, 2, function(cc) var(cc, na.rm = TRUE) == 0) # X is a matrix
  names(v)[v]
})))
cat("dropping zero-variance-in-some-cohort features:", 
    paste(sort(zv_features), collapse = ", "), "\n")
X <- X[, setdiff(colnames(X), zv_features), drop = FALSE]

cat("== WS-8 step 2: environment stability (rf + knockoff_equi, 200 bootstraps) ==\n")
t0 <- Sys.time()
env_rf <- stabl_environment_stability(
  x               = X,
  y               = y_age,
  environments    = env,
  lambda_grid     = "auto",
  base_learner    = "rf",
  artificial_type = "knockoff_equi",
  n_bootstraps    = 200L,
  sample_fraction = 1,
  replace         = TRUE,
  random_state    = 20261013
)
cat(sprintf("rf env stability done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("== WS-8 step 3: sensitivity (lasso + knockoff_equi) ==\n")
t0 <- Sys.time()
env_lasso <- stabl_environment_stability(
  x               = X,
  y               = y_age,
  environments    = env,
  lambda_grid     = "auto",
  base_learner    = "lasso",
  artificial_type = "knockoff_equi",
  n_bootstraps    = 200L,
  sample_fraction = 1,
  replace         = TRUE,
  random_state    = 20261013
)
cat(sprintf("lasso env stability done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Selection tables --------------------------------------------------------
sel_rf <- data.frame(
  feature    = names(get_min_scores(env_rf)),
  min_score  = as.numeric(get_min_scores(env_rf)),
  mean_score = as.numeric(colMeans(env_rf$env_scores))[names(get_min_scores(env_rf))],
  row.names  = NULL, stringsAsFactors = FALSE
)
sel_rf$selected_minagg <- sel_rf$min_score >= 0.5
sel_rf$selected_lasso  <- as.numeric(get_min_scores(env_lasso))[sel_rf$feature] >= 0.5

# Per-cohort directional read: Spearman correlation of feature vs age.
cohorts <- sort(unique(env))
cor_mat <- matrix(NA_real_, ncol(X), length(cohorts),
                  dimnames = list(colnames(X), cohorts))
for (cv in cohorts) {
  idx <- env == cv
  for (j in seq_len(ncol(X))) {
    cor_mat[j, cv] <- suppressWarnings(
      stats::cor(X[idx, j], y_age[idx], method = "spearman",
                 use = "complete.obs")
    )
  }
}
sel_cols <- paste0("rho_", cohorts)
sel_rf <- cbind(sel_rf, as.data.frame(cor_mat[sel_rf$feature, , drop = FALSE]))
rownames(sel_rf) <- NULL
names(sel_rf)[match(cohorts, names(sel_rf))] <- sel_cols

write.csv(sel_rf, file.path(out_dir, "ws8_results.csv"), row.names = FALSE)

# ---- LOCO age prediction -----------------------------------------------------
cat("== WS-8 step 4: LOCO age prediction ==\n")
suppressPackageStartupMessages(library(glmnet))
selected_rf <- sel_rf$feature[sel_rf$selected_minagg]
cat(sprintf("min-agg selected (rf): %d features: %s\n", length(selected_rf),
            paste(selected_rf, collapse = ", ")))

loco_mae <- function(Xt, yt, ev) {
  cohs <- sort(unique(ev))
  res <- do.call(rbind, lapply(cohs, function(h) {
    tr <- ev != h; te <- ev == h
    if (sum(te) < 5L) return(NULL)
    cvfit <- cv.glmnet(Xt[tr, , drop = FALSE], yt[tr], alpha = 0,
                       nfolds = 5)
    pr <- as.numeric(predict(cvfit, newx = Xt[te, , drop = FALSE],
                             s = "lambda.min"))
    data.frame(cohort = h, n = sum(te), mae = mean(abs(pr - yt[te])),
               slope = unname(coef(lm(pr ~ yt[te]))[2]),
               pred = I(pr), obs = I(yt[te]))
  }))
  res
}
loco_all  <- loco_mae(X, y_age, env)
loco_sel  <- loco_mae(X[, selected_rf, drop = FALSE], y_age, env)
cat("LOCO MAE (all features): ", sprintf("%.2f", mean(loco_all$mae)), "\n")
cat("LOCO MAE (selected):     ", sprintf("%.2f", mean(loco_sel$mae)), "\n")

# ---- Summary -----------------------------------------------------------------
paper_consistent <- c("CT:MAIT", "CT:CD4+ RTE", "CT:CD8+ RTE",
                      "CT:CD8+ T Naive", "CT:CD8+ Temra", "CT:CD4+ Temra",
                      "CT:CD8+ Tem GZMK+", "CT:CD8+ Tem GZMB+")
hit <- paper_consistent[paper_consistent %in% selected_rf]

summary_lines <- c(
  "# WS-8 — Aging-atlas application showcase summary",
  "",
  sprintf("- Data: %d donors, %d features (59 cell types + covariates), %d cohorts.",
          nrow(X), ncol(X), length(cohorts)),
  sprintf("- Primary (rf + knockoff_equi, 200 bootstraps): %d/%d features pass min-agg (>= 0.5).",
          length(selected_rf), ncol(X)),
  sprintf("- Paper-consistent features recovered by min-agg: %s.",
          if (length(hit)) paste(hit, collapse = ", ") else "none"),
  sprintf("- LOCO age MAE: all-feature ridge %.2f y; selected-feature model %.2f y (paper clock: 8.2 y full / 8.9 y reduced).",
          mean(loco_all$mae), mean(loco_sel$mae)),
  sprintf("- Lasso sensitivity min-agg selected: %d features.",
          sum(sel_rf$selected_lasso)),
  "",
  "## Per-cohort LOCO MAE",
  "",
  sprintf("| cohort | n | MAE all | MAE selected |"),
  sprintf("|---|---|---|---|"),
  sprintf("| %s | %d | %.2f | %.2f |",
          loco_all$cohort, loco_all$n, loco_all$mae, loco_sel$mae),
  "",
  "Application case study — NOT calibration evidence.",
  ""
)
writeLines(summary_lines, file.path(out_dir, "ws8_summary.md"))

# Save R objects for figure generation without refitting.
saveRDS(list(env_rf = env_rf, env_lasso = env_lasso, sel_rf = sel_rf,
             loco_all = loco_all, loco_sel = loco_sel, X = X, y_age = y_age,
             env = env, cor_mat = cor_mat),
        file.path(out_dir, "ws8_objects.rds"))
cat("WS-8 done.\n")
