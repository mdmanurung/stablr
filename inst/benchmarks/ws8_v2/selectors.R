# WS-8 v2 M3 — Selector matrix (M0–M5), training-only.
#
# Each selector is a function(X_train, y_train, cohort_train) returning a list:
#   ranked : character vector of feature names in rank order (deterministic)
#   native : character vector of the native selected panel (may be empty)
#   scores : named numeric training score used for ranking
# All operate on the 59 cell-type features of the TRAINING cohorts only.
#
# M4 (environment-aware STABL) is the primary method.  M3 (pooled STABL) is the
# controlled pair isolating the value of environment-aware aggregation: same
# learner (rf), artificial features (knockoff_equi), and bootstrap budget (200);
# only the cohort aggregation differs.  M0/M1/M2 are scientific comparators.
#
# Application case study — NOT calibration evidence.

# Deterministic ranking helper: score desc, tie-break name asc.
.rank_by <- function(scores) {
  o <- order(-scores, names(scores))
  names(scores)[o]
}

# Per-cohort Spearman rho and BH-adjusted association p (linear model
# feature ~ age is NOT used here; the primary endpoint uses Spearman rho, so the
# comparators are built on the same association scale).  For M1's vote we use a
# per-cohort Spearman correlation test with BH adjustment across features.
.per_cohort_assoc <- function(X, y, cohort) {
  feats <- colnames(X); cohs <- sort(unique(cohort))
  rho <- matrix(NA_real_, length(feats), length(cohs),
                dimnames = list(feats, cohs))
  p   <- rho
  for (cv in cohs) {
    idx <- cohort == cv
    for (f in feats) {
      ct <- suppressWarnings(stats::cor.test(X[idx, f], y[idx],
                                             method = "spearman",
                                             exact = FALSE))
      rho[f, cv] <- unname(ct$estimate)
      p[f, cv]   <- ct$p.value
    }
  }
  list(rho = rho, p = p)
}

# M0 — Pooled age-association meta-analysis (cohort-balanced mean rho).
sel_m0_meta <- function(X, y, cohort) {
  a <- .per_cohort_assoc(X, y, cohort)
  sc <- setNames(rowMeans(a$rho, na.rm = TRUE), colnames(X))
  list(ranked = .rank_by(sc),
       native = names(sc)[abs(sc) >= 0],  # native = all (no selection); see note
       scores = sc)
}

# M1 — Original-style cohort significance vote (>= 5/7 BH-adjusted).
sel_m1_vote <- function(X, y, cohort, vote_min = 5L) {
  a <- .per_cohort_assoc(X, y, cohort)
  padj <- apply(a$p, 2, p.adjust, method = "BH")
  votes <- rowSums(padj < 0.05, na.rm = TRUE)
  meta_rho <- rowMeans(a$rho, na.rm = TRUE)
  # rank by vote count desc, then |meta rho| desc, then name
  ord <- order(-votes, -abs(meta_rho), colnames(X))
  sc <- setNames(votes + abs(meta_rho) / 1000, colnames(X))  # for reproducibility
  list(ranked = colnames(X)[ord],
       native = colnames(X)[votes >= vote_min],
       scores = setNames(as.numeric(votes), colnames(X)))
}

# M2 — Pooled sparse multivariable selection (lasso stability, same budget).
sel_m2_lasso <- function(X, y, cohort, n_boot = 200L, seed = 20261013) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  freq <- setNames(numeric(p), colnames(X))
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    fit <- tryCatch(glmnet::cv.glmnet(X[idx, ], y[idx], alpha = 1, nfolds = 5),
                    error = function(e) NULL)
    if (is.null(fit)) next
    cf <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    freq <- freq + as.integer(cf != 0)
  }
  freq <- freq / n_boot
  list(ranked = .rank_by(freq),
       native = names(freq)[freq >= 0.5],
       scores = freq)
}

# M3 — Pooled STABL (ignores cohort partition).
sel_m3_pooled_stabl <- function(X, y, cohort, n_boot = 200L, seed = 20261013) {
  fit <- stabl_fit(x = X, y = y, lambda_grid = "auto", base_learner = "rf",
                   artificial_type = "knockoff_equi", n_bootstraps = n_boot,
                   sample_fraction = 1, replace = TRUE, random_state = seed)
  sc <- get_importances(fit)  # max-over-lambda stability score per feature
  if (is.null(names(sc))) names(sc) <- colnames(X)
  list(ranked = .rank_by(sc),
       native = names(sc)[sc >= 0.5],
       scores = sc)
}

# M4 — Environment-aware STABL (PRIMARY).  Per-cohort STABL, min-across-cohorts.
sel_m4_env_stabl <- function(X, y, cohort, n_boot = 200L, seed = 20261013) {
  es <- stabl_environment_stability(x = X, y = y, environments = cohort,
                                    lambda_grid = "auto", base_learner = "rf",
                                    artificial_type = "knockoff_equi",
                                    n_bootstraps = n_boot, sample_fraction = 1,
                                    replace = TRUE, random_state = seed)
  mins <- get_min_scores(es)
  means <- setNames(rowMeans(es$env_scores), rownames(es$env_scores))
  # rank by min desc, tie-break mean desc, then name (locked rule)
  ord <- order(-mins, -means[names(mins)], names(mins))
  list(ranked = names(mins)[ord],
       native = names(mins)[mins >= 0.5],
       scores = mins, env_scores = es$env_scores, mean_scores = means)
}

# M5 — Env-aware STABL, majority aggregation (secondary).  >= 5/7 cohorts >= 0.5.
sel_m5_env_majority <- function(X, y, cohort, n_boot = 200L, seed = 20261013,
                                vote_min = 5L) {
  es <- stabl_environment_stability(x = X, y = y, environments = cohort,
                                    lambda_grid = "auto", base_learner = "rf",
                                    artificial_type = "knockoff_equi",
                                    n_bootstraps = n_boot, sample_fraction = 1,
                                    replace = TRUE, random_state = seed)
  esc <- es$env_scores
  votes <- rowSums(esc >= 0.5)
  means <- rowMeans(esc)
  ord <- order(-votes, -means, rownames(esc))
  list(ranked = rownames(esc)[ord],
       native = rownames(esc)[votes >= vote_min],
       scores = setNames(as.numeric(votes), rownames(esc)),
       env_scores = esc)
}

# Registry: method id -> selector function.
ws8v2_selectors <- list(
  M0 = sel_m0_meta,
  M1 = sel_m1_vote,
  M2 = sel_m2_lasso,
  M3 = sel_m3_pooled_stabl,
  M4 = sel_m4_env_stabl,
  M5 = sel_m5_env_majority
)
