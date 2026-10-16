# WS-8 v2 M3 — Selector matrix (M0–M5), training-only.
#
# Each selector is a function(X_train, y_train, cohort_train, covar = NULL)
# returning a list:
#   ranked : character vector of feature names in rank order (deterministic)
#   native : character vector of the native selected panel (may be empty)
#   scores : named numeric training score used for ranking
#   detail : optional features x cohorts audit matrix (see detail_type)
# All operate on the 59 cell-type features of the TRAINING cohorts only.
# `covar` (sex, ethnicity; rows aligned to training donors) is used only by
# M1's adjusted-LM vote (protocol section 8); all other selectors ignore it.
#
# M4 (environment-aware STABL) is the primary method.  M3 (pooled STABL) is the
# controlled pair isolating the value of environment-aware aggregation: same
# learner (rf), artificial features (knockoff_equi), and bootstrap budget (200);
# only the cohort aggregation differs.  M0/M1/M2 are scientific comparators.
#
# Amendment 1 (2026-10-15): M0 is now a DerSimonian-Laird random-effects meta
# estimate (Fisher-z) ranked by magnitude, and M1 uses the locked adjusted-LM
# vote of protocol section 8; both were brought into conformance with the lock
# after the M7 adversarial review (findings 4 and 5).
#
# Application case study — NOT calibration evidence.

# Deterministic ranking helper: score desc, tie-break name asc.
.rank_by <- function(scores) {
  o <- order(-scores, names(scores))
  names(scores)[o]
}

# Per-cohort Spearman rho and BH-adjusted association p.  The primary endpoint
# uses Spearman rho, so the comparators' association scale is Spearman; the
# adjusted-LM p-values of the locked M1 rule come from .per_cohort_assoc_lm.
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

# DerSimonian-Laird random-effects meta-analysis of per-cohort Spearman rhos
# on the Fisher-z scale (protocol section 4, M0 row).  Returns pooled estimates
# back on the rho scale.  rho_mat: features x cohorts; n_vec: cohort sizes.
.re_meta_rho <- function(rho_mat, n_vec) {
  z <- atanh(pmin(pmax(rho_mat, -0.9999), 0.9999))
  v <- matrix(rep(1 / (n_vec - 3), each = nrow(z)), nrow = nrow(z))
  wf <- 1 / v
  mu_f <- rowSums(z * wf) / rowSums(wf)
  Q <- rowSums(wf * (z - mu_f)^2)
  df <- ncol(z) - 1
  Cc <- rowSums(wf) - rowSums(wf^2) / rowSums(wf)
  tau2 <- pmax(0, (Q - df) / Cc)
  wr <- 1 / (v + tau2)
  mu_re <- rowSums(z * wr) / rowSums(wr)
  setNames(tanh(mu_re), rownames(rho_mat))
}

# Per-cohort adjusted association: lm(feature ~ age + sex + ancestry), the
# locked conceptual comparator of protocol section 8 (the source paper's rule),
# adapted to 7 training cohorts.  Covariate terms constant within a cohort are
# dropped for that cohort (ancestry is constant in ABF300/CIMA/OneK1K/SPAC).
# Returns the features x cohorts matrix of age-coefficient p-values.
.per_cohort_assoc_lm <- function(X, y, cohort, covar) {
  feats <- colnames(X); cohs <- sort(unique(cohort))
  p <- matrix(NA_real_, length(feats), length(cohs),
              dimnames = list(feats, cohs))
  for (cv in cohs) {
    idx <- cohort == cv
    sex_cv <- droplevels(factor(covar$sex[idx]))
    eth_cv <- droplevels(factor(covar$ethnicity[idx]))
    rhs <- c("age", if (nlevels(sex_cv) > 1L) "sex",
             if (nlevels(eth_cv) > 1L) "eth")
    df <- data.frame(age = y[idx], sex = sex_cv, eth = eth_cv)
    df$Y <- X[idx, , drop = FALSE]  # matrix response -> multivariate lm
    fml <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")))
    fit <- stats::lm(fml, data = df)
    sms <- summary(fit)  # one coefficient table per response, in column order
    for (j in seq_along(feats)) {
      cf <- sms[[j]]$coefficients
      p[feats[j], cv] <- if ("age" %in% rownames(cf)) cf["age", 4] else NA_real_
    }
  }
  p
}

# M0 — Pooled age-association meta-analysis (random-effects, Fisher-z DL).
sel_m0_meta <- function(X, y, cohort, covar = NULL) {
  a <- .per_cohort_assoc(X, y, cohort)
  n_vec <- as.numeric(table(cohort)[colnames(a$rho)])
  sc <- .re_meta_rho(a$rho, n_vec)
  list(ranked = .rank_by(abs(sc)),   # magnitude ranking: direction-agnostic
       native = character(0),        # no native selection rule (ranking comparator)
       scores = sc,
       detail = a$rho, detail_type = "spearman_rho")
}

# M1 — Original-style cohort significance vote (>= 5/7 BH-adjusted), with the
# locked adjusted-LM per-cohort test of protocol section 8.
sel_m1_vote <- function(X, y, cohort, covar = NULL, vote_min = 5L) {
  stopifnot("M1 requires covar (sex, ethnicity) for the adjusted LM" = !is.null(covar))
  a <- .per_cohort_assoc(X, y, cohort)               # meta-evidence scale
  p_lm <- .per_cohort_assoc_lm(X, y, cohort, covar)  # locked adjusted rule
  padj <- apply(p_lm, 2, p.adjust, method = "BH")
  votes <- rowSums(padj < 0.05, na.rm = TRUE)
  n_vec <- as.numeric(table(cohort)[colnames(a$rho)])
  meta_rho <- .re_meta_rho(a$rho, n_vec)
  # rank by vote count desc, then |RE meta rho| desc, then name
  ord <- order(-votes, -abs(meta_rho), colnames(X))
  list(ranked = colnames(X)[ord],
       native = colnames(X)[votes >= vote_min],
       scores = setNames(as.numeric(votes), colnames(X)),
       detail = padj, detail_type = "lm_padj")
}

# M2 — Pooled sparse multivariable selection (lasso stability, same budget).
sel_m2_lasso <- function(X, y, cohort, covar = NULL, n_boot = 200L, seed = 20261013) {
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
sel_m3_pooled_stabl <- function(X, y, cohort, covar = NULL, n_boot = 200L, seed = 20261013) {
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
sel_m4_env_stabl <- function(X, y, cohort, covar = NULL, n_boot = 200L, seed = 20261013) {
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
       scores = mins, env_scores = es$env_scores, mean_scores = means,
       detail = es$env_scores, detail_type = "sel_freq")
}

# M5 — Env-aware STABL, majority aggregation (secondary).  >= 5/7 cohorts >= 0.5.
sel_m5_env_majority <- function(X, y, cohort, covar = NULL, n_boot = 200L, seed = 20261013,
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
       env_scores = esc,
       detail = esc, detail_type = "sel_freq")
}

# ---- S2 probe machinery (followup_s1_s2.md; NOT part of the locked matrix) ----
# Cross-cohort quantile aggregation with the locked M4 tie-break
# (aggregate desc, mean desc, name asc).  q = 0 reproduces min-aggregation.
# Quantile type 7 (R default) interpolates between order statistics.
aggregate_env_scores <- function(esc, q = 0) {
  agg <- setNames(apply(esc, 1, stats::quantile, probs = q, names = FALSE),
                  rownames(esc))
  means <- setNames(rowMeans(esc), rownames(esc))
  ord <- order(-agg, -means[names(agg)], names(agg))
  list(ranked = names(agg)[ord], native = names(agg)[agg >= 0.5],
       scores = agg, mean_scores = means)
}

# M4 variant factory: quantile q aggregation and/or alternative knockoff
# construction.  Used by probe selectors M4q10/M4q25/M4mvr below and by
# benchmarks/ws8_v2/s2b_aggregation_quantile.R (which re-aggregates the
# persisted M4 env scores without refitting — same seed, same fits).
make_m4_variant <- function(q = 0, artificial_type = "knockoff_equi") {
  force(q); force(artificial_type)
  function(X, y, cohort, covar = NULL, n_boot = 200L, seed = 20261013) {
    es <- stabl_environment_stability(x = X, y = y, environments = cohort,
                                      lambda_grid = "auto", base_learner = "rf",
                                      artificial_type = artificial_type,
                                      n_bootstraps = n_boot, sample_fraction = 1,
                                      replace = TRUE, random_state = seed)
    agg <- aggregate_env_scores(es$env_scores, q)
    c(agg, list(env_scores = es$env_scores,
                detail = es$env_scores, detail_type = "sel_freq"))
  }
}

# Registry: method id -> selector function.
ws8v2_selectors <- list(
  M0 = sel_m0_meta,
  M1 = sel_m1_vote,
  M2 = sel_m2_lasso,
  M3 = sel_m3_pooled_stabl,
  M4 = sel_m4_env_stabl,
  M5 = sel_m5_env_majority,
  # S2 probes — run only via explicit WS8_METHODS subsets
  M4q10 = make_m4_variant(0.10),
  M4q25 = make_m4_variant(0.25),
  M4mvr = make_m4_variant(0, "knockoff_mvr")
)
