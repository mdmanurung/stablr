# WS-8 v2 — metric and export helpers (M1 export repairs + shared metrics).
#
# Fixes the three confirmed v1 defects in versioned v2 code and provides the
# macro/micro aggregation used throughout the benchmark.  Pure functions; no
# side effects.  Sourced by run.R, null_probe.R, and figures.R.
#
# Application case study — NOT calibration evidence.

# ---- Macro / micro aggregation ------------------------------------------------
# macro: unweighted mean of per-cohort values (each cohort one generalization unit).
# micro: donor-weighted mean (each donor equal).  They differ when cohort sizes differ.
macro_mean <- function(values, groups) {
  # values: numeric; groups: per-value cohort label
  stopifnot(length(values) == length(groups))
  per <- tapply(values, groups, mean)
  mean(per)
}
micro_mean <- function(values, groups = NULL) {
  mean(values)
}

# Per-cohort MAE table -> macro and micro MAE from raw donor errors.
# errors: numeric |pred - obs|; cohort: per-donor label.
mae_macro_micro <- function(errors, cohort) {
  per <- tapply(errors, cohort, mean)
  list(macro = mean(per), micro = mean(errors), per_cohort = per)
}

# ---- Export repairs (v1 defects A1/A2) ----------------------------------------
# Correct per-feature mean score: rowMeans over the feature x cohort score matrix,
# with name-safe indexing.
safe_mean_scores <- function(env_scores) {
  # env_scores: features x cohorts matrix with feature rownames
  stopifnot(!is.null(rownames(env_scores)))
  setNames(rowMeans(env_scores, na.rm = TRUE), rownames(env_scores))
}

# Correct cross-object selection flag: match by feature name BEFORE numeric coercion.
# scores: named numeric (e.g. get_min_scores(env_lasso)); features: character vector.
safe_match_scores <- function(scores, features) {
  stopifnot(!is.null(names(scores)))
  out <- unname(scores[features])
  names(out) <- features
  out
}

# ---- Primary endpoint: signed-rho replication ---------------------------------
# features: character vector (the panel); train_sign: named {-1,0,+1};
# heldout_rho: named Spearman rho in the held-out cohort.
# Returns mean(d_j * rho_j,e) over the panel, plus per-feature signed values.
signed_rho_replication <- function(features, train_sign, heldout_rho) {
  stopifnot(all(features %in% names(train_sign)),
            all(features %in% names(heldout_rho)))
  signed <- unname(train_sign[features]) * unname(heldout_rho[features])
  names(signed) <- features
  list(R = mean(signed), per_feature = signed,
       sign_agreement = mean(signed > 0))
}
