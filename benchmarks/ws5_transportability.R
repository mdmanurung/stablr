#!/usr/bin/env Rscript
# ============================================================================
# WS-5 — Environment-stability transportability simulation (predeclared)
# ============================================================================
# E = 3 environments, n_env = 80, p = 30 (3 signal). Isolated mechanisms:
#   M1 location shift on 10 NULL features; M2 location shift on 2 SIGNAL
#   features; M3 outcome prevalence drift (binomial).
# Pooled single-fit STABL vs stabl_environment_stability() (min-aggregation).
# Predeclared question: does min-aggregation keep empirical FDP at or below
# the pooled fit's on shifted environments, and what TPR does it cost?
# Run: Rscript benchmarks/ws5_transportability.R   (~1-2 h)

.libPaths(c("/workspace/Rlib-shared", .libPaths()))
if (file.exists("DESCRIPTION")) pkgload::load_all(".", compile = TRUE, quiet = TRUE) else library(stablr)

seeds <- 1:10
n_env <- 80; p <- 30; E <- 3
n_bootstraps <- 100L
out_dir <- "benchmarks"

make_data <- function(mechanism, seed) {
  set.seed(seed)
  Sigma <- matrix(0.2, p, p); diag(Sigma) <- 1
  for (i in seq_len(p - 3L)) for (j in seq_len(p - 3L)) Sigma[3+i, 3+j] <- 0.5^abs(i-j)
  Sigma[1:3, 1:3] <- 0.3; diag(Sigma) <- 1
  L <- t(chol(Sigma))
  X <- matrix(rnorm(E * n_env * p), E * n_env, p) %*% L
  colnames(X) <- paste0("f", seq_len(p)); rownames(X) <- paste0("s", seq_len(E * n_env))
  env <- factor(rep(paste0("env", seq_len(E)), each = n_env))
  shift_size <- switch(mechanism, M1 = 0.6, M2 = 0.6, M3 = 0)
  if (mechanism %in% c("M1", "M2")) {
    idx <- if (mechanism == "M1") 4:13 else 1:2          # nulls vs signals
    for (e in seq_len(E)) {
      rows <- which(env == paste0("env", e))
      X[rows, idx] <- X[rows, idx, drop = FALSE] + shift_size * (e - 2)  # -0.6, 0, +0.6
    }
  }
  eta <- 0.9 * X[, 1] + 0.7 * X[, 2] - 0.8 * X[, 3]
  prev <- switch(mechanism, M3 = 1/(1 + exp(-eta + (as.integer(env) - 2) * 1.2)),
                 1/(1 + exp(-eta)))
  y <- rbinom(E * n_env, 1, prev)
  list(x = X, y = setNames(y, rownames(X)), environments = env,
       signals = paste0("f", seq_len(3L)), nulls = paste0("f", 4:p))
}

fdp_tpr <- function(selected, signals, nulls) {
  selected <- as.character(selected); n_sel <- length(selected)
  if (n_sel == 0L) return(data.frame(n_selected = 0L, fdp = 0, tpr = 0))
  data.frame(n_selected = n_sel,
             fdp = sum(selected %in% nulls) / n_sel,
             tpr = sum(selected %in% signals) / length(signals))
}

rows <- list()
for (mech in c("M1_null_shift", "M2_signal_shift", "M3_prevalence_drift")) {
  m <- sub("_.*", "", mech)
  for (seed in seeds) {
    d <- make_data(m, seed)
    # Pooled single fit
    pooled <- tryCatch(withCallingHandlers({
      fit <- stabl_fit(x = d$x, y = d$y, lambda_grid = "auto", family = "binomial",
                       artificial_type = "random_permutation",
                       n_bootstraps = n_bootstraps, random_state = seed)
      sel <- as.character(get_feature_names(fit))
      data.frame(mechanism = mech, method = "pooled", environment = "all", seed = seed,
                 fdp_tpr(sel, d$signals, d$nulls))
    }, warning = function(w) if (grepl("anti-conservative|not validated", conditionMessage(w)))
      invokeRestart("muffleWarning")),
      error = function(e) data.frame(mechanism = mech, method = "pooled", environment = "all", seed = seed,
                                     n_selected = NA, fdp = NA, tpr = NA))
    rows[[length(rows) + 1L]] <- pooled
    # Environment-stability workflow (min-aggregation) + per-env selections
    envres <- tryCatch(withCallingHandlers({
      res <- stabl_environment_stability(d$x, d$y, d$environments,
                                         lambda_grid = "auto", family = "binomial",
                                         artificial_type = "random_permutation",
                                         n_bootstraps = n_bootstraps,
                                         random_state = seed)
      mins <- get_min_scores(res)
      sel_min <- names(mins)[mins >= 0.5]  # documented min-score contract (print method; get_support 0.5)
      out_all <- cbind(data.frame(mechanism = mech, method = "env_minagg", environment = "all", seed = seed),
                       fdp_tpr(sel_min, d$signals, d$nulls))
      per_env <- lapply(levels(d$environments), function(lv) {
        sel_e <- as.character(get_feature_names(res$env_fits[[lv]]))  # FDP+ threshold per env
        cbind(data.frame(mechanism = mech, method = "env_per_env", environment = lv, seed = seed),
              fdp_tpr(sel_e, d$signals, d$nulls))
      })
      rbind(out_all, do.call(rbind, per_env))
    }, warning = function(w) if (grepl("anti-conservative|not validated", conditionMessage(w)))
      invokeRestart("muffleWarning")),
      error = function(e) data.frame(mechanism = mech, method = "env_minagg", environment = "error", seed = seed,
                                     n_selected = NA, fdp = NA, tpr = NA))
    rows[[length(rows) + 1L]] <- envres
  }
}
results <- do.call(rbind, rows)
write.csv(results, file.path(out_dir, "ws5_results.csv"), row.names = FALSE)

agg <- aggregate(cbind(n_selected, fdp, tpr) ~ mechanism + method + environment, results, mean, na.rm = TRUE)
sink(file.path(out_dir, "ws5_summary.md"))
cat("# WS-5 — Environment-stability transportability summary\n\n- Predeclared: MATURITY_PREDESIGN.md; E = 3, n_env = 80, p = 30 (3 signal); binomial; random_permutation; 10 seeds; 100 bootstraps.\n\n## Mean FDP / TPR by mechanism, method, environment\n\n")
print(agg, row.names = FALSE)
cat("\n## Predeclared reading\n\n- M1 (null shifts) is the calibrated null probe: min-aggregation should keep FDP at or below the pooled fit.\n",
    "- M2 (signal shifts) is the honest cost probe: min-aggregation trades TPR for transportability.\n",
    "- M3 (prevalence drift) tests the prior-shift regime.\n",
    "- FDP+ is a diagnostic minimizer, not a universal FDR guarantee; simulation-only.\n")
sink()
cat("WS-5 done.\n")
