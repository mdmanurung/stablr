#!/usr/bin/env Rscript
# ============================================================================
# B1 — study-disjoint nested CV benchmark (0.3.0 cross-study round)
# ============================================================================
# NOT a testthat test: a self-contained benchmark script, run via
#   Rscript benchmarks/b1_study_disjoint_cv.R [--smoke]
#
# Question: how optimistic is the sample-i.i.d. nested-CV estimate of
# held-out BER, relative to a study-disjoint estimate, when studies carry
# batch structure?
#
# Design: multi-study simulation, 6 studies x 25 samples (n = 150), p = 50
# features, binary outcome.  Three true signal features carry a fixed class
# effect.  Two batch configurations are compared (per seed, same data for
# both CV arms):
#   * location_shift  : ten null features carry a per-study location shift
#                       drawn once per study from N(0, tau); outcome
#                       prevalence is balanced within every study.
#   * prevalence_drift: per-study outcome prevalence drawn from
#                       Uniform(0.3, 0.7); no feature location shifts.
# For each configuration the SAME data are evaluated twice with
# stabl_multiomic_nested_cv():
#   * sample-iid      : study_labels = NULL (folds split samples)
#   * study-disjoint  : study_labels supplied (folds split studies)
# with matched outer_v = 3 and inner_v = 2, lasso base learner, and the
# nested-CV default artificial_type = "random_permutation".
#
# Expected direction: sample-iid BER is optimistic (lower) relative to
# study-disjoint BER when study-level structure is shared across the split.
# The location-shift-only configuration is a deliberate null probe: a
# class-independent location shift on null features under balanced
# prevalence gives sample-iid folds nothing to leak, so the gap should be
# approximately zero there; outcome-prevalence drift (prior shift) is the
# mechanism expected to produce a positive gap.
#
# Honesty caveats, recorded in the outputs:
# * Simulation-only evidence with two batch models; real cohorts may
#   combine location shift, prevalence drift, and covariance drift.
# * Both arms use the same fold counts; the comparison isolates the CV unit
#   (sample vs study), not the number of folds.
# * The study-disjoint arm trains on less effective data per fold (whole
#   studies held out), so part of any gap is reduced training size rather
#   than leakage alone; the direction of that confound is conservative for
#   the reliability claim (it can only overstate, not understate, the gap).
# * No reharmonization or batch correction is applied or claimed.
# * Results depend on the seed grid; the CSV ships all per-seed rows.
# ============================================================================

.libPaths(c("/workspace/Rlib-shared", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

n_studies    <- 6L
n_per_study  <- 25L
p            <- 50L
n_signal     <- 3L
n_batched    <- 10L
beta_signal  <- 0.8
tau_batch    <- 0.6
prev_range   <- c(0.3, 0.7)
smoke_seeds  <- 1:2
full_seeds   <- 1:10
smoke_boot   <- 20L
full_boot    <- 100L
outer_v      <- 3L
inner_v      <- 2L

configs <- list(
  location_shift   = list(prevalence_drift = FALSE),
  prevalence_drift = list(prevalence_drift = TRUE)
)

seeds        <- if (smoke) smoke_seeds else full_seeds
n_bootstraps <- if (smoke) smoke_boot else full_boot

# ---- Package loading (F5 pattern) -------------------------------------------
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", compile = TRUE, quiet = TRUE)
} else if (requireNamespace("stablr", quietly = TRUE)) {
  library(stablr)
} else {
  stop("Neither pkgload nor an installed stablr is available.")
}

run_one <- function(seed) {
  run_config <- function(config_name) {
    prevalence_drift <- configs[[config_name]]$prevalence_drift
    set.seed(seed)
    n <- n_studies * n_per_study
    study_labels <- setNames(
      rep(paste0("study", seq_len(n_studies)), each = n_per_study),
      paste0("s", seq_len(n))
    )

    if (isTRUE(prevalence_drift)) {
      # per-study outcome prevalence drawn once per study (prior shift)
      prev <- stats::runif(n_studies, prev_range[1L], prev_range[2L])
      y_chr <- unlist(lapply(seq_len(n_studies), function(s) {
        sample(c("A", "B"), n_per_study, replace = TRUE,
               prob = c(1 - prev[s], prev[s]))
      }))
      y <- setNames(factor(y_chr, levels = c("A", "B")), names(study_labels))
    } else {
      # balanced outcome within every study (up to odd sizes)
      y <- setNames(
        factor(unlist(replicate(n_studies, rep(c("A", "B"),
                                                length.out = n_per_study)))),
        names(study_labels)
      )
    }

    x <- matrix(rnorm(n * p, sd = 1), n, p,
                dimnames = list(names(study_labels), paste0("f", seq_len(p))))
    # true signal: fixed class effect on the first `n_signal` features
    x[, seq_len(n_signal)] <- x[, seq_len(n_signal)] +
      as.numeric(y == "B") * beta_signal
    if (!isTRUE(prevalence_drift)) {
      # batch: per-study location shift on a null-feature subset
      study_shifts <- matrix(rnorm(n_studies * n_batched, sd = tau_batch),
                             n_studies, n_batched,
                             dimnames = list(paste0("study",
                                                    seq_len(n_studies))))
      batch_cols <- n_signal + seq_len(n_batched)
      for (s in seq_len(n_studies)) {
        idx <- study_labels == paste0("study", s)
        x[idx, batch_cols] <- x[idx, batch_cols] +
          matrix(study_shifts[s, ], sum(idx), n_batched, byrow = TRUE)
      }
    }
    x_list <- list(rna = x)

    # B2 diagnostics: record the planted shift magnitude actually realisable
    diag <- stabl_multiomic_study_diagnostics(
      x_list = x_list, y = y, study_labels = study_labels
    )
    shift_row <- diag$feature_shift[diag$feature_shift$block == "rna", ]
    max_shift <- max(shift_row$smd_max_abs)

    run_arm <- function(scheme) {
      labels <- if (scheme == "study_disjoint") study_labels else NULL
      fit <- suppressWarnings(stabl_multiomic_nested_cv(
        x_list = x_list,
        y = y,
        lambda_grid = "auto",
        outer_v = outer_v,
        inner_v = inner_v,
        n_bootstraps = n_bootstraps,
        random_state = seed,
        study_labels = labels
      ))
      perf <- data.frame(
        seed = seed,
        config = config_name,
        scheme = scheme,
        outer_v = outer_v,
        inner_v = inner_v,
        n_bootstraps = n_bootstraps,
        ber = fit$performance$balanced_error_rate,
        accuracy = fit$performance$accuracy,
        max_abs_smd = max_shift,
        stringsAsFactors = FALSE
      )
      # Per-fold selection sets are persisted so the Round C
      # stability-of-stability analysis can run retroactively on these
      # artifacts without re-simulating.
      sel <- fit$selected_features
      sel$seed <- seed
      sel$config <- config_name
      sel$scheme <- scheme
      list(perf = perf, sel = sel)
    }
    arms <- list(run_arm("sample_iid"), run_arm("study_disjoint"))
    list(
      perf = do.call(rbind, lapply(arms, `[[`, "perf")),
      sel = do.call(rbind, lapply(arms, `[[`, "sel"))
    )
  }
  per_config <- lapply(names(configs), run_config)
  list(
    perf = do.call(rbind, lapply(per_config, `[[`, "perf")),
    sel = do.call(rbind, lapply(per_config, `[[`, "sel"))
  )
}

runs <- lapply(seeds, run_one)
results <- do.call(rbind, lapply(runs, `[[`, "perf"))
selections <- do.call(rbind, lapply(runs, `[[`, "sel"))

out_dir <- file.path("inst", "benchmarks")
write.csv(results, file.path(out_dir, "b1_results.csv"), row.names = FALSE)
write.csv(selections, file.path(out_dir, "b1_selected_features.csv"),
          row.names = FALSE)

summarise_cell <- function(config, scheme) {
  rows <- results[results$config == config & results$scheme == scheme, ]
  c(mean_ber = mean(rows$ber), sd_ber = stats::sd(rows$ber))
}

lines <- c(
  "# B1 — study-disjoint nested CV benchmark",
  "",
  sprintf("Configuration: %d studies x %d samples, p = %d (%d signal); outer_v = %d, inner_v = %d, n_bootstraps = %d; %d seeds; lasso base learner; artificial_type = \"random_permutation\".",
          n_studies, n_per_study, p, n_signal,
          outer_v, inner_v, n_bootstraps, length(seeds)),
  "",
  "Two batch configurations (isolated mechanisms): `location_shift` =",
  sprintf("per-study N(0, %s) location shift on %d null features, balanced prevalence;", format(tau_batch), n_batched),
  sprintf("`prevalence_drift` = per-study outcome prevalence ~ Uniform(%s, %s), no feature shifts.",
          format(prev_range[1L]), format(prev_range[2L])),
  "",
  "| Configuration | sample-iid BER (SD) | study-disjoint BER (SD) | gap |",
  "|---|---|---|---|"
)
for (config in names(configs)) {
  iid <- summarise_cell(config, "sample_iid")
  dis <- summarise_cell(config, "study_disjoint")
  gap <- dis[["mean_ber"]] - iid[["mean_ber"]]
  lines <- c(lines, sprintf(
    "| %s | %.4f (%.4f) | %.4f (%.4f) | %+.4f |",
    config, iid[["mean_ber"]], iid[["sd_ber"]],
    dis[["mean_ber"]], dis[["sd_ber"]], gap
  ))
}
lines <- c(
  lines,
  "",
  "Interpretation: the gap is sample-iid minus study-disjoint mean held-out",
  "BER. A positive gap means the sample-i.i.d. folds understate held-out",
  "error when studies carry the corresponding structure. The",
  "location-shift-only row is a deliberate null probe: class-independent",
  "location shifts on null features under balanced prevalence give",
  "sample-iid folds nothing to leak, so an approximately zero gap is the",
  "calibrated outcome there, and a positive gap is expected only under",
  "outcome-prevalence drift (prior shift). See the honesty caveats in the",
  "script header (simulation-only; two batch models; reduced effective",
  "training size in the study-disjoint arm; no reharmonization claim).",
  ""
)
writeLines(lines, file.path(out_dir, "b1_summary.md"))

for (config in names(configs)) {
  iid <- summarise_cell(config, "sample_iid")
  dis <- summarise_cell(config, "study_disjoint")
  cat(sprintf("B1 done [%s]: sample-iid BER %.4f vs study-disjoint BER %.4f (gap %+.4f)\n",
              config, iid[["mean_ber"]], dis[["mean_ber"]],
              dis[["mean_ber"]] - iid[["mean_ber"]]))
}
