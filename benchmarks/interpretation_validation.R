# IV — Interpretation & validated-workflow benchmark (PR6-B / Gate B).
# Predeclared design: interpretation_validation_predesign.md (2026-10-15),
# written BEFORE any run.  Do not edit the design after running; deviations
# are recorded as dated addenda in the predesign and in iv_summary.md.
#
# Honesty caveats (fixed): simulation-only evidence; FDP+ remains a
# learner x decoy x design diagnostic (no universal FDR claim); the
# leaky-vs-nested contrast (IV-5) is directional on average across seeds,
# not a per-replicate requirement.
#
# Run integrity: every planned design x seed x arm run is attempted and
# accounted for; failures are recorded with error class/message; the script
# exits non-zero if any planned run is missing.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)

out_dir <- Sys.getenv("IV_OUT_DIR", "benchmarks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Smoke mode (CI): IV_SEEDS="1,2" and IV_N_BOOT=10 shrink the run. The
# completeness gate stays self-consistent because planned = configured.
# Gate values at smoke size are NOT comparable to the predeclared 10-seed
# run; smoke asserts end-to-end execution and run accounting only.
seeds <- local({
  s <- Sys.getenv("IV_SEEDS", "")
  if (nzchar(s)) as.integer(strsplit(s, ",", fixed = TRUE)[[1]]) else 1:10
})
n_boot <- as.integer(Sys.getenv("IV_N_BOOT", "50"))
lambda_grid <- data.frame(lambda = c(0.3, 0.2, 0.1))
git_sha <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE),
  error = function(e) NA_character_
)
pkg_version <- tryCatch(
  as.character(utils::packageVersion("stablr")),
  error = function(e) NA_character_
)

# ── Data-generating processes ────────────────────────────────────────────────

# Binary DGP: eta = x %*% beta, centred to prevalence ~0.5, y ~ Bernoulli.
dgp_binary <- function(n, beta, seed, x = NULL) {
  set.seed(seed)
  p <- length(beta)
  if (is.null(x)) {
    x <- matrix(rnorm(n * p), n, p,
                dimnames = list(paste0("s", seq_len(n)), names(beta)))
  }
  eta <- as.vector(x %*% beta)
  eta <- eta - mean(eta)
  pr <- plogis(eta)
  y <- setNames(
    factor(ifelse(stats::rbinom(n, 1L, pr) == 1L, "B", "A"),
           levels = c("A", "B")),
    rownames(x)
  )
  list(x = x, y = y)
}

make_design_data <- function(design, seed) {
  switch(
    design,
    iv0_null = {
      beta <- setNames(rep(0, 10), paste0("f", 1:10))
      dgp_binary(120L, beta, seed)
    },
    iv1_sparse = {
      beta <- setNames(c(1.5, -1.2, 1.0, rep(0, 7)), paste0("f", 1:10))
      dgp_binary(120L, beta, seed)
    },
    iv2_substitutes = {
      set.seed(seed)
      n <- 120L
      f1a <- rnorm(n)
      f1b <- 0.85 * f1a + sqrt(1 - 0.85^2) * rnorm(n)
      f1c <- 0.85 * f1a + sqrt(1 - 0.85^2) * rnorm(n)
      rest <- matrix(rnorm(n * 7L), n, 7L)
      x <- cbind(f1a, f1b, f1c, rest)
      colnames(x) <- c("f1a", "f1b", "f1c", paste0("f", 2:8))
      rownames(x) <- paste0("s", seq_len(n))
      beta <- setNames(c(1.5, 0, 0, -1.2, rep(0, 6)), colnames(x))
      dgp_binary(n, beta, seed, x = x)
    },
    iv3_single = {
      beta <- setNames(c(1.5, rep(0, 4)), paste0("f", 1:5))
      dgp_binary(80L, beta, seed)
    },
    iv4_study_shift = {
      set.seed(seed)
      n_per <- 40L
      n <- 4L * n_per
      p <- 10L
      x <- matrix(rnorm(n * p), n, p,
                  dimnames = list(paste0("s", seq_len(n)),
                                  paste0("f", seq_len(p))))
      studies <- rep(paste0("study", 1:4), each = n_per)
      # Study 4 has a shifted feature scale.
      x[studies == "study4", ] <- x[studies == "study4", ] * 3
      beta <- setNames(c(1.5, -1.2, rep(0, 8)), paste0("f", seq_len(p)))
      eta <- as.vector(x %*% beta)
      shifts <- c(study1 = -1, study2 = -0.3, study3 = 0.3, study4 = 1)
      eta <- eta + unname(shifts[studies])
      eta <- eta - mean(eta)
      pr <- plogis(eta)
      y <- setNames(
        factor(ifelse(stats::rbinom(n, 1L, pr) == 1L, "B", "A"),
               levels = c("A", "B")),
        rownames(x)
      )
      list(x = x, y = y,
           study_labels = setNames(studies, rownames(x)))
    },
    iv5_leaky = {
      beta <- setNames(c(1.5, -1.2, 1.0, rep(0, 7)), paste0("f", 1:10))
      dgp_binary(120L, beta, seed)
    },
    stop("Unknown design: ", design, call. = FALSE)
  )
}

truth_sets <- list(
  iv0_null = character(),
  iv1_sparse = c("f1", "f2", "f3"),
  iv2_substitutes = c("f1a", "f2"),
  iv3_single = "f1",
  iv4_study_shift = c("f1", "f2"),
  iv5_leaky = c("f1", "f2", "f3")
)
true_signs <- list(
  iv1_sparse = c(f1 = 1.5, f2 = -1.2, f3 = 1.0),
  iv5_leaky = c(f1 = 1.5, f2 = -1.2, f3 = 1.0)
)

# ── Arms ─────────────────────────────────────────────────────────────────────

run_exploratory <- function(design, seed, dat) {
  ex <- suppressWarnings(stabl_analyze(
    dat$x, dat$y, family = "binomial", mode = "exploratory",
    lambda_grid = lambda_grid, n_bootstraps = n_boot,
    artificial_type = "random_permutation", random_state = seed,
    n_repeats = 10L
  ))
  ev <- ex$interpretation$evidence
  sel <- ev$feature[ev$selected]
  list(
    result = list(
      n_selected = length(sel),
      selected_features = paste(sel, collapse = "|")
    ),
    effects = data.frame(
      design = design, seed = seed, feature = ev$feature,
      selected = ev$selected, stability = ev$stability,
      coefficient = ev$coefficient, effect_direction = ev$effect_direction,
      correlation_group = ev$correlation_group,
      stringsAsFactors = FALSE
    )
  )
}

run_validated <- function(design, seed, dat) {
  study_labels <- if (identical(design, "iv4_study_shift")) {
    dat$study_labels
  } else {
    NULL
  }
  outer_v <- if (is.null(study_labels)) 3L else 4L
  inner_v <- if (is.null(study_labels)) 3L else 2L
  va <- suppressWarnings(stabl_analyze(
    dat$x, dat$y, family = "binomial", mode = "validated",
    preprocessing = c("center", "scale"),
    lambda_grid = lambda_grid, outer_v = outer_v, inner_v = inner_v,
    n_bootstraps = n_boot, artificial_type = "random_permutation",
    random_state = seed, n_repeats = 10L,
    study_labels = study_labels
  ))
  pooled <- va$performance$pooled
  preds <- va$predictions
  sel_counts <- table(va$selected_features$fold_id)
  fold_df <- do.call(rbind, lapply(split(preds, preds$fold_id), function(fp) {
    n_sel <- if (fp$fold_id[[1L]] %in% names(sel_counts)) {
      as.integer(sel_counts[[fp$fold_id[[1L]]]])
    } else {
      0L
    }
    data.frame(
      design = design, seed = seed, fold_id = fp$fold_id[[1L]],
      n_selected = n_sel,
      prob_degenerate = length(unique(fp$event_probability)) < 2L,
      stringsAsFactors = FALSE
    )
  }))
  imp <- va$heldout_importance$summary
  list(
    result = list(
      n_selected = as.integer(stats::median(fold_df$n_selected)),
      selected_features = paste(
        sort(unique(va$selected_features$feature)), collapse = "|"
      ),
      accuracy = pooled$accuracy, ber = pooled$balanced_error_rate,
      macro_f1 = pooled$macro_f1, mcc = pooled$mcc, brier = pooled$brier,
      log_loss = pooled$log_loss, auroc = pooled$auroc
    ),
    folds = fold_df,
    importance = data.frame(
      design = design, seed = seed, feature = imp$feature,
      mean_importance = imp$mean_importance, sd_importance = imp$sd_importance,
      n_folds_selected = imp$n_folds_selected,
      stringsAsFactors = FALSE
    )
  )
}

# Leaky arm (IV-5): selection on the FULL data, then fold-wise downstream
# prediction on the leaked panel over the SAME outer folds the nested arm
# uses.  Per-fold preprocessing is fitted on fold training data in both
# arms, so the contrast isolates selection leakage.
run_leaky <- function(design, seed, dat) {
  x <- dat$x
  y <- dat$y
  fit_full <- suppressWarnings(stabl_fit(
    x, y, lambda_grid = lambda_grid, family = "binomial",
    n_bootstraps = n_boot, artificial_type = "random_permutation",
    random_state = seed
  ))
  support <- get_support(fit_full)
  sel <- names(support)[support]

  strata_labels <- factor(y)
  names(strata_labels) <- names(y)
  folds <- .make_repeated_cv_folds(
    y = strata_labels, v = 3L, repeats = 1L, stratified = TRUE,
    random_state = seed
  )
  pred_rows <- list()
  for (fold in folds) {
    x_tr <- x[fold$train_ids, , drop = FALSE]
    x_va <- x[fold$valid_ids, , drop = FALSE]
    pp <- .fit_fold_preprocessing(c("center", "scale"), x_tr)
    x_tr <- stabl_preprocess_apply(pp, x_tr)
    x_va <- stabl_preprocess_apply(pp, x_va)
    pred <- .fit_selected_predictor(
      x_train = x_tr[, sel, drop = FALSE],
      y_train = y[fold$train_ids],
      x_valid = x_va[, sel, drop = FALSE],
      levels = levels(y),
      event_class = "B"
    )
    pred_rows[[fold$fold_id]] <- data.frame(
      truth = as.character(y[fold$valid_ids]),
      probability = pred$probability,
      stringsAsFactors = FALSE
    )
  }
  preds <- do.call(rbind, pred_rows)
  m <- .binary_probability_metrics(
    factor(preds$truth, levels = levels(y)), preds$probability, "B"
  )
  list(result = list(
    n_selected = length(sel),
    selected_features = paste(sel, collapse = "|"),
    brier = m$brier, log_loss = m$log_loss, auroc = m$auroc
  ))
}

# ── Run grid ─────────────────────────────────────────────────────────────────

planned <- rbind(
  expand.grid(design = c("iv0_null", "iv1_sparse", "iv2_substitutes",
                         "iv3_single"),
              arm = c("exploratory", "validated"), seed = seeds,
              stringsAsFactors = FALSE),
  expand.grid(design = "iv4_study_shift", arm = "validated", seed = seeds,
              stringsAsFactors = FALSE),
  expand.grid(design = "iv5_leaky", arm = c("leaky", "nested"),
              seed = seeds, stringsAsFactors = FALSE)
)
planned <- planned[order(planned$design, planned$arm, planned$seed), ]
rownames(planned) <- NULL
cat(sprintf("IV benchmark: %d planned runs (%s)\n", nrow(planned), git_sha))

results <- list()
effects_rows <- list()
fold_rows <- list()
importance_rows <- list()
failures <- list()

for (i in seq_len(nrow(planned))) {
  design <- planned$design[[i]]
  arm <- planned$arm[[i]]
  seed <- planned$seed[[i]]
  key <- paste(design, arm, seed, sep = "_")
  t0 <- Sys.time()
  out <- tryCatch(
    {
      dat <- make_design_data(design, seed)
      runner <- switch(
        arm,
        exploratory = run_exploratory,
        validated = run_validated,
        leaky = run_leaky,
        nested = run_validated
      )
      runner(design, seed, dat)
    },
    error = function(e) {
      list(error_class = paste(class(e), collapse = ","),
           error_message = conditionMessage(e))
    }
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (!is.null(out$error_class)) {
    failures[[key]] <- data.frame(
      design = design, arm = arm, seed = seed,
      error_class = out$error_class, error_message = out$error_message,
      stringsAsFactors = FALSE
    )
    cat(sprintf("%s ERROR [%s]: %s\n", key, out$error_class,
                out$error_message))
    next
  }

  res <- out$result
  results[[key]] <- data.frame(
    design = design, arm = arm, seed = seed,
    n_selected = res$n_selected,
    selected_features = res$selected_features,
    accuracy = res$accuracy %||% NA_real_,
    ber = res$ber %||% NA_real_,
    macro_f1 = res$macro_f1 %||% NA_real_,
    mcc = res$mcc %||% NA_real_,
    brier = res$brier %||% NA_real_,
    log_loss = res$log_loss %||% NA_real_,
    auroc = res$auroc %||% NA_real_,
    runtime_sec = runtime,
    git_sha = git_sha, pkg_version = pkg_version,
    stringsAsFactors = FALSE
  )
  if (!is.null(out$effects)) effects_rows[[key]] <- out$effects
  if (!is.null(out$folds)) fold_rows[[key]] <- out$folds
  if (!is.null(out$importance)) importance_rows[[key]] <- out$importance
  cat(sprintf("%s -> n_sel %s auroc %s (%.1fs)\n", key, res$n_selected,
              if (is.null(res$auroc)) "-" else sprintf("%.3f", res$auroc),
              runtime))
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL
failures_df <- if (length(failures)) do.call(rbind, failures) else {
  data.frame(design = character(), arm = character(), seed = integer(),
             error_class = character(), error_message = character(),
             stringsAsFactors = FALSE)
}
effects_df <- if (length(effects_rows)) do.call(rbind, effects_rows) else NULL
folds_df <- if (length(fold_rows)) do.call(rbind, fold_rows) else NULL
importance_df <- if (length(importance_rows)) {
  do.call(rbind, importance_rows)
} else {
  NULL
}

utils::write.csv(results_df, file.path(out_dir, "iv_results.csv"),
                 row.names = FALSE)
utils::write.csv(failures_df, file.path(out_dir, "iv_failures.csv"),
                 row.names = FALSE)
if (!is.null(effects_df)) {
  utils::write.csv(effects_df, file.path(out_dir, "iv_effects.csv"),
                   row.names = FALSE)
}
if (!is.null(folds_df)) {
  utils::write.csv(folds_df, file.path(out_dir, "iv_folds.csv"),
                   row.names = FALSE)
}
if (!is.null(importance_df)) {
  utils::write.csv(importance_df, file.path(out_dir, "iv_importance.csv"),
                   row.names = FALSE)
}

# ── Gate evaluation (criteria fixed in the predesign) ────────────────────────

n_planned <- nrow(planned)
n_completed <- nrow(results_df)
n_failed <- nrow(failures_df)

gates <- list()

# G-IV-2: effect-direction recovery on IV-1 exploratory.
iv1_eff <- effects_df[effects_df$design == "iv1_sparse" &
                        effects_df$selected & !is.na(effects_df$selected), ]
sig <- names(true_signs$iv1_sparse)
iv1_sig <- iv1_eff[iv1_eff$feature %in% sig, ]
dir_match <- mapply(
  function(feat, dir) {
    if (is.na(dir)) return(NA)
    (dir == "positive" && true_signs$iv1_sparse[[feat]] > 0) ||
      (dir == "negative" && true_signs$iv1_sparse[[feat]] < 0)
  },
  iv1_sig$feature, iv1_sig$effect_direction
)
g2_rate <- if (length(dir_match)) mean(dir_match, na.rm = TRUE) else NA_real_
gates[["G-IV-2 effect-direction recovery >= 0.90 (IV-1)"]] <- list(
  value = g2_rate, met = is.finite(g2_rate) && g2_rate >= 0.90
)

# G-IV-3: substitutes share one correlation group in >= 80% of IV-2 runs.
iv2_eff <- effects_df[effects_df$design == "iv2_substitutes", ]
g3_rates <- vapply(split(iv2_eff, iv2_eff$seed), function(df) {
  sub <- df[df$feature %in% c("f1a", "f1b", "f1c"), ]
  grp <- sub$correlation_group
  grp <- grp[!is.na(grp)]
  length(grp) == 3L && length(unique(grp)) == 1L
}, logical(1L))
g3_rate <- mean(g3_rates)
gates[["G-IV-3 substitutes in one group >= 0.80 (IV-2)"]] <- list(
  value = g3_rate, met = g3_rate >= 0.80
)

# G-IV-4: leaky mean AUROC >= nested mean AUROC (IV-5, directional).
iv5 <- results_df[results_df$design == "iv5_leaky", ]
leaky_auroc <- mean(iv5$auroc[iv5$arm == "leaky"], na.rm = TRUE)
nested_auroc <- mean(iv5$auroc[iv5$arm == "nested"], na.rm = TRUE)
gates[["G-IV-4 leaky AUROC >= nested AUROC (IV-5)"]] <- list(
  value = leaky_auroc - nested_auroc,
  met = is.finite(leaky_auroc) && is.finite(nested_auroc) &&
    leaky_auroc >= nested_auroc
)

# G-IV-5: null sanity (IV-0).
iv0 <- results_df[results_df$design == "iv0_null" & results_df$arm == "validated", ]
iv0_median_auroc <- stats::median(iv0$auroc, na.rm = TRUE)
gates[["G-IV-5 null median AUROC in [0.35, 0.65] (IV-0)"]] <- list(
  value = iv0_median_auroc,
  met = is.finite(iv0_median_auroc) && iv0_median_auroc >= 0.35 &&
    iv0_median_auroc <= 0.65 && nrow(iv0) == length(seeds)
)

# G-IV-6: one-feature path non-degenerate (IV-3).
iv3_folds <- folds_df[folds_df$design == "iv3_single", ]
iv3_single <- iv3_folds[iv3_folds$n_selected == 1L, ]
g6_ok <- nrow(iv3_folds) > 0L &&
  (nrow(iv3_single) == 0L || !any(iv3_single$prob_degenerate))
gates[["G-IV-6 one-feature probabilities non-degenerate (IV-3)"]] <- list(
  value = nrow(iv3_single), met = g6_ok
)

# ── Summary ──────────────────────────────────────────────────────────────────

# na.action = na.pass keeps arms that record only a subset of metrics
# (exploratory: selection only; leaky: probability metrics only); metrics
# an arm does not record are shown as NA, not NaN.
agg <- aggregate(
  cbind(n_selected, accuracy, ber, brier, log_loss, auroc) ~ design + arm,
  data = results_df,
  FUN = function(v) {
    m <- mean(v, na.rm = TRUE)
    if (is.nan(m)) NA_real_ else m
  },
  na.action = na.pass
)

lines <- c(
  "# IV — Interpretation & validated-workflow benchmark summary",
  "",
  sprintf("- Predeclared design: interpretation_validation_predesign.md (2026-10-15); %d seeds.", length(seeds)),
  sprintf("- Git SHA: %s; stablr version: %s.", git_sha, pkg_version),
  sprintf("- Run completeness: %d planned, %d completed, %d failed (iv_failures.csv).",
          n_planned, n_completed, n_failed),
  "- Simulation-only evidence; no application-data claim. FDP+ remains a",
  "  learner x decoy x design diagnostic.",
  "",
  "## Mean metrics per design x arm",
  "",
  capture.output(print(agg)),
  "",
  "## Gate criteria",
  ""
)
for (g in names(gates)) {
  lines <- c(lines, sprintf("- %s: value %s -> %s", g,
                            format(gates[[g]]$value, digits = 3L),
                            if (isTRUE(gates[[g]]$met)) "MET" else "NOT MET"))
}
lines <- c(
  lines, "",
  "Notes:",
  "- IV-5 contrasts leaky full-data selection with correctly nested",
  "  selection over identical outer folds; the criterion is directional",
  "  (pooled across seeds), not per-replicate.",
  "- Empirical selection truth uses the known simulated signal sets;",
  "  these are distinct from the fitted FDP+ estimate.",
  ""
)
writeLines(lines, file.path(out_dir, "iv_summary.md"))

cat("\nIV done. Gates:\n")
for (g in names(gates)) {
  cat(sprintf("  %s: %s\n", g, if (isTRUE(gates[[g]]$met)) "MET" else "NOT MET"))
}

# Completeness gate: fail the run if any planned combination is missing.
if (n_failed > 0L || n_completed != n_planned) {
  message(sprintf(
    "IV INCOMPLETE: %d of %d planned runs failed or are missing; see iv_failures.csv.",
    n_planned - n_completed, n_planned
  ))
  quit(save = "no", status = 1L)
}
