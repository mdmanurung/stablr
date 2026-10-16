# WS-9 — Ordinal learner graduation benchmark (corrected).
# Predeclared design: MATURITY_PREDESIGN.md, section WS-9 (written before the
# original run).  Proportional-odds-style simulations; graduation = ordinal +
# knockoff_equi mean FDP <= 0.15 AND mean TPR >= 0.8 on BOTH designs.
#
# Corrections relative to the historical script (see the dated erratum in
# MATURITY_PREDESIGN.md):
#   1. Truth set: `signals` is now the set of NONZERO coefficients of beta
#      (3 of 10 predictors), not `paste0("f", seq_along(beta))` (all 10).
#      The historical truth set made empirical FDP/TPR invalid.
#   2. Run integrity: every design x decoy x seed combination is attempted and
#      accounted for.  Failures are recorded with seed, error class, and
#      message instead of being silently dropped by `try(...); next`.
#   3. Completeness gate: the script fails (nonzero exit) if any planned run
#      is missing from the results, so a partial CSV can never be summarised
#      as complete.
#   4. Empirical FDP is kept distinct from the fitted FDP+ estimate; the
#      fitted FDP+ threshold is recorded per run for inspection.
#
# The historical DGP (Gaussian latent noise + data-derived cut() thresholds)
# is preserved unchanged so this remains a repair of the original experiment,
# not a new preregistered design.  A genuinely cumulative-logit design is a
# separate, separately predeclared benchmark (ws9_revalidation.R).

# Run against THIS checkout, not an installed copy (the ordinal routing and
# current RF adapter live here).
pkgload::load_all(".", compile = FALSE, quiet = TRUE)

out_dir <- Sys.getenv("WS9_OUT_DIR", "benchmarks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_seeds <- 10L
n_boot  <- 100L

make_design <- function(name, n, p, rho, beta) {
  list(name = name, n = n, p = p, rho = rho, beta = beta,
       # Truth set: nonzero coefficients only.
       signals = paste0("f", which(beta != 0)))
}

designs <- list(
  make_design("O1_linear_indep", 150, 10, 0,   c(1.2, 1.0, 0.8, rep(0, 7))),
  make_design("O2_equicorr03",   150, 10, 0.3, c(1.2, 1.0, 0.8, rep(0, 7)))
)

simulate_one <- function(dsg, seed, artificial_type) {
  set.seed(seed)
  if (dsg$rho > 0) {
    Sigma <- matrix(dsg$rho, dsg$p, dsg$p); diag(Sigma) <- 1
    X <- matrix(MASS::mvrnorm(dsg$n, rep(0, dsg$p), Sigma), dsg$n, dsg$p)
  } else {
    X <- matrix(rnorm(dsg$n * dsg$p), dsg$n, dsg$p)
  }
  colnames(X) <- paste0("f", seq_len(dsg$p))
  rownames(X) <- paste0("s", seq_len(dsg$n))
  eta  <- as.vector(X %*% dsg$beta) + rnorm(dsg$n, sd = 0.6)
  y    <- cut(eta, breaks = 4, labels = paste0("L", 1:4),
              ordered_result = TRUE)
  y    <- setNames(y, rownames(X))
  fit <- suppressWarnings(stabl_fit(
    x               = X,
    y               = y,
    lambda_grid     = "auto",
    base_learner    = "ordinal",
    family          = "ordinal",
    n_bootstraps    = n_boot,
    artificial_type = artificial_type,
    random_state    = seed
  ))
  sel   <- get_feature_names(fit)
  n_sel <- length(sel)
  false_pos <- setdiff(sel, dsg$signals)
  # Empirical FDP over the selected set; 0 when nothing is selected.
  fdp <- if (n_sel > 0) length(false_pos) / n_sel else 0
  tpr <- length(intersect(sel, dsg$signals)) / length(dsg$signals)
  list(design = dsg$name, artificial_type = artificial_type, seed = seed,
       n_selected = n_sel, fdp = fdp, tpr = tpr,
       selected_features = paste(sel, collapse = "|"),
       fdp_plus_threshold = fit$fdr_min_threshold_ %||% NA_real_)
}

# Full planned grid: every design x decoy x seed combination.
planned <- expand.grid(
  design = vapply(designs, `[[`, character(1L), "name"),
  artificial_type = c("knockoff_equi", "random_permutation"),
  seed = seq_len(n_seeds),
  stringsAsFactors = FALSE
)
n_planned <- nrow(planned)

rows <- list()
failures <- list()
for (dsg in designs) {
  for (art in c("knockoff_equi", "random_permutation")) {
    for (sd_ in seq_len(n_seeds)) {
      r <- tryCatch(
        simulate_one(dsg, sd_, art),
        error = function(e) e
      )
      if (inherits(r, "error")) {
        failures[[length(failures) + 1]] <- list(
          design = dsg$name, artificial_type = art, seed = sd_,
          error_class = paste(class(r), collapse = ","),
          error_message = conditionMessage(r)
        )
        cat(sprintf("%s %s seed %d ERROR [%s]: %s\n", dsg$name, art, sd_,
                    class(r)[1L], conditionMessage(r)))
        next
      }
      rows[[length(rows) + 1]] <- r
      cat(sprintf("%s %s seed %d -> n_sel %d fdp %.3f tpr %.3f\n",
                  dsg$name, art, sd_, r$n_selected, r$fdp, r$tpr))
    }
  }
}

res <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
write.csv(res, file.path(out_dir, "ws9_results.csv"), row.names = FALSE)

if (length(failures)) {
  fail_df <- do.call(rbind, lapply(failures, function(f) as.data.frame(f, stringsAsFactors = FALSE)))
  write.csv(fail_df, file.path(out_dir, "ws9_failures.csv"), row.names = FALSE)
} else {
  fail_df <- data.frame(
    design = character(), artificial_type = character(), seed = integer(),
    error_class = character(), error_message = character(),
    stringsAsFactors = FALSE
  )
  write.csv(fail_df, file.path(out_dir, "ws9_failures.csv"), row.names = FALSE)
}

# Completeness gate: planned vs completed.
n_completed <- nrow(res)
n_failed <- length(failures)
cat(sprintf("\nPlanned runs: %d | completed: %d | failed: %d\n",
            n_planned, n_completed, n_failed))
if (n_completed + n_failed != n_planned) {
  stop("WS-9 accounting error: completed + failed != planned.", call. = FALSE)
}

agg <- aggregate(cbind(n_selected, fdp, tpr) ~ design + artificial_type,
                 data = res, FUN = mean)
cat("\n== WS-9 means (corrected truth set) ==\n")
print(agg)

crit <- agg[agg$artificial_type == "knockoff_equi", ]
graduated <- nrow(crit) == 2L && all(crit$fdp <= 0.15) && all(crit$tpr >= 0.8)
cat(sprintf("\nWS-9 done. Graduated: %s\n", graduated))

lines <- c(
  "# WS-9 — Ordinal learner graduation summary (corrected)",
  "",
  sprintf("- Predeclared design: MATURITY_PREDESIGN.md (WS-9); %d seeds x 2 designs x 2 decoy types x %d bootstraps.",
          n_seeds, n_boot),
  "- Nominal FDP+ level: 0.1; graduation: knockoff_equi mean FDP <= 0.15 AND TPR >= 0.8 on BOTH designs.",
  "- Truth set: nonzero coefficients of beta (3 of 10 predictors). The",
  "  historical script used all 10 indices; see the dated erratum in",
  "  MATURITY_PREDESIGN.md.",
  sprintf("- Run completeness: %d planned, %d completed, %d failed (ws9_failures.csv).",
          n_planned, n_completed, n_failed),
  "",
  "## Mean empirical FDP / TPR per cell",
  "",
  capture.output(print(agg)),
  "",
  sprintf("Outcome: %s", if (graduated)
    "GRADUATED — ordinal promoted with knockoff_equi scope"
  else "NOT MET — ordinal stays experimental with documented scope"),
  "",
  "- Empirical FDP/TPR are over the corrected 3-predictor truth set; they are",
  "  distinct from the fitted FDP+ estimate.",
  "- Permutation rows are expected anti-conservative (characterization only).",
  "- Simulation-only evidence; no application-data claim.",
  ""
)
writeLines(lines, file.path(out_dir, "ws9_summary.md"))

# Fail the run if any planned combination is missing from the results.
if (n_failed > 0L) {
  message(sprintf(
    "WS-9 INCOMPLETE: %d of %d planned runs failed; see ws9_failures.csv. Graduation inference is not valid on a partial grid.",
    n_failed, n_planned
  ))
  quit(save = "no", status = 1L)
}
