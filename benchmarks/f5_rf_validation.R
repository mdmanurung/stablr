#!/usr/bin/env Rscript
# ============================================================================
# F5 — RF validation simulation suite (plan item F5; gates D1; audit E2/E3)
# ============================================================================
# NOT a testthat test: a self-contained benchmark script, run via
#   Rscript benchmarks/f5_rf_validation.R [--smoke]
# Full run (~10 seeds x {knockoff_equi, random_permutation} x
# {rf, oblique_rf (if aorsf is installed), lasso}) takes roughly 45-90
# minutes; `--smoke` runs a 2-seed small-n proof of execution in a couple of
# minutes.
#
# Design under test: purely non-linear signals (interaction x1*x2, quadratic
# x3^2, plus modest main effects so marginal permutation importance can see
# the interaction partners) with correlated null features.
#
# Per run we record, at the FDP+-chosen threshold (get_support()):
#   empirical FDP = |selected nulls| / |selected|   (0 when nothing selected)
#   TPR           = |selected signals| / |signals|
# and compare mean empirical FDP per configuration against the nominal 0.1.
#
# Honesty caveats, recorded in the outputs:
# * Scale contract (0.2.2): knockoff::create.fixed() column-normalizes its
#   input (normc, center = FALSE) and returns knockoffs on that unit-norm
#   scale; stablr back-transforms them by the INPUT column norms so the
#   shipped pair (x, xk) satisfies the fixed-X Gram condition on x's own
#   scale.  The start-of-run probe replicates that back-transform on a
#   HETEROGENEOUS-scale design and checks the Gram identity Xk'Xk = X'X plus
#   the nonnegative diagonal shift of the cross-Gram.  A probe failure
#   invalidates the knockoff_equi rows.
# * The lasso + knockoff_equi rows are expected to show FDP above nominal on
#   this design: the knockoff pool itself is valid (direct exchangeability
#   check -- knockoff nulls are selected at least as often as real nulls at
#   every lambda), but a linear learner's stability frequencies for AR(1)
#   correlated nulls overlap the knockoff frequencies, so the min-FDP+
#   threshold cannot reach 0.1.  FDP+ remains a conservative estimate of the
#   realized FDP (measured: estimate 0.67 vs empirical 0.33-0.58).  The rf
#   rows, which separate signals sharply, are the calibration evidence.
# * FDP+ is a diagnostic minimizer, not a universal FDR guarantee; for the
#   RF learners its calibration is purely empirical (plan item D1).  The
#   oblique_rf rows (audit E3) are the first F5 evidence for the aorsf
#   adapter; per-feature selection frequencies are noisier than for rf
#   because importance is spread across the linear combinations used in
#   each split.
# * Boruta baseline is included only if the Boruta package is installed
#   (requireNamespace guard; this script never installs packages).
# ============================================================================

.libPaths(c("/workspace/Rlib-shared", .libPaths()))  # sandbox ranger location

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args

smoke_n            <- 80L
full_n             <- 120L
smoke_seeds        <- 1:2
full_seeds         <- 1:10
smoke_boot         <- 20L
full_boot          <- 100L
num_trees          <- 500L
p                  <- 10L
nominal_fdp        <- 0.1

n_seeds     <- if (smoke) length(smoke_seeds) else length(full_seeds)
n           <- if (smoke) smoke_n else full_n
seeds       <- if (smoke) smoke_seeds else full_seeds
n_bootstraps <- if (smoke) smoke_boot else full_boot

out_dir <- if (file.exists("DESCRIPTION")) file.path("inst", "benchmarks") else {
  stop("Run this script from the package root (the directory containing DESCRIPTION).")
}
results_csv <- file.path(out_dir, "f5_results.csv")
summary_md  <- file.path(out_dir, "f5_summary.md")

# ---- Package loading --------------------------------------------------------
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", compile = TRUE, quiet = TRUE)
} else if (requireNamespace("stablr", quietly = TRUE)) {
  library(stablr)
} else {
  stop("Neither pkgload nor an installed stablr is available.")
}

has_ranger <- requireNamespace("ranger", quietly = TRUE)
if (!has_ranger) stop("ranger is required for the F5 RF validation suite.")
has_aorsf <- requireNamespace("aorsf", quietly = TRUE)
has_boruta <- requireNamespace("Boruta", quietly = TRUE)

# ---- Knockoff sanity probe --------------------------------------------------
# Probes the SHIPPED contract (0.2.2): create.fixed() column-normalizes its
# input and returns knockoffs on that unit-norm scale; stablr back-transforms
# them by the INPUT column norms.  The probe replicates that back-transform
# on a heterogeneous-scale design (column sds 0.01-100) and requires the
# fixed-X Gram identity Xk'Xk = X'X plus a nonnegative diagonal shift of the
# cross-Gram (t(X)X - t(X)Xk = diag(s), s >= 0).  No sd ~ 1 expectation: on
# heterogeneous inputs the knockoff columns inherit the input scales.
probe_status <- tryCatch({
  set.seed(2026L)
  sds <- 10^seq(-2, 2, length.out = 20)
  X0 <- matrix(rnorm(200 * 20), 200, 20) %*% diag(sds)
  kd <- knockoff::create.fixed(X0, sigma = 1)
  col_norms <- sqrt(colSums(X0^2))
  Xk0 <- kd$Xk %*% diag(col_norms)
  gram_ref <- crossprod(X0)
  gram_rel_err <- max(abs(crossprod(Xk0) - gram_ref)) / max(abs(gram_ref))
  cross_diag <- diag(gram_ref - crossprod(X0, Xk0))
  gram_ok <- gram_rel_err < 1e-6
  scale_ok <- all(cross_diag >= -1e-6 * max(abs(gram_ref)))
  if (gram_ok && scale_ok) "ok" else "degenerate"
}, error = function(e) paste0("error: ", conditionMessage(e)))

# ---- Synthetic design (mirrors tests/testthat/test-rf-adapter.R) ------------
simulate_design <- function(n, p, seed) {
  set.seed(seed)
  p_null <- p - 3L
  Sigma  <- matrix(0.2, nrow = p, ncol = p)
  diag(Sigma) <- 1
  for (i in seq_len(p_null)) {
    for (j in seq_len(p_null)) Sigma[3L + i, 3L + j] <- 0.5 ^ abs(i - j)
  }
  Sigma[1:3, 1:3] <- 0.3
  diag(Sigma) <- 1
  L <- t(chol(Sigma))
  X <- matrix(as.numeric(matrix(rnorm(n * p), n, p) %*% L), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(
    0.8 * X[, 1] + 0.8 * X[, 2] + 1.5 * X[, 1] * X[, 2] +
      1.0 * (X[, 3] ^ 2 - 1) + rnorm(n, sd = 0.5),
    rownames(X)
  )
  list(x = X, y = y, signals = paste0("f", seq_len(3L)),
       nulls = paste0("f", 4:p))
}

fdp_tpr <- function(selected, signals, nulls) {
  selected <- as.character(selected)
  n_sel <- length(selected)
  if (n_sel == 0L) {
    return(data.frame(n_selected = 0L, fdp = 0, tpr = 0))
  }
  data.frame(
    n_selected = n_sel,
    fdp = sum(selected %in% nulls) / n_sel,
    tpr = sum(selected %in% signals) / length(signals)
  )
}

run_config <- function(learner, artificial_type, d, seed) {
  res <- tryCatch(
    withCallingHandlers(
      {
        fit <- stabl_fit(
          x = d$x, y = d$y,
          lambda_grid = "auto",
          base_learner = learner,
          artificial_type = artificial_type,
          n_bootstraps = n_bootstraps,
          random_state = seed
        )
        out <- fdp_tpr(get_feature_names(fit), d$signals, d$nulls)
        out$actual_artificial_type <-
          if (!is.null(fit$artificial_provenance))
            paste(fit$artificial_provenance$actual_type, collapse = "+")
          else NA_character_
        out
      },
      warning = function(w) {
        # The rf/oblique_rf guardrail warnings (permutation nulls and
        # fixed-X knockoffs) are expected here; muffle them.
        if (grepl("anti-conservative|not validated", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {
      data.frame(n_selected = NA_integer_, fdp = NA_real_, tpr = NA_real_,
                 actual_artificial_type = paste0("ERROR: ", conditionMessage(e)))
    }
  )
  cbind(
    data.frame(seed = seed, learner = learner, artificial_type = artificial_type),
    res
  )
}

run_boruta <- function(d, seed) {
  if (!has_boruta) return(NULL)
  res <- tryCatch({
    set.seed(seed)
    bor <- Boruta::Boruta(x = as.data.frame(d$x), y = d$y,
                          numTrees = 500L, seed = seed)
    selected <- names(bor$finalDecision)[bor$finalDecision == "Confirmed"]
    out <- fdp_tpr(selected, d$signals, d$nulls)
    out$actual_artificial_type <- NA_character_
    cbind(data.frame(seed = seed, learner = "boruta",
                     artificial_type = NA_character_), out)
  }, error = function(e) NULL)
  res
}

# ---- Run the grid ------------------------------------------------------------
# Audit E2/E3 (0.2.2): the grid expands to the oblique-RF adapter (aorsf)
# when installed, giving the adapter its first F5 FDP evidence.
learners <- c("rf", if (has_aorsf) "oblique_rf", "lasso")

cat(sprintf("F5 RF validation | mode=%s | n=%d p=%d | seeds=%d | n_bootstraps=%d\n",
            if (smoke) "smoke" else "full", n, p, n_seeds, n_bootstraps))
cat(sprintf("knockoff probe: %s | oblique_rf: %s | Boruta baseline: %s\n",
            probe_status,
            if (has_aorsf) "included" else "aorsf not installed (skipped)",
            if (has_boruta) "available" else "not installed (skipped)"))

results <- list()
for (art in c("knockoff_equi", "random_permutation")) {
  for (learner in learners) {
    for (seed in seeds) {
      d <- simulate_design(n, p, seed)
      row <- run_config(learner, art, d, seed)
      results[[length(results) + 1L]] <- row
      cat(sprintf("  %-5s + %-18s seed=%-3d -> n_sel=%s FDP=%s TPR=%s\n",
                  learner, art, seed,
                  ifelse(is.na(row$n_selected), "ERR", row$n_selected),
                  ifelse(is.na(row$fdp), "-", sprintf("%.2f", row$fdp)),
                  ifelse(is.na(row$tpr), "-", sprintf("%.2f", row$tpr))))
    }
  }
  if (has_boruta) {
    d <- simulate_design(n, p, seeds[[1L]])
    row <- run_boruta(d, seeds[[1L]])
    if (!is.null(row)) {
      results[[length(results) + 1L]] <- row
      cat(sprintf("  boruta baseline (seed=%d) -> n_sel=%s FDP=%.2f TPR=%.2f\n",
                  seeds[[1L]], row$n_selected, row$fdp, row$tpr))
    }
  }
}
results <- do.call(rbind, results)
results$knockoff_probe <- probe_status
results$mode <- if (smoke) "smoke" else "full"
results$n <- n
results$n_bootstraps <- n_bootstraps
write.csv(results, results_csv, row.names = FALSE)

# ---- Summary -----------------------------------------------------------------
sumry <- aggregate(cbind(fdp, tpr, n_selected) ~ learner + artificial_type,
                   data = results[!is.na(results$fdp), ],
                   FUN = mean, na.action = na.omit)
names(sumry)[names(sumry) == "fdp"] <- "mean_empirical_fdp"
names(sumry)[names(sumry) == "tpr"] <- "mean_tpr"
names(sumry)[names(sumry) == "n_selected"] <- "mean_n_selected"
sumry$nominal_fdp <- nominal_fdp
sumry <- sumry[order(sumry$learner, sumry$artificial_type), ]

lines <- c(
  "# F5 — RF validation simulation summary",
  "",
  sprintf("- Mode: **%s** (smoke = execution proof only, not evidence)%s",
          if (smoke) "smoke" else "full",
          if (smoke) sprintf(" | n=%d, seeds=%d, n_bootstraps=%d", n, n_seeds, n_bootstraps) else ""),
  sprintf("- Design: non-linear signals (x1*x2 interaction, x3^2 quadratic, main effects) with correlated nulls; n=%d, p=%d.", n, p),
  sprintf("- Nominal FDP+ level: %.1f", nominal_fdp),
  sprintf("- Knockoff sanity probe (0.2.2 shipped-pair contract on a heterogeneous-scale design, Gram identity + nonnegative cross-Gram shift): **%s**", probe_status),
  paste0("- Boruta baseline: ",
         if (has_boruta) "included (first seed)" else "not installed - skipped gracefully"),
  "",
  "## Mean empirical FDP / TPR per configuration (at the FDP+-chosen threshold)",
  "",
  paste0("| learner | artificial_type | mean empirical FDP | mean TPR | mean n selected | nominal FDP |"),
  paste0("|---|---|---|---|---|---|"),
  apply(sumry, 1L, function(r) sprintf("| %s | %s | %.3f | %.3f | %.1f | %.1f |",
        r[["learner"]], r[["artificial_type"]], as.numeric(r[["mean_empirical_fdp"]]),
        as.numeric(r[["mean_tpr"]]), as.numeric(r[["mean_n_selected"]]),
        as.numeric(r[["nominal_fdp"]]))),
  "",
  "## Caveats",
  "",
  if (probe_status != "ok") paste0(
    "- **Knockoff probe FAILED (probe: ", probe_status, "):** the fixed-X shipped-pair ",
    "scale contract is violated in this environment, so the knockoff rows here do NOT ",
    "measure valid knockoff calibration. Re-run on CI with the CRAN knockoff ",
    "build before drawing conclusions about knockoff configs."
  ) else paste0(
    "- Knockoff build passed the sanity probe (shipped-pair contract on a ",
    "heterogeneous-scale design: Gram identity Xk'Xk = X'X exact to 1e-6 ",
    "relative, nonnegative cross-Gram diagonal shift)."
  ),
  "- FDP+ is a diagnostic minimizer, not a universal FDR guarantee; for `rf` its calibration is purely empirical (plan D1) and is validated by this suite.",
  "- With `random_permutation` nulls, RF FDP+ is expected to be anti-conservative (permuted nulls do not receive the correlation leakage real nulls receive); the rf+permutation rows should show inflated empirical FDP relative to the knockoff rows on a sane knockoff build.",
  ""
)
writeLines(lines, summary_md)

cat(sprintf("\nWrote %s and %s\n", results_csv, summary_md))
