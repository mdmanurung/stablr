# WS-9 — Ordinal learner graduation benchmark.
# Predeclared design: MATURITY_PREDESIGN.md, section WS-9 (written before
# this run).  Proportional-odds simulations; graduation = ordinal +
# knockoff_equi mean FDP <= 0.15 AND mean TPR >= 0.8 on BOTH designs.

suppressPackageStartupMessages(library(stablr))

out_dir <- Sys.getenv("WS9_OUT_DIR", "inst/benchmarks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_seeds <- 10L
n_boot  <- 100L

make_design <- function(name, n, p, rho, beta) {
  list(name = name, n = n, p = p, rho = rho, beta = beta,
       signals = paste0("f", seq_along(beta)))
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
  fdp <- if (n_sel > 0) length(false_pos) / n_sel else 0
  tpr <- length(intersect(sel, dsg$signals)) / length(dsg$signals)
  list(design = dsg$name, artificial_type = artificial_type, seed = seed,
       n_selected = n_sel, fdp = fdp, tpr = tpr)
}

rows <- list()
for (dsg in designs) {
  for (art in c("knockoff_equi", "random_permutation")) {
    for (sd_ in seq_len(n_seeds)) {
      r <- try(simulate_one(dsg, sd_, art), silent = TRUE)
      if (inherits(r, "try-error")) {
        cat(sprintf("%s %s seed %d ERROR: %s\n", dsg$name, art, sd_,
                    conditionMessage(attr(r, "condition"))))
        next
      }
      rows[[length(rows) + 1]] <- r
      cat(sprintf("%s %s seed %d -> n_sel %d fdp %.3f tpr %.3f\n",
                  dsg$name, art, sd_, r$n_selected, r$fdp, r$tpr))
    }
  }
}

res <- do.call(rbind, lapply(rows, as.data.frame))
write.csv(res, file.path(out_dir, "ws9_results.csv"), row.names = FALSE)

agg <- aggregate(cbind(n_selected, fdp, tpr) ~ design + artificial_type,
                 data = res, FUN = mean)
cat("\n== WS-9 means ==\n")
print(agg)

crit <- agg[agg$artificial_type == "knockoff_equi", ]
graduated <- nrow(crit) == 2L && all(crit$fdp <= 0.15) && all(crit$tpr >= 0.8)
cat(sprintf("\nWS-9 done. Graduated: %s\n", graduated))

lines <- c(
  "# WS-9 — Ordinal learner graduation summary",
  "",
  sprintf("- Predeclared design: MATURITY_PREDESIGN.md (WS-9); %d seeds x 2 designs x 2 decoy types x %d bootstraps.",
          n_seeds, n_boot),
  sprintf("- Nominal FDP+ level: 0.1; graduation: knockoff_equi mean FDP <= 0.15 AND TPR >= 0.8 on BOTH designs."),
  "",
  "## Mean empirical FDP / TPR per cell",
  "",
  capture.output(print(agg)),
  "",
  sprintf("Outcome: %s", if (graduated)
    "GRADUATED — ordinal promoted with knockoff_equi scope"
  else "NOT MET — ordinal stays experimental with documented scope"),
  "",
  "- Permutation rows are expected anti-conservative (characterization only).",
  "- Simulation-only evidence; no application-data claim.",
  ""
)
writeLines(lines, file.path(out_dir, "ws9_summary.md"))
