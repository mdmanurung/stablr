#!/usr/bin/env Rscript
# ============================================================================
# WS-3 — Oblique RF adapter: diagnosis + predeclared remedies (fix-or-fence)
# ============================================================================
# Part 1 (diagnosis): decoy selection-frequency noise proxy, rf vs oblique_rf
#   on the F5 non-linear design B: mean/max decoy importance across seeds.
# Part 2 (remedies, benchmark-level only, no runtime edits):
#   baseline | R1 column-standardize x | R2 n_bootstraps 200 | R3 = R1 + R2
#   on design B and a heterogeneous-scale variant (column SDs log-uniform
#   [0.1, 10]), oblique_rf + knockoff_equi, 10 seeds.
# Success (predeclared): mean empirical FDP <= 0.15 on BOTH designs.
# Run: Rscript benchmarks/ws3_oblique_remedy.R   (~1.5-2.5 h)

.libPaths(c("/workspace/Rlib-shared", .libPaths()))
if (!requireNamespace("aorsf", quietly = TRUE)) stop("aorsf required")
if (file.exists("DESCRIPTION")) pkgload::load_all(".", compile = TRUE, quiet = TRUE) else library(stablr)

seeds <- 1:10
out_dir <- "benchmarks"

make_design <- function(variant, seed) {
  set.seed(seed)
  n <- 120; p <- 10
  Sigma <- matrix(0.2, p, p); diag(Sigma) <- 1
  for (i in seq_len(p - 3L)) for (j in seq_len(p - 3L)) Sigma[3+i, 3+j] <- 0.5^abs(i-j)
  Sigma[1:3, 1:3] <- 0.3; diag(Sigma) <- 1
  L <- t(chol(Sigma))
  X <- matrix(rnorm(n*p), n, p) %*% L
  if (variant == "B_het") {
    sds <- 10^runif(p, -1, 1)   # log-uniform [0.1, 10]
    X <- X %*% diag(sds)
  }
  colnames(X) <- paste0("f", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
  y <- 0.8*X[,1] + 0.8*X[,2] + 1.5*X[,1]*X[,2] + 1.0*(X[,3]^2 - 1) + rnorm(n, sd = 0.5)
  list(x = X, y = setNames(y, rownames(X)),
       signals = paste0("f", seq_len(3L)), nulls = paste0("f", 4:p))
}

fdp_tpr <- function(selected, signals, nulls) {
  selected <- as.character(selected); n_sel <- length(selected)
  if (n_sel == 0L) return(data.frame(n_selected = 0L, fdp = 0, tpr = 0))
  data.frame(n_selected = n_sel,
             fdp = sum(selected %in% nulls) / n_sel,
             tpr = sum(selected %in% signals) / length(signals))
}

fit_one <- function(d, learner, decoy, boot, standardize, seed) {
  x <- d$x
  if (standardize) x <- scale(x)
  res <- tryCatch(
    withCallingHandlers({
      fit <- stabl_fit(x = x, y = d$y, lambda_grid = "auto",
                       base_learner = learner, artificial_type = decoy,
                       n_bootstraps = boot, random_state = seed)
      imp <- get_importances(fit)
      sel <- get_feature_names(fit)
      out <- fdp_tpr(sel, d$signals, d$nulls)
      out$decoy_mean_imp <- mean(imp[d$nulls])
      out$decoy_max_imp  <- max(imp[d$nulls])
      out
    },
    warning = function(w) if (grepl("anti-conservative|not validated", conditionMessage(w)))
      invokeRestart("muffleWarning")),
    error = function(e) data.frame(n_selected = NA_integer_, fdp = NA_real_, tpr = NA_real_,
                                   decoy_mean_imp = NA_real_, decoy_max_imp = NA_real_))
  cbind(res, data.frame(seed = seed))
}

# ---- Part 1: diagnosis (decoy-frequency noise, rf vs oblique_rf, design B) --
diag_rows <- list()
for (seed in seeds) {
  d <- make_design("B", seed)
  for (learner in c("rf", "oblique_rf")) {
    r <- fit_one(d, learner, "knockoff_equi", 100L, FALSE, seed)
    diag_rows[[length(diag_rows) + 1L]] <-
      cbind(data.frame(learner = learner), r)
  }
}
diag <- do.call(rbind, diag_rows)
write.csv(diag, file.path(out_dir, "ws3_diagnosis.csv"), row.names = FALSE)

# ---- Part 2: remedies --------------------------------------------------------
remedy_rows <- list()
for (variant in c("B", "B_het")) {
  for (remedy in c("baseline", "R1_standardize", "R2_boot200", "R3_both")) {
    boot <- if (remedy %in% c("R2_boot200", "R3_both")) 200L else 100L
    std  <- remedy %in% c("R1_standardize", "R3_both")
    for (seed in seeds) {
      d <- make_design(variant, seed)
      r <- fit_one(d, "oblique_rf", "knockoff_equi", boot, std, seed)
      remedy_rows[[length(remedy_rows) + 1L]] <-
        cbind(data.frame(variant = variant, remedy = remedy), r)
    }
  }
}
remedies <- do.call(rbind, remedy_rows)
write.csv(remedies, file.path(out_dir, "ws3_results.csv"), row.names = FALSE)

agg <- aggregate(cbind(n_selected, fdp, tpr, decoy_mean_imp, decoy_max_imp) ~ variant + remedy,
                 remedies, mean, na.rm = TRUE)
ok <- agg[agg$variant == "B" & agg$remedy %in% c("baseline","R1_standardize","R2_boot200","R3_both"), ]
ok_het <- agg[agg$variant == "B_het", ]
pass_B <- ok$fdp[ok$remedy != "baseline"] <= 0.15
pass_het <- ok_het$fdp[ok_het$remedy != "baseline"] <= 0.15
any_pass <- any(pass_B & pass_het, na.rm = TRUE)

sink(file.path(out_dir, "ws3_summary.md"))
cat("# WS-3 — Oblique RF fix-or-fence summary\n\n- Predeclared: MATURITY_PREDESIGN.md; nominal FDP+ 0.1; success criterion mean FDP <= 0.15 on BOTH B and B_het.\n\n")
cat("## Diagnosis: decoy selection-frequency noise proxy (design B, knockoff_equi, 10 seeds)\n\n")
print(aggregate(cbind(decoy_mean_imp, decoy_max_imp, fdp) ~ learner, diag, mean, na.rm = TRUE))
cat("\n## Remedies (oblique_rf + knockoff_equi)\n\n")
print(agg)
cat("\n## Outcome\n\n")
cat(if (any_pass) "A predeclared remedy met the criterion on both designs — candidate for package-level adoption (separate commit citing this file)." else "NO predeclared remedy met the criterion — adapter is FENCED: not-recommended scope goes into man pages/pkgdown; adverse result retained.", "\n\n")
cat("## Caveats\n\n- Remedies tested at benchmark level (input-level transforms), not runtime changes.\n",
    "- FDP+ is a diagnostic minimizer, not a universal FDR guarantee. Simulation-only.\n")
sink()
cat("WS-3 done. Any remedy passed:", any_pass, "\n")
