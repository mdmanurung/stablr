#!/usr/bin/env Rscript
# ============================================================================
# WS-2 — RF adapter graduation benchmark (predeclared: MATURITY_PREDESIGN.md)
# ============================================================================
# Designs A-D x {knockoff_equi, random_permutation} x 10 seeds, learner "rf".
# Graduation criterion (predeclared): mean empirical FDP <= 0.15 with
# TPR >= 0.8 on designs A-D for rf+knockoff_equi. Permutation rows
# characterize the guard; they do not graduate the adapter.
# Run: Rscript benchmarks/ws2_rf_graduation.R   (~1.5-2.5 h)

.libPaths(c("/workspace/Rlib-shared", .libPaths()))
if (!requireNamespace("ranger", quietly = TRUE)) stop("ranger required")
if (file.exists("DESCRIPTION")) pkgload::load_all(".", compile = TRUE, quiet = TRUE) else library(stablr)

nominal_fdp <- 0.1
seeds <- 1:10
n_bootstraps <- 100L
out_dir <- "benchmarks"

# ---- Knockoff sanity probe (F5 contract, heterogeneous scale) ---------------
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
  if (gram_rel_err < 1e-6 && all(cross_diag >= -1e-6 * max(abs(gram_ref)))) "ok" else "degenerate"
}, error = function(e) paste0("error: ", conditionMessage(e)))
cat("knockoff probe:", probe_status, "\n")

# ---- Designs (predeclared) ---------------------------------------------------
make_design <- function(design, seed) {
  set.seed(seed)
  if (design == "A_linear") {
    n <- 120; p <- 10
    Sigma <- matrix(0.2, p, p); diag(Sigma) <- 1
    for (i in seq_len(p - 3L)) for (j in seq_len(p - 3L)) Sigma[3+i, 3+j] <- 0.6^abs(i-j)
    L <- t(chol(Sigma))
    X <- matrix(rnorm(n*p), n, p) %*% L
    colnames(X) <- paste0("f", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
    y <- 0.8*X[,1] + 0.6*X[,2] - 0.7*X[,3] + rnorm(n, sd = 0.5)
  } else if (design == "B_nonlinear") {
    n <- 120; p <- 10
    Sigma <- matrix(0.2, p, p); diag(Sigma) <- 1
    for (i in seq_len(p - 3L)) for (j in seq_len(p - 3L)) Sigma[3+i, 3+j] <- 0.5^abs(i-j)
    Sigma[1:3, 1:3] <- 0.3; diag(Sigma) <- 1
    L <- t(chol(Sigma))
    X <- matrix(rnorm(n*p), n, p) %*% L
    colnames(X) <- paste0("f", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
    y <- 0.8*X[,1] + 0.8*X[,2] + 1.5*X[,1]*X[,2] + 1.0*(X[,3]^2 - 1) + rnorm(n, sd = 0.5)
  } else if (design == "C_highdim_nonlinear") {
    n <- 200; p <- 50
    Sigma <- matrix(0.15, p, p)  # 0.2 floor + AR(1) block not PD at p=50 (min eig -0.22); 0.15 is PD (+0.33), same structure
    diag(Sigma) <- 1
    for (i in seq_len(p - 3L)) for (j in seq_len(p - 3L)) Sigma[3+i, 3+j] <- 0.5^abs(i-j)
    Sigma[1:3, 1:3] <- 0.3; diag(Sigma) <- 1
    L <- t(chol(Sigma))
    X <- matrix(rnorm(n*p), n, p) %*% L
    colnames(X) <- paste0("f", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
    y <- 0.8*X[,1] + 0.8*X[,2] + 1.5*X[,1]*X[,2] + 1.0*(X[,3]^2 - 1) + rnorm(n, sd = 0.5)
  } else if (design == "D_nonlinear_binomial") {
    n <- 150; p <- 30
    Sigma <- matrix(0.2, p, p); diag(Sigma) <- 1
    for (i in seq_len(p - 3L)) for (j in seq_len(p - 3L)) Sigma[3+i, 3+j] <- 0.5^abs(i-j)
    Sigma[1:3, 1:3] <- 0.3; diag(Sigma) <- 1
    L <- t(chol(Sigma))
    X <- matrix(rnorm(n*p), n, p) %*% L
    colnames(X) <- paste0("f", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
    eta <- 1.2*X[,1]*X[,2] + 0.8*(X[,3]^2 - 1) - 0.6*X[,1]
    y <- rbinom(n, 1, 1/(1 + exp(-eta)))
  } else stop("unknown design")
  list(x = X, y = setNames(y, rownames(X)),
       signals = paste0("f", seq_len(3L)), nulls = paste0("f", 4:p),
       family = if (design == "D_nonlinear_binomial") "binomial" else "gaussian")
}

fdp_tpr <- function(selected, signals, nulls) {
  selected <- as.character(selected); n_sel <- length(selected)
  if (n_sel == 0L) return(data.frame(n_selected = 0L, fdp = 0, tpr = 0))
  data.frame(n_selected = n_sel,
             fdp = sum(selected %in% nulls) / n_sel,
             tpr = sum(selected %in% signals) / length(signals))
}

rows <- list()
for (design in c("A_linear", "B_nonlinear", "C_highdim_nonlinear", "D_nonlinear_binomial")) {
  for (decoy in c("knockoff_equi", "random_permutation")) {
    for (seed in seeds) {
      d <- make_design(design, seed)
      res <- tryCatch(
        withCallingHandlers({
          fit <- stabl_fit(x = d$x, y = d$y, lambda_grid = "auto",
                           base_learner = "rf", family = d$family,
                           artificial_type = decoy,
                           n_bootstraps = n_bootstraps, random_state = seed)
          out <- fdp_tpr(get_feature_names(fit), d$signals, d$nulls)
          out$actual_artificial_type <- if (!is.null(fit$artificial_provenance))
            paste(fit$artificial_provenance$actual_type, collapse = "+") else NA_character_
          out
        },
        warning = function(w) if (grepl("anti-conservative|not validated", conditionMessage(w)))
          invokeRestart("muffleWarning")),
        error = function(e) data.frame(n_selected = NA_integer_, fdp = NA_real_,
                                       tpr = NA_real_, actual_artificial_type = paste0("ERROR: ", conditionMessage(e))))
      rows[[length(rows) + 1L]] <- cbind(data.frame(design = design, seed = seed,
                                                    artificial_type = decoy), res)
      cat(design, decoy, "seed", seed, "-> n_sel", rows[[length(rows)]]$n_selected,
          "fdp", round(rows[[length(rows)]]$fdp, 3), "\n")
    }
  }
}
results <- do.call(rbind, rows)
write.csv(results, file.path(out_dir, "ws2_results.csv"), row.names = FALSE)

agg <- aggregate(cbind(n_selected, fdp, tpr) ~ design + artificial_type, results, mean, na.rm = TRUE)
grad_rows <- agg[agg$artificial_type == "knockoff_equi", ]
graduated <- all(grad_rows$fdp <= 0.15, na.rm = TRUE) && all(grad_rows$tpr >= 0.8, na.rm = TRUE)

sink(file.path(out_dir, "ws2_summary.md"))
cat("# WS-2 — RF adapter graduation summary\n\n- Predeclared design: MATURITY_PREDESIGN.md\n",
    "- Knockoff sanity probe:", probe_status, "\n",
    "- Nominal FDP+ level:", nominal_fdp, "\n\n## Mean empirical FDP / TPR per cell\n\n")
print(agg)
cat("\n## Predeclared graduation criterion (rf + knockoff_equi: mean FDP <= 0.15 AND mean TPR >= 0.8 on ALL designs)\n\n")
cat("Outcome:", if (isTRUE(graduated)) "MET — adapter qualifies for promotion with design scope A-D" else "NOT MET — adapter keeps a design-scope restriction", "\n\n")
cat("## Caveats\n\n- FDP+ is a diagnostic minimizer, not a universal FDR guarantee; rf calibration is purely empirical (plan D1).\n",
    "- Permutation rows are expected anti-conservative (permuted nulls lack correlation leakage); they characterize the guard, not graduation.\n",
    "- Simulation-only evidence; no application-data claim.\n")
sink()
cat("WS-2 done. Graduated:", isTRUE(graduated), "\n")
