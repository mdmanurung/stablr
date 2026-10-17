# Phase 7 mechanism probe: why does the min-agg native panel saturate under the null?
#
# Attribution target (findings.md lead, Phase 7): S0 null selects 59/59 features
# in 50/50 permutations. Hypothesis (code-level):
#   (1) auto_rf_threshold_grid floors the threshold path at the MEDIAN of the
#       full-data permutation-importance distribution (probs = seq(0.5, 0.999)).
#   (2) The per-cohort stability score is MAX-over-lambda of the per-bootstrap
#       selection frequency, so the densest (median-bar) point defines the score.
#   (3) Under a null, refit-time importances are noise whose scale exceeds the
#       fixed full-data median bar, so nearly every feature -- real AND decoy --
#       passes the densest threshold in nearly every bootstrap -> score ~ 1.
#   (4) The native panel rule (min_score >= 0.5) never consults the FDP+
#       machinery, so nothing corrects this leniency.
# Probe: pure-null stabl_fit with rf + knockoff_equi, auto grid vs a sparse-only
# grid (quantiles 0.9-0.999). Prediction: auto -> max-over-lambda scores ~ 1 for
# reals and decoys; sparse-only -> scores drop toward 0.
library(pkgload)
load_all("/workspace/integrate/stablr", quiet = TRUE)

set.seed(20261013)
n <- 80
p <- 20
x <- matrix(rnorm(n * p), nrow = n,
            dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
y <- setNames(rnorm(n), rownames(x))          # pure null, no signal
yp <- y[sample.int(n)]                        # within-sample permutation (S0 analog)

cat("== grids ==\n")
imp <- .rf_fit_permutation_importance(x, yp, 500L, 1L)
grid_auto <- auto_rf_threshold_grid(x, yp, n_grid = 6L, num_trees = 500L, seed = 1L)
probs_sparse <- seq(0.9, 0.999, length.out = 6L)
grid_sparse <- data.frame(
  lambda = sort(as.numeric(quantile(imp, probs = probs_sparse, names = FALSE)),
                decreasing = TRUE)
)
print(round(grid_auto$threshold, 4))
cat("auto grid: thresholds are quantiles", paste(round(seq(0.5, 0.999, length.out = 6), 3), collapse = ", "), "of full-data importance\n")
cat("median importance (densest auto bar):", round(min(grid_auto$threshold), 4),
    " max importance:", round(max(imp), 4), "\n")
print(round(grid_sparse$lambda, 4))

fit_auto <- stabl_fit(
  x = x, y = yp, lambda_grid = grid_auto, family = "gaussian",
  base_learner = "rf", n_bootstraps = 50L,
  artificial_type = "knockoff_equi", sample_fraction = 1, replace = TRUE,
  random_state = 2L
)
fit_sparse <- stabl_fit(
  x = x, y = yp, lambda_grid = grid_sparse, family = "gaussian",
  base_learner = "rf", n_bootstraps = 50L,
  artificial_type = "knockoff_equi", sample_fraction = 1, replace = TRUE,
  random_state = 2L
)

report <- function(fit, label) {
  sc  <- fit$stabl_scores_                   # real features x lambdas
  sca <- fit$stabl_scores_artificial_        # decoys x lambdas
  mol  <- apply(sc, 1L, max)                 # max-over-lambda per real feature
  mola <- apply(sca, 1L, max)
  cat("\n==", label, "==\n")
  cat("per-lambda mean selection frequency (reals): ",
      paste(round(colMeans(sc), 2), collapse = " "), "\n")
  cat("per-lambda mean selection frequency (decoys): ",
      paste(round(colMeans(sca), 2), collapse = " "), "\n")
  cat("max-over-lambda score, reals:  median %.3f  q1 %.3f  min %.3f  n>=0.5: %d/%d\n" |>
        sprintf(median(mol), quantile(mol, 0.25), min(mol), sum(mol >= 0.5), length(mol)))
  cat("max-over-lambda score, decoys: median %.3f  q1 %.3f  min %.3f  n>=0.5: %d/%d\n" |>
        sprintf(median(mola), quantile(mola, 0.25), min(mola), sum(mola >= 0.5), length(mola)))
}
report(fit_auto,   "AUTO grid (median floor, probs 0.5-0.999) under null")
report(fit_sparse, "SPARSE-only grid (probs 0.9-0.999) under null")
cat("\nDONE\n")
