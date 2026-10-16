# Phase 6 figure data extraction — same DGP/seeds as phase6_dgp_validation.R
# Captures per-seed FDP+ curves and max-over-lambda score distributions.
suppressMessages(library(stablr))
setwd("/workspace")
N <- 100L; P <- 150L; P_SIGNAL <- 10L; N_BOOT <- 100L

make_x <- function(seed, correlated = FALSE) {
  set.seed(seed)
  if (!correlated) return(matrix(rnorm(N * P), N, P,
    dimnames = list(paste0("s", seq_len(N)), paste0("f", seq_len(P)))))
  block <- 10L; nblk <- P %/% block
  x <- matrix(NA_real_, N, P)
  for (b in seq_len(nblk)) {
    idx <- (b - 1L) * block + seq_len(block)
    eps <- matrix(rnorm(N * block), N, block)
    x[, idx] <- t(apply(eps, 1L, function(row) {
      out <- numeric(block); out[1L] <- row[1L]
      for (j in 2L:block) out[j] <- 0.6 * out[j - 1L] + row[j]
      out }))
  }
  dimnames(x) <- list(paste0("s", seq_len(N)), paste0("f", seq_len(P)))
  x
}
make_y <- function(x, seed, signal = TRUE) {
  set.seed(seed + 10000L)
  eta <- if (signal) 0.8 * x[, seq_len(P_SIGNAL)] %*% rep(1 / sqrt(P_SIGNAL), P_SIGNAL) else rep(0, N)
  setNames(as.numeric(eta + rnorm(N, sd = 0.5)), rownames(x))
}

thr <- seq(0, 0.99, by = 0.01)
out <- list()
for (design in c("D1_sparse", "D2_null", "D3_correlated")) {
  correlated <- design == "D3_correlated"
  signal <- design != "D2_null"
  for (art in c("random_permutation", "knockoff_equi")) {
    for (seed in seq_len(20L)) {
      x <- make_x(seed, correlated); y <- make_y(x, seed, signal)
      grid <- auto_lambda_grid(x, y, family = "gaussian", n_lambda = 20L)
      fit <- stabl_fit(x, y, lambda_grid = grid, n_bootstraps = N_BOOT,
                       artificial_type = art, random_state = seed)
      out[[length(out) + 1L]] <- list(
        design = design, artificial = art, seed = seed,
        FDRs = fit$FDRs_, min_fdr = fit$min_fdr_,
        thr_star = fit$fdr_min_threshold_,
        scores_real = matrixStats::rowMaxs(get_stabl_scores(fit)),
        scores_art  = matrixStats::rowMaxs(fit$stabl_scores_artificial_))
    }
  }
  cat(design, "done\n")
}
saveRDS(out, "/workspace/phase6_out/curves.rds")
cat("extracted", length(out), "curve records\n")
