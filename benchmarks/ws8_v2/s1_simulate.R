# WS-8 v2 S1 — synthetic calibration simulations (followup_s1_s2.md).
#
# Generates donor_tables.rds per scenario x replicate, anchored to the real
# atlas: per-cohort mean logit-proportions, pooled within-cohort covariance,
# empirical per-cohort age / total-count / sex / ethnicity distributions.
# Donor generation: x_i ~ MVN(mu_c + beta_c * a_i, Sigma) on the logit scale
# with a_i = (age_i - 50)/10; counts ct_ij ~ Binomial(total_i, plogis(x_ij)).
# The analysis pipeline then re-applies the same logit((ct+0.5)/(total+0.5))
# transform as on real data — the generative and analysis transforms match.
#
# Scenarios (effect calibration targets population Spearman rho):
#   s1a_hom    10 true effects, consistent across cohorts (rho ~ U(0.20,0.35))
#   s1b_het    8 consistent + 2 sign-inconsistent (beta_jc ~ N(0, tau^2))
#   s1c_null   no true effects; real correlation structure (null probe)
#   s1d_sparse 3 strong consistent effects (rho = 0.50)
#
# True features: the 10 features at middle ranks (25-34) of real mean |rho|
# with age — moderate real association, not the atlas's top hits; signs from
# the real data.  Cohort names/sizes mirror the real design (incl. a 33-donor
# cohort, so the rank-deficient knockoff stress is preserved).
#
# Usage: Rscript benchmarks/ws8_v2/s1_simulate.R   # from repo root
#   env WS8_S1_OUT   output root (default /workspace/ws8_out/s1_inputs)
#   env WS8_S1_REPS  replicates per scenario (default 5)
#
# Inputs are regenerable from the committed manifest (seeds + MD5); the rds
# files themselves are session scratch and are not committed.

pkgload::load_all(".", compile = FALSE, quiet = TRUE)

out_root <- Sys.getenv("WS8_S1_OUT", "/workspace/ws8_out/s1_inputs")
n_reps   <- as.integer(Sys.getenv("WS8_S1_REPS", "5"))
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
base_seed <- 20261013L

# ---- Anchor statistics from the real atlas ------------------------------------
obj <- readRDS("benchmarks/ws8_v2/results/m1/donor_tables.rds")
prop <- (obj$ct_mat + 0.5) / (obj$total + 0.5)
X <- log(prop / (1 - prop))
meta <- obj$meta[match(rownames(X), obj$meta$donor), ]
cohort <- as.character(meta$dataset)
cohorts <- sort(unique(cohort))

mu_c <- sapply(cohorts, function(cc) colMeans(X[cohort == cc, , drop = FALSE]))
colnames(mu_c) <- cohorts
Xc <- X - t(mu_c[, cohort])           # within-cohort centred
Sigma <- stats::cov(Xc)
Sigma <- as.matrix(corpcor::make.positive.definite(Sigma))
Sig_chol <- chol(Sigma)               # upper U with Sigma = t(U) %*% U; x = mu + z %*% U

a_real <- (as.numeric(meta$age) - 50) / 10
sd_a <- stats::sd(a_real)

# Real per-cohort Spearman association with age (for feature choice/signs)
rho_real <- sapply(cohorts, function(cc) {
  idx <- cohort == cc
  apply(X[idx, ], 2, function(v)
    suppressWarnings(stats::cor(v, a_real[idx], method = "spearman")))
})
mean_abs <- rowMeans(abs(rho_real))
ord <- order(-mean_abs)
true10 <- rownames(rho_real)[ord[25:34]]        # middle ranks 25-34
sign10 <- sign(rowMeans(rho_real)[true10]); sign10[sign10 == 0] <- 1

beta_for_rho <- function(rho, j) {
  rho * sqrt(Sigma[j, j]) / (sd_a * sqrt(1 - rho^2))
}

scenarios <- c("s1a_hom", "s1b_het", "s1c_null", "s1d_sparse")

# ---- Generators ----------------------------------------------------------------
sim_cohort <- function(cc, beta_c, seed) {
  set.seed(seed)
  idx <- which(cohort == cc)
  n <- length(idx)
  ages <- sample(as.numeric(meta$age)[idx], n, replace = TRUE)
  tots <- sample(obj$total[idx], n, replace = TRUE)
  sx <- as.character(meta$sex)[idx]; sx <- sample(sx, n, replace = TRUE)
  et <- as.character(meta$ethnicity)[idx]; et <- sample(et, n, replace = TRUE)
  a <- (ages - 50) / 10
  z <- matrix(stats::rnorm(n * ncol(X)), n, ncol(X))
  x <- sweep(z %*% Sig_chol, 2, mu_c[, cc], `+`)
  if (length(beta_c))
    x[, names(beta_c)] <- x[, names(beta_c), drop = FALSE] +
      outer(a, beta_c)                    # a_i * beta_j added to feature j
  p <- plogis(pmin(pmax(x, -20), 20))
  ct <- matrix(stats::rbinom(n * ncol(X), size = rep(tots, ncol(X)),
                             prob = as.vector(p)), n, ncol(X),
               dimnames = list(NULL, colnames(X)))
  list(ct = ct, total = tots,
       meta = data.frame(age = ages, sex = sx, ethnicity = et,
                         stringsAsFactors = FALSE))
}

truth_rows <- list(); man_rows <- list(); beta_rows <- list()
for (si in seq_along(scenarios)) {
  sc <- scenarios[si]
  for (r in seq_len(n_reps)) {
    seed <- base_seed + si * 100000L + r
    set.seed(seed)
    # --- true effects for this replicate ---
    beta_c <- lapply(cohorts, function(cc) setNames(numeric(0), character(0)))
    names(beta_c) <- cohorts
    if (sc %in% c("s1a_hom", "s1b_het")) {
      rho_t <- stats::runif(10, 0.20, 0.35)
      b <- sign10 * mapply(beta_for_rho, rho_t, true10)
      for (cc in cohorts) beta_c[[cc]] <- setNames(b, true10)
      if (sc == "s1b_het") {
        incons <- true10[1:2]             # deterministic pair
        tau <- stats::median(abs(b))
        for (cc in cohorts) {
          set.seed(seed + match(cc, cohorts))
          beta_c[[cc]][incons] <- stats::rnorm(2, 0, tau)
        }
      }
      truth <- data.frame(scenario = sc, feature = true10,
                          consistent = !(true10 %in% true10[1:2]) | sc == "s1a_hom",
                          base_beta = as.numeric(b), stringsAsFactors = FALSE)
    } else if (sc == "s1d_sparse") {
      b <- sign10[1:3] * mapply(beta_for_rho, rep(0.5, 3), true10[1:3])
      for (cc in cohorts) beta_c[[cc]] <- setNames(b, true10[1:3])
      truth <- data.frame(scenario = sc, feature = true10[1:3],
                          consistent = TRUE, base_beta = as.numeric(b),
                          stringsAsFactors = FALSE)
    } else {
      truth <- data.frame(scenario = character(0), feature = character(0),
                          consistent = logical(0), base_beta = numeric(0),
                          stringsAsFactors = FALSE)
    }
    for (cc in cohorts)
      if (length(beta_c[[cc]]))
        beta_rows[[length(beta_rows) + 1]] <- data.frame(
          scenario = sc, rep = r, cohort = cc, feature = names(beta_c[[cc]]),
          beta = as.numeric(beta_c[[cc]]), stringsAsFactors = FALSE)

    # --- simulate cohorts ---
    ct_l <- list(); tot_l <- list(); meta_l <- list()
    for (cc in cohorts) {
      out <- sim_cohort(cc, beta_c[[cc]], seed + match(cc, cohorts) * 1000L)
      ids <- paste0("s1_", sc, "_r", r, "_", cc, "_", seq_len(nrow(out$ct)))
      rownames(out$ct) <- ids
      ct_l[[cc]] <- out$ct; tot_l[[cc]] <- out$total
      meta_l[[cc]] <- data.frame(donor = ids, dataset = cc, out$meta,
                                 stringsAsFactors = FALSE)
    }
    sim <- list(ct_mat = do.call(rbind, ct_l),
                total = unlist(tot_l),
                meta = do.call(rbind, meta_l))
    d <- file.path(out_root, paste0(sc, "_rep", r))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    saveRDS(sim, file.path(d, "donor_tables.rds"))
    man_rows[[length(man_rows) + 1]] <- data.frame(
      scenario = sc, rep = r, seed = seed,
      path = file.path(d, "donor_tables.rds"),
      md5 = unname(tools::md5sum(file.path(d, "donor_tables.rds"))),
      stringsAsFactors = FALSE)
    if (nrow(truth)) {
      truth$rep <- r
      truth_rows[[length(truth_rows) + 1]] <- truth
    }
    cat(sprintf("%s rep %d done (seed %d)\n", sc, r, seed))
  }
}
write.csv(do.call(rbind, man_rows),
          "benchmarks/ws8_v2/results/s1_input_manifest.csv", row.names = FALSE)
write.csv(do.call(rbind, truth_rows),
          "benchmarks/ws8_v2/results/s1_truth.csv", row.names = FALSE)
if (length(beta_rows))
  write.csv(do.call(rbind, beta_rows),
            "benchmarks/ws8_v2/results/s1_truth_by_cohort.csv", row.names = FALSE)
cat("S1 inputs written to", out_root, "\n")
