# Gaussian-outcome extension of the all-learner comparison: do the binomial
# rankings transfer beyond classification?
#
# Design (mirrors benchmarks/learner_comparison, binomial):
#   6 learners x 2 DGPs x seeds, native stabl_fit(), n=200 single cohort,
#   p=150, 10 signal features (beta = 1.2/sqrt(10)), family = "gaussian",
#   y = X*beta + N(0, 1.2^2)  ->  R^2 ~ 0.5.
#   splsda is classification-only (mixOmics sPLS-DA has no gaussian
#   variant) and is excluded by design; ordinal remains excluded.
#   DGPs:  independent  - i.i.d. N(0,1) features, signals scattered across
#                         the 15 contiguous groups of 10 (misaligned)
#          block        - 15 blocks of 10, AR(0.6) within block, signals
#                         at block starts (aligned)
#   Pairings: random_permutation for every learner; knockoff_equi sub-arm
#   for all. Seeds: 10 (5 for rf/oblique_rf). n_bootstraps = 100,
#   lambda_grid = "auto", random_state = seed.
#
# Checkpoints: /mnt/shared-workspace/shared/learner_bench_gauss/<dgp>_<learner>_<art>_seed%03d.rds
# Output:      /mnt/results/learner_comparison/{lb_gauss_per_seed.csv,lb_gauss_summary.csv}

suppressMessages({
  library(stablr); library(matrixStats)
})
CKPT <- "/mnt/shared-workspace/shared/learner_bench_gauss"
OUT  <- "/mnt/results/learner_comparison"
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT,  showWarnings = FALSE, recursive = TRUE)

N_SEEDS_FAST <- as.integer(Sys.getenv("N_SEEDS_FAST", "10"))
N_SEEDS_SLOW <- as.integer(Sys.getenv("N_SEEDS_SLOW", "5"))
if (Sys.getenv("PILOT", "0") == "1") { N_SEEDS_FAST <- 2L; N_SEEDS_SLOW <- 1L }

N <- 200L; P <- 150L; P_SIG <- 10L
BETA <- 1.2 / sqrt(P_SIG)
NOISE_SD <- 1.2          # signal var = 10 * (1.2^2/10) = 1.44 -> R^2 ~ 0.5
RHO <- 0.6; BLOCK <- 10L
N_BOOT <- 100L

gen_data <- function(seed, dgp = c("independent", "block")) {
  dgp <- match.arg(dgp)
  set.seed(seed + 70000L)
  if (dgp == "independent") {
    x <- matrix(rnorm(N * P), N, P,
                dimnames = list(paste0("s", seq_len(N)), paste0("f", seq_len(P))))
    sig <- paste0("f", c(3, 17, 24, 38, 45, 57, 66, 73, 89, 142))
  } else {
    x <- do.call(cbind, lapply(seq_len(P %/% BLOCK), function(b) {
      z <- matrix(rnorm(N * BLOCK), N, BLOCK)
      out <- z
      if (RHO > 0) for (j in 2:BLOCK) {
        out[, j] <- RHO * out[, j - 1L] + sqrt(1 - RHO^2) * z[, j]
      }
      out
    }))
    dimnames(x) <- list(paste0("s", seq_len(N)), paste0("f", seq_len(P)))
    sig <- paste0("f", seq(1L, by = BLOCK, length.out = P_SIG))
  }
  eta <- as.vector(x[, sig] %*% rep(BETA, P_SIG)) + rnorm(N, sd = NOISE_SD)
  y <- setNames(eta, rownames(x))
  list(x = x, y = y, signal = sig)
}

FAST <- c("lasso", "elastic_net", "adaptive_lasso", "sparse_group_lasso")
SLOW <- c("rf", "oblique_rf")

tasks <- list()
for (dgp in c("independent", "block")) {
  for (learner in c(FAST, SLOW)) {
    n_seeds <- if (learner %in% SLOW) N_SEEDS_SLOW else N_SEEDS_FAST
    for (art in c("random_permutation", "knockoff_equi")) {
      for (seed in seq_len(n_seeds)) {
        tasks <- c(tasks, list(list(dgp, learner, art, seed)))
      }
    }
  }
}

run_one <- function(dgp, learner, art, seed) {
  f <- file.path(CKPT, sprintf("%s_%s_%s_seed%03d.rds", dgp, learner, art, seed))
  if (file.exists(f)) return(as.integer(readRDS(f)$done))
  d <- gen_data(seed, dgp)
  args <- list(
    x = d$x, y = d$y, lambda_grid = "auto",
    base_learner = learner, family = "gaussian",
    n_bootstraps = N_BOOT, artificial_type = art,
    random_state = seed
  )
  if (learner == "sparse_group_lasso") {
    args$feature_groups <- factor(rep(seq_len(P %/% BLOCK), each = BLOCK))
  }
  if (learner == "elastic_net") args$l1_ratio <- 0.5
  t0 <- Sys.time()
  fit <- suppressMessages(suppressWarnings(do.call(stabl_fit, args)))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sel <- get_feature_names(fit)
  tp <- length(intersect(sel, d$signal))
  n_sel <- length(sel)
  rec <- list(
    done = 1L, dgp = dgp, learner = learner, art = art, seed = seed,
    n_sel = n_sel, tp = tp,
    fdp = if (n_sel == 0L) 0 else (n_sel - tp) / n_sel,
    tpr = tp / length(d$signal),
    reported = fit$min_fdr_, elapsed = elapsed, selected = sel
  )
  saveRDS(rec, f)
  1L
}

t_all <- Sys.time()
res <- parallel::mclapply(tasks, function(tk) {
  suppressMessages(suppressWarnings(
    run_one(tk[[1]], tk[[2]], tk[[3]], tk[[4]])
  ))
}, mc.cores = 8L, mc.preschedule = FALSE)
cat("total wall:", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1),
    "min; failed tasks:", sum(!vapply(res, is.integer, logical(1))), "\n")

# ---- Collect ----------------------------------------------------------------
recs <- list()
for (tk in tasks) {
  f <- file.path(CKPT, sprintf("%s_%s_%s_seed%03d.rds", tk[[1]], tk[[2]], tk[[3]], tk[[4]]))
  if (file.exists(f)) recs[[length(recs) + 1L]] <- readRDS(f)
}
per_seed <- do.call(rbind, lapply(recs, function(r) data.frame(
  dgp = r$dgp, learner = r$learner, art = r$art, seed = r$seed,
  n_sel = r$n_sel, tpr = r$tpr, fdp = r$fdp,
  reported = r$reported, elapsed = r$elapsed
)))
per_seed$calib <- ifelse(per_seed$reported > 0,
                         per_seed$fdp / pmax(per_seed$reported, 1e-12), NA_real_)
write.csv(per_seed, file.path(OUT, "lb_gauss_per_seed.csv"), row.names = FALSE)

jac <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}
summary_rows <- lapply(split(recs, list(vapply(recs, `[[`, "", "dgp"),
                                        vapply(recs, `[[`, "", "learner"),
                                        vapply(recs, `[[`, "", "art")),
                               drop = TRUE), function(rs) {
  sels <- lapply(rs, `[[`, "selected")
  jj <- NA_real_
  if (length(sels) > 1L) {
    pr <- utils::combn(length(sels), 2L)
    jj <- mean(vapply(seq_len(ncol(pr)), function(k)
      jac(sels[[pr[1, k]]], sels[[pr[2, k]]]), numeric(1)),
      na.rm = TRUE)
  }
  data.frame(
    dgp = rs[[1]]$dgp, learner = rs[[1]]$learner, art = rs[[1]]$art,
    n_seeds = length(rs),
    mean_tpr = mean(vapply(rs, `[[`, 0, "tpr")),
    mean_fdp = mean(vapply(rs, `[[`, 0, "fdp")),
    mean_n_sel = mean(vapply(rs, `[[`, 0, "n_sel")),
    median_reported = median(vapply(rs, `[[`, 0, "reported")),
    median_calib = median(vapply(rs, `[[`, 0, "fdp") /
                            pmax(vapply(rs, `[[`, 0, "reported"), 1e-12)),
    mean_jaccard = jj,
    mean_secs = mean(vapply(rs, `[[`, 0, "elapsed"))
  )
})
summ <- do.call(rbind, summary_rows)
summ <- summ[order(summ$dgp, summ$learner, summ$art), ]
write.csv(summ, file.path(OUT, "lb_gauss_summary.csv"), row.names = FALSE)
cat("== summary ==\n")
print(summ[, c("dgp", "learner", "art", "mean_tpr", "mean_fdp", "mean_n_sel",
               "median_reported", "median_calib", "mean_jaccard", "mean_secs")],
      digits = 3, row.names = FALSE)
