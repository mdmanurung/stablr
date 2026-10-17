# Workstream A: does a model-X knockoff pairing (knockoff_mvr) repair the
# knockoff_equi (KE) miscalibration seen in the MINT benchmark?
#
# Full factorial on the SAME multi-cohort DGP as benchmarks/mint_comparison:
#   base_learner in {lasso, splsda}  x  artificial_type in
#   {random_permutation, knockoff_equi, knockoff_mvr}, 20 seeds, native
#   stabl_fit() machinery for every arm (the sPLS-DA arms use the new
#   base_learner = "splsda" path, not the wrapper).
# Fixed-X knockoffs are impossible here (n=240 < 2p=300; create.fixed
# augments rows), so the model-X question is equi vs mvr.
#
# Checkpoints: /mnt/shared-workspace/shared/mx_test/<learner>_<art>_seed%03d.rds
# Output:      /mnt/results/knockoff_mx/{mx_per_seed.csv,mx_summary.csv}

suppressMessages({
  library(stablr); library(matrixStats)
})
CKPT <- "/mnt/shared-workspace/shared/mx_test"
OUT  <- "/mnt/results/knockoff_mx"
dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT,  showWarnings = FALSE, recursive = TRUE)

N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "20"))
if (Sys.getenv("PILOT", "0") == "1") N_SEEDS <- 2L

P <- 150L; P_SIG <- 10L; N_K <- 80L; K <- 3L
BETA <- 1.2 / sqrt(P_SIG)
BATCH_IDX <- 21:40; BATCH_SHIFT <- c(0, 0.3, -0.3)
PREV_SHIFT <- c(0, 0.2, -0.2)
N_BOOT <- 100L
THR <- seq(0, 0.99, by = 0.01)

gen_data <- function(seed) {
  set.seed(seed + 50000L)
  x <- do.call(rbind, lapply(seq_len(K), function(k) {
    xk <- matrix(rnorm(N_K * P), N_K, P,
                 dimnames = list(paste0("s", (k - 1L) * N_K + seq_len(N_K)),
                                 paste0("f", seq_len(P))))
    xk[, BATCH_IDX] <- xk[, BATCH_IDX] + BATCH_SHIFT[k]
    xk }))
  eta <- as.vector(x[, seq_len(P_SIG)] %*% rep(BETA, P_SIG)) +
    rep(PREV_SHIFT, each = N_K)
  y <- setNames(rbinom(K * N_K, 1L, plogis(eta)), rownames(x))
  list(x = x, y = y, signal = paste0("f", seq_len(P_SIG)))
}

ARMS <- expand.grid(
  learner = c("lasso", "splsda"),
  art     = c("random_permutation", "knockoff_equi", "knockoff_mvr"),
  stringsAsFactors = FALSE
)
ARMS$arm <- paste0(ARMS$learner, "+", sub("knockoff_", "", ARMS$art),
                   ifelse(ARMS$art == "random_permutation", "", ""))

run_one <- function(learner, art, seed) {
  f <- file.path(CKPT, sprintf("%s_%s_seed%03d.rds", learner, art, seed))
  if (file.exists(f)) return(as.integer(readRDS(f)$done))
  d <- gen_data(seed)
  t0 <- Sys.time()
  fit <- stabl_fit(
    x               = d$x,
    y               = d$y,
    lambda_grid     = "auto",
    base_learner    = learner,
    family          = "binomial",
    n_bootstraps    = N_BOOT,
    artificial_type = art,
    random_state    = seed,
    splsda_ncomp    = 2L
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sel <- get_feature_names(fit)
  sig <- d$signal
  tp  <- length(intersect(sel, sig))
  n_sel <- length(sel)
  fdp <- if (n_sel == 0L) 0 else (n_sel - tp) / n_sel
  rec <- list(
    done = 1L, learner = learner, art = art, seed = seed,
    n_sel = n_sel, tp = tp, fdp = fdp,
    tpr = tp / length(sig),
    reported = fit$min_fdr_,
    elapsed = elapsed,
    selected = sel
  )
  saveRDS(rec, f)
  1L
}

tasks <- list()
for (i in seq_len(nrow(ARMS))) for (seed in seq_len(N_SEEDS)) {
  tasks <- c(tasks, list(list(ARMS$learner[i], ARMS$art[i], seed)))
}

t_all <- Sys.time()
res <- parallel::mclapply(tasks, function(tk) {
  suppressMessages(suppressWarnings(
    run_one(tk[[1]], tk[[2]], tk[[3]])
  ))
}, mc.cores = 8L, mc.preschedule = FALSE)
cat("total wall:", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1),
    "min; failed tasks:", sum(!vapply(res, is.integer, logical(1))), "\n")

# ---- Collect ----------------------------------------------------------------
recs <- list()
for (i in seq_len(nrow(ARMS))) for (seed in seq_len(N_SEEDS)) {
  f <- file.path(CKPT, sprintf("%s_%s_seed%03d.rds",
                               ARMS$learner[i], ARMS$art[i], seed))
  if (file.exists(f)) recs[[length(recs) + 1L]] <- readRDS(f)
}
per_seed <- do.call(rbind, lapply(recs, function(r) data.frame(
  learner = r$learner, art = r$art, seed = r$seed,
  n_sel = r$n_sel, tpr = r$tpr, fdp = r$fdp,
  reported = r$reported, elapsed = r$elapsed
)))
per_seed$calib <- ifelse(per_seed$reported > 0, per_seed$fdp / per_seed$reported, NA_real_)
per_seed$arm <- paste0(per_seed$learner, "+", sub("knockoff_", "", per_seed$art))
write.csv(per_seed, file.path(OUT, "mx_per_seed.csv"), row.names = FALSE)

jac <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}
summary_rows <- lapply(split(recs, list(vapply(recs, `[[`, "", "learner"),
                                      vapply(recs, `[[`, "", "art")),
                               drop = TRUE), function(rs) {
  sels <- lapply(rs, `[[`, "selected")
  jj <- NULL
  if (length(sels) > 1L) {
    pr <- utils::combn(length(sels), 2L)
    jj <- mean(vapply(seq_len(ncol(pr)), function(k)
      jac(sels[[pr[1, k]]], sels[[pr[2, k]]]), numeric(1)),
      na.rm = TRUE)
  }
  data.frame(
    learner = rs[[1]]$learner, art = rs[[1]]$art,
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
summ$arm <- paste0(summ$learner, "+", sub("knockoff_", "", summ$art))
summ <- summ[order(summ$learner, summ$art), ]
write.csv(summ, file.path(OUT, "mx_summary.csv"), row.names = FALSE)
cat("== summary ==\n")
print(summ[, c("arm", "mean_tpr", "mean_fdp", "mean_n_sel",
               "median_reported", "median_calib", "mean_jaccard", "mean_secs")],
      digits = 3, row.names = FALSE)
