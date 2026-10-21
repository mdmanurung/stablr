# Benchmark: stablr vs MINT sPLS-DA (mixOmics) on a multi-cohort classification DGP
#
# Design: K=3 cohorts x n=80, p=150, 10 shared true features (beta=0.8/sqrt(10)),
# cohort prevalence offsets (+0.3/-0.3), batch location shifts on 20 nuisance
# features (cohorts 2,3) that are marginally y-associated through cohort.
#
# Part A (selection, pooled n=240, 20 seeds):
#   1. stablr lasso+RP          (package default argmin-threshold rule)
#   2. stablr lasso+KE
#   3. splsda-stabl-RP          (wrapper: bootstrap splsda + RP columns + FDP+ rule)
#   4. splsda-stabl-KE          (wrapper with fixed-X equi knockoffs)
#   5. mint.splsda              (study-aware, tuned keepX/ncomp; no FDP claim)
#   6. splsda pooled            (study-blind control)
# Metrics: TPR, realized FDP, |selected|, reported min_fdr_ (calibration),
#          pairwise-Jaccard stability, runtime.
# Part B (prediction, LOCO over cohorts, 10 seeds): BER on held-out cohort for
#   mint.splsda, splsda pooled, glmnet-on-stablr-RP-selected, splsda-on-wrapper-selected.
#
# Checkpoints: one RDS per seed in CKPT_DIR; reruns skip completed seeds.
# Usage: N_SEEDS_A=20 N_SEEDS_B=10 Rscript mint_benchmark.R   (PILOT=1 for 2/1)

suppressMessages({
  library(stablr); library(mixOmics); library(glmnet); library(matrixStats)
})
set.seed <- base::set.seed
CKPT_DIR <- "/mnt/shared-workspace/shared/mint_bench"
OUT_DIR  <- "/mnt/results/mint_benchmark"
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)

N_SEEDS_A <- as.integer(Sys.getenv("N_SEEDS_A", "20"))
N_SEEDS_B <- as.integer(Sys.getenv("N_SEEDS_B", "10"))
if (Sys.getenv("PILOT", "0") == "1") { N_SEEDS_A <- 2L; N_SEEDS_B <- 1L }

P <- 150L; P_SIG <- 10L; N_K <- 80L; K <- 3L
BETA <- 1.2 / sqrt(P_SIG)
BATCH_IDX <- 21:40; BATCH_SHIFT <- c(0, 0.3, -0.3)
PREV_SHIFT <- c(0, 0.2, -0.2)
KEEPX_GRID <- c(2, 5, 8, 12, 18, 25, 35, 50)
NCOMP <- 2L; N_BOOT_WRAP <- 100L
THR <- seq(0, 0.99, by = 0.01)
LOAD_TOL <- 1e-10

# ---------------------------------------------------------------- DGP --------
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
  list(x = x, y = y, cohort = factor(rep(paste0("C", seq_len(K)), each = N_K)),
       signal = seq_len(P_SIG))
}
ber <- function(truth, pred) {
  truth <- truth == 1L; pred <- pred == 1L
  (mean(pred[!truth]) + mean(!pred[truth])) / 2
}
jaccard <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# ------------------------------------------------- splsda-stabl wrapper ------
# Mirrors stablr's scheme: artificial columns generated once (seeded), the
# augmented matrix is resampled by bootstrap, per (bootstrap x keepX) fit the
# learner, accumulate selection frequencies, max-over-grid score, FDP+ curve
# (inv_prop=1): FDRs = (n_art + 1) / max(n_real, 1), strict exceedance,
# argmin threshold (first index wins), support = score > thr_star.
splsda_stabl_fit <- function(x, y01, artificial = c("RP", "KE"), seed,
                             keepx_grid = KEEPX_GRID, ncomp = NCOMP,
                             n_boot = N_BOOT_WRAP) {
  artificial <- match.arg(artificial)
  n <- nrow(x); p <- ncol(x)
  # Use stablr's own artificial-feature machinery (same generators the
  # package uses internally; internal since Phase 14 PR-A):
  # RP = column-permuted competitors, KE = equi knockoffs.
  art <- stablr:::make_artificial_features(
    x, p, if (artificial == "RP") "random_permutation" else "knockoff_equi",
    random_state = seed + 900L, knockoff_fallback = "stop")
  x_aug <- art$x_augmented
  art_pos <- p + seq_len(p)   # artificial block = last p columns (stablr layout);
  p_art <- p                  # noise_col_indices is source-feature provenance only
  cn <- colnames(x_aug); cn[art_pos] <- paste0("art.", seq_len(p_art))
  colnames(x_aug) <- cn   # mixOmics requires unique column ids (glmnet does not)
  y <- factor(y01, levels = c("0", "1"))
  # bootstrap stream: mirror stablr defaults read from its formals
  sf  <- eval(formals(stabl_fit)$sample_fraction)
  rep_flag <- isTRUE(eval(formals(stabl_fit)$replace))
  set.seed(seed + 901L)
  boots <- lapply(seq_len(n_boot), function(b) {
    if (rep_flag) sample.int(n, n, replace = TRUE)
    else sample.int(n, max(2L, round(sf * n)))
  })
  nl <- length(keepx_grid)
  real_cols <- seq_len(p)
  cnt_r <- matrix(0, p, nl); cnt_a <- matrix(0, p, nl); n_ok <- 0L
  for (b in boots) {
    xw <- x_aug[b, , drop = FALSE]; yw <- y[b]
    rownames(xw) <- paste0("r", seq_len(nrow(xw))); names(yw) <- rownames(xw)
    if (length(unique(yw)) < 2L) next
    for (l in seq_len(nl)) {
      fit <- tryCatch(
        mixOmics::splsda(xw, yw, keepX = rep(keepx_grid[l], ncomp), ncomp = ncomp),
        error = function(e) NULL)
      if (is.null(fit)) next
      load <- fit$loadings$X
      if (is.null(dim(load))) next
      mask <- rowMaxs(abs(as.matrix(load))) > LOAD_TOL
      cnt_r[, l] <- cnt_r[, l] + mask[real_cols]
      cnt_a[, l] <- cnt_a[, l] + mask[art_pos]
      n_ok <- n_ok + 1L
    }
  }
  if (n_ok == 0L) stop("all splsda fits failed")
  sreal <- rowMaxs(cnt_r / n_ok); sart <- rowMaxs(cnt_a / n_ok)
  n_real_t <- vapply(THR, function(t) sum(sreal > t), integer(1))
  n_art_t  <- vapply(THR, function(t) sum(sart  > t), integer(1))
  fdrs <- (n_art_t + 1) / pmax(1, n_real_t)
  i <- which.min(fdrs)                        # first index wins (stablr ties)
  list(scores_real = sreal, scores_art = sart, FDRs = fdrs,
       min_fdr = fdrs[i], thr_star = THR[i],
       selected = which(sreal > THR[i]))
}

# ------------------------------------------------------ MINT / sPLS-DA -------
tune_splsda_family <- function(x, y01, study = NULL, grid = KEEPX_GRID,
                               nrepeat = 3L) {
  y <- factor(y01, levels = c("0", "1"))
  best <- list(ber = Inf, ncomp = 1L, keepx = 5)
  for (nc in 1:3) {
    t1 <- tryCatch(
      mixOmics::tune.splsda(x, y, study = study, ncomp = nc,
                            test.keepX = grid, validation = "Mfold",
                            folds = 5L, nrepeat = nrepeat,
                            measure = "BER", dist = "max.dist"),
      error = function(e) tryCatch(
        mixOmics::tune.splsda(x, y, ncomp = nc, test.keepX = grid,
                              validation = "Mfold", folds = 5L,
                              nrepeat = nrepeat, measure = "BER",
                              dist = "max.dist"),
        error = function(e2) NULL))
    if (is.null(t1)) next
    b <- min(t1$error.rate[, 1])
    if (b < best$ber) {
      best <- list(ber = b, ncomp = nc,
                   keepx = as.integer(t1$choice.keepX[[nc]]))
    }
  }
  if (!is.finite(best$ber)) best <- list(ber = NA_real_, ncomp = 1L, keepx = 5L)
  best
}
fit_predict_family <- function(x, y01, study, x_new, study_new,
                               tuned, mint = TRUE) {
  y <- factor(y01, levels = c("0", "1"))
  if (mint) { study <- droplevels(study); study_new <- droplevels(study_new) }
  fit <- tryCatch(
    if (mint) mixOmics::mint.splsda(x, y, study = study,
                                    keepX = rep(tuned$keepx, tuned$ncomp),
                                    ncomp = tuned$ncomp)
    else mixOmics::splsda(x, y, keepX = rep(tuned$keepx, tuned$ncomp),
                          ncomp = tuned$ncomp),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  load <- as.matrix(fit$loadings$X)
  selected <- which(rowMaxs(abs(load)) > LOAD_TOL)
  pr <- tryCatch({
    if (mint) predict(fit, newdata = x_new,
                      study = study_new)$class$max.dist[, tuned$ncomp]
    else predict(fit, newdata = x_new)$class$max.dist[, tuned$ncomp]
  }, error = function(e) tryCatch(
    predict(fit, newdata = x_new)$class$max.dist[, tuned$ncomp],
    error = function(e2) NULL))
  if (is.null(pr)) return(list(selected = selected, pred = NULL))
  list(selected = selected,
       pred = as.integer(as.character(pr)))  # "0"/"1" (char or factor) -> 0/1
}

# ================================================================ PART A =====
resA <- list()
for (seed in seq_len(N_SEEDS_A)) {
  ck <- file.path(CKPT_DIR, sprintf("seedA_%03d.rds", seed))
  if (file.exists(ck)) { resA[[seed]] <- readRDS(ck); next }
  d <- gen_data(seed); x <- d$x; y <- d$y; sig <- d$signal
  rec <- list(seed = seed, signal = sig)
  t0 <- Sys.time()
  grid <- auto_lambda_grid(x, y, family = "binomial", n_lambda = 20L)
  f_rp <- stabl_fit(x, y, lambda_grid = grid, n_bootstraps = 100L,
                    artificial_type = "random_permutation", random_state = seed)
  rec$stablr_rp <- list(sel = which(rowMaxs(get_stabl_scores(f_rp)) > f_rp$fdr_min_threshold_),
                        min_fdr = f_rp$min_fdr_)
  f_ke <- stabl_fit(x, y, lambda_grid = grid, n_bootstraps = 100L,
                    artificial_type = "knockoff_equi", random_state = seed)
  rec$stablr_ke <- list(sel = which(rowMaxs(get_stabl_scores(f_ke)) > f_ke$fdr_min_threshold_),
                        min_fdr = f_ke$min_fdr_)
  rec$stablr_time <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  w_rp <- splsda_stabl_fit(x, y, "RP", seed)
  w_ke <- splsda_stabl_fit(x, y, "KE", seed)
  rec$wrap_rp <- list(sel = w_rp$selected, min_fdr = w_rp$min_fdr)
  rec$wrap_ke <- list(sel = w_ke$selected, min_fdr = w_ke$min_fdr)
  rec$wrap_time <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  tuned_m <- tune_splsda_family(x, y, study = d$cohort)
  fm <- fit_predict_family(x, y, d$cohort, x, d$cohort, tuned_m, mint = TRUE)
  tuned_p <- tune_splsda_family(x, y, study = NULL)
  fp <- fit_predict_family(x, y, NULL, x, NULL, tuned_p, mint = FALSE)
  rec$mint  <- list(sel = fm$selected, tuned = tuned_m)
  rec$splsda<- list(sel = fp$selected, tuned = tuned_p)
  rec$mint_time <- as.numeric(Sys.time() - t0, units = "secs")
  saveRDS(rec, ck); resA[[seed]] <- rec
  cat(sprintf("seed %d/%d done (stablr %.0fs, wrap %.0fs, mint %.0fs)\n",
              seed, N_SEEDS_A, rec$stablr_time, rec$wrap_time, rec$mint_time))
}

# ================================================================ PART B =====
resB <- list()
for (seed in seq_len(N_SEEDS_B)) {
  ck <- file.path(CKPT_DIR, sprintf("seedB_%03d.rds", seed))
  if (file.exists(ck)) { resB[[seed]] <- readRDS(ck); next }
  d <- gen_data(seed); rec <- list(seed = seed)
  for (k in seq_len(K)) {
    tr <- d$cohort != levels(d$cohort)[k]; te <- !tr
    xtr <- d$x[tr, ]; ytr <- d$y[tr]; xte <- d$x[te, ]; yte <- d$y[te]
    ctr <- droplevels(d$cohort[tr]); cte <- droplevels(d$cohort[te])
    maj <- as.integer(names(which.max(table(ytr))))
    # stablr-RP selection + glmnet on selected
    grid <- auto_lambda_grid(xtr, ytr, family = "binomial", n_lambda = 20L)
    f <- stabl_fit(xtr, ytr, lambda_grid = grid, n_bootstraps = 100L,
                   artificial_type = "random_permutation", random_state = seed)
    sel <- which(rowMaxs(get_stabl_scores(f)) > f$fdr_min_threshold_)
    pred <- if (length(sel) >= 1L) {
      cv <- cv.glmnet(xtr[, sel, drop = FALSE], ytr, family = "binomial",
                      alpha = 1, nfolds = 5)
      as.integer(as.character(predict(cv, xte[, sel, drop = FALSE], type = "class",
                                      s = "lambda.min")))
    } else rep(maj, sum(te))
    rec[[paste0("stablr_rp_C", k)]] <- ber(yte, pred)
    # wrapper-RP selection + splsda on selected
    w <- splsda_stabl_fit(xtr, ytr, "RP", seed)
    pred <- if (length(w$selected) >= 2L) {
      yf <- factor(ytr, levels = c("0", "1"))
      fw <- mixOmics::splsda(xtr[, w$selected, drop = FALSE], yf,
                             keepX = rep(min(5L, length(w$selected)), NCOMP),
                             ncomp = NCOMP)
      pr <- predict(fw, newdata = xte[, w$selected, drop = FALSE])$class$max.dist[, NCOMP]
      as.integer(as.character(pr))
    } else rep(maj, sum(te))
    rec[[paste0("wrap_rp_C", k)]] <- ber(yte, pred)
    # mint.splsda
    tm <- tune_splsda_family(xtr, ytr, study = ctr, nrepeat = 2L)
    fm <- fit_predict_family(xtr, ytr, ctr, xte, cte, tm, mint = TRUE)
    rec[[paste0("mint_C", k)]] <- if (is.null(fm$pred)) NA_real_ else ber(yte, fm$pred)
    # pooled splsda
    tp <- tune_splsda_family(xtr, ytr, study = NULL, nrepeat = 2L)
    fp <- fit_predict_family(xtr, ytr, NULL, xte, NULL, tp, mint = FALSE)
    rec[[paste0("splsda_C", k)]] <- if (is.null(fp$pred)) NA_real_ else ber(yte, fp$pred)
  }
  saveRDS(rec, ck); resB[[seed]] <- rec
  cat(sprintf("LOCO seed %d/%d done\n", seed, N_SEEDS_B))
}

saveRDS(list(partA = resA, partB = resB),
        file.path(CKPT_DIR, "all_results.rds"))
cat("benchmark complete\n")
