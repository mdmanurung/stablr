# WS-8 v3 — Candidate selectors for cross-cohort reproducibility (M6–M9).
#
# New arms benchmarked against the locked WS-8 v2 matrix (M0–M5):
#   M6  meta_homog     — DL random-effects meta-analysis of per-cohort Spearman
#                        associations with an explicit homogeneity filter
#                        (Candidate 1: effect-consistency selection).
#   M7  mint_stabl     — mixOmics::mint.spls (regression, study = cohort) inside
#                        the STABL bootstrap + artificial-feature machinery
#                        (Candidate 2: integrative joint modeling).
#   M7c mint_tuned     — tuned mint.spls without stability (control for M7).
#   M8  contrast_stabl — common + cohort-deviation (contrast-penalized) lasso
#                        inside STABL stability selection
#                        (Candidate 3: explicit homogeneity penalization).
#   M9  combat_stabl   — frozen-hyperparameter ComBat harmonization + pooled
#                        STABL (mandatory batch-correction baseline).
#
# All selectors follow the locked WS-8 v2 contract:
#   fn(X_train, y_train, cohort_train, covar = NULL)
#     -> list(ranked, native, scores, detail = NULL, detail_type, ...)
# Training data only.  Seed 20261013, n_boot = 200 (locked).
#
# Application case study — NOT calibration evidence.

# v2 helpers (.per_cohort_assoc, .re_meta_rho, .rank_by) come from
# engine/selectors.R, sourced by the driver before this file.

# ======================================================================
# M6 — RE meta-analysis with homogeneity filter
# ======================================================================

# DL random-effects meta of per-cohort Spearman rhos (Fisher-z), with pooled
# inference, Cochran's Q, tau2 and I2 per feature.
.re_meta_full <- function(rho_mat, n_vec) {
  z <- atanh(pmin(pmax(rho_mat, -0.9999), 0.9999))
  v <- matrix(rep(1 / (n_vec - 3), each = nrow(z)), nrow = nrow(z))
  wf <- 1 / v
  mu_f <- rowSums(z * wf) / rowSums(wf)
  Q <- rowSums(wf * (z - mu_f)^2)
  df <- ncol(z) - 1
  Cc <- rowSums(wf) - rowSums(wf^2) / rowSums(wf)
  tau2 <- pmax(0, (Q - df) / Cc)
  wr <- 1 / (v + tau2)
  mu_re <- rowSums(z * wr) / rowSums(wr)
  se_re <- sqrt(1 / rowSums(wr))
  p_pool <- 2 * pnorm(-abs(mu_re / se_re))
  p_Q <- pchisq(Q, df = df, lower.tail = FALSE)
  I2 <- pmax(0, (Q - df) / pmax(Q, 1e-12))
  list(rho_re = setNames(tanh(mu_re), rownames(rho_mat)),
       p_pool = setNames(p_pool, rownames(rho_mat)),
       p_Q = setNames(p_Q, rownames(rho_mat)),
       tau2 = setNames(tau2, rownames(rho_mat)),
       I2 = setNames(I2, rownames(rho_mat)))
}

sel_m6_meta_homog <- function(X, y, cohort, covar = NULL,
                              q_thresh = 0.10, bh_thresh = 0.05) {
  a <- .per_cohort_assoc(X, y, cohort)
  n_vec <- as.numeric(table(cohort)[colnames(a$rho)])
  m <- .re_meta_full(a$rho, n_vec)
  padj <- p.adjust(m$p_pool, method = "BH")
  homog_pass <- m$p_Q >= q_thresh
  native <- names(which(padj < bh_thresh & homog_pass))
  list(ranked = .rank_by(abs(m$rho_re)),   # same magnitude ranking as M0
       native = native,
       scores = m$rho_re,
       detail = a$rho, detail_type = "spearman_rho",
       extra = list(pooled_padj = padj, q_p = m$p_Q, tau2 = m$tau2,
                    I2 = m$I2, homog_pass = homog_pass))
}

# ======================================================================
# M7 — MINT-STABL (mint.spls regression inside STABL machinery)
# ======================================================================

# keepX grid = lambda analog (locked).  ncomp = 2 (locked).
M7_KEEPX <- c(2L, 5L, 8L, 12L, 18L, 25L, 35L, 50L)
M7_NCOMP <- 2L
M7_THR <- seq(0, 0.99, by = 0.01)
M7_LOAD_TOL <- 1e-10

.mint_fit_mask <- function(x_aug, yv, study, keepx, ncomp = M7_NCOMP) {
  fit <- tryCatch(
    mixOmics::mint.spls(x_aug, Y = yv, study = study, ncomp = ncomp,
                        keepX = rep(keepx, ncomp), scale = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  load <- fit$loadings$X
  if (is.null(dim(load))) return(NULL)
  matrixStats::rowMaxs(abs(as.matrix(load))) > M7_LOAD_TOL
}

sel_m7_mint_stabl <- function(X, y, cohort, covar = NULL, n_boot = 200L,
                              seed = 20261013) {
  n <- nrow(X); p <- ncol(X)
  art <- stablr:::make_artificial_features(
    X, p, "knockoff_equi", random_state = seed + 900L,
    knockoff_fallback = "stop")
  x_aug <- art$x_augmented
  art_pos <- p + seq_len(p)
  cn <- colnames(x_aug); cn[art_pos] <- paste0("art.", seq_len(p))
  colnames(x_aug) <- make.unique(cn)   # mixOmics requires unique column ids
  study <- factor(cohort)
  set.seed(seed + 901L)
  boots <- lapply(seq_len(n_boot), function(b) sample.int(n, n, replace = TRUE))
  nl <- length(M7_KEEPX)
  cnt_r <- matrix(0, p, nl); cnt_a <- matrix(0, p, nl); n_ok <- 0L
  n_fail <- 0L
  for (b in boots) {
    xw <- x_aug[b, , drop = FALSE]; yw <- y[b]; sw <- study[b]
    rownames(xw) <- paste0("r", seq_len(nrow(xw))); names(yw) <- rownames(xw)
    if (length(unique(sw)) < 2L) next
    for (l in seq_len(nl)) {
      mask <- .mint_fit_mask(xw, yw, sw, M7_KEEPX[l])
      if (is.null(mask)) { n_fail <- n_fail + 1L; next }
      cnt_r[, l] <- cnt_r[, l] + mask[seq_len(p)]
      cnt_a[, l] <- cnt_a[, l] + mask[art_pos]
      n_ok <- n_ok + 1L
    }
  }
  if (n_ok == 0L) stop("all mint.spls fits failed")
  sreal <- matrixStats::rowMaxs(cnt_r / n_ok)
  sart <- matrixStats::rowMaxs(cnt_a / n_ok)
  n_real_t <- vapply(M7_THR, function(t) sum(sreal > t), integer(1))
  n_art_t <- vapply(M7_THR, function(t) sum(sart > t), integer(1))
  fdrs <- (n_art_t + 1) / pmax(1, n_real_t)
  i <- which.min(fdrs)                       # first index wins (stablr rule)
  thr <- M7_THR[i]
  sc <- setNames(sreal, colnames(X))
  native <- names(sc)[sc > thr]
  list(ranked = .rank_by(sc), native = native, scores = sc,
       detail = NULL, detail_type = NULL,
       extra = list(min_fdr = setNames(fdrs[i], "min_fdr"),
                    thr_star = setNames(thr, "thr_star"),
                    frac_fit_fail = setNames(n_fail / (n_fail + n_ok),
                                             "frac_fit_fail"),
                    fdp_curve = fdrs))
}

# M7c — tuned mint.spls control (no stability, no artificial features).
sel_m7c_mint_tuned <- function(X, y, cohort, covar = NULL, seed = 20261013) {
  n <- nrow(X); p <- ncol(X)
  study <- factor(cohort)
  set.seed(seed + 902L)
  foldid <- sample(rep(seq_len(5L), length.out = n))
  mse <- sapply(M7_KEEPX, function(kx) {
    pr <- sapply(seq_len(5L), function(f) {
      tr <- foldid != f; te <- !tr
      fit <- tryCatch(
        mixOmics::mint.spls(X[tr, ], Y = y[tr], study = study[tr],
                            ncomp = M7_NCOMP, keepX = rep(kx, M7_NCOMP),
                            scale = TRUE),
        error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      pd <- tryCatch(predict(fit, newdata = X[te, , drop = FALSE],
                             study = factor(droplevels(study[te]))),
                     error = function(e) NULL)
      if (is.null(pd)) return(NA_real_)
      yh <- pd$predict[, , M7_NCOMP]
      mean((y[te] - as.numeric(yh))^2)
    })
    mean(pr, na.rm = TRUE)
  })
  best_k <- M7_KEEPX[which.min(replace(mse, is.na(mse), Inf))]
  fit <- mixOmics::mint.spls(X, Y = y, study = study, ncomp = M7_NCOMP,
                             keepX = rep(best_k, M7_NCOMP), scale = TRUE)
  load <- as.matrix(fit$loadings$X)
  mx <- matrixStats::rowMaxs(abs(load))
  sc <- setNames(mx, colnames(X))
  list(ranked = .rank_by(sc), native = names(sc)[sc > M7_LOAD_TOL],
       scores = sc, detail = NULL, detail_type = NULL,
       extra = list(keepX_cv = setNames(best_k, "keepX"),
                    cv_mse = setNames(mse, paste0("keepX_", M7_KEEPX))))
}

# ======================================================================
# M8 — Contrast-penalized (common + cohort-deviation) lasso in STABL
# ======================================================================

# Expanded design: per feature j, gamma_j (common) + delta_jk (deviation for
# cohort k, centered indicator).  penalty.factor: 1 on gamma, omega on delta.
.expand_contrast <- function(X, cohort) {
  cohs <- sort(unique(cohort))
  cols <- list(X)
  cn <- colnames(X)
  for (k in cohs) {
    w <- as.numeric(cohort == k); w <- w - mean(w)
    cols[[k]] <- sweep(X, 1L, w, `*`)
    cn <- c(cn, paste0(colnames(X), "|dev:", k))
  }
  out <- do.call(cbind, cols)
  colnames(out) <- cn
  out
}

sel_m8_contrast_stabl <- function(X, y, cohort, covar = NULL, n_boot = 200L,
                                  seed = 20261013, omega = 2) {
  n <- nrow(X); p <- ncol(X)
  cohs <- sort(unique(cohort)); K <- length(cohs)
  art <- stablr:::make_artificial_features(
    X, p, "knockoff_equi", random_state = seed + 900L,
    knockoff_fallback = "stop")
  x_aug <- art$x_augmented
  xe <- .expand_contrast(x_aug, cohort)
  p2 <- 2L * p
  gamma_idx <- seq_len(p2)                       # first block = common terms
  pf <- rep(c(1, omega), c(p2, p2 * K))
  set.seed(seed + 901L)
  boots <- lapply(seq_len(n_boot), function(b) sample.int(n, n, replace = TRUE))
  # lambda path from a full-data fit (deterministic given training data)
  path <- glmnet::glmnet(xe, y, alpha = 1, penalty.factor = pf, nlambda = 50)
  lambdas <- path$lambda
  nl <- length(lambdas)
  cnt_r <- matrix(0, p, nl); cnt_a <- matrix(0, p, nl); n_ok <- 0L
  for (b in boots) {
    fit <- tryCatch(
      glmnet::glmnet(xe[b, , drop = FALSE], y[b], alpha = 1,
                     penalty.factor = pf, lambda = lambdas),
      error = function(e) NULL)
    if (is.null(fit)) next
    bg <- as.matrix(fit$beta)[gamma_idx, , drop = FALSE] != 0
    cnt_r <- cnt_r + bg[seq_len(p), , drop = FALSE]
    cnt_a <- cnt_a + bg[p + seq_len(p), , drop = FALSE]
    n_ok <- n_ok + 1L
  }
  if (n_ok == 0L) stop("all contrast-lasso fits failed")
  sreal <- matrixStats::rowMaxs(cnt_r / n_ok)
  sart <- matrixStats::rowMaxs(cnt_a / n_ok)
  n_real_t <- vapply(M7_THR, function(t) sum(sreal > t), integer(1))
  n_art_t <- vapply(M7_THR, function(t) sum(sart > t), integer(1))
  fdrs <- (n_art_t + 1) / pmax(1, n_real_t)
  i <- which.min(fdrs)
  thr <- M7_THR[i]
  sc <- setNames(sreal, colnames(X))
  list(ranked = .rank_by(sc), native = names(sc)[sc > thr], scores = sc,
       detail = NULL, detail_type = NULL,
       extra = list(min_fdr = setNames(fdrs[i], "min_fdr"),
                    thr_star = setNames(thr, "thr_star"),
                    omega = setNames(omega, "omega"),
                    fdp_curve = fdrs))
}

# ======================================================================
# M9 — Frozen-hyperparameter ComBat + pooled STABL
# ======================================================================

# Parametric ComBat location/scale with empirical Bayes (Johnson et al. 2007),
# reimplemented to expose the fitted parameters for frozen application to a
# new batch.  Mirrors sva::ComBat (parametric, mean.only = FALSE) numerics,
# including the iterative it.sol posterior (MAP) updates.
.combat_postmean <- function(g_hat, g_bar, n, d_star, t2)
  (n * t2 * g_hat + d_star * g_bar) / (n * t2 + d_star)
.combat_postvar <- function(sum2, n, a, b)
  (0.5 * sum2 + b) / (0.5 * n + a - 1)

# Iterative EB posterior for one batch, given that batch's hyperprior scalars
# (g_bar, t2, a, b).  In ComBat the EB prior is PER BATCH, with moments taken
# over features within the batch (Johnson et al. 2007; sva:::it.sol).
# Z_sub: n_k x p standardized matrix for the batch.
.combat_itsol <- function(Z_sub, g_bar, t2, a, b, conv = 1e-4) {
  nk <- nrow(Z_sub)
  g_hat <- colMeans(Z_sub)
  g_old <- g_hat
  d_old <- pmax(apply(Z_sub, 2L, stats::var), 1e-12)
  change <- 1; iter <- 0L
  # Mirror sva:::it.sol exactly, including its signed-denominator stopping rule
  # (guard only against exact-zero division, which never binds in practice).
  while (change > conv && iter < 10000L) {
    g_new <- .combat_postmean(g_hat, g_bar, nk, d_old, t2)
    sum2 <- colSums(sweep(Z_sub, 2L, g_new, `-`)^2)
    d_new <- .combat_postvar(sum2, nk, a, b)
    change <- max(c(abs(g_new - g_old) / ifelse(g_old == 0, 1e-300, g_old),
                    abs(d_new - d_old) / d_old))
    g_old <- g_new; d_old <- d_new
    iter <- iter + 1L
  }
  list(gamma_star = g_new, delta_star = d_new)
}

# Per-batch hyperprior moments over features (sva aprior/bprior on the
# batch's gamma_hat / delta_hat rows).
.combat_priors <- function(gamma_hat, delta_hat) {
  # gamma_hat, delta_hat: features x batches
  g_bar <- colMeans(gamma_hat)
  t2 <- apply(gamma_hat, 2L, stats::var)
  t2 <- pmax(t2, 1e-12)
  m <- colMeans(delta_hat); s2 <- pmax(apply(delta_hat, 2L, stats::var), 1e-12)
  a <- (2 * s2 + m^2) / s2
  b <- (m * s2 + m^3) / s2
  list(gamma_bar = g_bar, tau2 = t2, lambda = a, theta = b)
}

.combat_fit <- function(dat, batch, mod) {
  # dat: n x p matrix; batch: coercible to factor; mod: model matrix WITHOUT
  # intercept (all-ones columns are dropped, matching sva).  Returns the
  # adjusted training matrix plus every parameter needed for frozen apply.
  batch <- factor(batch)
  batches <- levels(batch)
  n_b <- as.numeric(table(batch)[batches])
  n <- nrow(dat)
  batchmod <- model.matrix(~ 0 + batch)
  colnames(batchmod) <- batches
  design <- cbind(batchmod, mod)
  keep <- !apply(design, 2L, function(cl) all(cl == 1))
  design <- design[, keep, drop = FALSE]
  B_hat <- solve(crossprod(design), crossprod(design, dat))
  K <- length(batches)
  grand_mean <- as.numeric(crossprod(n_b / n, B_hat[seq_len(K), , drop = FALSE]))
  tmp <- design; tmp[, seq_len(K)] <- 0
  stand_mean <- sweep(tmp %*% B_hat, 2L, grand_mean, `+`)  # n x p
  resid <- dat - design %*% B_hat
  var_pooled <- pmax(colSums(resid^2) / n, 1e-12)
  Z <- sweep(dat - stand_mean, 2L, sqrt(var_pooled), `/`)
  gamma_hat <- sapply(batches, function(k) colMeans(Z[batch == k, , drop = FALSE]))
  delta_hat <- sapply(batches, function(k)
    apply(Z[batch == k, , drop = FALSE], 2L, stats::var))
  gamma_hat <- matrix(gamma_hat, nrow = ncol(dat),
                      dimnames = list(colnames(dat), batches))
  delta_hat <- pmax(matrix(delta_hat, nrow = ncol(dat),
                           dimnames = list(colnames(dat), batches)), 1e-12)
  pri <- .combat_priors(gamma_hat, delta_hat)   # per-batch scalars
  gamma_star <- matrix(NA_real_, ncol(dat), K, dimnames = dimnames(gamma_hat))
  delta_star <- gamma_star
  for (k in seq_len(K)) {
    post <- .combat_itsol(Z[batch == batches[k], , drop = FALSE],
                          pri$gamma_bar[k], pri$tau2[k],
                          pri$lambda[k], pri$theta[k])
    gamma_star[, k] <- post$gamma_star
    delta_star[, k] <- post$delta_star
  }
  adj <- Z
  for (k in batches) {
    idx <- batch == k
    adj[idx, ] <- sweep(adj[idx, , drop = FALSE], 2L, gamma_star[, k], `-`)
    adj[idx, ] <- sweep(adj[idx, , drop = FALSE], 2L, sqrt(delta_star[, k]), `/`)
  }
  adj <- sweep(adj, 2L, sqrt(var_pooled), `*`) + stand_mean
  list(adjusted = adj, stand_mean = stand_mean, grand_mean = grand_mean,
       var_pooled = var_pooled,
       priors = pri,
       gamma_hat = gamma_hat, delta_hat = delta_hat,
       gamma_star = gamma_star, delta_star = delta_star,
       batches = batches, B_hat = B_hat, design = design)
}

# Apply a frozen ComBat fit to a NEW batch.  The new-data stand_mean is the
# training GRAND MEAN only: the covariate (age) term of stand_mean requires
# labels, so it cannot be included for unlabeled new data.  Consequence of
# this label-free constraint: in eb_newbatch mode the new batch's location
# estimate absorbs any between-cohort mean-age shift (documented limitation
# of correct-then-predict; within-cohort age structure is preserved).
# mode = "eb_newbatch": treat the new data as one new batch and run the full
# ComBat EB machinery for it — estimate its location/scale moments AND its
# own per-batch hyperpriors from its (label-free) feature values, exactly as
# ComBat does per training batch.  Only grand_mean / var_pooled are frozen
# from training.  With a covariate-free fit, applying to a training batch
# reproduces the training adjustment to machine precision.
# mode = "standardize_only": strict training-only transform (no new-batch
# moment estimation).
.combat_apply <- function(fit, new_dat,
                          mode = c("eb_newbatch", "standardize_only")) {
  mode <- match.arg(mode)
  Z <- sweep(new_dat - rep(fit$grand_mean, each = nrow(new_dat)),
             2L, sqrt(fit$var_pooled), `/`)
  if (mode == "eb_newbatch") {
    g_hat <- colMeans(Z)
    d_hat <- pmax(apply(Z, 2L, stats::var), 1e-12)
    pri <- .combat_priors(matrix(g_hat, ncol = 1L), matrix(d_hat, ncol = 1L))
    post <- .combat_itsol(Z, pri$gamma_bar[1], pri$tau2[1],
                          pri$lambda[1], pri$theta[1])
    Z <- sweep(Z, 2L, post$gamma_star, `-`)
    Z <- sweep(Z, 2L, sqrt(pmax(post$delta_star, 1e-12)), `/`)
  }
  sweep(Z, 2L, sqrt(fit$var_pooled), `*`) +
    rep(fit$grand_mean, each = nrow(new_dat))
}

sel_m9_combat_stabl <- function(X, y, cohort, covar = NULL, n_boot = 200L,
                                seed = 20261013,
                                apply_mode = c("eb_newbatch", "standardize_only"),
                                protect_age = TRUE) {
  apply_mode <- match.arg(apply_mode)
  if (is.null(names(y)) && !is.null(rownames(X))) names(y) <- rownames(X)
  # protect_age = TRUE (locked M9): covariate-protected ComBat (~ age).
  # protect_age = FALSE (M9b sensitivity): mod = NULL, so the training
  # stand_mean IS the grand mean and the label-free test transform is exactly
  # consistent with the training adjustment (no covariate add-back mismatch).
  mod <- if (protect_age) model.matrix(~ y)[, -1, drop = FALSE] else NULL
  cf <- .combat_fit(X, cohort, mod)
  Xa <- cf$adjusted
  colnames(Xa) <- colnames(X)
  fit <- stablr::stabl_fit(x = Xa, y = y, lambda_grid = "auto",
                           base_learner = "rf", artificial_type = "knockoff_equi",
                           n_bootstraps = n_boot, sample_fraction = 1,
                           replace = TRUE, random_state = seed)
  sc <- stablr::get_importances(fit)
  if (is.null(names(sc))) names(sc) <- colnames(X)
  test_transform <- function(Xnew) {
    out <- .combat_apply(cf, Xnew, mode = apply_mode)
    colnames(out) <- colnames(Xnew)
    out
  }
  list(ranked = .rank_by(sc), native = names(sc)[sc >= 0.5], scores = sc,
       detail = NULL, detail_type = NULL,
       train_transformed = Xa, test_transform = test_transform,
       extra = list(apply_mode = setNames(apply_mode, "apply_mode"),
                    protect_age = setNames(as.numeric(protect_age),
                                           "protect_age")))
}

# ---- Registry -----------------------------------------------------------------
ws8v3_selectors <- list(
  M6  = sel_m6_meta_homog,
  M7  = sel_m7_mint_stabl,
  M7c = sel_m7c_mint_tuned,
  M8  = sel_m8_contrast_stabl,
  M9  = sel_m9_combat_stabl,
  M9s = function(X, y, cohort, covar = NULL)
    sel_m9_combat_stabl(X, y, cohort, covar, apply_mode = "standardize_only"),
  M9b = function(X, y, cohort, covar = NULL)
    sel_m9_combat_stabl(X, y, cohort, covar, protect_age = FALSE)
)
