# PR-2 (Phase 13) RNG-hygiene regression tests: F-CV-01, F-CV-02,
# F-CV-03/P4-F28.
#
# Fix summary:
#   F-CV-02        stabl_fit saves/restores the global `.Random.seed`.
#   F-CV-01        nested-CV outer folds get one seed each; the NULL case is
#                  drawn in the parent before forking, so results no longer
#                  depend on `cv_workers`.
#   F-CV-03/P4-F28 one independent stream per consumer (per-omic stabl fits,
#                  early/late fusion, cooperative foldid, downstream
#                  cv.glmnet predictors) instead of one reused seed.

# Deterministic fixtures (no RNG at construction time, so tests control the
# ambient stream exactly).
.rng_hygiene_x <- function(n = 40L) {
  x <- outer(seq_len(n), seq_len(6L), function(i, j) sin(i * j) * cos(i + j))
  dimnames(x) <- list(paste0("s", seq_len(n)), paste0("f", seq_len(6L)))
  x
}
.rng_hygiene_y_gaussian <- function(x) {
  setNames(x[, 1] + 0.5 * x[, 2], rownames(x))
}

test_that("F-CV-02: stabl_fit leaves the global .Random.seed untouched", {
  set.seed(20261)
  old_seed <- get(".Random.seed", envir = globalenv())

  suppressWarnings(stabl_fit(
    x = .rng_hygiene_x(),
    y = .rng_hygiene_y_gaussian(.rng_hygiene_x()),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    n_bootstraps = 5L,
    random_state = 7L
  ))

  expect_identical(get(".Random.seed", envir = globalenv()), old_seed)
})

test_that("F-CV-02: stabl_fit restores the no-seed state when none existed", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(list = ".Random.seed", envir = globalenv())
  }
  on.exit(set.seed(20261), add = TRUE)

  suppressWarnings(stabl_fit(
    x = .rng_hygiene_x(),
    y = .rng_hygiene_y_gaussian(.rng_hygiene_x()),
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    n_bootstraps = 5L,
    random_state = 7L
  ))

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("F-CV-02: consecutive NULL-seed stabl_fit calls are ambient-deterministic", {
  # Before the fix the second call continued the advanced stream and the two
  # fits differed; after the fix each call starts from the same ambient state.
  x <- .rng_hygiene_x()
  y <- .rng_hygiene_y_gaussian(x)
  lg <- data.frame(lambda = c(0.2, 0.1))

  set.seed(20262)
  fit1 <- suppressWarnings(stabl_fit(x = x, y = y, lambda_grid = lg,
                                     n_bootstraps = 5L))
  fit2 <- suppressWarnings(stabl_fit(x = x, y = y, lambda_grid = lg,
                                     n_bootstraps = 5L))

  expect_identical(fit1$stabl_scores_, fit2$stabl_scores_)
})

test_that("F-CV-03/P4-F28: .derive_stream_seeds is deterministic and distinct", {
  seeds1 <- .derive_stream_seeds(42L, 5L)
  seeds2 <- .derive_stream_seeds(42L, 5L)

  expect_length(seeds1, 5L)
  expect_identical(seeds1, seeds2)
  expect_false(any(duplicated(seeds1)))

  set.seed(20263)
  old_seed <- get(".Random.seed", envir = globalenv())
  invisible(.derive_stream_seeds(42L, 5L))
  expect_identical(get(".Random.seed", envir = globalenv()), old_seed)
})

test_that("F-CV-03/P4-F28: NULL derivation draws once from the ambient stream", {
  set.seed(20264)
  old_seed <- get(".Random.seed", envir = globalenv())

  seeds <- .derive_stream_seeds(NULL, 4L)

  expect_length(seeds, 4L)
  expect_false(any(duplicated(seeds)))
  expect_false(identical(get(".Random.seed", envir = globalenv()), old_seed))
})

test_that("F-CV-03/P4-F28: identical omics get independent bootstrap streams", {
  x <- .rng_hygiene_x()
  y <- .rng_hygiene_y_gaussian(x)
  lg <- data.frame(lambda = c(0.2, 0.1))
  fit_params <- list(
    base_learner = "lasso", family = "gaussian", n_bootstraps = 7L,
    artificial_type = NULL, hard_threshold = 1e-9, sample_groups = NULL,
    stratify_bootstrap = FALSE, bootstrap_strata = NULL, l1_ratio = NULL,
    # Deliberately different from the per-omic seeds: the override must win.
    random_state = 999L
  )

  out <- .fit_multiomic_per_omic(
    x_train_list = list(omic1 = x, omic2 = x),
    y_train = y,
    lambda_by_omic = list(omic1 = lg, omic2 = lg),
    x_valid_list = NULL,
    omic_names = c("omic1", "omic2"),
    fit_params = fit_params,
    omic_seeds = c(omic1 = 101L, omic2 = 202L)
  )

  expect_false(identical(
    out$fits$omic1$stabl_scores_,
    out$fits$omic2$stabl_scores_
  ))
})

test_that("F-CV-01: nested CV with random_state=NULL is worker-count-independent", {
  skip_on_cran()
  n <- 40L
  x <- .rng_hygiene_x()
  y_cont <- .rng_hygiene_y_gaussian(x)
  y_bin <- setNames(as.factor(y_cont > median(y_cont)), names(y_cont))
  xl <- list(mrna = x, mirna = x[, 6:1, drop = FALSE])
  lg <- data.frame(lambda = c(0.2, 0.1))
  run_ncv <- function(workers) {
    suppressWarnings(stabl_multiomic_nested_cv(
      x_list = xl, y = y_bin,
      candidates = list(list(name = "mrna", blocks = "mrna")),
      lambda_grid = lg, outer_v = 2L, outer_repeats = 1L, inner_v = 2L,
      metric = "ber", family = "binomial", artificial_type = NULL,
      hard_threshold = 1e-9, n_bootstraps = 3L, sample_fraction = 1.0,
      random_state = NULL, cv_workers = workers
    ))
  }

  set.seed(20265)
  seq_fit <- run_ncv(1L)
  set.seed(20265)
  par_fit <- run_ncv(2L)

  expect_identical(seq_fit$outer_predictions, par_fit$outer_predictions)
})

test_that("F-CV-03/P4-F28: train_validate is hermetic given random_state", {
  # Two runs with the same random_state but different ambient streams must
  # agree: covers the per-omic fan-out, the fusion seeds, and the seeded
  # downstream cv.glmnet predictors.
  x <- .rng_hygiene_x()
  y <- .rng_hygiene_y_gaussian(x)
  xl <- list(mrna = x, mirna = x[, 6:1, drop = FALSE])
  lg <- data.frame(lambda = c(0.2, 0.1))
  run_tv <- function() {
    suppressWarnings(stabl_multiomic_train_validate(
      x_train_list = xl, y_train = y, lambda_grid = lg,
      n_bootstraps = 5L, artificial_type = NULL, hard_threshold = 1e-9,
      random_state = 42L
    ))
  }

  set.seed(111)
  tv1 <- run_tv()
  set.seed(99999)
  tv2 <- run_tv()

  expect_identical(tv1$selected_features, tv2$selected_features)
  expect_identical(
    tv1$fits$mrna$stabl_scores_, tv2$fits$mrna$stabl_scores_
  )
  expect_identical(
    tv1$fits$mirna$stabl_scores_, tv2$fits$mirna$stabl_scores_
  )
})

test_that("F-CV-03/P4-F28: multiclass late fusion is hermetic given random_state", {
  # The downstream per-omic predictor is a cv.glmnet fit whose internal
  # foldid previously came from the ambient RNG; two runs with the same
  # random_state but different ambient streams must now agree.
  n <- 42L
  x <- .rng_hygiene_x(n)
  y_mc <- setNames(
    factor(rep(c("A", "B", "C"), length.out = n), levels = c("A", "B", "C")),
    rownames(x)
  )
  xl <- list(mrna = x, mirna = x[, 6:1, drop = FALSE])
  lg <- data.frame(lambda = c(0.2, 0.1))
  run_lf <- function(late_fusion_training) {
    suppressWarnings(stabl_multiomic_train_validate(
      x_train_list = xl, y_train = y_mc, lambda_grid = lg,
      family = "multinomial",
      n_bootstraps = 5L, artificial_type = NULL, hard_threshold = 1e-9,
      late_fusion = TRUE, late_fusion_training = late_fusion_training,
      n_iter_lf = 50L, random_state = 42L
    ))
  }

  set.seed(222)
  lf1 <- run_lf("oof")
  set.seed(888888)
  lf2 <- run_lf("oof")
  expect_identical(
    lf1$late_fusion$train_predictions,
    lf2$late_fusion$train_predictions
  )

  # 1.0.0 fencing: the legacy algorithm is internal-only; hermeticity is
  # still characterized by calling the retained unexported path directly.
  run_legacy <- function() {
    .late_fusion_multiomic_fit(
      selected_train = xl,
      selected_valid = NULL,
      y_train = y_mc, y_valid = NULL,
      omic_names = names(xl),
      family = "multinomial",
      n_iter_lf = 50L, random_state = 42L
    )
  }
  set.seed(333)
  lf3 <- run_legacy()
  set.seed(777777)
  lf4 <- run_legacy()
  expect_identical(
    lf3$train_predictions,
    lf4$train_predictions
  )
})

test_that("F-CV-05: nested-seed offset streams are mutually collision-free", {
  # .derive_nested_seed(base, i, off) = (base + off + i*7919) mod INT_MAX.
  # Two streams sharing a base collide iff (off1 - off2) is a multiple of
  # 7919 (an index shift then aligns them).  Production structure (verified
  # 2026-10-18, F-CV-05 audit):
  #   * the constant offsets below are pairwise collision-free for ALL
  #     indices (no pairwise difference is divisible by 7919);
  #   * the dynamic candidate family (offset k*100 in
  #     `.evaluate_stabl_candidates_inner`) is only ever used under a base
  #     pre-derived with offset 2000, so its values exceed every in-fold
  #     constant stream (offsets <= 4000) by construction:
  #     2000 + 100*k + 7919*fold_i >= 2000 + 100 + 7919 > 4000;
  #   * within the candidate family, collisions need k >= 7919 candidates
  #     (gcd(100, 7919) = 1).
  base <- 12345L
  offsets <- c(10L, 1000L, 2000L, 3000L, 4000L, 7000L, 29001L, 31001L,
               39001L, 41000L, 42000L, 47001L)
  diffs <- as.numeric(dist(offsets))
  expect_false(any(diffs %% 7919 == 0))

  # Concrete seeds over realistic index ranges are all distinct.
  idx <- seq_len(50L)
  seeds <- unlist(lapply(offsets, function(off) {
    vapply(idx, function(i) .derive_nested_seed(base, i, off), integer(1L))
  }))
  expect_identical(anyDuplicated(seeds), 0L)

  # Candidate family: internally distinct for any plausible candidate count.
  dyn <- unlist(lapply(seq_len(100L), function(k) {
    vapply(idx, function(i) .derive_nested_seed(base, i, k * 100L), integer(1L))
  }))
  expect_identical(anyDuplicated(dyn), 0L)

  # Production invariant: under the offset-2000 pre-derived base, every
  # candidate-stream value strictly exceeds every in-fold constant stream
  # (offsets 1000/2000/3000/4000) at the same outer index.
  outer_i <- 3L
  base2 <- .derive_nested_seed(base, outer_i, 2000L)
  cand_vals <- as.numeric(base2) + 100 * seq_len(100L) +
    7919 * rep(idx, each = 100L)
  const_vals <- as.numeric(base) + c(1000, 2000, 3000, 4000) + 7919 * outer_i
  expect_true(min(cand_vals) > max(const_vals))
})

test_that("F-019: environment-stability seeds come from the hashed stream derivation", {
  # The per-environment seeds must match .derive_stream_seeds() output
  # (decorrelated), not the historical consecutive random_state + i.
  seeds <- .derive_stream_seeds(11L, 3L)
  expect_false(any(diff(as.numeric(sort(seeds))) == 1))
  expect_identical(seeds, .derive_stream_seeds(11L, 3L))

  d <- local({
    set.seed(41L)
    n <- 36L
    x <- matrix(rnorm(n * 6L), n, 6L,
                dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(6L))))
    y <- setNames(0.8 * x[, 1] + rnorm(n, sd = 0.5), rownames(x))
    list(x = x, y = y, env = factor(rep(c("A", "B", "C"), each = 12L)))
  })
  res <- stabl_environment_stability(
    d$x, d$y, d$env,
    lambda_grid = data.frame(lambda = c(0.3, 0.15)),
    n_bootstraps = 4L, artificial_type = NULL, hard_threshold = 0.5,
    random_state = 11L
  )
  # Seeded with the hashed streams, the per-environment fits equal fits
  # seeded explicitly with each derived stream seed.
  for (i in seq_along(res$environments)) {
    lv <- res$environments[[i]]
    keep <- which(d$env == lv)
    ref <- stabl_fit(
      x = d$x[keep, ], y = d$y[keep],
      lambda_grid = data.frame(lambda = c(0.3, 0.15)),
      n_bootstraps = 4L, artificial_type = NULL, hard_threshold = 0.5,
      random_state = seeds[[i]]
    )
    expect_identical(res$env_fits[[lv]]$stabl_scores_, ref$stabl_scores_)
  }
})
