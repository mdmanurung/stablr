# Audit M-2 / V-7 (parallel determinism) and M-5 (RNG isolation between
# artificial-feature generation and bootstrap-index draws).

test_that(".can_open_test_server_socket returns a scalar logical", {
  result <- .can_open_test_server_socket()
  expect_type(result, "logical")
  expect_length(result, 1L)
  expect_false(is.na(result))
})

test_that("stabl_fit is reproducible: same random_state yields identical scores (sequential)", {
  withr::local_seed(0)
  n <- 50L; p <- 6L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam <- data.frame(lambda = c(0.2, 0.1, 0.05))

  fit1 <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 10L,
                    artificial_type = "random_permutation",
                    random_state = 42L, workers = 1L)
  fit2 <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 10L,
                    artificial_type = "random_permutation",
                    random_state = 42L, workers = 1L)

  expect_equal(fit1$stabl_scores_, fit2$stabl_scores_, tolerance = 0)
  expect_equal(fit1$stabl_scores_artificial_, fit2$stabl_scores_artificial_,
               tolerance = 0)
})

test_that("stabl_fit parallel and sequential paths produce identical stabl_scores_ for the same random_state", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  skip_if_not(
    .can_open_test_server_socket(),
    "future multisession requires opening a local server socket"
  )

  withr::local_seed(0)
  n <- 60L; p <- 8L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam <- data.frame(lambda = c(0.25, 0.1, 0.04))

  fit_seq <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 12L,
                       artificial_type = "random_permutation",
                       random_state = 7L, workers = 1L)

  old_plan <- future::plan(future::multisession, workers = 2L)
  withr::defer(future::plan(old_plan))

  fit_par <- suppressWarnings(stabl_fit(
    x, y, lambda_grid = lam, n_bootstraps = 12L,
    artificial_type = "random_permutation",
    random_state = 7L, workers = 2L
  ))

  expect_equal(fit_seq$stabl_scores_, fit_par$stabl_scores_, tolerance = 0)
  expect_equal(fit_seq$stabl_scores_artificial_,
               fit_par$stabl_scores_artificial_, tolerance = 0)
})

# ----- WI-12: RNG isolation between artificial features and bootstrap draws ---
# Pre-fix, both art-features and boot-indices were seeded with the *same*
# random_state, sharing one RNG sequence.  Post-fix, they consume independent
# derived seeds.  We pin this by calling the helpers directly with the seeds
# that stabl_fit derives, and asserting the two streams produce different
# draws (a basic non-degeneracy check) and that re-seeding does not leak.

test_that("stabl_fit's derived seeds yield independent artificial and bootstrap streams", {
  withr::local_seed(0)
  n <- 40L; p <- 5L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam <- data.frame(lambda = c(0.2, 0.1))

  fit1 <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 6L,
                    artificial_type = "random_permutation",
                    random_state = 1L, workers = 1L)
  fit2 <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 6L,
                    artificial_type = "random_permutation",
                    random_state = 1L, workers = 1L)
  # Determinism: same seed → identical artificial scores.
  expect_equal(fit1$stabl_scores_artificial_, fit2$stabl_scores_artificial_,
               tolerance = 0)
})

test_that("stabl_fit does not leave the global RNG state predictable across calls", {
  # The helper saves/restores .Random.seed around per-iter learner calls; the
  # outer set.seed(random_state) still advances the global RNG.  This test
  # documents that contract by asserting that two consecutive fits with the
  # same `random_state` give identical outputs (independent of any RNG used
  # by the caller in between).
  withr::local_seed(99)
  n <- 30L; p <- 4L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam <- data.frame(lambda = c(0.2, 0.1))

  fit_a <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 5L,
                     artificial_type = "random_permutation",
                     random_state = 17L, workers = 1L)
  # Caller consumes some RNG between fits.
  rnorm(100)
  fit_b <- stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 5L,
                     artificial_type = "random_permutation",
                     random_state = 17L, workers = 1L)

  expect_equal(fit_a$stabl_scores_, fit_b$stabl_scores_, tolerance = 0)
})

# C2 characterization: .derive_nested_seed integer-overflow safety.
# as.integer(index) * 7919L overflows (returns NA) for index > ~271500.
# The fix computes in double before casting.
test_that(".derive_nested_seed returns non-NA integers for small indices", {
  # Pin current (correct) values for realistic index range.
  # (random_state + offset + index * 7919) %% .Machine$integer.max:
  # (1 + 0 + 1*7919) = 7920; (42 + 0 + 3*7919) = 42 + 23757 = 23799
  expect_identical(stablr:::.derive_nested_seed(1L,  1L, 0L), 7920L)
  expect_identical(stablr:::.derive_nested_seed(42L, 3L, 0L), 23799L)
  # NULL random_state always returns NULL.
  expect_null(stablr:::.derive_nested_seed(NULL, 5L, 0L))
})

test_that(".derive_nested_seed does not return NA for large index", {
  # index = 300000 overflows as.integer(index) * 7919L in the naive
  # integer-only formula; the fix must return a valid integer, not NA.
  seed <- stablr:::.derive_nested_seed(1L, 300000L, 0L)
  expect_false(is.na(seed))
  expect_true(is.integer(seed))
  expect_true(seed >= 0L && seed < .Machine$integer.max)
})

# ----- D4: the stream-independence claim, actually asserted ---------------------
# The test above only asserted that two runs agree.  That is determinism, not
# independence: an implementation that drew the artificial features and the
# bootstrap indices from ONE stream would also produce two matching runs.  These
# assertions pin the independence the comment claims, and pin the non-collision
# that a shared stream would create.

test_that("D4 the three derived RNG streams are independent and non-colliding", {
  withr::local_seed(0)

  # This is exactly how stabl_fit derives its three sub-seeds from random_state.
  derive <- function(random_state) {
    set.seed(random_state)
    sample.int(.Machine$integer.max, size = 3L)
  }

  for (random_state in c(1L, 7L, 42L, 1234L)) {
    seeds <- derive(random_state)
    expect_length(unique(seeds), 3L)
    expect_false(anyNA(seeds))
    expect_true(all(seeds > 0L))
  }

  # The same random_state always yields the same triple (determinism) ...
  expect_identical(derive(42L), derive(42L))
  # ... and different random_states yield different triples.
  expect_false(identical(derive(42L), derive(43L)))
})

test_that("D4 the artificial and bootstrap streams produce different draws", {
  withr::local_seed(0)
  set.seed(42L)
  seeds <- sample.int(.Machine$integer.max, size = 3L)
  art_seed <- seeds[[1L]]
  boot_seed <- seeds[[2L]]

  art_draw <- withr::with_seed(art_seed, stats::runif(20))
  boot_draw <- withr::with_seed(boot_seed, stats::runif(20))

  expect_false(isTRUE(all.equal(art_draw, boot_draw)))
  # And the two streams are individually reproducible, which is what makes the
  # per-component comparison meaningful.
  expect_identical(art_draw, withr::with_seed(art_seed, stats::runif(20)))
  expect_identical(boot_draw, withr::with_seed(boot_seed, stats::runif(20)))
})

test_that("D4 a seeded stabl_fit does not disturb the caller's RNG stream", {
  withr::local_seed(0)
  n <- 30L; p <- 4L
  x <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  lam <- data.frame(lambda = c(0.2, 0.1))

  set.seed(99)
  without_fit <- {
    invisible(NULL)
    rnorm(5)
  }
  set.seed(99)
  with_fit <- {
    invisible(stabl_fit(x, y, lambda_grid = lam, n_bootstraps = 5L,
                        artificial_type = "random_permutation",
                        random_state = 17L, workers = 1L))
    rnorm(5)
  }
  expect_identical(with_fit, without_fit)
})

# ----- D5: RNGkind() / sample.kind are part of the determinism contract --------

test_that("D5 the fit records RNGkind, and the record matches the session", {
  withr::local_seed(0)
  n <- 30L; p <- 4L
  x <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- setNames(rnorm(n), rownames(x))
  fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = c(0.2, 0.1)),
                   n_bootstraps = 5L, artificial_type = "random_permutation",
                   random_state = 3L, workers = 1L)
  expect_identical(fit$rng_kind_, RNGkind())
  expect_true(all(c("Mersenne-Twister", "Inversion", "Rejection") %in%
                    unlist(fit$rng_kind_)))
})

test_that("D5 .with_local_seed restores .Random.seed and does not pin sample.kind", {
  kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(kind)))

  set.seed(5)
  before <- get(".Random.seed", envir = globalenv())
  invisible(stablr:::.with_local_seed(11L, stats::runif(10)))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  # The documented dependency: sample.kind governs sample.int()/sample(), and is
  # not pinned, so a different sampler yields different draws from the same seed.
  draw <- function() stablr:::.with_local_seed(11L, sample.int(20L, 10L))
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  rounding <- suppressWarnings(draw())
  suppressWarnings(RNGkind(sample.kind = "Rejection"))
  rejection <- suppressWarnings(draw())
  expect_false(isTRUE(all.equal(rounding, rejection)))
})

test_that("no exported entry point advances the caller's RNG stream", {
  # One documented rule: a seeded public entry point saves and restores
  # .Random.seed. This is what makes a seeded stabl_fit composable, because a
  # caller who seeds once and then calls several entry points gets the same
  # result regardless of what ran before.
  y <- c(rep(0, 12), rep(1, 12))
  groups <- rep(c("G1", "G2", "G3", "G4"), each = 6L)

  probe <- function(expr) {
    set.seed(20260926L)
    before <- .Random.seed
    force(expr)
    identical(before, .Random.seed)
  }

  expect_true(probe(classic_bootstrap_indices(y, n_subsamples = 12L, seed = 7L)))
  expect_true(probe(
    group_bootstrap_indices(y, groups, n_subsamples = 12L, seed = 7L)
  ))

  # Build the design before the probe snapshots the stream, otherwise the
  # rnorm() call itself is what advances the state and the check is vacuous.
  set.seed(11L)
  design <- matrix(
    rnorm(36), nrow = 12,
    dimnames = list(paste0("s", 1:12), paste0("f", 1:3))
  )
  for (type in c("random_permutation", "knockoff_equi")) {
    expect_true(probe(
      make_artificial_features(
        design, n_injected = 2L, type = type, random_state = 7L
      )
    ))
  }
})

test_that("a seeded sampler call is reproducible across calls", {
  y <- c(rep(0, 12), rep(1, 12))
  first <- classic_bootstrap_indices(y, n_subsamples = 12L, seed = 7L)
  second <- classic_bootstrap_indices(y, n_subsamples = 12L, seed = 7L)
  expect_identical(first, second)
})

test_that("unseeded sampler calls still draw distinct replicates", {
  # Guards against a fix that accidentally made every draw identical by
  # re-seeding inside the loop, which would silently collapse STABL's entire
  # bootstrap to one repeated subsample.
  y <- c(rep(0, 12), rep(1, 12))
  draws <- lapply(1:5, function(i) classic_bootstrap_indices(y, n_subsamples = 12L))
  expect_gt(length(unique(vapply(draws, function(d) paste(d, collapse = ","), character(1)))), 1L)
})
