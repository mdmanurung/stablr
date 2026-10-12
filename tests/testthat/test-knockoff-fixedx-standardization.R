# Regression tests for improvement item A1: the fixed-X knockoff path must
# standardize the design before calling knockoff::create.fixed() (whose
# documented contract requires a normalized n-by-p matrix) and must record
# the standardization in the returned provenance metadata.
#
# Scale-contract note: create.fixed() column-normalizes its input to unit
# Euclidean norm and returns knockoffs on that internal scale (entries
# ~ 1/sqrt(n)).  The package rescales $Xk by the standardized column norms
# before shipping, so the knockoffs satisfy t(Xk) Xk = t(X) X on the scale
# actually paired with the design.  Quality and calibration assertions
# probe this shipped contract and run wherever the knockoff package is
# available.  Standardization-contract assertions run unconditionally (the
# fixed-X branch is exercised through a mocked create.fixed where needed).

.fixedx_probe_sd_tolerance <- 0.5

# Probes the SHIPPED fixed-X contract end to end: create.fixed() returns
# knockoffs on its internal unit-norm column scale (entries ~ 1/sqrt(n));
# the package rescales them by the standardized column norms before use.
# The probe replicates that rescale and requires column sds ~ 1 plus the
# defining Gram identity t(Xk) Xk = t(X) X.  Returns FALSE when the
# knockoff package is unavailable.
.fixedx_knockoff_probe_ok <- function(seed = 20261012L, n = 120L, p = 8L,
                                      tolerance = .fixedx_probe_sd_tolerance) {
  if (!requireNamespace("knockoff", quietly = TRUE)) {
    return(FALSE)
  }
  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x <- scale(x)
  kd <- suppressWarnings(knockoff::create.fixed(x, sigma = 1))
  col_norms <- sqrt(colSums(x^2))
  xk <- kd$Xk %*% diag(col_norms)
  sds <- apply(xk, 2L, sd)
  gram_ok <- max(abs(crossprod(xk) - crossprod(x))) < 1e-6
  all(abs(sds - 1) < tolerance) && gram_ok
}

.fixedx_probe_ok <- .fixedx_knockoff_probe_ok()

# Heterogeneous-scale design: column sds span 0.01 to 100.
.fixedx_heterogeneous_design <- function(n, p = 5L, seed = 20261013L) {
  set.seed(seed)
  scales <- exp(seq(log(0.01), log(100), length.out = p))
  x <- matrix(rnorm(n * p), nrow = n) %*% diag(scales)
  colnames(x) <- paste0("f", seq_len(p))
  x
}

# ── Standardization contract (unconditional) ─────────────────────────────────

test_that(".standardize_columns handles heterogeneous scales", {
  set.seed(20261014L)
  n <- 500L
  scales <- c(0.01, 0.1, 1, 10, 100)
  x <- matrix(rnorm(n * length(scales)), nrow = n) %*% diag(scales)
  colnames(x) <- paste0("f", seq_along(scales))

  out <- stablr:::.standardize_columns(x)

  expect_identical(dim(out$x), dim(x))
  expect_identical(colnames(out$x), colnames(x))
  expect_equal(unname(colMeans(out$x)), rep(0, length(scales)), tolerance = 1e-12)
  expect_equal(unname(apply(out$x, 2L, sd)), rep(1, length(scales)), tolerance = 1e-12)
  expect_equal(out$means, colMeans(x))
  expect_equal(out$sds, apply(x, 2L, sd))
  # Affine per column: (x - mean) / sd
  expect_equal(out$x[, 1L], (x[, 1L] - mean(x[, 1L])) / sd(x[, 1L]))
})

test_that(".standardize_columns leaves zero-variance columns at scale 1", {
  set.seed(20261015L)
  n <- 50L
  x <- cbind(matrix(rnorm(n * 2L), nrow = n), rep(3.5, n))

  out <- stablr:::.standardize_columns(x)

  expect_identical(out$sds[[3L]], 1)          # scale left at 1 ...
  expect_equal(out$means[[3L]], 3.5)          # ... but still centered
  expect_true(all(out$x[, 3L] == 0))
  expect_equal(apply(out$x[, 1L:2L, drop = FALSE], 2L, sd), c(1, 1),
               tolerance = 1e-12)
})

test_that("fixed-X branch passes a standardized matrix to create.fixed", {
  skip_if_not_installed("knockoff")
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    create.fixed = function(x, sigma = NULL, ...) {
      captured$x <- x
      captured$sigma <- sigma
      list(Xk = x + 1)  # deterministic, dimension-preserving stand-in
    },
    .package = "knockoff"
  )

  n <- 60L
  x <- .fixedx_heterogeneous_design(n, p = 5L, seed = 20261016L)
  x <- cbind(x, constant = rep(-7, n))

  result <- make_artificial_features(
    x = x,
    n_injected = 3L,
    type = "knockoff",
    random_state = 20261016L
  )

  # The matrix handed to create.fixed is standardized: centered columns,
  # unit sd for non-constant columns, zero column for the constant one.
  expect_identical(captured$sigma, 1)
  expect_equal(unname(colMeans(captured$x)), rep(0, ncol(x)), tolerance = 1e-12)
  expect_equal(unname(apply(captured$x[, 1L:5L, drop = FALSE], 2L, sd)),
               rep(1, 5L), tolerance = 1e-12)
  expect_true(all(captured$x[, "constant"] == 0))

  # Provenance records the standardization actually applied.
  prov <- result$artificial_provenance
  expect_equal(prov$actual_type, "knockoff")
  expect_false(prov$chunks$fallback)
  expect_true(isTRUE(prov$standardization$standardized))
  expect_equal(prov$standardization$means, colMeans(x))
  sds_applied <- apply(x, 2L, sd)
  sds_applied[["constant"]] <- 1  # zero-variance column left at scale 1
  expect_equal(unname(prov$standardization$sds), unname(sds_applied))
  expect_equal(prov$standardization$sds[["constant"]], 1)
  expect_equal(sd(x[, "constant"]), 0)

  # The injected knockoff block is the mocked constructor's output.
  xk_block <- result$x_augmented[, ncol(x) + seq_len(3L), drop = FALSE]
  # The shipped block is the mocked constructor's output RESCALED to the
  # standardized scale (create.fixed's unit-norm-scale return times the
  # standardized column norms).
  mock_xk <- (captured$x + 1) %*% diag(sqrt(colSums(captured$x^2)))
  matched <- apply(xk_block, 2L, function(col) {
    any(vapply(seq_len(ncol(mock_xk)), function(j) {
      isTRUE(all.equal(col, mock_xk[, j]))
    }, logical(1L)))
  })
  expect_true(all(matched))
})

test_that("provenance standardization field reflects the path actually taken", {
  skip_if_not_installed("knockoff")
  n <- 120L
  x <- .fixedx_heterogeneous_design(n, p = 5L, seed = 20261017L)

  result <- withCallingHandlers(
    make_artificial_features(
      x = x,
      n_injected = 3L,
      type = "knockoff",
      random_state = 20261017L
    ),
    warning = function(w) {
      if (grepl("falling back to random permutation", conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  prov <- result$artificial_provenance
  if (identical(prov$actual_type, "knockoff")) {
    expect_true(isTRUE(prov$standardization$standardized))
    expect_equal(prov$standardization$means, colMeans(x))
    expect_equal(prov$standardization$sds, apply(x, 2L, sd))
  } else {
    # Fallback path (e.g. constructor failure): no standardization applied.
    expect_identical(prov$actual_type, "random_permutation")
    expect_false(isTRUE(prov$standardization$standardized))
    expect_null(prov$standardization$means)
    expect_null(prov$standardization$sds)
  }
})

# ── Knockoff quality and calibration (probe-gated; validated on CI) ──────────

test_that("fixed-X knockoffs reproduce the standardized design's Gram matrix", {
  skip_if_not_installed("knockoff")
  skip_if_not(.fixedx_probe_ok,
              paste0("local knockoff build produces degenerate knockoffs; ",
                     "quality assertions are validated on CI with CRAN knockoff"))

  n <- 120L
  x <- .fixedx_heterogeneous_design(n, p = 8L, seed = 20261018L)

  result <- make_artificial_features(
    x = x,
    n_injected = ncol(x),
    type = "knockoff",
    random_state = 20261018L
  )
  expect_identical(result$artificial_provenance$actual_type, "knockoff")

  standardized <- stablr:::.standardize_columns(x)
  xk <- result$x_augmented[, ncol(x) + seq_len(ncol(x)), drop = FALSE]

  # Fixed-X knockoff defining identities on the standardized design:
  # t(Xk)Xk = t(X)X exactly, and D = t(X)X - t(X)Xk is a nonnegative diagonal
  # with s_j <= 2 * G_jj such that the joint Gram matrix stays PSD.
  gram_x  <- crossprod(standardized$x) / n
  gram_k  <- crossprod(xk) / n
  gram_xk <- crossprod(standardized$x, xk) / n
  expect_equal(gram_k, gram_x, tolerance = 1e-8)

  s_mat <- gram_x - gram_xk
  expect_equal(unname(s_mat - diag(diag(s_mat))),
               matrix(0, ncol(x), ncol(x)), tolerance = 1e-8)
  expect_gte(min(diag(s_mat)), 0)
  expect_lte(max(diag(s_mat) - 2 * diag(gram_x)), 1e-8)
  joint <- rbind(cbind(gram_x, gram_xk), cbind(t(gram_xk), gram_k))
  expect_gte(min(eigen(joint, symmetric = TRUE, only.values = TRUE)$values),
             -1e-8)
})

test_that("model-X knockoffs calibrate FDP+ end to end on a heterogeneous-scale design", {
  skip_if_not_installed("knockoff")

  # End-to-end calibration is asserted in the bootstrap-valid regime:
  # model-X knockoff rows are i.i.d. Gaussian draws and survive the row
  # subsampling of the bootstrap loop.  Fixed-X knockoffs are construction-
  # tested above (Gram identities); under sample_fraction < 1 their
  # exchangeability does not survive resampling, so they are not the right
  # subject for an end-to-end calibration assertion.
  n <- 150L
  p <- 12L
  x <- .fixedx_heterogeneous_design(n, p = p, seed = 20261019L)
  rownames(x) <- paste0("s", seq_len(n))
  # Plant signals on the large-scale features (f10-f12): after
  # standardization their standardized effect (beta * sd) is ~2-3.4 sigma,
  # so the test has power.  (On f1-f3, sd ~ 0.01-0.03 makes the planted
  # signal invisible and the test uninformative.)
  signal_idx <- c(10L, 11L, 12L)
  signals <- paste0("f", signal_idx)
  nulls <- setdiff(colnames(x), signals)
  beta <- numeric(p)
  beta[signal_idx] <- c(2.0, -1.5, 1.2)

  empirical_fdp <- numeric(3L)
  empirical_tpr <- numeric(3L)
  for (s in seq_len(3L)) {
    set.seed(2700L + s)
    y <- setNames(as.numeric(x %*% beta + rnorm(n, 0, 0.3)), rownames(x))
    fit <- stabl_fit(
      x = x,
      y = y,
      lambda_grid = data.frame(lambda = c(0.6, 0.4, 0.25)),
      n_bootstraps = 25L,
      artificial_type = "knockoff_equi",
      random_state = 2700L + s
    )
    selected <- get_feature_names_out(fit)
    empirical_fdp[[s]] <- if (length(selected)) {
      sum(selected %in% nulls) / length(selected)
    } else {
      0
    }
    empirical_tpr[[s]] <- sum(signals %in% selected) / length(signals)
  }

  # Acceptance: empirical FDP at the FDP+-chosen threshold at or below
  # nominal (0.1) across seeds, with the planted signals recovered.
  expect_lte(mean(empirical_fdp), 0.1)
  expect_gte(mean(empirical_tpr), 2 / 3)
})
