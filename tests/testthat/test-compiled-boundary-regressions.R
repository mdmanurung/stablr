# Phase 13 PR-1 regression tests: compiled-boundary fail-closed validation
# in the vendored cooperative engine (R/cooperative-multiview.path.R,
# R/cooperative-get_start.R).
#
# References:
# - docs/audit/2026-10-17-adversarial/remediation-matrix.md (P0/P1 table)
# - docs/audit/2026-10-17-adversarial/phase12-synthesis.md sections 4 and 7
# - findings.md F-011 (P0), F-012 (P1), F-013 (P1)
#
# These tests exercise the boundary between R and the vendored glmnetpp-based
# C++ engine (wls_exp / spwls_exp) — the coverage gap identified by the
# Phase 9 test-suite audit.

test_that("F-011: sparse x without xm/xs attributes fails closed, not OOB", {
  # Regression: the sparse branch passed as.double(attr(x, "xm")) straight
  # to spwls_exp; a missing attribute became double(0) and the engine read
  # xm[0..ni-1] out of bounds (SIGABRT observed; silent heap OOB under
  # -DNDEBUG). Must now stop with an informative error.
  set.seed(11)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] * 1.5 + rnorm(n)
  w <- rep(1, n)
  xs_raw <- Matrix::Matrix(x, sparse = TRUE)
  expect_null(attr(xs_raw, "xm"))
  expect_null(attr(xs_raw, "xs"))
  expect_error(
    stablr:::elnet.fit(xs_raw, y, w, lambda = 0.05, intercept = FALSE),
    regexp = "standardization attributes"
  )
})

test_that("F-011: sparse x with wrong-length attributes fails closed", {
  set.seed(11)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] * 1.5 + rnorm(n)
  w <- rep(1, n)
  bad <- Matrix::Matrix(x, sparse = TRUE)
  attr(bad, "xm") <- colMeans(x)[seq_len(p - 1L)]
  attr(bad, "xs") <- apply(x, 2, sd)
  expect_error(
    stablr:::elnet.fit(bad, y, w, lambda = 0.05, intercept = FALSE),
    regexp = "standardization attributes"
  )
})

test_that("F-011: sparse x with non-positive xs fails closed", {
  # xs = 0 would divide by zero inside the engine; legitimate callers
  # (multiview preparation) give constant columns xs = 1 and exclude them.
  set.seed(11)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] * 1.5 + rnorm(n)
  w <- rep(1, n)
  bad <- Matrix::Matrix(x, sparse = TRUE)
  attr(bad, "xm") <- colMeans(x)
  attr(bad, "xs") <- c(rep(1, p - 1L), 0)
  expect_error(
    stablr:::elnet.fit(bad, y, w, lambda = 0.05, intercept = FALSE),
    regexp = "finite and positive"
  )
})

test_that("F-011: well-formed sparse x matches the dense standardized fit", {
  # Guard the fix against over-blocking: a properly attributed sparse
  # matrix must still fit, and must agree with the dense branch on the
  # standardized scale (agreement is at machine precision, ~1e-17).
  set.seed(11)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] * 1.5 + rnorm(n)
  w <- rep(1, n)
  xm <- colMeans(x)
  xs <- apply(x, 2, sd)
  sp <- Matrix::Matrix(x, sparse = TRUE)
  attr(sp, "xm") <- xm
  attr(sp, "xs") <- xs
  fit_sp <- stablr:::elnet.fit(sp, y, w, lambda = 0.05, intercept = FALSE)
  fit_dn <- stablr:::elnet.fit(scale(x, center = xm, scale = xs), y, w,
                               lambda = 0.05, intercept = FALSE)
  expect_identical(fit_sp$jerr, 0L)
  expect_identical(fit_dn$jerr, 0L)
  expect_true(all(is.finite(as.numeric(fit_sp$beta))))
  expect_lt(max(abs(as.numeric(fit_sp$beta) - as.numeric(fit_dn$beta))), 1e-10)
})

test_that("F-012: NA in dense x fails closed (no silent column drop)", {
  set.seed(12)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] + rnorm(n)
  w <- rep(1, n)
  x_na <- x
  x_na[1L, 1L] <- NA
  expect_error(
    stablr:::elnet.fit(x_na, y, w, lambda = 0.05, intercept = FALSE),
    regexp = "'x' contains NA"
  )
})

test_that("F-012: NA in sparse x fails closed", {
  set.seed(12)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] + rnorm(n)
  w <- rep(1, n)
  sp <- Matrix::Matrix(x, sparse = TRUE)
  attr(sp, "xm") <- colMeans(x)
  attr(sp, "xs") <- apply(x, 2, sd)
  sp[1L, 1L] <- NA
  expect_error(
    stablr:::elnet.fit(sp, y, w, lambda = 0.05, intercept = FALSE),
    regexp = "'x' contains NA"
  )
})

test_that("F-012: NA in y and in weights fails closed", {
  set.seed(12)
  n <- 30L; p <- 5L
  x <- matrix(rnorm(n * p), n, p)
  y <- x[, 1L] + rnorm(n)
  w <- rep(1, n)
  y_na <- y
  y_na[5L] <- NA
  w_na <- w
  w_na[7L] <- NA
  expect_error(
    stablr:::elnet.fit(x, y_na, w, lambda = 0.05, intercept = FALSE),
    regexp = "'y' contains NA"
  )
  expect_error(
    stablr:::elnet.fit(x, y, w_na, lambda = 0.05, intercept = FALSE),
    regexp = "'weights' contains NA"
  )
})

test_that("F-013: maxit non-convergence warns and returns the jerr code", {
  # Regression: elnet.fit returned list(jerr = -1) silently on maxit
  # non-convergence; every caller propagated it bare (no warning anywhere
  # on the mvglmnet.fit / cooperative paths). Must now warn at the
  # boundary with a code-specific message and keep the return contract.
  set.seed(13)
  n <- 50L; p <- 4L
  x <- matrix(rnorm(n * p), n, p)
  y <- x %*% c(2, -1.5, 1, 0.5) + rnorm(n)
  w <- rep(1, n)
  expect_warning(
    fit <- stablr:::elnet.fit(x, y, w, lambda = 0.001, maxit = 1L,
                              intercept = FALSE),
    regexp = "Convergence not reached"
  )
  expect_identical(fit$jerr, -1L)
  expect_identical(names(fit), "jerr")
})

test_that("F-013: jerr mappers give distinct messages per code range", {
  # The old mappers collapsed every negative code into "convergence" and
  # every positive code below 7777 into "memory allocation error", which
  # is wrong for 7777/8888/9999/10000 (glmnetpp util::exceptions.hpp).
  msg_for <- function(code) stablr:::jerr.glmnetfit(code, 100000)$msg

  # non-fatal ranges are distinct
  expect_match(msg_for(-1L), "Convergence not reached", fixed = TRUE)
  expect_match(msg_for(-2L), "Convergence not reached", fixed = TRUE)
  expect_match(msg_for(-10001L), "active-set", fixed = TRUE)
  expect_match(msg_for(-20001L), "variance", fixed = TRUE)

  # fatal codes are distinct and no longer mislabeled as memory errors
  expect_match(msg_for(1L), "Memory allocation", fixed = TRUE)
  expect_match(msg_for(7777L), "excluded", fixed = TRUE)
  expect_match(msg_for(8888L), "Negative response", fixed = TRUE)
  expect_match(msg_for(9999L), "weight", fixed = TRUE)
  expect_match(msg_for(10000L), "penalty", fixed = TRUE)
  expect_match(msg_for(10001L), "Internal engine error", fixed = TRUE)
  expect_match(msg_for(12345L), "Unknown engine error code", fixed = TRUE)

  # fatal/non-fatal classification
  expect_false(stablr:::jerr.glmnetfit(-1L, 100000)$fatal)
  expect_true(stablr:::jerr.glmnetfit(7777L, 100000)$fatal)

  # multiview mapper adds lambda context on non-fatal codes
  mv <- stablr:::jerr.multiview(-2L, 100000, k = 2)
  expect_match(mv$msg, "Lambda 2", fixed = TRUE)
  expect_match(mv$msg, "Solutions for larger lambdas returned", fixed = TRUE)
  expect_false(mv$fatal)
})
