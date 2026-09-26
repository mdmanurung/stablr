# Numeric-contract regressions for the MVR S-matrix solver. Each block targets
# one numerical defect that previously produced a plausible-looking but wrong
# S matrix, or silently swallowed an infeasibility.

.mvr_body <- function(name) {
  paste(deparse(get(name, envir = asNamespace("stablr"))), collapse = " ")
}

.mvr_count_calls <- function(body_text, pattern) {
  regmatches(body_text, gregexpr(pattern, body_text))[[1]]
}


# ── A1: the PSD rescale must not return the zero matrix ──────────────────────

test_that("PSD rescale aborts typed instead of returning the zero matrix", {
  # 2 * Sigma has minimum eigenvalue 2e-7, below the 1e-5 feasibility
  # tolerance, so no gamma in (0, 1] is feasible and the binary search never
  # advances `lo` past its initial value of 0.
  Sigma <- diag(c(1e-7, 1, 1, 1), 4)
  S <- diag(c(1e-7, 1, 1, 1), 4)

  condition <- expect_error(
    stablr:::.scale_until_psd(Sigma, S, tol = 1e-5, num_iter = 10L),
    class = "stablr_mvr_infeasible"
  )
  expect_s3_class(condition, "stablr_numerical_infeasibility")
})

test_that("PSD rescale returns a strictly positive diagonal on the feasible path", {
  Sigma <- toeplitz(0.55^(0:7))
  S <- diag(0.4, 8)

  rescaled <- stablr:::.scale_until_psd(Sigma, S, tol = 1e-5, num_iter = 10L)

  expect_true(all(is.finite(diag(rescaled$S))))
  expect_true(all(diag(rescaled$S) > 0))
  expect_gte(
    min(eigen(
      2 * Sigma - rescaled$gamma * S,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    1e-5
  )
})

test_that("solve_mvr never reports a zero S-matrix as a solution", {
  # A1 end-to-end: this Sigma makes every coordinate infeasible, so the old
  # binary search returned S = 0 and knockoff::create.gaussian was handed
  # diag_s = 0 ("Knockoffs will have no power").
  Sigma <- diag(c(1e-7, 1, 1, 1, 1, 1), 6)

  expect_error(
    stablr:::solve_mvr(Sigma, num_iter = 50L, random_state = 7101L),
    class = "stablr_mvr_infeasible"
  )
})


# ── A2: Cholesky downdate radicand ───────────────────────────────────────────

test_that("Cholesky downdate aborts typed on a materially infeasible radicand", {
  condition <- expect_error(
    stablr:::.mvr_cholupdate(chol(diag(3)), c(2, 0, 0), add = FALSE),
    class = "stablr_mvr_infeasible"
  )
  expect_s3_class(condition, "stablr_numerical_infeasibility")
  expect_match(conditionMessage(condition), "downdate")
})

test_that("Cholesky downdate radicand guard clamps round-off at the boundary", {
  # One ulp past feasibility: the radicand is negative only at double
  # precision, which is a boundary case rather than a magnitude violation.
  one_ulp_past <- 1 + .Machine$double.eps
  expect_gt(one_ulp_past, 1)
  expect_lt(one_ulp_past^2, 1 + 4 * .Machine$double.eps)
  boundary <- expect_error(
    stablr:::.mvr_cholupdate(chol(diag(3)), c(one_ulp_past, 0, 0),
                             add = FALSE),
    class = "stablr_mvr_infeasible"
  )
  expect_match(conditionMessage(boundary), "boundary")

  expect_error(
    stablr:::.mvr_cholupdate(chol(diag(3)), c(1.5, 0, 0), add = FALSE),
    class = "stablr_mvr_infeasible"
  )
})

test_that("downdate radicand helper returns the exact radius when feasible", {
  expect_equal(stablr:::.mvr_cholupdate_radicand(5, 3, add = FALSE, 1L), 4)
  expect_equal(stablr:::.mvr_cholupdate_radicand(3, 4, add = TRUE, 1L), 5)
  expect_equal(stablr:::.mvr_cholupdate_radicand(1, 0, add = FALSE, 1L), 1)
  expect_error(
    stablr:::.mvr_cholupdate_radicand(0, 0, add = FALSE, 1L),
    class = "stablr_mvr_infeasible"
  )
  expect_error(
    stablr:::.mvr_cholupdate_radicand(1, NaN, add = TRUE, 1L),
    class = "stablr_mvr_infeasible"
  )
})


# ── A3: divergence must not be reported as convergence ───────────────────────

test_that("non-finite MVR loss aborts typed rather than ending the descent", {
  # min(diag(S)) == 0 passes the .mvr_loss sign test but makes 1 / diag(S)
  # infinite, so the loss is Inf without any eigenvalue being negative.
  expect_error(
    stablr:::.mvr_loss_checked(diag(2), diag(c(0, 1))),
    class = "stablr_mvr_infeasible"
  )
  expect_error(
    stablr:::.mvr_loss_checked(diag(2), diag(c(-1, 1))),
    class = "stablr_mvr_infeasible"
  )
  expect_equal(stablr:::.mvr_loss_checked(diag(2), diag(c(1, 1))), 4)
})

test_that("solve_mvr aborts typed when the loss diverges", {
  testthat::local_mocked_bindings(
    .mvr_loss = function(...) Inf,
    .package = "stablr"
  )
  expect_error(
    stablr:::solve_mvr(toeplitz(0.5^(0:5)), random_state = 7102L),
    class = "stablr_mvr_infeasible"
  )
})


# ── A6 / A7: dead structure must not come back ───────────────────────────────

test_that("PSD rescale has no dead error handler around a logical", {
  body_text <- .mvr_body(".scale_until_psd")
  expect_false(grepl("tryCatch", body_text, fixed = TRUE))
  expect_false(grepl("stablr_numerical_infeasibility", body_text, fixed = TRUE))
})

test_that("the MVR quadratic solves its roots exactly once", {
  calls <- c(
    .mvr_count_calls(.mvr_body(".solve_mvr_quadratic"), "polyroot\\("),
    .mvr_count_calls(.mvr_body(".mvr_real_roots"), "polyroot\\(")
  )
  expect_length(calls, 1L)
})


# ── A8: the feasible interval may invert under the fixed slack ───────────────

test_that("an inverted feasible interval aborts typed instead of degrading", {
  slack <- .Machine$double.eps * 1e4
  cd <- 1e13
  sj <- 1e-13
  expect_true(-sj + slack > 1 / cd - slack)

  condition <- expect_error(
    stablr:::.solve_mvr_quadratic(cn = -1e-26, cd = cd, sj = sj),
    class = "stablr_mvr_infeasible"
  )
  expect_s3_class(condition, "stablr_numerical_infeasibility")
  expect_match(conditionMessage(condition), "inverted")
})

test_that("a representable feasible interval is not mistaken for an inversion", {
  expect_silent(
    stablr:::.solve_mvr_quadratic(cn = -0.5, cd = 1, sj = 0.3)
  )
  expect_error(
    stablr:::.solve_mvr_quadratic(cn = -0.5, cd = 0, sj = 0.3),
    class = "stablr_mvr_infeasible"
  )
})


# ── A9: the real-root filter must be relative to the root magnitude ──────────

test_that("the real-root filter rejects genuinely complex roots", {
  # z^2 + 1e-30 has roots 0 +/- 1e-15i. A filter with an absolute 1e-14 floor
  # admits both as real roots.
  expect_length(stablr:::.mvr_real_roots(1e-30, 0, 1), 0L)
  expect_length(stablr:::.mvr_real_roots(1, 0, 1), 0L)
})

test_that("the real-root filter keeps real roots carrying round-off residue", {
  # z^2 + 2.3 z - 0.5 has two real roots; polyroot returns them with a
  # ~1e-23 imaginary residue that must not be read as complexity.
  roots <- stablr:::.mvr_real_roots(-0.5, 2.3, -0.5)
  expect_length(roots, 2L)
  expect_equal(sort(roots), sort(Re(polyroot(c(-0.5, 2.3, -0.5)))))
})


# ── A12: singular solves must not return finite garbage ──────────────────────

test_that("a singular Sigma is rejected typed instead of silently shifted", {
  expect_error(
    stablr:::solve_mvr(diag(6) * 0),
    class = "stablr_mvr_infeasible"
  )
  expect_error(
    stablr:::solve_mvr(matrix(0, 4, 4)),
    class = "stablr_mvr_infeasible"
  )
  expect_error(
    stablr:::solve_mvr(matrix(NA_real_, 3, 3)),
    class = "stablr_mvr_infeasible"
  )
})

test_that("coordinate solves recover the defining quantities of a PD factor", {
  set.seed(7103L)
  M <- crossprod(matrix(rnorm(36), 6, 6)) + diag(6)
  R_upper <- chol(M)
  M_inverse <- solve(M)

  for (j in seq_len(6L)) {
    solved <- stablr:::.mvr_coordinate_quadratic(R_upper, j)
    expect_gt(solved$cd, 0)
    expect_lt(solved$cn, 0)
    expect_equal(solved$cd, M_inverse[j, j], tolerance = 1e-10)
  }
})

test_that("coordinate solves abort typed on a singular or NaN factor", {
  good <- chol(diag(3) + 0.5)
  zero_pivot <- good
  zero_pivot[2L, 2L] <- 0
  expect_error(
    stablr:::.mvr_coordinate_quadratic(zero_pivot, 1L),
    class = "stablr_mvr_infeasible"
  )
  expect_error(
    stablr:::.mvr_coordinate_quadratic(good * Inf, 1L),
    class = "stablr_mvr_infeasible"
  )
  expect_error(
    stablr:::.mvr_coordinate_quadratic(matrix(0, 3, 3), 1L),
    class = "stablr_mvr_infeasible"
  )
})

test_that("the Cholesky factor is verified against the matrix it claims", {
  M <- crossprod(matrix(c(1, 2, 0, 0, 3, 1, 1, 1, 4), 3, 3)) + diag(3)

  expect_silent(stablr:::.mvr_verify_cholesky_factor(chol(M), M))
  expect_error(
    stablr:::.mvr_verify_cholesky_factor(2 * chol(M), M),
    class = "stablr_mvr_infeasible"
  )
  singular <- chol(M)
  singular[3L, 3L] <- 0
  expect_error(
    stablr:::.mvr_verify_cholesky_factor(singular, M),
    class = "stablr_mvr_infeasible"
  )
})

test_that("a drifted Cholesky factor aborts the solve instead of yielding S", {
  # base::forwardsolve / base::backsolve raise only on an exactly-zero pivot,
  # so a factor that has drifted from 2 * Sigma - S still returns finite cd and
  # cn. Nothing in the old solver noticed.
  testthat::local_mocked_bindings(
    .mvr_cholupdate = function(R, x, add) R * 1.5,
    .package = "stablr"
  )
  expect_error(
    stablr:::solve_mvr(toeplitz(0.55^(0:7)), num_iter = 3L,
                       random_state = 7104L),
    class = "stablr_mvr_infeasible"
  )
})
