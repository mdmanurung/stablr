test_that("MVR solver returns a feasible diagonal S matrix", {
  Sigma <- toeplitz(0.55^(0:7))

  s_diag <- stablr:::solve_mvr(
    Sigma,
    num_iter = 50L,
    random_state = 4101L
  )
  S <- diag(s_diag)

  expect_length(s_diag, nrow(Sigma))
  expect_true(all(is.finite(s_diag)))
  expect_true(all(s_diag > 0))
  expect_gte(
    min(eigen(2 * Sigma - S, symmetric = TRUE, only.values = TRUE)$values),
    0
  )
})

test_that("MVR solver does not worsen its initialization objective", {
  Sigma <- toeplitz(0.7^(0:9))
  initial_S <- diag(
    min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values),
    nrow(Sigma)
  )

  s_diag <- stablr:::solve_mvr(
    Sigma,
    num_iter = 50L,
    random_state = 4102L
  )
  solved_S <- diag(s_diag)

  expect_lte(
    stablr:::.mvr_loss(Sigma, solved_S),
    stablr:::.mvr_loss(Sigma, initial_S) + 1e-10
  )
})

test_that("MVR solver is deterministic for a fixed random state", {
  Sigma <- toeplitz(0.6^(0:6))

  first <- stablr:::solve_mvr(
    Sigma,
    num_iter = 30L,
    random_state = 4103L
  )
  second <- stablr:::solve_mvr(
    Sigma,
    num_iter = 30L,
    random_state = 4103L
  )

  expect_identical(first, second)
})

test_that("seeded MVR solving preserves the caller RNG state", {
  Sigma <- toeplitz(0.5^(0:5))
  set.seed(4104L)
  caller_seed <- .Random.seed

  invisible(stablr:::solve_mvr(
    Sigma,
    num_iter = 20L,
    random_state = 991L
  ))

  expect_identical(.Random.seed, caller_seed)
})

# ── A4: the per-coordinate degradation handler must be scoped ────────────────

test_that("the infeasible-quadratic condition is not reachable as a warning", {
  condition <- tryCatch(
    stablr:::.mvr_abort_infeasible_quadratic(cn = -1, cd = 2, sj = 1),
    error = identity
  )

  expect_s3_class(condition, "stablr_mvr_infeasible_quadratic")
  expect_s3_class(condition, "stablr_numerical_infeasibility")
  expect_s3_class(condition, "error")
  expect_false(inherits(condition, "warning"))
  expect_identical(condition$cn, -1)
  expect_identical(condition$sj, 1)
})

test_that("no feasible root raises only the dedicated condition", {
  expect_error(
    stablr:::.solve_mvr_quadratic(cn = -1e-26, cd = 1e8, sj = 1e-8),
    class = "stablr_mvr_infeasible_quadratic"
  )
})

test_that("an unrelated warning inside the quadratic reaches the caller", {
  # The old handler was keyed on "warning", so it caught anything raised
  # anywhere inside .solve_mvr_quadratic and silently replaced the result with
  # a 1%-of-bound heuristic.
  real_roots <- stablr:::.mvr_real_roots
  root_calls <- 0L
  testthat::local_mocked_bindings(
    .mvr_real_roots = function(...) {
      root_calls <<- root_calls + 1L
      warning("unrelated upstream warning")
      real_roots(...)
    },
    .package = "stablr"
  )

  raised <- character()
  solved <- withCallingHandlers(
    stablr:::solve_mvr(toeplitz(0.55^(0:7)), random_state = 7105L),
    warning = function(w) {
      raised <<- c(raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_gt(root_calls, 0L)
  # Every warning raised inside the quadratic reached the caller; none was
  # consumed and replaced by a substitute delta.
  expect_identical(
    sum(grepl("unrelated upstream warning", raised, fixed = TRUE)),
    root_calls
  )
  expect_false(any(grepl("not the MVR minimiser", raised, fixed = TRUE)))
  expect_true(all(is.finite(solved)))
})

test_that("a degraded MVR coordinate is reported to the user", {
  testthat::local_mocked_bindings(
    .solve_mvr_quadratic = function(...) {
      stablr:::.mvr_abort_infeasible_quadratic(cn = -1, cd = 2, sj = 1)
    },
    .package = "stablr"
  )

  expect_warning(
    solved <- stablr:::solve_mvr(
      toeplitz(0.55^(0:7)),
      num_iter = 2L,
      random_state = 7106L
    ),
    "not the MVR minimiser"
  )
  expect_true(all(is.finite(solved)))
})
