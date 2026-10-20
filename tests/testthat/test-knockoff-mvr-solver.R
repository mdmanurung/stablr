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

test_that(".scale_until_psd warns when no feasible gamma exists (L9 regression)", {
  expect_warning(
    res <- stablr:::.scale_until_psd(diag(2), diag(1e6, 2),
                                     tol = 1e-5, num_iter = 10L),
    "degenerate"
  )
  expect_identical(res$gamma, 0)
})
