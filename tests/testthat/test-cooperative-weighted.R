# F3 — weighted cooperative-backend coverage (pairs with B1/B2 fixes).
#
# stablr's public API never exposes case weights or offsets to the vendored
# multiview engine, so these tests call the internal backend directly
# (`:::` is intentional). Two regression guards:
#   * a weighted fit must complete cleanly (no NaN, monotone dev.ratio);
#   * constant weights c must give the same fit as unit weights — this is
#     what breaks if the warm-start residual mishandles weights (B2).

# Helper: run expr, muffling and collecting warnings (binomial initialize
# warns about non-integer successes when weights are fractional).
.muffled <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = msgs)
}

test_that("cooperative backend accepts non-unit case weights (binomial, rho > 0)", {
  set.seed(2026)
  n <- 80L; p <- 15L
  x_a <- matrix(rnorm(n * p), nrow = n, dimnames = list(NULL, paste0("a", seq_len(p))))
  x_b <- matrix(rnorm(n * p), nrow = n, dimnames = list(NULL, paste0("b", seq_len(p))))
  eta <- 0.8 * x_a[, 1] - 0.6 * x_b[, 2] + 0.4 * x_a[, 5]
  y <- rbinom(n, 1, stats::plogis(eta))
  # fixed non-unit weight pattern: heavy heads, light tails
  w <- rep(1, n); w[1:10] <- 20; w[71:80] <- 0.25

  out <- .muffled(
    stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y,
                                      family = "binomial", rho = 0.3, weights = w)
  )
  fit <- out$value
  expect_s3_class(fit, "multiview")
  beta <- as.matrix(fit$beta)
  expect_false(anyNA(beta))
  expect_false(anyNA(fit$dev.ratio))
  # objective proxy: explained deviance must be non-decreasing along the path
  expect_true(all(diff(fit$dev.ratio) >= -1e-8))
  expect_true(sum(beta[, ncol(beta)] != 0) > 0)

  # stress: extreme weight contrast + aggressive lambdas must not produce NaN
  w_extreme <- rep(1, n); w_extreme[1:5] <- 5000; w_extreme[76:80] <- 1e-4
  out2 <- .muffled(
    stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y,
                                      family = "binomial", rho = 0.3,
                                      weights = w_extreme,
                                      lambda = c(0.05, 0.01, 0.002, 0.0005))
  )
  expect_s3_class(out2$value, "multiview")
  expect_false(anyNA(as.matrix(out2$value$beta)))
})

test_that("constant case weights reproduce the unit-weight fit (B2 regression guard)", {
  set.seed(2026)
  n <- 80L; p <- 15L
  x_a <- matrix(rnorm(n * p), nrow = n, dimnames = list(NULL, paste0("a", seq_len(p))))
  x_b <- matrix(rnorm(n * p), nrow = n, dimnames = list(NULL, paste0("b", seq_len(p))))
  eta <- 0.8 * x_a[, 1] - 0.6 * x_b[, 2] + 0.4 * x_a[, 5]
  y <- rbinom(n, 1, stats::plogis(eta))
  y_g <- 0.9 * x_a[, 1] - 0.7 * x_b[, 2] + rnorm(n, sd = 0.3)

  # binomial, rho > 0: exercises multiview.fit / multiview.path (the B2 code)
  fit_u <- .muffled(stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y,
                      family = "binomial", rho = 0.3))$value
  fit_c <- .muffled(stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y,
                      family = "binomial", rho = 0.3, weights = rep(2.5, n)))$value
  expect_identical(dim(as.matrix(fit_u$beta)), dim(as.matrix(fit_c$beta)))
  b_u <- as.matrix(fit_u$beta); b_c <- as.matrix(fit_c$beta)
  expect_true(max(abs(b_u - b_c)) < 1e-3)
  # same selected set at the densest end of the path
  expect_equal(which(b_u[, ncol(b_u)] != 0), which(b_c[, ncol(b_c)] != 0))

  # gaussian, rho = 0 (early fusion, direct glmnet branch): must be exact
  fit_ug <- stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y_g,
                                              family = "gaussian", rho = 0)
  fit_cg <- stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y_g,
                                              family = "gaussian", rho = 0,
                                              weights = rep(2.5, n))
  expect_equal(as.matrix(fit_cg$beta), as.matrix(fit_ug$beta), tolerance = 1e-7)
})

test_that("weighted fit with a non-NULL offset completes cleanly (B2 residual path)", {
  set.seed(2026)
  n <- 80L; p <- 15L
  x_a <- matrix(rnorm(n * p), nrow = n)
  x_b <- matrix(rnorm(n * p), nrow = n)
  eta <- 0.8 * x_a[, 1] - 0.6 * x_b[, 2]
  y <- rbinom(n, 1, stats::plogis(eta))
  w <- rep(1, n); w[1:10] <- 5; w[71:80] <- 0.5
  offset <- rnorm(n)

  out <- .muffled(
    stablr:::.cooperative_backend_fit(list(a = x_a, b = x_b), y,
                                      family = "binomial", rho = 0.3,
                                      weights = w, offset = offset)
  )
  fit <- out$value
  expect_s3_class(fit, "multiview")
  expect_false(anyNA(as.matrix(fit$beta)))
  expect_false(anyNA(fit$dev.ratio))
  expect_true(all(diff(fit$dev.ratio) >= -1e-8))
})

test_that("build_predmat subsets per-fold offsets correctly (B1 regression guard)", {
  set.seed(2026)
  n <- 80L; p <- 15L
  x_a <- matrix(rnorm(n * p), nrow = n)
  x_b <- matrix(rnorm(n * p), nrow = n)
  y <- 0.9 * x_a[, 1] - 0.7 * x_b[, 2] + rnorm(n, sd = 0.3)
  offset <- rnorm(n)
  foldid <- rep(1:4, length.out = n)

  # fold fits trained with the fold-specific offset, as cv.multiview does
  fits <- lapply(1:4, function(i) {
    k <- foldid != i
    stablr:::.mv_multiview(list(a = x_a[k, , drop = FALSE], b = x_b[k, , drop = FALSE]),
                           y[k], rho = 0.3, offset = offset[k])
  })
  lambda <- fits[[1]]$mvlambda

  # direct call: must not error (pre-B1 it referenced undefined `offset_sub`)
  predmat <- stablr:::build_predmat(fits, lambda, list(a = x_a, b = x_b), offset, foldid,
                                    alignment = "lambda", y = y,
                                    weights = rep(1, n), grouped = TRUE,
                                    type.measure = "mse", family = stats::gaussian())
  expect_false(anyNA(predmat))

  # per-fold predictions must match a manual predict() with the fold offset
  for (i in 1:4) {
    k <- foldid == i
    manual <- predict(fits[[i]], list(a = x_a[k, , drop = FALSE], b = x_b[k, , drop = FALSE]),
                      s = lambda, newoffset = offset[k])
    expect_equal(unname(predmat[k, seq(ncol(manual)), drop = FALSE]),
                 unname(as.matrix(manual)),
                 info = paste("fold", i))
  }

  # end-to-end: the cv.multiview caller with an offset must also be clean
  cv_fit <- .muffled(
    stablr:::.cooperative_backend_cv(list(a = x_a, b = x_b), y, family = "gaussian",
                                     rho = 0.3, foldid = foldid, offset = offset)
  )
  expect_s3_class(cv_fit$value, "cv.multiview")
  expect_false(anyNA(cv_fit$value$cvm))
})
