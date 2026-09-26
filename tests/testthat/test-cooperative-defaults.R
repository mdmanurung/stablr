# G1: `cooperative_fusion = TRUE` must not silently degenerate to early fusion.
#
# The defect: `rho = NULL` resolved to `0`, `unique(as.numeric(0))` left a
# single rho candidate, so `.cooperative_select_cv()` fitted `cv.multiview` once
# and returned it, `.mv_multiview()` took the plain-`glmnet` early-fusion branch,
# and the fit still reported `cooperative_rho = 0` with
# `cooperative_selection = "cv"`.  A plain `cooperative_fusion = TRUE` call
# therefore ran no cooperative learning at all while claiming to.
#
# Fix: the `NULL` default now resolves to a strictly positive multi-candidate
# ladder, so the documented default can no longer be a no-agreement-penalty
# sweep.  The resolved ladder is what the returned fit reports as `rho_grid`,
# and it is what `get_cooperative_diagnostics()` enumerates, so the agreement
# penalties actually swept are visible in the returned diagnostics.

.cf_fixture <- function(n = 24L, seed = 71L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  list(
    x_a = matrix(stats::rnorm(n * 4L), n, 4L,
                 dimnames = list(ids, paste0("a", 1:4))),
    x_b = matrix(stats::rnorm(n * 4L), n, 4L,
                 dimnames = list(ids, paste0("b", 1:4))),
    y = stats::setNames(stats::rnorm(n), ids)
  )
}

.cf_fit <- function(d, ...) {
  stabl_multiomic_train_validate(
    x_train_list = list(omic_a = d$x_a, omic_b = d$x_b),
    y_train = d$y,
    lambda_grid = data.frame(lambda = c(0.2, 0.1)),
    artificial_type = NULL,
    hard_threshold = 0.3,
    n_bootstraps = 2L,
    family = "gaussian",
    random_state = 9L,
    cooperative_fusion = TRUE,
    cooperation_selection = "cv",
    cooperation_nfolds = 3L,
    ...
  )
}

test_that("G1 the NULL rho default resolves to a strictly positive, multi-candidate ladder", {
  d <- .cf_fixture()

  args <- stablr:::.normalize_cooperative_multiomic_args(
    x_list = list(omic_a = d$x_a, omic_b = d$x_b),
    family = "gaussian",
    rho = NULL,
    cooperative_selection = "cv",
    cooperation_nfolds = 3L
  )

  expect_gt(length(args$rho), 1L)
  expect_true(all(args$rho > 0))
  expect_identical(args$rho, sort(unique(args$rho)))
  expect_true(isTRUE(args$rho_defaulted))
})

test_that("G1 a caller-supplied rho grid is passed through untouched", {
  d <- .cf_fixture()
  base <- list(
    x_list = list(omic_a = d$x_a, omic_b = d$x_b),
    family = "gaussian",
    cooperative_selection = "cv",
    cooperation_nfolds = 3L
  )

  supplied <- do.call(
    stablr:::.normalize_cooperative_multiomic_args,
    c(list(rho = c(0.3, 0.1, 0.3)), base)
  )
  expect_identical(supplied$rho, c(0.3, 0.1))
  expect_false(isTRUE(supplied$rho_defaulted))

  # rho = 0 stays a legal, deliberate request for the early-fusion baseline.
  zero <- do.call(
    stablr:::.normalize_cooperative_multiomic_args,
    c(list(rho = 0), base)
  )
  expect_identical(zero$rho, 0)
  expect_false(isTRUE(zero$rho_defaulted))
})

test_that("G1 cooperative_fusion = TRUE with no rho reports a real agreement penalty", {
  d <- .cf_fixture()

  fit <- .cf_fit(d)
  cf <- fit$cooperative_fusion

  expect_gt(length(cf$rho_grid), 1L)
  expect_true(all(cf$rho_grid > 0))
  expect_true(cf$rho > 0)
  expect_identical(cf$rho, cf$rho_grid[[which(cf$diagnostics$selected)]])

  # The returned diagnostics enumerate exactly the ladder that was swept, so the
  # agreement penalties behind the reported fit are readable by a consumer.
  diagnostics <- get_cooperative_diagnostics(fit)
  expect_setequal(unique(diagnostics$rho), cf$rho_grid)
  expect_equal(sum(diagnostics$selected), 1L)
  expect_true(all(diagnostics$rho > 0))
})

test_that("G1 the default cooperative fit is not the plain early-fusion branch", {
  d <- .cf_fixture()

  fit <- .cf_fit(d, early_fusion = TRUE)
  cf <- fit$cooperative_fusion

  expect_true(inherits(cf$fit, "cv.multiview"))
  expect_gt(length(cf$selected_features), 0L)
  # Early fusion is the rho = 0 branch; the default sweep must not reproduce it.
  early_union <- sort(unique(unlist(fit$early_fusion$selected_features, use.names = FALSE)))
  cooperative_union <- sort(unique(unlist(cf$selected_features, use.names = FALSE)))
  expect_false(
    identical(cooperative_union, early_union) && length(early_union) > 0L,
    info = "the default cooperative ladder must not collapse onto early fusion"
  )
})

test_that("G1 an explicit zero-only ladder still reports rho = 0 rather than claiming cooperation", {
  d <- .cf_fixture()

  fit <- .cf_fit(d, rho = 0)
  cf <- fit$cooperative_fusion

  # Deliberate early-fusion baseline: the report must not dress it up.
  expect_identical(cf$rho, 0)
  expect_identical(cf$rho_grid, 0)
  expect_true(all(cf$diagnostics$rho == 0))
})
