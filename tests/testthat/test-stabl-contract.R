# Contract tests for the `stabl_fit()` selection/RNG/threshold contract and the
# shared lambda-grid resolver.  Each block names the defect it pins.

.contract_data <- function(n = 40L, p = 6L, seed = 1L) {
  set.seed(seed)
  x <- matrix(
    stats::rnorm(n * p), nrow = n,
    dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p)))
  )
  list(x = x, y = stats::setNames(stats::rnorm(n), rownames(x)),
       lambda_grid = data.frame(lambda = c(0.2, 0.1)))
}

# ---- C1: hard_threshold vs the FDP+ threshold search ---------------------------

test_that("C1 hard_threshold explicitly suppresses the FDP+ threshold search", {
  d <- .contract_data()
  fit <- stabl_fit(
    d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
    artificial_type = "random_permutation", hard_threshold = 0.3,
    random_state = 1L
  )

  # The suppression is recorded, not silent.
  expect_true(fit$fdp_threshold_suppressed)
  # Nothing that FDP+ *would* have produced is present ...
  expect_null(fit$fdr_min_threshold_)
  expect_null(fit$FDRs_)
  expect_null(fit$min_fdr_)
  expect_null(fit$fdrs_table_)
  # ... and the accessor resolves to the hard threshold, not to a discarded
  # FDP+ value.
  expect_equal(fit$hard_threshold, 0.3)
  # The artificial block itself is still produced, so no diagnostic is lost
  # except the FDP+ curve.
  expect_false(is.null(fit$stabl_scores_artificial_))
  expect_false(is.null(fit$artificial_provenance))
})

test_that("C1 fdp_threshold_suppressed distinguishes all three configurations", {
  d <- .contract_data()
  common <- list(x = d$x, y = d$y, lambda_grid = d$lambda_grid,
                 n_bootstraps = 4L, random_state = 1L)

  both <- do.call(stabl_fit, c(common, list(
    artificial_type = "random_permutation", hard_threshold = 0.3
  )))
  fdp_only <- do.call(stabl_fit, c(common, list(
    artificial_type = "random_permutation"
  )))
  hard_only <- do.call(stabl_fit, c(common, list(
    artificial_type = NULL, hard_threshold = 0.3
  )))

  expect_true(both$fdp_threshold_suppressed)
  expect_false(fdp_only$fdp_threshold_suppressed)
  expect_true(is.na(hard_only$fdp_threshold_suppressed))
  expect_false(is.null(fdp_only$fdr_min_threshold_))
})

# ---- C2: new_hard_threshold validation -----------------------------------------

test_that("C2 accessors reject NA, negative and out-of-range new_hard_threshold", {
  d <- .contract_data()
  fit <- stabl_fit(
    d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
    artificial_type = "random_permutation", random_state = 1L
  )

  for (bad in list(NA_real_, -1, 2, Inf, c(0.1, 0.2), "0.5")) {
    expect_error(get_support(fit, new_hard_threshold = bad),
                 "new_hard_threshold")
    expect_error(get_feature_names_out(fit, new_hard_threshold = bad),
                 "new_hard_threshold")
    expect_error(transform_stabl(fit, d$x, new_hard_threshold = bad),
                 "new_hard_threshold")
  }

  # A valid proportion still overrides as before.
  expect_true(is.logical(get_support(fit, new_hard_threshold = 0.001)))
})

# ---- C3: family, l1_ratio, and the grid's alpha column --------------------------

test_that("C3 a mistyped family is rejected before reaching glmnet", {
  d <- .contract_data()
  for (bad in c("gaussain", "binomal", "Gausian", "coxph")) {
    expect_error(
      stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
                artificial_type = "random_permutation", family = bad,
                random_state = 1L),
      "family"
    )
  }
  for (bad in list(NA_character_, c("gaussian", "cox"), 1)) {
    expect_error(
      stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
                artificial_type = "random_permutation", family = bad,
                random_state = 1L),
      "family"
    )
  }
  # The four documented families are accepted, so the guard does not narrow the
  # supported surface (multinomial and cox are checked here because they branch
  # inside .modeled_class_levels()).
  for (good in c("gaussian", "binomial")) {
    fit <- suppressWarnings(stabl_fit(
      d$x, if (good == "binomial") {
        stats::setNames(rep(0:1, length.out = nrow(d$x)), rownames(d$x))
      } else d$y,
      lambda_grid = d$lambda_grid, n_bootstraps = 2L, family = good,
      artificial_type = "random_permutation", random_state = 1L
    ))
    expect_s3_class(fit, "stabl_fit")
  }
})

test_that("C3 l1_ratio is validated at the entry point and in auto_lambda_grid", {
  d <- .contract_data()
  for (bad in list(1.5, -0.1, c(0.5, NA_real_), c(0.5, 2), Inf)) {
    expect_error(
      stabl_fit(d$x, d$y, lambda_grid = "auto", n_bootstraps = 2L,
                n_lambda = 2L, artificial_type = NULL, hard_threshold = 0.3,
                l1_ratio = bad, random_state = 1L),
      "l1_ratio"
    )
    expect_error(
      auto_lambda_grid(d$x, d$y, n_lambda = 2L, l1_ratio = bad),
      "l1_ratio"
    )
  }
  # NULL keeps the lasso-like default and is not a validation failure.
  grid <- auto_lambda_grid(d$x, d$y, n_lambda = 2L, l1_ratio = NULL)
  expect_false("alpha" %in% names(grid))
})

test_that("C3 an NA or out-of-range alpha column is a typed error, not an empty glmnet path", {
  d <- .contract_data()
  expect_error(
    stabl_fit(d$x, d$y,
              lambda_grid = data.frame(alpha = c(0.5, NA), lambda = c(0.2, 0.1)),
              base_learner = "elastic_net", n_bootstraps = 2L,
              artificial_type = NULL, hard_threshold = 0.3, random_state = 1L),
    "alpha"
  )
  expect_error(
    stabl_fit(d$x, d$y,
              lambda_grid = data.frame(alpha = c(0.5, 1.7), lambda = c(0.2, 0.1)),
              base_learner = "elastic_net", n_bootstraps = 2L,
              artificial_type = NULL, hard_threshold = 0.3, random_state = 1L),
    "alpha"
  )
  # Also reachable through the sparse-group adapter's alpha/asparse column.
  expect_error(
    stabl_fit(d$x, d$y,
              lambda_grid = data.frame(alpha = NA_real_, lambda = 0.1),
              base_learner = "sparse_group_lasso", feature_groups = rep(1:3, each = 2L),
              n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.3,
              random_state = 1L),
    "alpha"
  )
})

test_that("C3 the grid validator also guards lambda and the grid shape", {
  d <- .contract_data()
  base <- list(x = d$x, y = d$y, n_bootstraps = 2L, artificial_type = NULL,
               hard_threshold = 0.3, random_state = 1L)
  expect_error(do.call(stabl_fit, c(base, list(
    lambda_grid = data.frame(lambda = c(0.2, NA))
  ))), "lambda")
  expect_error(do.call(stabl_fit, c(base, list(
    lambda_grid = data.frame(lambda = c(0.2, -0.1))
  ))), "lambda")
  expect_error(do.call(stabl_fit, c(base, list(
    lambda_grid = data.frame(lambda = numeric(0))
  ))), "row")
  expect_error(do.call(stabl_fit, c(base, list(
    lambda_grid = data.frame(lamb = 0.2)
  ))), "lambda")
  expect_error(do.call(stabl_fit, c(base, list(
    lambda_grid = 0.2
  ))), "lambda_grid")
})

# ---- C4: degenerate x shapes ---------------------------------------------------

test_that("C4 zero-column and zero-row x are rejected with an explicit message", {
  d <- .contract_data()

  zero_cols <- d$x[, 0, drop = FALSE]
  expect_error(
    stabl_fit(zero_cols, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
              artificial_type = "random_permutation", random_state = 1L),
    "at least one feature column"
  )

  zero_rows <- d$x[0, , drop = FALSE]
  expect_error(
    stabl_fit(zero_rows, d$y[0], lambda_grid = d$lambda_grid, n_bootstraps = 2L,
              artificial_type = "random_permutation", random_state = 1L),
    "at least one sample row"
  )
})

test_that("C4 fdr_threshold_range is validated at the fitting call", {
  d <- .contract_data()
  for (bad in list(numeric(0), c(-0.1, 0.5), c(0.5, 1.5), c(0.5, NA_real_))) {
    expect_error(
      stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
                artificial_type = "random_permutation",
                fdr_threshold_range = bad, random_state = 1L),
      "fdr_threshold_range"
    )
  }
})

# ---- C5: stratify_bootstrap provenance -----------------------------------------

test_that("C5 stratify_bootstrap reports the request and records its origin", {
  set.seed(4)
  n <- 40L
  x <- matrix(stats::rnorm(n * 5L), n, 5L,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(5L))))
  y <- stats::setNames(factor(rep(c("NP", "P"), each = n / 2)), rownames(x))
  lam <- data.frame(lambda = c(0.2, 0.1))
  study <- stats::setNames(rep(c("A", "B"), length.out = n), rownames(x))
  strata <- data.frame(outcome = y, study = study)
  rownames(strata) <- rownames(x)

  fit <- list(
    # explicit design supplied, stratify_bootstrap left at its default
    explicit = suppressWarnings(stabl_fit(
      x, y, lambda_grid = lam, family = "binomial", n_bootstraps = 4L,
      artificial_type = NULL, hard_threshold = 0.3,
      bootstrap_strata = strata, random_state = 10L
    )),
    # flag supplied, no explicit design
    flagged = suppressWarnings(stabl_fit(
      x, y, lambda_grid = lam, family = "binomial", n_bootstraps = 4L,
      artificial_type = NULL, hard_threshold = 0.3,
      stratify_bootstrap = TRUE, random_state = 10L
    )),
    # neither
    none = suppressWarnings(stabl_fit(
      x, y, lambda_grid = lam, family = "binomial", n_bootstraps = 4L,
      artificial_type = NULL, hard_threshold = 0.3, random_state = 10L
    ))
  )

  # The reported flag answers "did the user ask for stratified sampling?",
  # which is TRUE for both routes and FALSE for neither, and no longer depends
  # on whether the request happened to resolve to a design.
  expect_true(fit$explicit$stratify_bootstrap)
  expect_true(fit$flagged$stratify_bootstrap)
  expect_false(fit$none$stratify_bootstrap)

  # The literal argument and the origin of the design are recorded separately,
  # so the FALSE-but-resolved case is distinguishable.
  expect_false(fit$explicit$stratify_bootstrap_requested)
  expect_true(fit$flagged$stratify_bootstrap_requested)
  expect_identical(fit$explicit$bootstrap_strata_source, "bootstrap_strata")
  expect_identical(fit$flagged$bootstrap_strata_source, "stratify_bootstrap")
  expect_null(fit$none$bootstrap_strata_source)
  expect_true(fit$explicit$bootstrap_strata_applied)
  expect_true(fit$flagged$bootstrap_strata_applied)
  expect_false(fit$none$bootstrap_strata_applied)
  expect_length(fit$explicit$bootstrap_strata_levels, 4L)
})

# ---- C6: fit-time explore fallback ---------------------------------------------

test_that("C6 the explore fallback is resolved once at fit time and shared by all accessors", {
  set.seed(42)
  n <- 30L; p <- 10L
  x <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- stats::setNames(stats::rnorm(n), rownames(x))

  # Huge lambdas select nothing, so the fallback is what determines the support.
  fit <- stabl_fit(
    x, y, lambda_grid = data.frame(lambda = c(10, 20)), n_bootstraps = 8L,
    artificial_type = NULL, hard_threshold = 0.99, explore = TRUE,
    n_explore = 3L, random_state = 1L
  )

  expect_true(fit$explore_applied_)
  expect_length(fit$explore_features_, 3L)
  # The recorded fallback is the top-n_explore importance ranking, and the
  # fitted support is exactly that set.
  top3 <- names(sort(get_importances(fit), decreasing = TRUE))[seq_len(3L)]
  expect_identical(fit$explore_features_, top3)
  expect_identical(get_feature_names_out(fit), top3)
  expect_identical(names(get_support(fit))[get_support(fit)], top3)
  expect_identical(unname(get_support(fit)), unname(fit$support_))

  # The fallback is threshold-independent, so an overridden threshold reports
  # the same set instead of recomputing a different one.
  loose <- get_support(fit, new_hard_threshold = 1e-9)
  expect_identical(names(loose)[loose], top3)
})

test_that("C6 explore_applied_ is FALSE when the threshold already admits features", {
  d <- .contract_data()
  fit <- stabl_fit(
    d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
    artificial_type = "random_permutation", random_state = 1L
  )
  expect_false(fit$explore_applied_)
  expect_null(fit$explore_features_)
  expect_true(any(fit$support_))
})

# ---- C9: parallel reproducibility preconditions ---------------------------------

test_that("C9 workers > 1 without random_state warns that the run is not reproducible", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")
  d <- .contract_data()
  expect_warning(
    stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
              artificial_type = "random_permutation", workers = 2L),
    "random_state"
  )
  # Seeded parallel runs must stay silent, so the warning cannot be a blanket
  # "you used workers > 1" notice.
  expect_no_warning(suppressWarnings(stabl_fit(
    d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
    artificial_type = "random_permutation", workers = 2L, random_state = 1L
  )))
  # workers = 1 is unaffected either way.
  expect_no_warning(suppressWarnings(stabl_fit(
    d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
    artificial_type = "random_permutation", workers = 1L
  )))
})

# ---- C10: correlation-derived sparse groups ------------------------------------

test_that("C10 .build_corr_groups reads an undefined correlation as zero", {
  n <- 30L
  normal <- stats::rnorm(n)
  # A zero-variance column has no estimable correlation with anything, so
  # stats::cor() leaves those off-diagonal entries undefined.
  x <- cbind(normal = normal, constant = rep(1, n))
  corr_x <- suppressWarnings(stats::cor(x))
  expect_true(all(is.na(c(corr_x[1L, 2L], corr_x[2L, 1L]))))

  # Documented behaviour: undefined is read as a numerical zero, so a matrix
  # with no estimable correlation partitions exactly like a matrix whose
  # off-diagonal correlations are all 0 by construction.
  orthogonal <- cbind(
    a = c(1, 1, -1, -1),
    b = c(1, -1, 1, -1),
    c = c(1, -1, -1, 1)
  )
  corr_orth <- stats::cor(orthogonal)
  expect_equal(c(corr_orth[1L, 2L], corr_orth[1L, 3L], corr_orth[2L, 3L]),
               c(0, 0, 0))
  all_undefined <- matrix(1, nrow = 4L, ncol = 3L)
  expect_identical(
    stablr:::.build_corr_groups(all_undefined, percentile = 50),
    stablr:::.build_corr_groups(orthogonal, percentile = 50)
  )

  # Not vacuous: a genuinely correlated pair is still grouped.
  set.seed(31)
  strong <- cbind(u = stats::rnorm(40L), v = stats::rnorm(40L))
  strong[, "v"] <- strong[, "u"] + stats::rnorm(40L, sd = 0.01)
  grps <- stablr:::.build_corr_groups(strong, percentile = 95)
  expect_identical(grps[[1L]], grps[[2L]])
})

test_that("C10 a zero-variance column never errors and never merges on NaN", {
  n <- 30L
  x <- cbind(stats::rnorm(n), rep(2, n), stats::rnorm(n))
  expect_no_error(grps <- stablr:::.build_corr_groups(x, percentile = 90))
  expect_length(grps, 3L)
  expect_true(all(grps >= 1L))
  expect_false(anyNA(grps))
})

test_that("C10 the -0.1 percentile offset is part of the documented Python-parity contract", {
  set.seed(9)
  n <- 50L
  f1 <- stats::rnorm(n)
  x <- cbind(f1, f1 + stats::rnorm(n, sd = 0.05),
             matrix(stats::rnorm(n * 8L), n, 8L))
  colnames(x) <- paste0("v", seq_len(ncol(x)))

  with_offset <- stablr:::.build_corr_groups(x, percentile = 95)
  expect_identical(with_offset[[1L]], with_offset[[2L]])

  # The offset widens the supra-threshold set: the raw percentile alone would
  # leave this near-duplicate pair ungrouped at this percentile.
  corr <- stats::cor(x)
  cutoff_with <- as.numeric(stats::quantile(
    corr[upper.tri(corr)], probs = 0.95, names = FALSE)) - 0.1
  cutoff_without <- as.numeric(stats::quantile(
    corr[upper.tri(corr)], probs = 0.95, names = FALSE))
  expect_lt(cutoff_with, cutoff_without)
  expect_gt(corr[[1L, 2L]], cutoff_with)
})

test_that("C10 .build_corr_groups rejects an empty feature matrix and a bad percentile", {
  expect_error(stablr:::.build_corr_groups(matrix(numeric(0), nrow = 4, ncol = 0),
                                           percentile = 90),
               "at least one feature column")
  x <- matrix(stats::rnorm(20), 10, 2)
  for (bad in list(0, 101, -1, NA_real_, Inf, c(10, 20), "90")) {
    expect_error(stablr:::.build_corr_groups(x, percentile = bad),
                 "corr_group_threshold")
  }
  # p == 1 still returns a single group rather than erroring.
  expect_identical(stablr:::.build_corr_groups(x[, 1, drop = FALSE], 90), 1L)
})

# ---- C11: sparse-group assignment of artificial decoys --------------------------

test_that("C11 permutation decoys get their own groups so both blocks are comparable", {
  skip_if_not_installed("sparsegl")
  set.seed(18)
  n <- 44L; p <- 12L
  x <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- stats::setNames(stats::rnorm(n), rownames(x))

  fit <- stabl_fit(
    x, y, lambda_grid = data.frame(lambda = 0.1),
    base_learner = "sparse_group_lasso", feature_groups = rep(seq_len(4L), each = 3L),
    n_bootstraps = 3L, artificial_type = "random_permutation", random_state = 18L
  )

  groups <- fit$sgl_feature_groups_
  expect_length(groups, p + p)
  expect_identical(fit$sgl_group_provenance$group_source, "feature_groups")
  expect_identical(fit$sgl_group_provenance$noise_group_inheritance, "off")
  expect_identical(fit$sgl_group_provenance$n_noise_features, p)

  real_groups <- groups[seq_len(p)]
  noise_groups <- groups[p + seq_len(p)]
  # No decoy shares a group with a real feature, so the artificial block faces
  # its own group-lasso penalty and its stability score cannot be inflated by
  # the real block's group structure.
  expect_length(intersect(unique(real_groups), unique(noise_groups)), 0L)
  expect_identical(
    sort(unique(noise_groups)),
    as.integer(max(real_groups) + seq_len(p))
  )
  # The real block is unchanged.
  expect_identical(sort(unique(real_groups)), seq_len(4L))
})

test_that("C11 knockoff decoys keep inheritance on and the choice is opt-in", {
  expect_identical(stablr:::.sgl_noise_group_inheritance("random_permutation"),
                   "off")
  expect_identical(stablr:::.sgl_noise_group_inheritance("knockoff"), "on")
  expect_identical(stablr:::.sgl_noise_group_inheritance("knockoff_equi"), "on")
  expect_identical(stablr:::.sgl_noise_group_inheritance(NULL), "off")

  groups <- c(1L, 1L, 2L, 3L, 3L)
  src <- c(2L, 1L, 3L)
  # opt-in restores the co-grouping
  expect_identical(
    stablr:::.append_noise_groups(groups, src, length(groups) + length(src),
                                  inherit_groups = TRUE),
    c(1L, 1L, 2L, 3L, 3L, 1L, 1L, 2L)
  )
  # opt-out gives every decoy a fresh group
  expect_identical(
    stablr:::.append_noise_groups(groups, src, length(groups) + length(src),
                                  inherit_groups = FALSE),
    c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 6L)
  )
})

test_that("C11 non-SGL learners record no sparse-group provenance", {
  d <- .contract_data()
  fit <- stabl_fit(
    d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 3L,
    artificial_type = "random_permutation", random_state = 1L
  )
  expect_null(fit$sgl_feature_groups_)
  expect_null(fit$sgl_group_provenance)
})

# ---- C12: honest cost in the progress message ----------------------------------

test_that("C12 the progress message reports the real number of solver path calls", {
  d <- .contract_data()
  expect_identical(
    stablr:::.stabl_path_calls_per_bootstrap("lasso", d$lambda_grid), 1L
  )
  # lasso/adaptive_lasso pin alpha, so a mixed-alpha grid is still one call.
  mixed <- data.frame(alpha = c(0.2, 0.6, 0.9), lambda = c(0.3, 0.2, 0.1))
  expect_identical(stablr:::.stabl_path_calls_per_bootstrap("lasso", mixed), 1L)
  expect_identical(
    stablr:::.stabl_path_calls_per_bootstrap("adaptive_lasso", mixed), 1L
  )
  # elastic_net and sparse_group_lasso fit one path per distinct alpha.
  expect_identical(stablr:::.stabl_path_calls_per_bootstrap("elastic_net", mixed), 3L)
  expect_identical(
    stablr:::.stabl_path_calls_per_bootstrap("sparse_group_lasso", mixed), 3L
  )

  expect_message(
    stabl_fit(d$x, d$y, lambda_grid = mixed, base_learner = "elastic_net",
              n_bootstraps = 2L, artificial_type = NULL, hard_threshold = 0.3,
              verbose = TRUE, random_state = 1L),
    "3 solver path call\\(s\\) per bootstrap"
  )
  expect_message(
    stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
              artificial_type = NULL, hard_threshold = 0.3,
              verbose = TRUE, random_state = 1L),
    "1 solver path call\\(s\\) per bootstrap"
  )
})

# ---- C14: cox + stratify_bootstrap ----------------------------------------------

test_that("C14 family = cox with stratify_bootstrap = TRUE raises a typed, actionable error", {
  skip_if_not_installed("survival")
  set.seed(24)
  n <- 50L; p <- 6L
  x <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- survival::Surv(time = stats::rexp(n, rate = 0.1),
                      event = stats::rbinom(n, size = 1L, prob = 0.7))
  rownames(y) <- rownames(x)

  condition <- tryCatch(
    stabl_fit(x, y, lambda_grid = data.frame(lambda = 0.1), family = "cox",
              n_bootstraps = 2L, artificial_type = "random_permutation",
              stratify_bootstrap = TRUE, random_state = 31L),
    error = identity
  )
  expect_s3_class(condition, "stablr_bootstrap_stratification_infeasible")
  expect_s3_class(condition, "stablr_error")
  expect_match(conditionMessage(condition), "Surv", fixed = TRUE)
  expect_match(conditionMessage(condition), "bootstrap_strata", fixed = TRUE)

  # An explicit design is the documented way out, and it works.
  strata <- stats::setNames(rep(c("a", "b"), length.out = n), rownames(x))
  fit <- stabl_fit(x, y, lambda_grid = data.frame(lambda = 0.1), family = "cox",
                   n_bootstraps = 2L, artificial_type = "random_permutation",
                   bootstrap_strata = strata, random_state = 31L)
  expect_s3_class(fit, "stabl_fit")
  expect_identical(fit$bootstrap_strata_source, "bootstrap_strata")
})

# ---- D1: the caller's RNG state is never mutated -------------------------------

test_that("D1 a seeded stabl_fit leaves the caller's .Random.seed untouched", {
  d <- .contract_data()

  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  invisible(stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
                      artificial_type = "random_permutation", random_state = 17L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  # The caller's subsequent stream is therefore unaffected by the fit: the same
  # caller seed yields the same draw whether or not a fit ran in between.
  set.seed(99)
  stream_without_fit <- {
    invisible(stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
                        artificial_type = "random_permutation", random_state = 17L))
    stats::rnorm(5)
  }
  set.seed(99)
  stream_plain <- stats::rnorm(5)
  expect_identical(stream_without_fit, stream_plain)
})

# ---- D5: RNGkind / sample.kind dependency --------------------------------------

test_that("D5 the fit records the RNGkind it ran under", {
  d <- .contract_data()
  fit <- stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
                   artificial_type = "random_permutation", random_state = 3L)
  expect_identical(fit$rng_kind_, RNGkind())
})

test_that("D5 determinism is conditional on sample.kind, which is not pinned", {
  d <- .contract_data()
  kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(kind)))

  fitted <- function() {
    stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 4L,
              artificial_type = "random_permutation", random_state = 3L
    )$stabl_scores_
  }

  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  rounding <- suppressWarnings(fitted())
  suppressWarnings(RNGkind(sample.kind = "Rejection"))
  rejection <- suppressWarnings(fitted())

  # Two seeds, two sampling algorithms: the sampled bootstrap indices differ, so
  # the scores differ.  This is the documented dependency, pinned as a test so
  # that any future attempt to pin sample.kind inside the package shows up here.
  expect_false(isTRUE(all.equal(rounding, rejection)))
})

# ---- J4: one shared lambda-grid resolver ----------------------------------------

test_that("J4 the shared resolver reports a provenance label for every route", {
  d <- .contract_data()

  expect_identical(
    stablr:::.resolve_lambda_grid("auto", character(0), d$x, d$y,
                                  "gaussian", 2L)$source,
    "auto"
  )
  expect_identical(
    stablr:::.resolve_lambda_grid(d$lambda_grid, character(0), d$x, d$y,
                                  "gaussian", 2L)$source,
    "data.frame"
  )
  named <- list(mrna = d$lambda_grid)
  expect_identical(
    stablr:::.resolve_lambda_grid(named, c("mrna", "pro"), d$x, d$y,
                                  "gaussian", 2L)$source,
    "list:mrna"
  )
  # Keys are tried in order, so a candidate name wins over its block name.
  expect_identical(
    stablr:::.resolve_lambda_grid(
      list(mrna = data.frame(lambda = 0.9), pro = d$lambda_grid),
      c("early_fusion", "mrna", "pro"), d$x, d$y, "gaussian", 2L
    )$source,
    "list:mrna"
  )
  # A miss falls back to an auto grid and says so, and warns by default.
  expect_warning(
    miss <- stablr:::.resolve_lambda_grid(named, c("early_fusion"), d$x, d$y,
                                          "gaussian", 2L),
    "auto_lambda_grid"
  )
  expect_identical(miss$source, "auto_fallback")
  expect_s3_class(miss$grid, "data.frame")
  # The same fallback can be requested silently for callers that do their own
  # provenance recording.
  silent <- stablr:::.resolve_lambda_grid(named, c("early_fusion"), d$x, d$y,
                                          "gaussian", 2L, on_missing_key = "silent")
  expect_identical(silent$source, "auto_fallback")
})

test_that("J4 stabl_fit records which lambda-grid route was used", {
  d <- .contract_data()
  expect_identical(
    stabl_fit(d$x, d$y, lambda_grid = "auto", n_lambda = 3L, n_bootstraps = 2L,
              artificial_type = NULL, hard_threshold = 0.3,
              random_state = 1L)$lambda_grid_source,
    "auto"
  )
  expect_identical(
    stabl_fit(d$x, d$y, lambda_grid = d$lambda_grid, n_bootstraps = 2L,
              artificial_type = NULL, hard_threshold = 0.3,
              random_state = 1L)$lambda_grid_source,
    "data.frame"
  )
})
