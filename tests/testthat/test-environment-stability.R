# D3: environment-stability workflow (stabl_environment_stability).
# Synthetic multi-environment design: n = 150, p = 10, 3 environments of 50
# samples each. f1 and f2 carry the SAME effect in all environments; f3 is
# planted in batchA only; f4-f10 are null.

.make_env_data <- function() {
  set.seed(2026L)
  n <- 150L
  p <- 10L
  environments <- factor(rep(c("batchA", "batchB", "batchC"), each = 50L))
  x <- matrix(
    rnorm(n * p),
    nrow = n,
    dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p)))
  )
  y <- setNames(
    0.9 * x[, "f1"] + 0.8 * x[, "f2"] +
      1.0 * x[, "f3"] * (environments == "batchA") +
      rnorm(n, sd = 0.5),
    rownames(x)
  )
  list(x = x, y = y, environments = environments)
}

.env_args <- list(
  lambda_grid     = data.frame(lambda = c(0.6, 0.4, 0.25)),
  n_bootstraps    = 15L,
  artificial_type = NULL,
  hard_threshold  = 0.5
)

.env_run <- function(d, random_state = 11L) {
  do.call(
    stabl_environment_stability,
    c(
      list(x = d$x, y = d$y, environments = d$environments,
           random_state = random_state),
      .env_args
    )
  )
}

test_that("D3: returns the documented structure", {
  d <- .make_env_data()
  res <- .env_run(d)

  expect_s3_class(res, "stabl_env_stability")
  expect_setequal(
    names(res),
    c("env_fits", "env_scores", "min_scores", "mean_scores",
      "environments", "pairwise_adjusted_similarity", "call")
  )
  expect_identical(res$environments, c("batchA", "batchB", "batchC"))
  expect_identical(names(res$env_fits), res$environments)
  expect_true(all(vapply(res$env_fits, inherits, logical(1L), "stabl_fit")))

  expect_true(is.numeric(res$env_scores))
  expect_identical(dim(res$env_scores), c(10L, 3L))
  expect_identical(colnames(res$env_scores), res$environments)
  expect_identical(rownames(res$env_scores), colnames(d$x))

  expect_identical(res$min_scores, apply(res$env_scores, 1L, min))
  expect_identical(res$mean_scores, rowMeans(res$env_scores))
  expect_true(is.call(res$call))

  sim <- res$pairwise_adjusted_similarity
  expect_identical(dim(sim), c(3L, 3L))
  expect_identical(dimnames(sim), list(res$environments, res$environments))
  expect_identical(unname(diag(sim)), rep(1, 3L))
  expect_true(all(sim[lower.tri(sim)] >= 0 & sim[lower.tri(sim)] <= 1))
})

test_that("D3: min aggregation separates stable from environment-specific features", {
  d <- .make_env_data()
  res <- .env_run(d)
  ms <- res$min_scores

  # f1/f2: same effect everywhere -> high in every environment.
  expect_true(all(ms[c("f1", "f2")] >= 0.5))
  # f3: planted in batchA only -> null-level in batchB/batchC, so its
  # weakest-environment score is low even though batchA detects it.
  expect_true(res$env_scores["f3", "batchA"] >= 0.5)
  expect_true(all(res$env_scores["f3", c("batchB", "batchC")] < 0.5))
  expect_true(ms[["f3"]] < 0.5)
  # Null features must not look transportable either.
  expect_true(all(ms[paste0("f", 4:10)] < 0.6))
})

test_that("D3: per-environment fits use deterministic derived seeds", {
  d <- .make_env_data()
  r1 <- .env_run(d, random_state = 11L)
  r2 <- .env_run(d, random_state = 11L)

  expect_identical(r1$env_scores, r2$env_scores)
  expect_identical(r1$min_scores, r2$min_scores)
  for (lv in r1$environments) {
    expect_identical(
      r1$env_fits[[lv]]$stabl_scores_,
      r2$env_fits[[lv]]$stabl_scores_
    )
  }
})

test_that("D3: character environments give the same result as factors", {
  d <- .make_env_data()
  res_f <- .env_run(d)
  res_c <- do.call(
    stabl_environment_stability,
    c(list(x = d$x, y = d$y, environments = as.character(d$environments),
           random_state = 11L),
      .env_args)
  )
  expect_identical(res_f$env_scores, res_c$env_scores)
  expect_identical(res_f$environments, res_c$environments)
})

test_that("D3: validation fails closed on NA, empty levels, and bad shapes", {
  d <- .make_env_data()

  env_na <- d$environments
  env_na[3L] <- NA
  expect_error(
    .env_run(list(x = d$x, y = d$y, environments = env_na)),
    regexp = "NA",
    class = "stablr_error"
  )

  env_empty <- factor(as.character(d$environments),
                      levels = c("batchA", "batchB", "batchC", "ghost"))
  expect_error(
    .env_run(list(x = d$x, y = d$y, environments = env_empty)),
    regexp = "ghost"
  )

  expect_error(
    stabl_environment_stability(
      d$x, d$y, d$environments[-1L],
      lambda_grid = data.frame(lambda = 0.4), n_bootstraps = 5L
    ),
    regexp = "one label per sample"
  )

  expect_error(
    stabl_environment_stability(
      d$x, d$y, factor(rep("batchA", 150L)),
      lambda_grid = data.frame(lambda = 0.4), n_bootstraps = 5L
    ),
    regexp = "At least two"
  )

  expect_error(
    stabl_environment_stability(
      d$x, d$y, as.integer(d$environments),
      lambda_grid = data.frame(lambda = 0.4), n_bootstraps = 5L
    ),
    regexp = "factor or character"
  )
})

test_that("D3: print method and get_min_scores accessor", {
  d <- .make_env_data()
  res <- .env_run(d)

  expect_output(print(res), "<stabl_env_stability>")
  expect_output(print(res), "batchA, batchB, batchC")
  expect_invisible(print(res))

  expect_identical(get_min_scores(res), res$min_scores)
  # get_min_scores dispatches on "stabl_env_stability" only.
  expect_error(get_min_scores(res$env_fits[[1L]]), regexp = "no applicable method")
})

test_that("D3: runs with the default artificial-feature (FDP+) path", {
  d <- .make_env_data()
  res <- stabl_environment_stability(
    d$x, d$y, d$environments,
    lambda_grid    = data.frame(lambda = c(0.6, 0.4, 0.25)),
    n_bootstraps   = 15L,
    random_state   = 7L
  )
  expect_s3_class(res, "stabl_env_stability")
  expect_true(is.matrix(res$pairwise_adjusted_similarity))
  expect_true(all(is.finite(res$pairwise_adjusted_similarity)))
})
