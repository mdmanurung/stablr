# Tests for the 0.3.0 two-tier API:
#   * stabl() — the beginner one-stop-shop wrapper (data.frame in,
#     plain-language selection out, safe defaults)
#   * the stabl S3 methods (print/summary/plot/support/feature_names)
#   * the 0.3.0 rename shims (warn once, delegate, no behavior change)
#   * stabl_fit(sample_groups =) with the deprecated `groups` spelling
#
# All seeds are explicit so the expectations are deterministic.

# Helper: capture all warning messages from an expression.
.catch_warnings <- function(expr) {
  w <- list()
  withCallingHandlers(
    force(expr),
    warning = function(cm) {
      w[[length(w) + 1L]] <<- conditionMessage(cm)
      invokeRestart("muffleWarning")
    }
  )
  w
}

.make_df <- function(n = 60L) {
  set.seed(1L)
  df <- data.frame(
    m1 = rnorm(n), m2 = rnorm(n), m3 = rnorm(n),
    m4 = rnorm(n), m5 = rnorm(n), n1 = rnorm(n)
  )
  df$y <- 0.9 * df$m1 + 0.9 * df$m2 + rnorm(n, sd = 0.4)
  df
}

# ---- stabl(): input handling -------------------------------------------------

test_that("stabl() accepts a data.frame with y as a column name and removes it from predictors", {
  df <- .make_df()
  fit <- stabl(df, y = "y", n_bootstraps = 30L, random_state = 1L)
  expect_s3_class(fit, "stabl")
  expect_s3_class(fit, "stabl_fit")
  # the outcome column must not reappear as a predictor
  expect_false("y" %in% fit$feature_names)
  expect_identical(fit$n_features_in_, 6L)
})

test_that("stabl() accepts a matrix plus a vector outcome", {
  df <- .make_df()
  xm <- as.matrix(df[, setdiff(names(df), "y")])
  fit <- stabl(xm, df$y, n_bootstraps = 30L, random_state = 1L)
  expect_s3_class(fit, "stabl")
  expect_identical(fit$n_features_in_, 6L)
})

test_that("stabl() aligns a named outcome vector by name", {
  df <- .make_df()
  xm <- as.matrix(df[, setdiff(names(df), "y")])
  yv <- setNames(df$y, paste0("s", seq_len(nrow(xm))))
  fit <- stabl(xm, yv, n_bootstraps = 30L, random_state = 1L)
  expect_s3_class(fit, "stabl")
})

test_that("stabl() rejects non-numeric data.frame columns with an actionable error", {
  df <- .make_df()
  df$label <- paste0("p", seq_len(nrow(df)))
  expect_error(
    stabl(df, y = "y", n_bootstraps = 10L),
    "non-numeric", fixed = FALSE
  )
})

test_that("stabl() rejects a y column name that is not in x", {
  df <- .make_df()
  expect_error(
    stabl(df, y = "nope", n_bootstraps = 10L),
    "non-numeric|must be"
  )
})

# ---- stabl(): safe defaults --------------------------------------------------

test_that("stabl() defaults to model-X knockoff artificial features and FDP+ calibration", {
  df <- .make_df()
  fit <- stabl(df, y = "y", n_bootstraps = 30L, random_state = 1L)
  expect_identical(fit$artificial_type, "knockoff_equi")
  # FDP+ mode: no hard threshold supplied
  expect_null(fit$hard_threshold)
})

test_that("stabl() forwards expert arguments through ...", {
  df <- .make_df()
  fit <- stabl(df, y = "y", n_bootstraps = 30L, random_state = 1L,
               artificial_proportion = 0.5)
  expect_identical(fit$artificial_proportion, 0.5)
})

# ---- stabl() S3 methods ------------------------------------------------------

test_that("print.stabl gives a plain-language answer and returns invisibly", {
  df <- .make_df()
  fit <- stabl(df, y = "y", n_bootstraps = 30L, random_state = 1L)
  out <- capture.output(print(fit))
  expect_match(out[1L], "<stabl stability selection>", fixed = TRUE)
  expect_true(any(grepl("feature(s) selected", out, fixed = TRUE)))
  expect_invisible(print(fit))
})

test_that("summary.stabl prints the plain-language block and returns a summary object", {
  df <- .make_df()
  fit <- stabl(df, y = "y", n_bootstraps = 30L, random_state = 1L)
  out <- capture.output(s <- summary(fit))
  expect_match(out[1L], "<stabl stability selection summary>", fixed = TRUE)
  expect_s3_class(s, "summary.stabl_fit")
  expect_identical(s$n_features_in, 6L)
})

test_that("plot.stabl delegates to plot.stabl_fit and returns a ggplot", {
  skip_if_not_installed("ggplot2")
  df <- .make_df()
  fit <- stabl(df, y = "y", n_bootstraps = 30L, random_state = 1L)
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})

# ---- argument renames on canonical functions (retained) -----------------

test_that("pearson_similarity_pairwise/summary: deprecated `d` warns and matches nb_total_elements", {
  sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
  w <- .catch_warnings(old <- pearson_similarity_pairwise(sets, NULL, d = 10L))
  expect_true(any(grepl("`d` is deprecated", w, fixed = TRUE)))
  expect_identical(old, pearson_similarity_pairwise(sets, 10L))
  w2 <- .catch_warnings(old2 <- pearson_similarity_summary(sets, NULL, d = 10L))
  expect_true(any(grepl("`d` is deprecated", w2, fixed = TRUE)))
  expect_identical(old2, pearson_similarity_summary(sets, 10L))
})

# ---- stabl_fit(sample_groups =) ----------------------------------------------

test_that("stabl_fit(groups =) warns once and matches sample_groups =", {
  df <- .make_df(60L)
  xm <- as.matrix(df[, setdiff(names(df), "y")])
  rownames(xm) <- paste0("s", seq_len(nrow(xm)))
  yv <- setNames(df$y, rownames(xm))
  grp <- setNames(rep(c("A", "B"), length.out = nrow(xm)), rownames(xm))

  w <- .catch_warnings(
    old <- stabl_fit(xm, yv, lambda_grid = "auto", n_bootstraps = 20L,
                     artificial_type = "knockoff_equi", groups = grp,
                     random_state = 1L)
  )
  expect_length(w, 1L)
  expect_match(w[[1L]], "`groups` is deprecated", fixed = TRUE)

  new <- stabl_fit(xm, yv, lambda_grid = "auto", n_bootstraps = 20L,
                   artificial_type = "knockoff_equi", sample_groups = grp,
                   random_state = 1L)
  expect_identical(old$stabl_scores_, new$stabl_scores_)
  expect_identical(old$min_fdr_, new$min_fdr_)

  # the new spelling alone does not warn
  w2 <- .catch_warnings(
    stabl_fit(xm, yv, lambda_grid = "auto", n_bootstraps = 20L,
              artificial_type = "knockoff_equi", sample_groups = grp,
              random_state = 1L)
  )
  expect_length(w2, 0L)
})

test_that("a named length-1 character y naming a column errors with guidance", {
  set.seed(1)
  df <- data.frame(a = rnorm(30), b = rnorm(30), out = rnorm(30))
  rownames(df) <- paste0("s", 1:30)
  expect_error(stabl(df, c(mycol = "out")), "named string")
})
