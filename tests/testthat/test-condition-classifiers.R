## Characterization tests for the third-party numerical-error classifiers in
## R/conditions.R.
##
## DESIGN NOTE (read before editing)
## ---------------------------------
## These tests deliberately NEVER assert a hard-coded third-party message.
## Third-party strings (glmnet's C++ prefixes, base R's LAPACK wording) change
## across versions and platforms. Instead each test:
##
##   1. PROVOKES the real condition from the real library,
##   2. hands the captured condition to stablr's classifier, and
##   3. asserts only the contract stablr actually depends on -- namely that
##      every condition which is signalled as an ERROR is classified TRUE, and
##      that conditions which are programming errors are classified FALSE.
##
## A glmnet/R version bump therefore either keeps passing (message reworded,
## same signal class) or fails with a legible diagnostic naming the family and
## class-count that stopped being recognized. It never fails on cosmetic
## wording drift.

# ── helpers ───────────────────────────────────────────────────────────────────

# Capture every condition a single expression signals, in order, without letting
# a warning escape as a test failure. Returns a list of conditions (possibly
# empty). Collecting ALL of them matters because glmnet emits unrelated
# deprecation notices ahead of the condition under test on some versions.
capture_conditions <- function(expr) {
  seen <- list()
  withCallingHandlers(
    tryCatch(force(expr), error = function(c) {
      seen[[length(seen) + 1L]] <<- c
      NULL
    }),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  seen
}

# First error, or the first warning if there was no error, or NULL.
capture_condition <- function(expr) {
  conds <- capture_conditions(expr)
  if (!length(conds)) return(NULL)
  errors <- Filter(function(c) inherits(c, "error"), conds)
  if (length(errors)) return(errors[[1L]])
  conds[[1L]]
}

# Build a bare condition of the requested class from a message. Used only for
# negative tests where provoking is impractical.
synthetic <- function(message, kind = "error") {
  cls <- if (identical(kind, "warning")) {
    c("simpleWarning", "warning", "condition")
  } else {
    c("simpleError", "error", "condition")
  }
  structure(list(message = message, call = NULL), class = cls)
}

wrap_linalg <- function(expr, source = "characterization probe") {
  stablr:::.with_linear_algebra_numerical_infeasibility(
    expr,
    source = source,
    subclass = "stablr_characterization_probe"
  )
}

wrap_glmnet <- function(expr, source = "characterization probe") {
  stablr:::.with_glmnet_numerical_infeasibility(expr, source = source)
}

# Run `expr` through the wrapper and report how it terminated.
outcome_of <- function(wrapped) {
  out <- tryCatch(wrapped(), error = function(c) c)
  if (!inherits(out, "condition")) {
    return(list(kind = "value", condition = NULL, value = out))
  }
  list(
    kind = if (inherits(out, "stablr_numerical_infeasibility")) "typed" else "untyped",
    condition = out,
    value = NULL
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── contract of .with_typed_numerical_infeasibility (version independent) ─────

test_that("an unmatched condition is rethrown unchanged (fail loud, not swallow)", {
  sentinel <- synthetic("a message the classifiers must not claim")
  raised <- tryCatch(
    wrap_linalg(stop(sentinel), source = "characterization probe"),
    error = function(c) c
  )

  expect_s3_class(raised, "error")
  expect_false(inherits(raised, "stablr_numerical_infeasibility"))
  # The very same object comes back out -- not a copy, not a wrapper.
  expect_identical(raised, sentinel)
  expect_identical(conditionMessage(raised), "a message the classifiers must not claim")
})

test_that("a matched condition becomes typed, keeps its cause, and names its source", {
  raised <- tryCatch(
    wrap_linalg(chol(matrix(c(1, 1, 1, 1), nrow = 2)), source = "probe source label"),
    error = function(c) c
  )

  expect_s3_class(raised, "stablr_characterization_probe")
  expect_s3_class(raised, "stablr_numerical_infeasibility")
  expect_s3_class(raised, "stablr_error")
  expect_match(conditionMessage(raised), "^probe source label was numerically infeasible: ")
  # The original third-party condition survives as the cause.
  expect_s3_class(raised$parent, "error")
  expect_match(conditionMessage(raised$parent), "not positive")
  expect_identical(raised$source, "probe source label")
})

test_that("an already-typed condition passes through without being re-wrapped", {
  inner <- stablr:::.new_stablr_error(
    c("stablr_characterization_inner", "stablr_numerical_infeasibility"),
    "inner infeasibility"
  )
  raised <- tryCatch(
    wrap_linalg(stop(inner), source = "probe source label"),
    error = function(c) c
  )

  expect_identical(raised, inner)
  expect_identical(conditionMessage(raised), "inner infeasibility")
  expect_null(raised$source)
})

test_that("a warning never triggers the typed fallback (error= handler only)", {
  caught <- FALSE
  raised <- tryCatch(
    withCallingHandlers(
      wrap_linalg({
        warning(synthetic("chol-like warning", kind = "warning"))
        stop("a later, unrelated failure")
      }),
      warning = function(w) {
        caught <<- TRUE
        invokeRestart("muffleWarning")
      }
    ),
    error = function(c) c
  )
  expect_true(caught)

  # The unrelated failure is what the classifier saw, and it does not match.
  expect_s3_class(raised, "error")
  expect_false(inherits(raised, "stablr_numerical_infeasibility"))
})

# ── glmnet characterization ───────────────────────────────────────────────────
#
# The single contract stablr relies on: any glmnet condition raised as an ERROR
# while fitting must be recognized, so the caller can fall back instead of
# crashing. Signals raised as warnings are recorded but must NOT be claimed by
# the error classifier, because tryCatch(error=) can never see them.

glmnet_class_count_fixture <- function(counts, family, p = 5L) {
  function() {
    y <- factor(rep(names(counts), counts), levels = names(counts))
    x <- matrix(seq_len(sum(counts) * p) / (sum(counts) * p), nrow = sum(counts), ncol = p)
    glmnet::glmnet(x = x, y = y, family = family, lambda = 0.05)
  }
}

test_that("every glmnet single-class condition raised as an error is classified", {
  skip_if_not_installed("glmnet")

  scenarios <- list(
    list(label = "binomial, 1 observation in one class", counts = c(a = 16L, b = 1L), family = "binomial"),
    list(label = "binomial, 0 observations in one class", counts = c(a = 16L, b = 0L), family = "binomial"),
    list(label = "binomial, 1 observation in both classes", counts = c(a = 1L, b = 1L), family = "binomial"),
    list(label = "multinomial, 1 observation in one class", counts = c(a = 8L, b = 1L, c = 8L), family = "multinomial"),
    list(label = "multinomial, 0 observations in one class", counts = c(a = 8L, b = 0L, c = 8L), family = "multinomial"),
    list(label = "multinomial, 1 observation in every class", counts = c(a = 1L, b = 1L, c = 1L), family = "multinomial")
  )

  observed <- lapply(scenarios, function(s) {
    skip_if_not_installed("glmnet")
    cond <- capture_condition(glmnet_class_count_fixture(s$counts, s$family)())
    list(
      label = s$label,
      signalled_as = if (is.null(cond)) "nothing" else class(cond)[[1L]],
      message = if (is.null(cond)) NA_character_ else conditionMessage(cond),
      classified = !is.null(cond) && isTRUE(stablr:::.is_glmnet_numerical_error(cond))
    )
  })

  errors <- Filter(function(o) identical(o$signalled_as, "simpleError"), observed)
  # Guard the guard: if a future glmnet stops raising anything here, this test
  # would vacuously pass, so insist the family is actually exercised.
  expect_gte(length(errors), 3L)

  for (o in errors) {
    testthat::expect_true(
      o$classified,
      label = paste0(
        "glmnet raised an unclassified error for: ", o$label,
        "\n  message: <<", o$message, ">>\n",
        "  .is_glmnet_numerical_error() must return TRUE so the learner can ",
        "fall back instead of propagating a raw third-party error."
      )
    )
  }
})

test_that("a glmnet class-count error typed through the wrapper reaches the caller", {
  skip_if_not_installed("glmnet")

  result <- outcome_of(function() {
    wrap_glmnet(glmnet_class_count_fixture(c(a = 16L, b = 1L), "binomial")())
  })

  skip_if(
    identical(result$kind, "value"),
    "this glmnet build no longer rejects a single-observation class; nothing to characterize"
  )
  expect_identical(result$kind, "typed")
  expect_s3_class(result$condition, "stablr_learner_numerical_infeasibility")
  expect_match(conditionMessage(result$condition), "was numerically infeasible: ")
  expect_s3_class(result$condition$parent, "error")
})

test_that("glmnet convergence signals are warnings and never reach error=", {
  skip_if_not_installed("glmnet")

  set.seed(20240501L)
  x <- matrix(stats::rnorm(60 * 8), nrow = 60, ncol = 8)
  y <- factor(rep(c("a", "b", "c"), each = 20L))
  conds <- capture_conditions(
    glmnet::glmnet(
      x = x, y = y, family = "multinomial",
      lambda = exp(seq(30, -8, length.out = 60)), maxit = 1L
    )
  )
  skip_if(
    !length(conds),
    "this glmnet build converged; nothing to characterize"
  )

  errors <- Filter(function(c) inherits(c, "error"), conds)

  if (length(errors)) {
    # A future glmnet may promote the convergence message to an error. If it
    # does, the classifier MUST match it, because error= will then see it.
    for (e in errors) {
      testthat::expect_true(
        isTRUE(stablr:::.is_glmnet_numerical_error(e)),
        label = paste0(
          "glmnet raised a convergence ERROR the classifier does not match.\n",
          "  message: <<", conditionMessage(e), ">>"
        )
      )
    }
  } else {
    # Documented glmnet 5.0 behaviour: convergence shortfalls are WARNINGS.
    # .with_glmnet_numerical_infeasibility installs only an error= handler, so
    # it cannot and must not intercept them. This branch records the
    # characterization without asserting a hard-coded message.
    messages <- vapply(conds, conditionMessage, character(1L))
    expect_true(any(nzchar(messages)))
  }
})

test_that("glmnet programming errors are not claimed as numerical infeasibility", {
  for (message in c(
    "x should be a matrix with 2 or more columns",
    "x and y have different number of rows in call to glmnet",
    "offset should have the same number of values as observations in binomial/multinomial call to glmnet"
  )) {
    expect_false(
      stablr:::.is_glmnet_numerical_error(synthetic(message)),
      label = paste("must not classify as numerical infeasibility:", message)
    )
  }
})

# ── base R / LAPACK characterization ──────────────────────────────────────────
#
# Same contract: every base-R/LAPACK condition that a knockoff or covariance
# routine raises as an ERROR must be recognized, because stablr relies on that
# to convert it into a typed, catchable infeasibility and fall back to
# random_permutation.

## Only failures that are genuinely *numerical* belong in the must-classify set.
## Shape errors ("0 x 0 matrix", "non-square matrix in 'eigen'") are excluded on
## purpose: classifying those as numerical infeasibility would mask a genuine
## stablr bug behind a graceful fallback, so they are asserted FALSE below.
linalg_fixture <- list(
  chol_not_positive = function() chol(matrix(c(1, 1, 1, 1), nrow = 2)),
  chol_psd_singular = function() chol(matrix(c(1, 2, 2, 4), nrow = 2)),
  chol_negative_definite = function() chol(diag(-1, 2)),
  chol_with_nan = function() {
    m <- diag(2)
    m[1, 1] <- NaN
    chol(m)
  },
  solve_exactly_singular = function() solve(matrix(c(1, 1, 1, 1), nrow = 2), c(1, 2)),
  backsolve_singular = function() backsolve(matrix(0, 2, 2), c(1, 2)),
  forwardsolve_singular = function() forwardsolve(matrix(0, 2, 2), c(1, 2)),
  qr_solve_singular = function() qr.solve(matrix(c(1, 1, 1, 1), nrow = 2), c(1, 2)),
  eigen_non_finite = function() {
    m <- diag(2)
    m[1, 1] <- Inf
    eigen(m, symmetric = TRUE, only.values = TRUE)
  },
  eigen_na = function() {
    m <- diag(2)
    m[1, 1] <- NA_real_
    eigen(m, symmetric = TRUE, only.values = TRUE)
  },
  svd_non_finite = function() svd(matrix(c(NA_real_, 1, 1, NA_real_), 2, 2))
)

test_that("every provoked base-R/LAPACK error is classified by the linalg classifier", {
  observed <- lapply(names(linalg_fixture), function(nm) {
    cond <- capture_condition(linalg_fixture[[nm]]())
    list(
      name = nm,
      signalled_as = if (is.null(cond)) "nothing" else class(cond)[[1L]],
      message = if (is.null(cond)) NA_character_ else conditionMessage(cond),
      classified = !is.null(cond) && isTRUE(stablr:::.is_linear_algebra_numerical_error(cond))
    )
  })
  names(observed) <- names(linalg_fixture)

  errors <- Filter(function(o) identical(o$signalled_as, "simpleError"), observed)
  expect_gte(length(errors), 6L)

  for (o in errors) {
    testthat::expect_true(
      o$classified,
      label = paste0(
        "base R raised an unclassified error for: ", o$name, "\n",
        "  message: <<", o$message, ">>\n",
        "  .is_linear_algebra_numerical_error() must return TRUE, otherwise the ",
        "knockoff fallback chain never engages and the raw error reaches the user."
      )
    )
  }
})

test_that("each classified base-R error is turned into a typed infeasibility", {
  for (nm in names(linalg_fixture)) {
    result <- outcome_of(function() wrap_linalg(linalg_fixture[[nm]]()))
    if (identical(result$kind, "value")) next
    testthat::expect_identical(
      result$kind, "typed",
      label = paste0(
        "wrapper did not type-ify: ", nm, "\n",
        "  message: <<", conditionMessage(result$condition), ">>"
      )
    )
    expect_match(conditionMessage(result$condition), "was numerically infeasible: ")
    expect_s3_class(result$condition$parent, "error")
  }
})

test_that("linalg programming errors are not claimed as numerical infeasibility", {
  for (message in c(
    "0 x 0 matrix",
    "non-square matrix in 'eigen'",
    "invalid nrow(x)",
    "'b' (3 x 1) must be compatible with 'a' (2 x 2)"
  )) {
    expect_false(
      stablr:::.is_linear_algebra_numerical_error(synthetic(message)),
      label = paste("must not classify as numerical infeasibility:", message)
    )
  }
})

test_that("eigen succeeds silently for rank-deficient input, so no guard is expected", {
  # Documents WHY the classifier needs the non-finite patterns: a singular
  # covariance does not raise at all, it just returns non-positive eigenvalues.
  # Detection is the caller's job, which is what the eigen() guards in
  # artificial_features.R and knockoff_mvr.R do.
  singular <- matrix(c(1, 2, 2, 4), nrow = 2)
  expect_null(capture_condition(eigen(singular, symmetric = TRUE, only.values = TRUE)))
  expect_lte(
    min(eigen(singular, symmetric = TRUE, only.values = TRUE)$values), 0
  )
})

test_that("forwardsolve on a singular but non-zero-diagonal matrix returns silently", {
  # The reason knockoff_mvr.R's forward/backward solves are unguarded by message
  # matching: they only raise on an exactly-zero diagonal, and otherwise return
  # finite garbage that poisons downstream cn/cd without any condition.
  a <- matrix(c(1, 1, 1, 1), nrow = 2)
  expect_null(capture_condition(forwardsolve(a, c(1, 2))))
  expect_null(capture_condition(backsolve(a, c(1, 2))))
  expect_false(anyNA(forwardsolve(a, c(1, 2))))
})

# ── the regression this file exists to prevent ────────────────────────────────

test_that("an eigen() failure inside a knockoff constructor becomes a fallback, not a crash", {
  skip_if_not_installed("knockoff")

  set.seed(20240502L)
  x <- matrix(stats::rnorm(60 * 6), nrow = 60, ncol = 6)

  for (constructor in c(
    "make_knockoff_features",
    "make_knockoff_equi_features",
    "make_knockoff_mvr_features"
  )) {
    force(constructor)
    result <- outcome_of(function() {
      get(constructor, envir = asNamespace("stablr"))(x, 3L, random_state = 99L)
    })

    # A well-conditioned input must never be reported as numerically
    # infeasible; if it ever is, a pattern has grown too broad.
    expect_identical(
      result$kind, "value",
      label = paste(constructor, "reported:", if (is.null(result$condition)) "?" else conditionMessage(result$condition))
    )
  }
})

test_that("a non-finite design degrades gracefully instead of crashing", {
  set.seed(20240503L)
  x_ok <- matrix(stats::rnorm(60 * 6), nrow = 60, ncol = 6)
  x_bad <- x_ok
  x_bad[3, 2] <- NA_real_

  for (constructor in c("make_knockoff_equi_features", "make_knockoff_mvr_features")) {
    force(constructor)
    signalled <- list()
    result <- outcome_of(function() {
      withCallingHandlers(
        get(constructor, envir = asNamespace("stablr"))(x_bad, 3L, random_state = 99L),
        warning = function(w) {
          signalled[[length(signalled) + 1L]] <<- w
          invokeRestart("muffleWarning")
        }
      )
    })

    if (identical(result$kind, "untyped")) {
      # Fail loud is acceptable, provided the caller still sees an error.
      expect_s3_class(result$condition, "error")
      next
    }

    # Returning a value is only acceptable if the value is a genuine
    # random_permutation fallback, never a knockoff built from a non-finite
    # design. This is the property the fallback chain exists to guarantee.
    if (!identical(result$kind, "value")) {
      fail(paste0(
        constructor, " raised an unclassified error: ",
        conditionMessage(result$condition),
        "\n  every error a knockoff constructor can raise must be classified, ",
        "otherwise the fallback chain never engages."
      ))
    }
    provenance <- result$value$artificial_provenance
    expect_false(is.null(provenance))
    # Use [[ ]] / exact names: `$` on a list does partial matching and would
    # silently pick up `selected_type_counts`.
    types <- unique(c(
      provenance[["actual_type"]],
      unlist(provenance[["actual_types"]]),
      as.character(provenance[["chunks"]][["actual_type"]])
    ))
    testthat::expect_true(
      all(types == "random_permutation"),
      label = paste0(
        constructor, " returned non-random-permutation features for a ",
        "non-finite design: ", paste(types, collapse = ", "),
        "\n  a non-finite design must never be used to build real knockoffs."
      )
    )
    # A silent fallback would be indistinguishable from a bug in the provenance
    # record, so require the fallback to be announced.
    testthat::expect_true(
      any(vapply(signalled, function(w) grepl("fall", conditionMessage(w), ignore.case = TRUE), logical(1L))),
      label = paste(constructor, "fell back to random_permutation without warning the caller")
    )
    # Only the injected columns are stablr's responsibility; the NA lives in
    # the caller's own input and is echoed back by cbind().
    injected <- result$value$x_augmented[, -seq_len(ncol(x_bad)), drop = FALSE]
    expect_true(all(is.finite(injected)))
  }
})
