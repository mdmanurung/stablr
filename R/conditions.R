.new_stablr_error <- function(subclass, message, call = NULL, ...) {
  structure(
    c(
      list(message = message, call = call),
      list(...)
    ),
    class = c(subclass, "stablr_error", "error", "condition")
  )
}

.abort_stablr <- function(subclass, message, call = NULL, ...) {
  stop(.new_stablr_error(subclass, message, call = call, ...))
}

.abort_numerical_infeasibility <- function(subclass, message, parent = NULL,
                                           call = NULL, ...) {
  .abort_stablr(
    c(subclass, "stablr_numerical_infeasibility"),
    message,
    call = call,
    parent = parent,
    ...
  )
}

.is_glmnet_numerical_error <- function(condition) {
  message <- conditionMessage(condition)
  patterns <- c(
    # glmnet signals a negative C++ error code as a non-fatal warning and a
    # positive code as a fatal error, so the sign must not constrain the
    # pattern. R prepends "(converted from warning) " when options(warn = 2)
    # promotes a warning, which defeats a leading anchor.
    "^(?:\\(converted from warning\\) )?from glmnet C\\+\\+ code \\(error code -?[0-9]+\\); (?:Numerical error|Convergence for)",
    # Reachable fatal codes that no other pattern covers. Cohort-level codes
    # (20000 initialization, 8888 all-censored, 9999 no positive weights) are
    # deliberately excluded: they describe the data, not numerical
    # conditioning, so a knockoff fallback must not absorb them.
    "^(?:\\(converted from warning\\) )?from glmnet C\\+\\+ code \\(error code -?[0-9]+\\); (?:Memory allocation error|All used predictors have zero variance|All penalty factors are <= 0)",
    # glmnet 5.0 emits the literal "1 or 0" for both arities. The live sparse
    # class message is "fewer than 8  observations" with a repeated space, so
    # match any count and tolerate the whitespace.
    "^one multinomial or binomial class has (?:1 or 0|fewer than [0-9]+ +) observations",
    "^an empty model has been returned"
  )
  any(vapply(patterns, grepl, logical(1L), x = message, perl = TRUE))
}

.is_linear_algebra_numerical_error <- function(condition) {
  message <- conditionMessage(condition)
  patterns <- c(
    "eigenvalue\\(s\\) converged",
    "did not converge",
    "not positive definite",
    "not positive semidefinite",
    "computationally singular",
    "system is exactly singular",
    "leading minor.*not positive",
    # base::forwardsolve, base::backsolve and base::qr.solve wording, which
    # differs from solve()'s LAPACK message.
    "singular matrix.*(?:forwardsolve|backsolve|solve)",
    # base::eigen and base::svd reject non-finite input up front rather than
    # reporting a failed decomposition.
    "infinite or missing values in 'x'",
    "NA/NaN/Inf in foreign function call"
  )
  any(vapply(
    patterns,
    grepl,
    logical(1L),
    x = message,
    ignore.case = TRUE,
    perl = TRUE
  ))
}

.with_typed_numerical_infeasibility <- function(expr, subclass, source,
                                                classifier) {
  tryCatch(
    expr,
    error = function(condition) {
      if (inherits(condition, "stablr_numerical_infeasibility")) {
        stop(condition)
      }
      if (!isTRUE(classifier(condition))) {
        stop(condition)
      }
      .abort_numerical_infeasibility(
        subclass,
        paste0(source, " was numerically infeasible: ",
               conditionMessage(condition)),
        parent = condition,
        source = source
      )
    }
  )
}

.with_glmnet_numerical_infeasibility <- function(expr, source) {
  .with_typed_numerical_infeasibility(
    expr,
    subclass = "stablr_learner_numerical_infeasibility",
    source = source,
    classifier = .is_glmnet_numerical_error
  )
}

.with_linear_algebra_numerical_infeasibility <- function(expr, source,
                                                         subclass) {
  .with_typed_numerical_infeasibility(
    expr,
    subclass = subclass,
    source = source,
    classifier = .is_linear_algebra_numerical_error
  )
}
