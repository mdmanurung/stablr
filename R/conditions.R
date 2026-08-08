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
    "^from glmnet C\\+\\+ code \\(error code -[0-9]+\\); Numerical error",
    "^from glmnet C\\+\\+ code \\(error code -[0-9]+\\); Convergence for",
    "^one multinomial or binomial class has (1 or 0|fewer than 2) observations"
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
    "leading minor.*not positive"
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
