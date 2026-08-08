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
