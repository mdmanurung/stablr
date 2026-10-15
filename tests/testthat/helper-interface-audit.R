# Shared helpers for tests that read the public-interface audit registry.
# The registry lives in the package source tree (inst/release/); in an
# R CMD check tree the source layout is not available, so callers skip.

.interface_source_root <- function() {
  candidates <- unique(c(
    normalizePath(testthat::test_path("..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  ))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "NAMESPACE")) &&
        dir.exists(file.path(candidate, "inst", "release"))) {
      return(candidate)
    }
  }
  testthat::skip("package source interface audit is not available")
}

.read_interface_audit <- function() {
  utils::read.csv(
    file.path(
      .interface_source_root(), "inst", "release",
      "public-interface-audit.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character()
  )
}
