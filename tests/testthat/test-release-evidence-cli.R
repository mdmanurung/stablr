.release_cli_source_path <- function() {
  source <- testthat::test_path(
    "..", "..", "inst", "analysis", "release_evidence_cli.R"
  )
  if (file.exists(source)) return(source)
  installed <- system.file(
    "analysis", "release_evidence_cli.R", package = "stablr"
  )
  if (nzchar(installed)) return(installed)
  stop("Could not locate the release-evidence CLI.", call. = FALSE)
}

test_that("RE-09 CLI exposes only prepare, run, and assemble", {
  rscript <- file.path(R.home("bin"), "Rscript")
  output <- system2(
    rscript,
    c("--vanilla", shQuote(.release_cli_source_path()), "--help"),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"))
  expect_match(paste(output, collapse = "\n"), " prepare ", fixed = TRUE)
  expect_match(paste(output, collapse = "\n"), " run ", fixed = TRUE)
  expect_match(paste(output, collapse = "\n"), " assemble ", fixed = TRUE)
})

test_that("RE-09 CLI maps invalid invocations to nonzero status", {
  rscript <- file.path(R.home("bin"), "Rscript")
  output <- suppressWarnings(system2(
    rscript,
    c("--vanilla", shQuote(.release_cli_source_path()), "unknown"),
    stdout = TRUE,
    stderr = TRUE
  ))
  expect_identical(attr(output, "status"), 70L)
  expect_match(paste(output, collapse = "\n"), "Unknown command")
})

test_that("RE-09 CLI contains no scientific gate implementation", {
  text <- readLines(.release_cli_source_path(), warn = FALSE)
  expect_false(any(grepl("run_methodology_validation", text, fixed = TRUE)))
  expect_false(any(grepl("run_late_fusion_validation", text, fixed = TRUE)))
  expect_false(any(grepl("Scientific Gate", text, fixed = TRUE)))
})
