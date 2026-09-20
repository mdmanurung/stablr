release_source_root <- function() {
  candidates <- unique(c(
    normalizePath(testthat::test_path("..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "DESCRIPTION")) &&
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
  }

  testthat::skip("package source root is not available")
}

release_source_lines <- function(path) {
  file <- file.path(release_source_root(), path)
  if (!file.exists(file)) {
    testthat::skip(sprintf("source file '%s' is not available", path))
  }
  readLines(file, warn = FALSE)
}

test_that(".Rbuildignore explicitly excludes local release-check byproducts", {
  lines <- release_source_lines(".Rbuildignore")

  expect_true(any(grepl("CONTEXT[.]md", lines, fixed = TRUE)))
  expect_true(any(grepl("Rcheck", lines, fixed = TRUE)))
  expect_true(any(grepl("tar[.]gz", lines, fixed = TRUE)) ||
              any(grepl("tar[.]xz", lines, fixed = TRUE)))
  expect_true(any(grepl("src/.*[.](o|so|dll)$", lines, fixed = TRUE)))
  expect_true(any(grepl("test-release-hygiene[.]R", lines, fixed = TRUE)))
  expect_true(any(grepl("test-vignette-readiness[.]R", lines, fixed = TRUE)))
  expect_true(any(vapply(
    lines,
    function(pattern) isTRUE(grepl(pattern, ".claude/reviews/evidence")),
    logical(1)
  )))
})

test_that("manual release-candidate workflow binds the artifact to main", {
  text <- paste(
    release_source_lines(".github/workflows/release-candidate.yaml"),
    collapse = "\n"
  )

  expect_match(text, "workflow_dispatch:", fixed = TRUE)
  expect_match(text, "candidate_sha:", fixed = TRUE)
  expect_match(text, "git rev-parse origin/main", fixed = TRUE)
  expect_match(text, "stablr_0.1.1.tar.gz.sha256", fixed = TRUE)
  expect_match(text, "release-candidate-build.json", fixed = TRUE)
  expect_match(text, "actions/upload-artifact@v4", fixed = TRUE)
  expect_match(text, "actions/download-artifact@v5", fixed = TRUE)
  expect_match(text, "sha256sum -c stablr_0.1.1.tar.gz.sha256", fixed = TRUE)
  expect_match(text, "stablr.exact-candidate-receipt/v1", fixed = TRUE)
  for (job_id in c(
    "linux-release", "linux-devel", "linux-oldrel", "macos-release",
    "windows-release", "r-4.1", "clean-install"
  )) {
    expect_match(text, job_id, fixed = TRUE)
  }
  build_calls <- gregexpr("R CMD build", text, fixed = TRUE)[[1L]]
  expect_length(build_calls[build_calls > 0L], 1L)
  expect_match(text, "needs: build-candidate", fixed = TRUE)
  expect_match(text, "steps.exact_check.outcome", fixed = TRUE)
  expect_match(text, "steps.install_smoke.outcome", fixed = TRUE)
})

test_that("SCI review packet preserves v1 and declares the corrected matrix", {
  v1 <- release_source_lines(
    "inst/analysis/contracts/SCI-04-independent-review-packet-v1.md"
  )
  v2 <- paste(release_source_lines(
    "inst/analysis/contracts/SCI-04-independent-review-packet-v2.md"
  ), collapse = "\n")
  expect_true(any(grepl("Packet ID: `SCI-04-independent-review-packet-v1`", v1,
                        fixed = TRUE)))
  expect_match(v2, "96 declared cells", fixed = TRUE)
  expect_match(v2, "Eight fixed-X/high-dimensional cells", fixed = TRUE)
  expect_match(v2, "88 applicable cells", fixed = TRUE)
  expect_match(v2, "44 are global-null cells", fixed = TRUE)
  expect_match(v2, "No formula, replicate count", fixed = TRUE)
})

test_that("repository license contains the complete GPL version 2 text", {
  lines <- release_source_lines("LICENSE.md")
  text <- paste(lines, collapse = "\n")

  expect_gt(length(lines), 300L)
  expect_true(any(grepl("Version 2, June 1991", lines, fixed = TRUE)))
  expect_true(any(grepl("TERMS AND CONDITIONS FOR COPYING", lines,
                        fixed = TRUE)))
  expect_match(text,
               "either version 2 of the License or, at your[[:space:]]+option")
})

test_that("bundled OOL provenance declares source, terms, transformations, and hashes", {
  text <- paste(release_source_lines("data-raw/README.md"), collapse = "\n")

  expect_match(text, "10.5061/dryad.stqjq2c7d", fixed = TRUE)
  expect_match(text, "10.1126/scitranslmed.abd9898", fixed = TRUE)
  expect_match(text, "CC0", fixed = TRUE)
  expect_match(text, "first 100", fixed = TRUE)
  expect_match(text, "SHA-256", fixed = TRUE)
  expect_match(text, "no additional normalization or imputation", fixed = TRUE)
})

test_that("cran-comments remain explicitly provisional before external checks", {
  text <- paste(release_source_lines("cran-comments.md"), collapse = "\n")

  expect_match(text, "DRAFT", fixed = TRUE)
  expect_match(text, "first CRAN submission", fixed = TRUE)
  expect_match(text, "stablr_0.1.1.tar.gz", fixed = TRUE)
  expect_false(grepl("stablr_0.1.0.tar.gz", text, fixed = TRUE))
  expect_false(grepl("0 errors \\| 0 warnings", text))
  expect_false(grepl("--no-manual --ignore-vignettes --no-build-vignettes", text))
})

test_that("cran-comments CRAN-skip notes match source skip_on_cran guards", {
  root <- release_source_root()
  test_files <- list.files(
    file.path(root, "tests", "testthat"),
    pattern = "\\.R$",
    full.names = TRUE
  )
  has_skip <- vapply(
    test_files,
    function(file) any(grepl("skip_on_cran\\(", readLines(file, warn = FALSE))),
    logical(1)
  )
  skip_files <- basename(test_files[has_skip])
  text <- paste(release_source_lines("cran-comments.md"), collapse = "\n")

  expect_match(text, sprintf("Tests skipped on CRAN \\(%d\\)", length(skip_files)))
  for (file in skip_files) {
    expect_match(text, file, fixed = TRUE)
  }
})

test_that("release notes do not retain stale no-vignette check metadata", {
  news_lines <- release_source_lines("NEWS.md")
  historical <- which(news_lines == "# stablr 0.1.0")
  expect_length(historical, 1L)
  news <- paste(news_lines[seq_len(historical - 1L)], collapse = "\n")

  expect_false(grepl("no-vignette", news, ignore.case = TRUE))
  expect_false(grepl("no-manual/no-vignette", news, ignore.case = TRUE))
  expect_false(grepl("PASS 1611", news, fixed = TRUE))
})

test_that("vendored cooperative S3 note matches the hand-maintained namespace", {
  cran <- paste(release_source_lines("cran-comments.md"), collapse = "\n")
  namespace <- paste(release_source_lines("NAMESPACE"), collapse = "\n")

  expect_match(namespace, "S3method(plot,multiview)", fixed = TRUE)
  expect_match(cran, "plot.multiview", fixed = TRUE)
  expect_false(grepl("Package-level startup messages", cran, fixed = TRUE))
  expect_false(grepl("coef_ordered.multiview", cran, fixed = TRUE))
})
