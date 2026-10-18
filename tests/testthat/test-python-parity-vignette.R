# Contract tests for the 1.0.0 vignette bundle (Phase 13 PR-B):
# inst/extdata/python_parity_vignette/ + the stablr-python-parity and
# stablr-reading-plots articles.
#
# Gates encoded here (bundle INTEGRATION.md):
#   1. every bundled asset is installed and readable;
#   2. each fixture matches its SHA-256 manifest entry;
#   3. each fixture is byte-identical to the canonical frozen fixture in
#      tests/testthat/fixtures/python_parity/ (copies, not substitutes);
#   4. both vignettes render without Python, using only installed-package
#      system.file() resources.

.parity_vignette_dir <- function() {
  d <- system.file("extdata", "python_parity_vignette", package = "stablr")
  if (!nzchar(d)) {
    # Source-tree fallback for devtools::test() without pkgload shims.
    d <- file.path("..", "..", "inst", "extdata", "python_parity_vignette")
  }
  normalizePath(d, mustWork = FALSE)
}

.parity_vignette_sha256 <- function(path) {
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    output <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
    return(tolower(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]))
  }
  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE)
    return(tolower(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]))
  }
  NA_character_
}

.parity_vignette_files <- c(
  "fdp_plus_by_lambda.csv",
  "fdp_plus_curve.csv",
  "provenance.json",
  "solver_cases.csv",
  "stability_scores.csv",
  "stabl_fit_reference.csv"
)

test_that("python_parity_vignette assets are installed and readable", {
  d <- .parity_vignette_dir()
  expect_true(dir.exists(d))
  for (f in c(.parity_vignette_files, "manifest.sha256", "README.md")) {
    p <- file.path(d, f)
    expect_true(file.exists(p), info = f)
    expect_gt(file.size(p), 0L)
  }
})

test_that("python_parity_vignette fixtures match their SHA-256 manifest", {
  d <- .parity_vignette_dir()
  manifest <- readLines(file.path(d, "manifest.sha256"), warn = FALSE)
  manifest <- manifest[nzchar(trimws(manifest))]
  expect_length(manifest, length(.parity_vignette_files))
  for (line in manifest) {
    parts <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    expect_gte(length(parts), 2L)
    expected <- tolower(parts[[1L]])
    fname <- sub("^\\*", "", parts[[2L]])
    observed <- .parity_vignette_sha256(file.path(d, fname))
    skip_if(is.na(observed), "no SHA-256 tool available")
    expect_identical(observed, expected, info = fname)
  }
})

test_that("python_parity_vignette fixtures are byte-identical to the canonical frozen fixtures", {
  d <- .parity_vignette_dir()
  canonical <- test_path("fixtures", "python_parity")
  skip_if_not(dir.exists(canonical), "canonical fixture tree not available")
  for (f in .parity_vignette_files) {
    a <- readBin(file.path(d, f), "raw", n = file.size(file.path(d, f)))
    b <- readBin(file.path(canonical, f), "raw",
                 n = file.size(file.path(canonical, f)))
    expect_identical(a, b, info = f)
  }
})

test_that("both bundle vignettes render without Python from installed resources", {
  skip_on_cran()
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")
  skip_if_not(dir.exists(.parity_vignette_dir()),
              "python_parity_vignette extdata not installed")

  vignette_src <- c(
    system.file("doc", "stablr-python-parity.Rmd", package = "stablr"),
    file.path("..", "..", "vignettes", "stablr-python-parity.Rmd")
  )
  vignette_src <- vignette_src[file.exists(vignette_src)]
  plots_src <- c(
    system.file("doc", "stablr-reading-plots.Rmd", package = "stablr"),
    file.path("..", "..", "vignettes", "stablr-reading-plots.Rmd")
  )
  plots_src <- plots_src[file.exists(plots_src)]
  expect_gt(length(vignette_src), 0L)
  expect_gt(length(plots_src), 0L)

  # Python must be absent and unnecessary: refuse to render if reticulate
  # is loaded, and render in a clean temp library path context.
  expect_false("reticulate" %in% loadedNamespaces(),
               info = "reticulate must not be needed for these articles")

  for (src in c(vignette_src[[1L]], plots_src[[1L]])) {
    out <- file.path(tempdir(), paste0(tools::file_path_sans_ext(basename(src)),
                                       ".html"))
    render_err <- tryCatch(
      {
        suppressWarnings(rmarkdown::render(
          src,
          output_format = "html_document",
          output_file = out,
          quiet = TRUE,
          envir = new.env(parent = globalenv())
        ))
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    expect_null(render_err, info = basename(src))
    expect_true(file.exists(out))
    expect_gt(file.size(out), 10000L)
  }
})
