source_root <- function() {
  candidates <- unique(c(
    normalizePath(testthat::test_path("..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "DESCRIPTION")) &&
        dir.exists(file.path(candidate, "vignettes"))) {
      return(candidate)
    }
  }

  testthat::skip("package source vignettes are not available")
}

vignette_lines <- function(file) {
  readLines(file.path(source_root(), "vignettes", file), warn = FALSE)
}

test_that("vignette references to package-local analysis files resolve", {
  root <- source_root()
  vignette_files <- list.files(
    file.path(root, "vignettes"),
    pattern = "\\.Rmd$",
    full.names = TRUE
  )
  text <- paste(unlist(lapply(vignette_files, readLines, warn = FALSE)), collapse = "\n")
  matches <- gregexpr("inst/analysis/[[:alnum:]_.-]+", text, perl = TRUE)
  refs <- unique(unlist(regmatches(text, matches), use.names = FALSE))
  refs <- refs[nzchar(refs)]

  missing <- refs[!file.exists(file.path(root, refs))]
  expect_equal(missing, character())
})

test_that("knockoff artificial-feature taxonomy is scientifically precise", {
  lines <- vignette_lines("stablr-advanced.Rmd")
  bad_model_x_claim <- grep(
    "Model-X knockoffs.*artificial_type = \"knockoff\"",
    lines,
    value = TRUE
  )

  expect_equal(bad_model_x_claim, character())
})

test_that("multiomic vignette uses the same matrix for lambda grid and fit", {
  lines <- vignette_lines("stablr-multiomic.Rmd")
  text <- paste(lines, collapse = "\n")
  direct_scaled_fit <- grep("x\\s*=\\s*scale\\(", lines, value = TRUE)

  expect_equal(direct_scaled_fit, character())
  expect_match(
    text,
    "lambda_prot <- auto_lambda_grid\\(\\s*proteomics_train_scaled",
    perl = TRUE
  )
  expect_match(
    text,
    "fit_prot <- stabl_fit\\([\\s\\S]*x\\s*=\\s*proteomics_train_scaled",
    perl = TRUE
  )
})

test_that("nested TCGA vignette does not access external artifacts by default", {
  text <- paste(vignette_lines("stablr-tcga-nestedcv.Rmd"), collapse = "\n")

  expect_match(text, "STABLR_USE_EXTERNAL_TCGA_CACHE")
  expect_match(text, "use_external_cache && file\\.exists\\(analysis_script\\)")
  expect_match(text, "<external TCGA cache disabled>", fixed = TRUE)
})

test_that("fusion interpretation prose avoids leakage and causal overclaiming", {
  text <- paste(
    vignette_lines("stablr-multiomic.Rmd"),
    vignette_lines("stablr-cooperative.Rmd"),
    collapse = "\n"
  )

  expect_false(grepl("weighted combination of predictions on the validation set", text))
  expect_false(grepl("independent predictive contribution", text))
  expect_false(grepl("actively helped each other", text))
  expect_false(grepl("not benefiting from cooperation", text))
})

test_that("bounded OOL tutorial discloses its exact non-evidence subset", {
  text <- paste(vignette_lines("stablr-multiomic.Rmd"), collapse = "\n")

  expect_match(text, "150 training samples", fixed = TRUE)
  expect_match(text, "21 aligned validation samples", fixed = TRUE)
  expect_match(text, "first 100 features from each", fixed = TRUE)
  expect_match(text, "not publication evidence", ignore.case = TRUE)
})

test_that("TCGA pages are explicitly non-executable research protocols", {
  for (file in c("stablr-tcga.Rmd", "stablr-tcga-nestedcv.Rmd")) {
    text <- paste(vignette_lines(file), collapse = "\n")
    expect_match(text, "Research protocol", ignore.case = TRUE, info = file)
    expect_match(text, "Non-executable research protocol", fixed = TRUE, info = file)
    expect_match(text, "No TCGA", fixed = TRUE, info = file)
  }
})

test_that("explicit eval-false package calls name their execution coverage", {
  expectations <- list(
    "stablr-advanced.Rmd" = "test-phase7.R",
    "stablr-cooperative.Rmd" = "test-multiomic-workflows.R",
    "stablr-intro.Rmd" = "test-phase7.R"
  )
  for (file in names(expectations)) {
    text <- paste(vignette_lines(file), collapse = "\n")
    expect_match(text, expectations[[file]], fixed = TRUE, info = file)
  }
})
