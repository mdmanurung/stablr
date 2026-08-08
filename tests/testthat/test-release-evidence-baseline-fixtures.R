.re01_fixture_dir <- function() {
  testthat::test_path("fixtures", "release_evidence_baseline", "v1")
}

.re01_fixture_sha256 <- function(path) {
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

test_that("RE-01 bounded fixtures have a complete SHA-256 inventory", {
  fixture_dir <- .re01_fixture_dir()
  golden <- c(
    "late_fusion_manifest.txt",
    "late_fusion_replicates.csv",
    "late_fusion_summary.csv",
    "late_fusion_warnings.csv",
    "methodology_validation_gates.csv",
    "methodology_validation_manifest.txt",
    "methodology_validation_replicates.csv",
    "methodology_validation_summary.csv",
    "methodology_validation_warnings.csv",
    "python_metrics_parity.csv"
  )
  expected <- sort(c(golden, "README.md", "manifest.sha256"))
  expect_equal(sort(list.files(fixture_dir, all.files = FALSE, no.. = TRUE)), expected)

  manifest <- readLines(file.path(fixture_dir, "manifest.sha256"), warn = FALSE)
  expect_equal(length(manifest), length(golden) + 1L)
  expect_true(all(grepl("^[0-9a-f]{64}  [A-Za-z0-9_.-]+$", manifest)))
  expected_hash <- substr(manifest, 1L, 64L)
  filenames <- substring(manifest, 67L)
  expect_equal(anyDuplicated(filenames), 0L)
  expect_setequal(filenames, c(golden, "README.md"))

  actual_hash <- vapply(
    file.path(fixture_dir, filenames),
    .re01_fixture_sha256,
    character(1L)
  )
  if (all(is.na(actual_hash))) {
    testthat::succeed("No platform SHA-256 utility is available; manifest structure was checked.")
  } else {
    expect_equal(unname(actual_hash), expected_hash)
  }
})

test_that("RE-01 normalization removes only declared volatile fields", {
  fixture_dir <- .re01_fixture_dir()
  replicate_names <- names(utils::read.csv(
    file.path(fixture_dir, "methodology_validation_replicates.csv"),
    nrows = 0L,
    check.names = FALSE
  ))
  summary_names <- names(utils::read.csv(
    file.path(fixture_dir, "methodology_validation_summary.csv"),
    nrows = 0L,
    check.names = FALSE
  ))
  expect_false("elapsed_sec" %in% replicate_names)
  expect_false("mean_elapsed_sec" %in% summary_names)
  expect_true(all(c("fit_seed", "status", "selected_features", "error") %in% replicate_names))
  expect_true(all(c("mean_empirical_fdp", "mean_tpr") %in% summary_names))

  methodology_manifest <- readLines(
    file.path(fixture_dir, "methodology_validation_manifest.txt"),
    warn = FALSE
  )
  expect_true("validation_profile: bounded" %in% methodology_manifest)
  expect_true("seed: 270627" %in% methodology_manifest)
  expect_true("target_fdp: 0.1" %in% methodology_manifest)
  expect_true("workers: 1" %in% methodology_manifest)
  expect_true(any(methodology_manifest == "created_at: <normalized>"))
  expect_true(any(methodology_manifest == "source_git_commit: <normalized>"))

  late_manifest <- readLines(
    file.path(fixture_dir, "late_fusion_manifest.txt"),
    warn = FALSE
  )
  expect_true("replicates_per_cell 1" %in% late_manifest)
  expect_true("n_bootstraps 2" %in% late_manifest)
  expect_true("n_iter 2" %in% late_manifest)
  expect_true("source_git_commit <normalized>" %in% late_manifest)
})

test_that("RE-01 fixtures are explicitly non-promotable", {
  readme <- paste(
    readLines(file.path(.re01_fixture_dir(), "README.md"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(readme, "c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1", fixed = TRUE)
  expect_match(readme, "04-r4-toolchain", fixed = TRUE)
  expect_match(readme, "not promotable release evidence", fixed = TRUE)
  expect_match(readme, "Late-fusion release gates failed; release must stop.", fixed = TRUE)
  expect_match(readme, "SCI-01 through SCI-10", fixed = TRUE)
})
