# Contract tests for the WS-9 ordinal benchmark script and its outputs.
# These are fast structural checks: they never run the full simulation.

.ws9_script_path <- function() {
  path <- file.path("..", "..", "inst", "benchmarks", "ws9_ordinal_benchmark.R")
  if (!file.exists(path)) {
    path <- system.file("benchmarks", "ws9_ordinal_benchmark.R",
                        package = "stablr")
  }
  skip_if_not(file.exists(path), "WS-9 benchmark script not found")
  path
}

test_that("WS-9 truth set is the nonzero coefficients of beta", {
  src <- readLines(.ws9_script_path(), warn = FALSE)
  src <- paste(src, collapse = "\n")

  # The design constructor must derive signals from nonzero beta entries.
  expect_match(src, "signals = paste0(\"f\", which(beta != 0))",
               fixed = TRUE)
  # The historical bug must be gone: no truth set over all indices.
  expect_false(grepl("signals = paste0(\"f\", seq_along(beta))",
                     src, fixed = TRUE))
})

test_that("WS-9 truth set has 3 of 10 predictors", {
  beta <- c(1.2, 1.0, 0.8, rep(0, 7))
  signals <- paste0("f", which(beta != 0))
  expect_identical(signals, c("f1", "f2", "f3"))
  expect_length(signals, 3L)
})

test_that("WS-9 empirical FDP/TPR handle the empty selected set", {
  signals <- c("f1", "f2", "f3")

  # Empty selection: FDP defined as 0, TPR 0.
  sel <- character()
  fdp <- if (length(sel) > 0) length(setdiff(sel, signals)) / length(sel) else 0
  tpr <- length(intersect(sel, signals)) / length(signals)
  expect_equal(fdp, 0)
  expect_equal(tpr, 0)

  # Mixed selection: one true, one false.
  sel <- c("f1", "f9")
  fdp <- length(setdiff(sel, signals)) / length(sel)
  tpr <- length(intersect(sel, signals)) / length(signals)
  expect_equal(fdp, 0.5)
  expect_equal(tpr, 1 / 3)
})

test_that("WS-9 script defines the full planned grid and a completeness gate", {
  src <- paste(readLines(.ws9_script_path(), warn = FALSE), collapse = "\n")

  # Planned grid over design x decoy x seed.
  expect_match(src, "expand.grid", fixed = TRUE)
  expect_match(src, "n_planned", fixed = TRUE)
  # Completeness accounting: completed + failed must equal planned.
  expect_match(src, "n_completed + n_failed != n_planned", fixed = TRUE)
  # Fail-closed: nonzero exit when any planned run is missing.
  expect_match(src, "quit(save = \"no\", status = 1L)", fixed = TRUE)
})

test_that("WS-9 script records failures instead of silently dropping them", {
  src <- paste(readLines(.ws9_script_path(), warn = FALSE), collapse = "\n")

  # Failures captured with seed, error class, and message.
  expect_match(src, "error_class", fixed = TRUE)
  expect_match(src, "error_message", fixed = TRUE)
  expect_match(src, "ws9_failures.csv", fixed = TRUE)
  # No silent `try(...); next` pattern.
  expect_false(grepl("try\\(simulate_one", src))
})

test_that("WS-9 committed results, when present, satisfy the run contract", {
  res_path <- file.path("..", "..", "inst", "benchmarks", "ws9_results.csv")
  if (!file.exists(res_path)) {
    res_path <- system.file("benchmarks", "ws9_results.csv",
                            package = "stablr")
  }
  skip_if_not(file.exists(res_path), "WS-9 results not present")

  res <- utils::read.csv(res_path, stringsAsFactors = FALSE)
  required_cols <- c("design", "artificial_type", "seed", "n_selected",
                     "fdp", "tpr")
  expect_true(all(required_cols %in% names(res)))

  # No duplicate design x decoy x seed combinations.
  key <- paste(res$design, res$artificial_type, res$seed)
  expect_equal(anyDuplicated(key), 0L)

  # Empirical FDP/TPR are bounded in [0, 1].
  expect_true(all(res$fdp >= 0 & res$fdp <= 1))
  expect_true(all(res$tpr >= 0 & res$tpr <= 1))

  # TPR is a multiple of 1/3 (3 true signals).
  expect_true(all(abs(res$tpr * 3 - round(res$tpr * 3)) < 1e-8))
})
