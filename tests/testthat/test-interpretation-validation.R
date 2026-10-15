# Contract tests for the IV interpretation/validated-workflow benchmark
# (PR6-B) and its committed outputs.  These are fast structural checks:
# they never run the full simulation.

.iv_benchmark_path <- function(file) {
  path <- file.path("..", "..", "inst", "benchmarks", file)
  if (!file.exists(path)) {
    path <- system.file("benchmarks", file, package = "stablr")
  }
  skip_if_not(file.exists(path), paste("IV benchmark file not found:", file))
  path
}

test_that("IV predesign exists, predates the run, and fixes the gate criteria", {
  src <- paste(readLines(.iv_benchmark_path(
    "interpretation_validation_predesign.md"), warn = FALSE), collapse = "\n")

  # Written-before-running declaration and fixed date.
  expect_match(src, "Written BEFORE any run", fixed = TRUE)
  expect_match(src, "Date: 2026-10-15", fixed = TRUE)
  # All six gate criteria are predeclared.
  for (g in paste0("G-IV-", 1:6)) {
    expect_match(src, g, fixed = TRUE)
  }
  # Honesty caveats are part of the frozen design.
  expect_match(src, "Simulation-only evidence", fixed = TRUE)
  expect_match(src, "Adverse results are committed, not suppressed",
               fixed = TRUE)
})

test_that("IV script defines the full planned grid and a completeness gate", {
  src <- paste(readLines(.iv_benchmark_path("interpretation_validation.R"),
                         warn = FALSE), collapse = "\n")

  # Planned grid: 4 designs x 2 arms + IV-4 validated + IV-5 leaky/nested.
  expect_match(src, "expand.grid", fixed = TRUE)
  expect_match(src, "n_planned", fixed = TRUE)
  # Fail-closed: nonzero exit when any planned run is missing.
  expect_match(src, "quit(save = \"no\", status = 1L)", fixed = TRUE)
  # Failures are recorded with class and message, never silently dropped.
  expect_match(src, "error_class", fixed = TRUE)
  expect_match(src, "error_message", fixed = TRUE)
  expect_match(src, "iv_failures.csv", fixed = TRUE)
})

test_that("IV script evaluates exactly the predeclared gate criteria", {
  src <- paste(readLines(.iv_benchmark_path("interpretation_validation.R"),
                         warn = FALSE), collapse = "\n")

  for (g in paste0("G-IV-", 2:6)) {
    expect_match(src, g, fixed = TRUE)
  }
  # Gate thresholds match the predesign.
  expect_match(src, "g2_rate >= 0.90", fixed = TRUE)
  expect_match(src, "g3_rate >= 0.80", fixed = TRUE)
  expect_match(src, "iv0_median_auroc >= 0.35", fixed = TRUE)
  expect_match(src, "iv0_median_auroc <= 0.65", fixed = TRUE)
})

test_that("IV committed results satisfy the run contract", {
  res <- utils::read.csv(.iv_benchmark_path("iv_results.csv"),
                         stringsAsFactors = FALSE)

  required_cols <- c("design", "arm", "seed", "n_selected",
                     "selected_features", "accuracy", "ber", "macro_f1",
                     "mcc", "brier", "log_loss", "auroc", "runtime_sec",
                     "git_sha", "pkg_version")
  expect_true(all(required_cols %in% names(res)))

  # Complete grid: 110 planned runs, no duplicates, no failures file rows.
  expect_equal(nrow(res), 110L)
  key <- paste(res$design, res$arm, res$seed)
  expect_equal(anyDuplicated(key), 0L)
  fail <- utils::read.csv(.iv_benchmark_path("iv_failures.csv"),
                          stringsAsFactors = FALSE)
  expect_equal(nrow(fail), 0L)

  # Arm coverage per design matches the predesign.
  arms_for <- function(d) sort(unique(res$arm[res$design == d]))
  for (d in c("iv0_null", "iv1_sparse", "iv2_substitutes", "iv3_single")) {
    expect_identical(arms_for(d), c("exploratory", "validated"))
  }
  expect_identical(arms_for("iv4_study_shift"), "validated")
  expect_identical(arms_for("iv5_leaky"), c("leaky", "nested"))

  # Probability metrics are bounded where recorded (the leaky arm records
  # probability metrics but no hard-class metrics).
  prob <- res[!is.na(res$auroc), ]
  expect_true(all(prob$auroc >= 0 & prob$auroc <= 1))
  expect_true(all(prob$brier >= 0 & prob$brier <= 1))
  expect_true(all(prob$log_loss >= 0))
  hard <- res[!is.na(res$accuracy), ]
  expect_true(all(hard$accuracy >= 0 & hard$accuracy <= 1))
  expect_true(all(hard$ber >= 0 & hard$ber <= 1))
  expect_true(all(hard$mcc >= -1 & hard$mcc <= 1))
})

test_that("IV committed detail tables are internally consistent", {
  res <- utils::read.csv(.iv_benchmark_path("iv_results.csv"),
                         stringsAsFactors = FALSE)
  eff <- utils::read.csv(.iv_benchmark_path("iv_effects.csv"),
                         stringsAsFactors = FALSE)
  folds <- utils::read.csv(.iv_benchmark_path("iv_folds.csv"),
                           stringsAsFactors = FALSE)
  imp <- utils::read.csv(.iv_benchmark_path("iv_importance.csv"),
                         stringsAsFactors = FALSE)

  # Effects: exploratory arms only, one row per feature per run.
  expect_true(all(unique(eff$design) %in%
                    c("iv0_null", "iv1_sparse", "iv2_substitutes",
                      "iv3_single")))
  expect_true(all(eff$effect_direction %in%
                    c("positive", "negative", "zero", NA)))
  expect_true(all(eff$stability >= 0 & eff$stability <= 1))

  # Folds: one row per outer fold of every validated/nested run.
  n_val <- sum(res$arm == "validated" & res$design != "iv4_study_shift")
  expect_equal(nrow(folds), n_val * 3L + 10L * 4L + 10L * 3L)
  expect_true(all(folds$prob_degenerate %in% c(TRUE, FALSE)))

  # Importance: one row per feature per validated/nested run.
  expect_true(all(imp$n_folds_selected >= 0L))
  expect_true(all(!is.na(imp$feature)))
})
