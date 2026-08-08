.methodology_runner_path <- function() {
  source <- testthat::test_path(
    "..", "..", "inst", "analysis", "run_methodology_validation.R"
  )
  if (file.exists(source)) return(source)
  installed <- system.file(
    "analysis", "run_methodology_validation.R", package = "stablr"
  )
  if (nzchar(installed)) return(installed)
  stop("Could not locate the methodology validation runner.", call. = FALSE)
}

.late_fusion_runner_path <- function() {
  source <- testthat::test_path(
    "..", "..", "inst", "analysis", "run_late_fusion_validation.R"
  )
  if (file.exists(source)) return(source)
  installed <- system.file(
    "analysis", "run_late_fusion_validation.R", package = "stablr"
  )
  if (nzchar(installed)) return(installed)
  stop("Could not locate the late-fusion validation runner.", call. = FALSE)
}

test_that("methodology validation runner writes bounded artifact schema", {
  skip_on_cran()

  env <- new.env(parent = globalenv())
  sys.source(
    .methodology_runner_path(),
    envir = env
  )

  out <- tempfile("stablr-methodology-validation-")
  dir.create(out)

  artifacts <- env$run_methodology_validation(
    out = out,
    replicates = 1L,
    n_bootstraps = 2L,
    n_lambda = 2L,
    artificial_types = "random_permutation",
    scenario_ids = "null_independent",
    seed = 270627L
  )

  expect_true(file.exists(artifacts$replicates))
  expect_true(file.exists(artifacts$summary))
  expect_true(file.exists(artifacts$warnings))
  expect_true(file.exists(artifacts$parity))
  expect_true(file.exists(artifacts$manifest))

  manifest <- readLines(artifacts$manifest, warn = FALSE)
  expect_true(any(startsWith(manifest, "source_git_commit: ")))
  expect_true(any(startsWith(manifest, "end_git_commit: ")))
  expect_true(any(startsWith(manifest, "source_git_stable: ")))

  replicate_rows <- utils::read.csv(artifacts$replicates, stringsAsFactors = FALSE)
  expect_named(
    replicate_rows,
    c(
      "family", "scenario", "regime", "correlation", "profile", "replicate",
      "artificial_type", "data_seed", "fit_seed", "status", "n", "p", "n_signal",
      "n_bootstraps", "n_lambda", "fallback_random_permutation_warnings",
      "fallback_equi_warnings", "warning_count", "elapsed_sec",
      "n_selected", "true_positives", "false_positives", "empirical_fdp",
      "tpr", "min_fdp_plus", "fdp_threshold", "mean_max_real_score",
      "mean_max_artificial_score", "max_artificial_score",
      "selected_features", "error"
    )
  )
  expect_equal(nrow(replicate_rows), 1L)
  expect_true(all(replicate_rows$status == "ok"))
  expect_equal(length(unique(replicate_rows$data_seed)), 1L)
  expect_equal(length(unique(replicate_rows$fit_seed)), 1L)
  expect_identical(replicate_rows$family, "gaussian")

  summary_rows <- utils::read.csv(artifacts$summary, stringsAsFactors = FALSE)
  expect_true(all(c(
    "sd_selected", "se_selected", "sd_empirical_fdp", "se_empirical_fdp",
    "sd_tpr", "se_tpr", "sd_empirical_fdp_exceedance",
    "se_empirical_fdp_exceedance_rate",
    "se_fallback_random_permutation_rate", "se_fallback_equi_rate"
  ) %in% names(summary_rows)))

  parity_rows <- utils::read.csv(artifacts$parity, stringsAsFactors = FALSE)
  expect_true(all(c("metric", "reference", "observed", "abs_error", "status") %in% names(parity_rows)))

  gate_rows <- utils::read.csv(artifacts$gates, stringsAsFactors = FALSE)
  expect_true(all(c(
    "family", "scenario", "artificial_type", "cell_status", "gate", "pass"
  ) %in% names(gate_rows)))
})

test_that("release provenance detects dirty and moving source trees", {
  env <- new.env(parent = globalenv())
  sys.source(.methodology_runner_path(), envir = env)
  clean <- list(commit = "abc", tree = "tree-a", dirty = FALSE)
  dirty <- list(commit = "abc", tree = "tree-a", dirty = TRUE)
  moved <- list(commit = "def", tree = "tree-b", dirty = FALSE)

  expect_true(env$.git_provenance_is_clean(clean))
  expect_false(env$.git_provenance_is_clean(dirty))
  expect_true(env$.git_provenance_is_stable(clean, clean))
  expect_false(env$.git_provenance_is_stable(clean, moved))

  env$.git_provenance <- function(root) dirty
  expect_error(
    env$run_methodology_validation(
      out = tempfile("method-dirty-release-"),
      profile = "release"
    ),
    "requires a clean Git source tree"
  )

  late_env <- new.env(parent = globalenv())
  sys.source(.late_fusion_runner_path(), envir = late_env)
  expect_true(late_env$.git_provenance_is_available(clean))
  expect_false(late_env$.git_provenance_is_available(list(
    commit = NA_character_, tree = NA_character_, dirty = NA
  )))
  late_env$.source_provenance_start <- dirty
  expect_error(
    late_env$run_late_fusion_validation(
      out = tempfile("late-fusion-dirty-release-"),
      replicates = 1L,
      n_bootstraps = 2L,
      n_iter = 2L
    ),
    "requires a clean Git source tree"
  )
})

test_that("locked release methodology profile covers the advertised matrix", {
  env <- new.env(parent = globalenv())
  sys.source(.methodology_runner_path(), envir = env)

  settings <- env$.release_profile_settings()
  expect_identical(settings$replicates, 100L)
  expect_identical(settings$n_bootstraps, 1000L)
  expect_identical(settings$n_lambda, 30L)
  expect_setequal(
    settings$families,
    c("gaussian", "binomial", "multinomial", "cox")
  )
  expect_setequal(
    settings$artificial_types,
    c("random_permutation", "knockoff", "knockoff_equi", "knockoff_mvr")
  )

  scenarios <- env$.default_scenarios()
  expect_true(all(c("null", "signal") %in% scenarios$regime))
  expect_true(all(c(0, 0.7) %in% scenarios$correlation))
  expect_true(any(scenarios$profile == "high_dim" & scenarios$regime == "null"))
  expect_true(any(scenarios$profile == "high_dim" & scenarios$regime == "signal"))
})

test_that("methodology scenarios generate valid outcomes for every family", {
  env <- new.env(parent = globalenv())
  sys.source(.methodology_runner_path(), envir = env)
  scenario <- env$.default_scenarios()[3L, , drop = FALSE]

  generated <- lapply(
    c("gaussian", "binomial", "multinomial", "cox"),
    function(family) env$.simulate_scenario(scenario, family, seed = 101L)
  )
  expect_true(is.numeric(generated[[1L]]$y))
  expect_setequal(unique(generated[[2L]]$y), c(0, 1))
  expect_s3_class(generated[[3L]]$y, "factor")
  expect_setequal(levels(generated[[3L]]$y), c("A", "B", "C"))
  expect_s3_class(generated[[4L]]$y, "Surv")
  expect_true(all(vapply(generated, function(x) length(x$signal) == 5L, logical(1L))))
  expect_true(all(vapply(generated, function(x) {
    identical(rownames(x$x), if (is.matrix(x$y)) rownames(x$y) else names(x$y))
  }, logical(1L))))
})

test_that("release gates are per cell and use replicate-level observations", {
  env <- new.env(parent = globalenv())
  sys.source(.methodology_runner_path(), envir = env)

  scenarios <- data.frame(
    scenario = c("null", "signal"),
    regime = c("null", "signal"),
    profile = "low_dim",
    correlation = 0,
    n = 50L,
    p = 20L,
    n_signal = c(0L, 5L),
    stringsAsFactors = FALSE
  )
  make_cell <- function(scenario, artificial_type, selected, fdp, tpr,
                        status = "ok") {
    n <- length(selected)
    data.frame(
      family = "gaussian",
      scenario = scenario,
      regime = if (scenario == "null") "null" else "signal",
      profile = "low_dim",
      correlation = 0,
      replicate = seq_len(n),
      artificial_type = artificial_type,
      status = status,
      p = 20L,
      n_selected = selected,
      empirical_fdp = fdp,
      tpr = tpr,
      stringsAsFactors = FALSE
    )
  }
  rows <- rbind(
    make_cell("null", "good", rep(0L, 100L), rep(0, 100L), rep(NA, 100L)),
    make_cell("signal", "good", rep(4L, 100L), rep(0.05, 100L), rep(0.8, 100L)),
    make_cell("null", "bad", rep(1L, 100L), rep(1, 100L), rep(NA, 100L)),
    make_cell("signal", "bad", rep(10L, 100L), rep(0.5, 100L), rep(0.2, 100L))
  )

  gates <- env$.release_gate_table(
    rows,
    scenarios = scenarios,
    families = "gaussian",
    artificial_types = c("good", "bad"),
    expected_replicates = 100L
  )
  expect_true(all(gates$pass[
    gates$artificial_type == "good" &
      gates$gate != "null_selected_fraction"
  ]))
  expect_false(gates$pass[
    gates$artificial_type == "good" &
      gates$gate == "null_selected_fraction"
  ])
  expect_false(all(gates$pass[gates$artificial_type == "bad"]))
  expect_false(gates$pass[
    gates$artificial_type == "bad" & gates$gate == "null_select_any"
  ])
  expect_false(gates$pass[
    gates$artificial_type == "bad" & gates$gate == "signal_mean_fdp"
  ])
  expect_true(all(gates$cell_status == "complete"))

  selected_fraction <- gates$bound[
    gates$artificial_type == "bad" &
      gates$gate == "null_selected_fraction"
  ]
  selected_fraction_row <- gates[
    gates$artificial_type == "bad" &
      gates$gate == "null_selected_fraction",
    , drop = FALSE
  ]
  expect_true(is.na(selected_fraction))
  expect_identical(
    selected_fraction_row$gate_version,
    "SCI-04-null-selected-fraction-v2-draft"
  )
  expect_identical(selected_fraction_row$decision_status, "proposed_unaccepted")
  expect_false(selected_fraction_row$pass)
})

test_that("SCI-04 statistical decision record is versioned and unaccepted", {
  source_record <- testthat::test_path(
    "..", "..", "inst", "analysis", "contracts",
    "SCI-04-null-selected-fraction-v2.md"
  )
  record <- if (file.exists(source_record)) source_record else system.file(
    "analysis", "contracts", "SCI-04-null-selected-fraction-v2.md",
    package = "stablr"
  )
  expect_true(nzchar(record))
  text <- readLines(record, warn = FALSE)
  document <- paste(text, collapse = " ")
  expect_true(any(text == "Status: proposed, not accepted"))
  expect_true(grepl("Hoeffding", document, fixed = TRUE))
  expect_true(grepl("must not be implemented", document, fixed = TRUE))
})

test_that("missing and errored release cells fail explicitly", {
  env <- new.env(parent = globalenv())
  sys.source(.methodology_runner_path(), envir = env)
  scenarios <- data.frame(
    scenario = "null",
    regime = "null",
    profile = "low_dim",
    correlation = 0,
    n = 50L,
    p = 20L,
    n_signal = 0L,
    stringsAsFactors = FALSE
  )
  rows <- data.frame(
    family = "gaussian",
    scenario = "null",
    regime = "null",
    profile = "low_dim",
    correlation = 0,
    replicate = 1:2,
    artificial_type = "random_permutation",
    status = c("ok", "error"),
    p = 20L,
    n_selected = 0L,
    empirical_fdp = 0,
    tpr = NA_real_,
    stringsAsFactors = FALSE
  )
  gates <- env$.release_gate_table(
    rows,
    scenarios = scenarios,
    families = c("gaussian", "binomial"),
    artificial_types = "random_permutation",
    expected_replicates = 2L
  )

  expect_true(all(gates$cell_status[gates$family == "gaussian"] == "fit_not_ok"))
  expect_true(all(gates$cell_status[gates$family == "binomial"] == "missing_cell"))
  expect_false(any(gates$pass))
})

test_that("methodology validation rejects duplicate strategy entries", {
  skip_on_cran()

  env <- new.env(parent = globalenv())
  sys.source(
    .methodology_runner_path(),
    envir = env
  )

  expect_error(
    env$run_methodology_validation(
      out = tempfile("stablr-methodology-validation-"),
      artificial_types = c("random_permutation", "random_permutation")
    ),
    "must not contain duplicates"
  )
})

test_that("methodology validation pairs distinct strategies by data and fit seed", {
  skip_on_cran()
  skip_if_not_installed("knockoff")

  env <- new.env(parent = globalenv())
  sys.source(
    .methodology_runner_path(),
    envir = env
  )
  out <- tempfile("stablr-methodology-validation-paired-")
  artifacts <- env$run_methodology_validation(
    out = out,
    replicates = 1L,
    n_bootstraps = 2L,
    n_lambda = 2L,
    artificial_types = c("random_permutation", "knockoff_equi"),
    scenario_ids = "null_independent",
    seed = 270627L
  )
  rows <- utils::read.csv(artifacts$replicates, stringsAsFactors = FALSE)

  expect_equal(nrow(rows), 2L)
  expect_setequal(rows$artificial_type,
                  c("random_permutation", "knockoff_equi"))
  expect_equal(length(unique(rows$data_seed)), 1L)
  expect_equal(length(unique(rows$fit_seed)), 1L)
})

test_that("methodology replicate workers preserve deterministic scientific results", {
  skip_on_cran()
  skip_if(.Platform$OS.type == "windows")
  env <- new.env(parent = globalenv())
  sys.source(.methodology_runner_path(), envir = env)
  common <- list(
    replicates = 2L, n_bootstraps = 2L, n_lambda = 2L,
    families = "gaussian", artificial_types = "random_permutation",
    scenario_ids = "null_independent", seed = 991L
  )
  sequential <- do.call(env$run_methodology_validation,
                        c(list(out = tempfile("method-seq-"), workers = 1L), common))
  parallel <- do.call(env$run_methodology_validation,
                      c(list(out = tempfile("method-par-"), workers = 2L), common))
  a <- utils::read.csv(sequential$replicates, stringsAsFactors = FALSE)
  b <- utils::read.csv(parallel$replicates, stringsAsFactors = FALSE)
  a$elapsed_sec <- NULL
  b$elapsed_sec <- NULL
  expect_identical(a, b)
})

test_that("late-fusion release simulations predeclare well-posed class counts", {
  env <- new.env(parent = globalenv())
  sys.source(
    .late_fusion_runner_path(),
    envir = env
  )
  for (family in c("binomial", "multinomial")) {
    for (seed in 1:10) {
      simulated <- env$.simulate_well_posed(family, "signal", seed)
      expect_gte(min(table(simulated$train_y)), 10L)
      expect_gte(min(table(simulated$test_y)), 10L)
    }
  }
})
