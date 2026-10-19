# Tests for stabl_null_check() (1.0.0 spec, A2 — outcome-null calibration
# workbench). Covers: exchangeability-unit design legality (refusals, never
# silent global permutation), null-generator correctness per design,
# replayable seeds, bounded failure recording, decoy-realization
# sensitivity structure, synthetic-vs-observational truth separation, and
# recovery of known permutation-null behavior on a tiny problem.

# ---- helpers ----------------------------------------------------------------

make_nv_data <- function(n = 48L, p = 6L, signal = TRUE) {
  set.seed(11)
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  y <- if (signal) {
    setNames(2 * x[, 1] - 1.5 * x[, 2] + rnorm(n, sd = 0.5), rownames(x))
  } else {
    setNames(rnorm(n), rownames(x))
  }
  list(x = x, y = y)
}

fast_spec <- list(n_bootstraps = 10, lambda_grid = "auto",
                  artificial_type = "random_permutation")

# ---- design legality --------------------------------------------------------

test_that("unsupported or contradicted designs are refused", {
  d <- make_nv_data()
  blk <- rep(letters[1:4], each = 12L)
  clu <- rep(letters[1:6], each = 8L)

  # the key safety refusal: global permutation with declared structure
  expect_error(
    stabl_null_check(d$x, d$y, design = "independent", clusters = clu,
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "global permutation|GLOBAL permutation"
  )
  expect_error(
    stabl_null_check(d$x, d$y, design = "independent", blocks = blk,
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "global permutation|GLOBAL permutation"
  )
  # missing required structure
  expect_error(
    stabl_null_check(d$x, d$y, design = "within_block",
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "blocks"
  )
  expect_error(
    stabl_null_check(d$x, d$y, design = "whole_cluster",
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "clusters"
  )
  # degenerate structures
  expect_error(
    stabl_null_check(d$x, d$y, design = "within_block",
                     blocks = rep("a", nrow(d$x)),
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "at least two blocks"
  )
  expect_error(
    stabl_null_check(d$x, d$y, design = "whole_cluster",
                     clusters = rep(letters[1:2], each = 24L),
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "at least three clusters"
  )
  # unequal cluster sizes cannot take whole-vector permutation
  expect_error(
    stabl_null_check(d$x, d$y, design = "whole_cluster",
                     clusters = c(rep("a", 20L), rep("b", 16L), rep("c", 12L)),
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "equal-sized clusters"
  )
  # simulated requires an explicit simulator
  expect_error(
    stabl_null_check(d$x, d$y, design = "simulated",
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "null_simulator"
  )
  # fit_spec may not hijack data or seeds; unused dots are an error
  expect_error(
    stabl_null_check(d$x, d$y, fit_spec = list(random_state = 1),
                     n_reps = 2, seed = 1),
    "random_state"
  )
  expect_error(
    stabl_null_check(d$x, d$y, fit_spec = fast_spec, n_reps = 2,
                     seed = 1, bogus = 1),
    "Unused arguments"
  )
  # misaligned structure vectors
  expect_error(
    stabl_null_check(d$x, d$y, design = "within_block",
                     blocks = blk[1:10],
                     fit_spec = fast_spec, n_reps = 2, seed = 1),
    "nrow"
  )
})

# ---- null generator correctness ----------------------------------------------

test_that("within_block permutation preserves within-block multisets", {
  d <- make_nv_data()
  blk <- factor(rep(letters[1:4], each = 12L))
  out <- stablr:::.generate_null_outcome(
    d$y, d$x, "within_block", blocks = blk, clusters = NULL,
    null_simulator = NULL, rep_id = 1L, seed = 99L
  )
  for (b in levels(blk)) {
    expect_identical(sort(unname(out[blk == b])), sort(unname(d$y[blk == b])))
  }
  # but the null is not the identity permutation (with this seed, safely)
  expect_false(identical(unname(out), unname(d$y)))
})

test_that("whole_cluster permutation moves whole outcome vectors", {
  d <- make_nv_data()
  clu <- factor(rep(letters[1:6], each = 8L))
  out <- stablr:::.generate_null_outcome(
    d$y, d$x, "whole_cluster", blocks = NULL, clusters = clu,
    null_simulator = NULL, rep_id = 1L, seed = 5L
  )
  # every cluster's outcome vector after permutation equals some original
  # cluster's outcome vector (within-cluster order preserved)
  original <- lapply(levels(clu), function(k) unname(d$y[clu == k]))
  for (k in levels(clu)) {
    moved <- unname(out[clu == k])
    expect_true(any(vapply(original, identical, logical(1L), moved)))
  }
})

test_that("independent permutation is a global permutation of y", {
  d <- make_nv_data()
  out <- stablr:::.generate_null_outcome(
    d$y, d$x, "independent", blocks = NULL, clusters = NULL,
    null_simulator = NULL, rep_id = 1L, seed = 3L
  )
  expect_identical(sort(unname(out)), sort(unname(d$y)))
})

# ---- end-to-end behavior -----------------------------------------------------

test_that("null check recovers known permutation-null behavior", {
  d <- make_nv_data()
  nv <- stabl_null_check(d$x, d$y, design = "independent",
                         fit_spec = fast_spec, n_reps = 4, seed = 21)
  expect_s3_class(nv, "stabl_null_validation")
  expect_identical(nv$truth_status, "induced_global_null_permutation")
  expect_identical(nv$n_successful, 4L)
  expect_identical(nv$n_failures, 0L)
  expect_length(nv$panel_sizes, 4L)
  expect_true(all(nv$panel_sizes >= 0L))
  expect_true(nv$p_any_selection >= 0 && nv$p_any_selection <= 1)
  # per-replicate provenance and replay seeds are stored
  seeds <- vapply(nv$replicates, function(r) r$perm_seed, integer(1L))
  expect_identical(anyDuplicated(seeds), 0L)
  expect_true(all(vapply(nv$replicates, function(r) {
    r$decoy_requested == "random_permutation"
  }, logical(1L))))
  # the signal features are not recovered under the null the way the
  # observed fit recovers them (stronger observed fit as reference)
  fit_obs <- stabl_fit(d$x, d$y, n_bootstraps = 50, lambda_grid = "auto",
                       artificial_type = "random_permutation",
                       random_state = 21)
  expect_true("f1" %in% names(which(get_support(fit_obs))))
  freq <- nv$feature_frequency
  if ("f1" %in% names(freq)) {
    expect_lt(unname(freq[["f1"]]), 1)
  }
})

test_that("replicates are exactly replayable from the recorded seed", {
  d <- make_nv_data()
  a <- stabl_null_check(d$x, d$y, design = "independent",
                        fit_spec = fast_spec, n_reps = 3, seed = 77)
  b <- stabl_null_check(d$x, d$y, design = "independent",
                        fit_spec = fast_spec, n_reps = 3, seed = 77)
  expect_identical(a$panel_sizes, b$panel_sizes)
  expect_identical(
    vapply(a$replicates, function(r) r$perm_seed, integer(1L)),
    vapply(b$replicates, function(r) r$perm_seed, integer(1L))
  )
  c <- stabl_null_check(d$x, d$y, design = "independent",
                        fit_spec = fast_spec, n_reps = 3, seed = 78)
  expect_false(identical(
    vapply(a$replicates, function(r) r$perm_seed, integer(1L)),
    vapply(c$replicates, function(r) r$perm_seed, integer(1L))
  ))
})

test_that("simulated design marks synthetic truth and uses the simulator", {
  d <- make_nv_data()
  calls <- integer(0)
  sim <- function(y, x, rep_id, seed) {
    calls <<- c(calls, rep_id)
    setNames(rnorm(length(y)), names(y))
  }
  nv <- stabl_null_check(d$x, d$y, design = "simulated", null_simulator = sim,
                         fit_spec = fast_spec, n_reps = 3, seed = 9)
  expect_identical(nv$truth_status, "synthetic_global_null")
  expect_identical(as.integer(calls), 1:3)
  expect_identical(nv$n_successful, 3L)
})

test_that("failures are recorded bounded and fail-closed", {
  d <- make_nv_data()
  # every fit fails -> abort once failures exceed max_failures
  expect_error(
    stabl_null_check(d$x, d$y, design = "independent",
                     fit_spec = list(base_learner = "nonexistent_learner"),
                     n_reps = 4, seed = 1, max_failures = 2),
    "max_failures"
  )
  # one null-generation failure is recorded, not dropped
  sim <- function(y, x, rep_id, seed) {
    if (rep_id == 2L) stop("boom")
    setNames(rnorm(length(y)), names(y))
  }
  nv <- stabl_null_check(d$x, d$y, design = "simulated", null_simulator = sim,
                         fit_spec = fast_spec, n_reps = 4, seed = 9,
                         max_failures = 3)
  expect_identical(nv$n_failures, 1L)
  expect_identical(nv$n_successful, 3L)
  expect_identical(nv$failures$stage, "null_generation")
  expect_match(nv$failures$message, "boom")
  expect_identical(nv$failures$rep_id, 2L)
})

test_that("decoy_reps separates decoy-realization sensitivity", {
  d <- make_nv_data()
  nv <- stabl_null_check(d$x, d$y, design = "independent",
                         fit_spec = fast_spec, n_reps = 2, decoy_reps = 2,
                         seed = 31)
  expect_identical(nv$n_successful, 4L)
  # same permutation seed within a rep, distinct fit seeds
  r <- nv$replicates
  expect_identical(r[[1]]$perm_seed, r[[2]]$perm_seed)
  expect_false(identical(r[[1]]$fit_seed, r[[2]]$fit_seed))
  expect_identical(r[[3]]$perm_seed, r[[4]]$perm_seed)
})

test_that("print and plot methods work", {
  skip_if_not_installed("ggplot2")
  d <- make_nv_data()
  nv <- stabl_null_check(d$x, d$y, design = "independent",
                         fit_spec = fast_spec, n_reps = 3, seed = 41)
  out <- capture.output(print(nv))
  expect_true(any(grepl("outcome-null validation", out)))
  expect_true(any(grepl("P(any selection", out, fixed = TRUE)))
  expect_true(any(grepl("permutation induces a global null", out)))
  plots <- plot(nv, which = c("panel_sizes", "reported_vs_realized",
                              "saturation"))
  expect_named(plots, c("panel_sizes", "reported_vs_realized", "saturation"))
  expect_true(all(vapply(plots, inherits, logical(1L), "ggplot")))
  # simulated design does not print the observational caveat
  sim <- function(y, x, rep_id, seed) setNames(rnorm(length(y)), names(y))
  nv2 <- stabl_null_check(d$x, d$y, design = "simulated",
                          null_simulator = sim, fit_spec = fast_spec,
                          n_reps = 2, seed = 43)
  out2 <- capture.output(print(nv2))
  expect_false(any(grepl("permutation induces a global null", out2)))
})

test_that("independent-design null survives stabl_fit name alignment (C1 regression)", {
  # Regression: sample(y) carried each value's name along, and stabl_fit's
  # name-based outcome alignment (y[rownames(x)]) then exactly undid the
  # permutation — every "null" replicate fitted the observed data.
  d <- make_nv_data()
  ids <- rownames(d$x)
  out <- stablr:::.generate_null_outcome(
    d$y, d$x, "independent", blocks = NULL, clusters = NULL,
    null_simulator = NULL, rep_id = 1L, seed = 3L
  )
  aligned <- out[ids]  # what stabl_fit does with the null outcome
  # Still a global permutation of the values...
  expect_identical(sort(unname(aligned)), sort(unname(d$y)))
  # ...with names at their ORIGINAL positions...
  expect_identical(names(aligned), ids)
  # ...so name alignment does NOT restore the observed pairing.
  expect_false(identical(unname(aligned), unname(d$y)))

  # End-to-end: with the same seed and data, a fit on the aligned null must
  # differ from a fit on the observed outcome (previously bit-identical).
  fit_obs <- stabl_fit(d$x, d$y, lambda_grid = data.frame(lambda = c(0.2, 0.1)),
                       n_bootstraps = 3L, artificial_type = NULL,
                       hard_threshold = 0.5, random_state = 8L)
  fit_null <- stabl_fit(d$x, aligned, lambda_grid = data.frame(lambda = c(0.2, 0.1)),
                        n_bootstraps = 3L, artificial_type = NULL,
                        hard_threshold = 0.5, random_state = 8L)
  expect_false(identical(fit_obs$stabl_scores_, fit_null$stabl_scores_))
})

test_that("systematic null-generation failures abort fail-closed (C2 regression)", {
  # Regression: generation failures bypassed the max_failures abort via
  # `next`, so a systematically failing null_simulator returned
  # n_successful = 0 instead of aborting.
  d <- make_nv_data()
  expect_error(
    stabl_null_check(d$x, d$y, design = "simulated",
                     null_simulator = function(...) stop("always fails"),
                     fit_spec = fast_spec, n_reps = 3L, seed = 1L,
                     max_failures = 0L),
    "max_failures"
  )
})
