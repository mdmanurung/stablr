# Tests for stabl_replication_profile() (1.0.0 spec A3; Phase 13 PR-E).
#
# Coverage: input validation (fail-closed), pooled_fit and effects
# validation, per feature-by-cohort eligibility reasons, hand-worked
# assembly math (within/LOSO frequencies, retention, sensitivity),
# inverse-variance direction pooling, held-out downstream-refit rows
# (gaussian MSE, binomial BER, empty-support reason), replayability,
# print smoke, and the read-only contract.

# ---- shared fixtures --------------------------------------------------------

.make_rep_data <- function() {
  set.seed(1)
  ids <- paste0("s", 1:60)
  x <- matrix(rnorm(60 * 8), 60, 8,
              dimnames = list(ids, paste0("f", 1:8)))
  y <- setNames(x[, 1] * 3 + rnorm(60, sd = 0.8), ids)
  studies <- setNames(rep(paste0("study", 1:3), each = 20L), ids)
  list(x = x, y = y, studies = studies)
}

fast_spec <- list(n_bootstraps = 50,
                  artificial_type = "random_permutation")

# ---- input validation -------------------------------------------------------

test_that("input validation is fail-closed", {
  d <- .make_rep_data()

  expect_error(
    stabl_replication_profile(unname(d$x), d$y, d$studies,
                              fit_spec = fast_spec, seed = 5),
    "`x` must be a matrix", fixed = TRUE
  )
  x_bad <- d$x
  rownames(x_bad) <- NULL
  expect_error(
    stabl_replication_profile(x_bad, d$y, d$studies,
                              fit_spec = fast_spec, seed = 5),
    "`x` must be a matrix", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, unname(d$y), d$studies,
                              fit_spec = fast_spec, seed = 5),
    "`y` must be named", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, unname(d$studies),
                              fit_spec = fast_spec, seed = 5),
    "`study_labels` must be named", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies[1:59],
                              fit_spec = fast_spec, seed = 5),
    "must have length nrow(x)", fixed = TRUE
  )
  two_studies <- setNames(rep(c("a", "b"), each = 30L), rownames(d$x))
  expect_error(
    stabl_replication_profile(d$x, d$y, two_studies,
                              fit_spec = fast_spec, seed = 5),
    "at least three studies", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies,
                              fit_spec = "notalist", seed = 5),
    "must be a named list", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies,
                              fit_spec = list(10), seed = 5),
    "fully named list", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies,
                              fit_spec = list(random_state = 1), seed = 5),
    "must not set", fixed = TRUE
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies,
                              fit_spec = fast_spec, seed = 5, bogus = 1),
    "Unused arguments", fixed = TRUE
  )
})

test_that("pooled_fit is validated and reused when valid", {
  d <- .make_rep_data()
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies, fit_spec = fast_spec,
                              pooled_fit = list(), seed = 5),
    "stabl_scores_", fixed = TRUE
  )
  pf_bad <- stabl_fit(
    d$x[, 1:4], d$y, lambda_grid = "auto", n_bootstraps = 10,
    random_state = 1
  )
  expect_error(
    stabl_replication_profile(d$x, d$y, d$studies, fit_spec = fast_spec,
                              pooled_fit = pf_bad, seed = 5),
    "feature names do not match", fixed = TRUE
  )
  pf <- stabl_fit(d$x, d$y, lambda_grid = "auto", n_bootstraps = 50,
                  random_state = 99)
  rp <- stabl_replication_profile(d$x, d$y, d$studies,
                                  fit_spec = fast_spec, pooled_fit = pf,
                                  seed = 5)
  expect_identical(rp$pooled_support, names(which(get_support(pf))))
})

test_that("effects validation is fail-closed", {
  d <- .make_rep_data()
  base <- data.frame(
    feature = c("f1", "f1", "f1"),
    study = c("study1", "study2", "study3"),
    estimate = c(1, 2, 3), se = c(0.5, 0.5, 1)
  )
  run <- function(eff) {
    stabl_replication_profile(d$x, d$y, d$studies, fit_spec = fast_spec,
                              effects = eff, seed = 5)
  }
  expect_error(run(base[, c("feature", "study", "estimate")]),
               "must be a data frame with columns", fixed = TRUE)
  bad <- base
  bad$estimate[1] <- Inf
  expect_error(run(bad), "finite estimates", fixed = TRUE)
  bad <- base
  bad$se[2] <- 0
  expect_error(run(bad), "strictly positive SEs", fixed = TRUE)
  bad <- base
  bad$feature[1] <- "nope"
  expect_error(run(bad), "unknown features or studies", fixed = TRUE)
  bad <- base
  bad$study[1] <- "study9"
  expect_error(run(bad), "unknown features or studies", fixed = TRUE)
  bad <- rbind(base, base[1, ])
  expect_error(run(bad), "duplicated (feature, study)", fixed = TRUE)
  rp <- run(base)
  expect_true(all(c("sign_concordance", "pooled_estimate", "het_Q") %in%
                    names(rp$features)))
})

# ---- eligibility ------------------------------------------------------------

test_that("eligibility reasons are exact per feature x cohort", {
  ids <- paste0("s", 1:12)
  x <- matrix(rnorm(12 * 3), 12, 3,
              dimnames = list(ids, c("f1", "f2", "f3")))
  studies <- setNames(rep(c("A", "B", "C"), each = 4L), ids)
  x[studies == "A", "f1"] <- NA_real_          # all missing in A
  x[studies == "B", "f2"] <- 7                 # constant in B
  elig <- stablr:::.replication_eligibility(x, as.factor(studies), c("A", "B", "C"))
  expect_equal(nrow(elig), 9L)
  f1 <- elig[elig$feature == "f1", ]
  expect_identical(f1$reason[f1$study == "A"], "all_missing")
  expect_true(all(f1$eligible[f1$study != "A"]))
  f2 <- elig[elig$feature == "f2", ]
  expect_identical(f2$reason[f2$study == "B"], "constant_in_cohort")
  expect_equal(sum(elig$eligible), 7L)
})

# ---- hand-worked assembly math ----------------------------------------------

.hand_worked <- function() {
  feature_names <- c("f1", "f2", "f3")
  studies <- c("A", "B", "C")
  eligibility <- data.frame(
    feature = rep(feature_names, each = 3L),
    study = rep(studies, 3L),
    eligible = c(TRUE, TRUE, TRUE,
                 TRUE, FALSE, TRUE,
                 TRUE, TRUE, TRUE),
    reason = c(NA, NA, NA, NA, "constant_in_cohort", NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )
  list(
    feature_names = feature_names, studies = studies,
    eligibility = eligibility,
    within_supports = list(A = c("f1", "f2"), B = "f1", C = c("f1", "f3")),
    loso_supports = list(A = "f1", B = c("f1", "f2"), C = character(0)),
    pooled_support = "f1"
  )
}

test_that("assembly math matches hand-worked values", {
  hw <- .hand_worked()
  out <- stablr:::.assemble_replication_features(
    hw$feature_names, hw$studies, hw$eligibility, hw$within_supports,
    hw$loso_supports, hw$pooled_support, effects = NULL
  )
  f1 <- out[out$feature == "f1", ]
  expect_equal(f1$n_measured_cohorts, 3L)
  expect_equal(f1$n_eligible_cohorts, 3L)
  expect_equal(f1$within_study_frequency, 1)
  expect_equal(f1$training_replication, 2 / 3)
  expect_equal(f1$fold_frequency, 2 / 3)
  expect_true(f1$in_pooled_support)
  expect_equal(f1$loso_retention, 2 / 3)
  expect_identical(f1$dropped_on_removal_of, "C")

  # f2 is ineligible in cohort B: within-frequency conditions on {A, C}
  f2 <- out[out$feature == "f2", ]
  expect_equal(f2$n_measured_cohorts, 3L)
  expect_equal(f2$n_eligible_cohorts, 2L)
  expect_equal(f2$within_study_frequency, 0.5)
  expect_equal(f2$training_replication, 1 / 3)
  expect_false(f2$in_pooled_support)
  expect_true(is.na(f2$loso_retention))
  expect_true(is.na(f2$dropped_on_removal_of))

  f3 <- out[out$feature == "f3", ]
  expect_equal(f3$within_study_frequency, 1 / 3)
  expect_equal(f3$training_replication, 0)
})

# ---- direction/effect replication math --------------------------------------

test_that("inverse-variance pooling and heterogeneity are exact", {
  effects <- data.frame(
    feature = c("f1", "f1", "f1", "f2"),
    study = c("A", "B", "C", "A"),
    estimate = c(1.0, 2.0, 3.0, 0.5),
    se = c(0.5, 0.5, 1.0, 0.25),
    stringsAsFactors = FALSE
  )
  r1 <- stablr:::.replication_direction_row("f1", effects)
  expect_equal(r1$n_estimates, 3L)
  expect_equal(r1$pooled_estimate, 15 / 9)
  expect_equal(r1$pooled_se, 1 / 3)
  expect_equal(r1$pooled_ci_lo, 15 / 9 - 1.96 / 3)
  expect_equal(r1$pooled_ci_hi, 15 / 9 + 1.96 / 3)
  expect_equal(r1$het_Q, 4.0)
  expect_equal(r1$het_df, 2L)
  expect_equal(r1$het_p, stats::pchisq(4, df = 2, lower.tail = FALSE))
  expect_equal(r1$sign_concordance, 1)

  r2 <- stablr:::.replication_direction_row("f2", effects)
  expect_equal(r2$n_estimates, 1L)
  expect_true(all(is.na(c(r2$pooled_estimate, r2$het_Q,
                          r2$sign_concordance))))

  discordant <- data.frame(
    feature = c("f3", "f3"), study = c("A", "B"),
    estimate = c(-1, 2), se = c(0.5, 0.5)
  )
  r3 <- stablr:::.replication_direction_row("f3", discordant)
  expect_equal(r3$sign_concordance, 0.5)
})

# ---- end-to-end structure ---------------------------------------------------

test_that("end-to-end gaussian profile has the pinned structure", {
  d <- .make_rep_data()
  rp <- stabl_replication_profile(d$x, d$y, d$studies,
                                  fit_spec = fast_spec, seed = 5)
  expect_s3_class(rp, "stabl_replication_profile")
  expect_equal(nrow(rp$features), 8L)
  expect_true(all(c(
    "feature", "n_measured_cohorts", "n_eligible_cohorts",
    "within_study_frequency", "training_replication", "fold_frequency",
    "in_pooled_support", "loso_retention", "dropped_on_removal_of"
  ) %in% names(rp$features)))
  expect_length(rp$estimands, 4L)
  expect_type(rp$pooled_support, "character")
  # deterministic for this seed: strong-signal f1 replicates perfectly
  expect_identical(rp$pooled_support,
                   c("f1", "f2", "f3", "f5", "f6", "f7", "f8"))
  expect_equal(rp$features$training_replication[1], 1)
  expect_equal(rp$features$loso_retention[1], 1)
  expect_true(all(rp$features$training_replication >= 0 &
                    rp$features$training_replication <= 1))
  expect_equal(nrow(rp$eligibility), 24L)
  expect_null(rp$heldout)
})

# ---- held-out performance ---------------------------------------------------

test_that("held-out gaussian rows carry MSE; empty support is explicit", {
  d <- .make_rep_data()
  rp <- stabl_replication_profile(d$x, d$y, d$studies,
                                  fit_spec = fast_spec,
                                  heldout_performance = TRUE, seed = 5)
  expect_equal(nrow(rp$heldout), 3L)
  expect_true(all(rp$heldout$metric_name == "mse"))
  expect_true(all(is.finite(rp$heldout$metric)))
  expect_true(all(rp$heldout$metric >= 0))
  expect_true(all(rp$heldout$n_panel >= 1L))

  rp0 <- stabl_replication_profile(
    d$x, d$y, d$studies,
    fit_spec = list(n_bootstraps = 10, hard_threshold = 1),
    heldout_performance = TRUE, seed = 5
  )
  expect_true(all(is.na(rp0$heldout$metric)))
  expect_true(all(rp0$heldout$reason == "empty_support"))
  expect_true(all(rp0$heldout$n_panel == 0L))
  expect_true(all(!rp0$features$in_pooled_support))
})

test_that("held-out binomial rows carry balanced error rate", {
  d <- .make_rep_data()
  yb <- setNames(as.integer(d$x[, 1] > 0), rownames(d$x))
  rp <- suppressWarnings(stabl_replication_profile(
    d$x, yb, d$studies,
    fit_spec = list(n_bootstraps = 10, family = "binomial"),
    heldout_performance = TRUE, seed = 5
  ))
  expect_equal(nrow(rp$heldout), 3L)
  expect_true(all(rp$heldout$metric_name == "balanced_error_rate"))
  expect_true(all(rp$heldout$metric >= 0 & rp$heldout$metric <= 1))
})

# ---- replayability ----------------------------------------------------------

test_that("same seed replays exactly; seeds are recorded", {
  d <- .make_rep_data()
  rp1 <- stabl_replication_profile(d$x, d$y, d$studies,
                                   fit_spec = fast_spec, seed = 5)
  rp2 <- stabl_replication_profile(d$x, d$y, d$studies,
                                   fit_spec = fast_spec, seed = 5)
  expect_identical(rp1$features, rp2$features)
  expect_identical(rp1$pooled_support, rp2$pooled_support)
  expect_identical(rp1$seed, 5L)
  expect_length(rp1$fit_seeds, 7L)  # 2 * 3 studies + pooled
})

# ---- print ------------------------------------------------------------------

test_that("print method is readable and invisible", {
  d <- .make_rep_data()
  rp <- stabl_replication_profile(d$x, d$y, d$studies,
                                  fit_spec = fast_spec,
                                  heldout_performance = TRUE, seed = 5)
  out <- capture.output(v <- withVisible(print(rp)))
  expect_false(v$visible)
  expect_identical(v$value, rp)
  expect_true(any(grepl("cohort replication profile", out, fixed = TRUE)))
  expect_true(any(grepl("LOSO retention", out, fixed = TRUE)))
  expect_true(any(grepl("held-out performance (mse)", out, fixed = TRUE)))

  fake <- structure(
    list(
      features = data.frame(
        feature = "f1", training_replication = NA_real_,
        in_pooled_support = FALSE, loso_retention = NA_real_,
        dropped_on_removal_of = NA_character_
      ),
      pooled_support = character(0), heldout = NULL,
      study_labels = factor(rep(c("A", "B", "C"), each = 2L))
    ),
    class = "stabl_replication_profile"
  )
  expect_output(print(fake), "no eligible feature-fold pairs",
                fixed = TRUE)
})

# ---- read-only contract -----------------------------------------------------

test_that("inputs are not modified", {
  d <- .make_rep_data()
  x0 <- d$x
  y0 <- d$y
  invisible(stabl_replication_profile(d$x, d$y, d$studies,
                                      fit_spec = fast_spec, seed = 5))
  expect_identical(d$x, x0)
  expect_identical(d$y, y0)
})
