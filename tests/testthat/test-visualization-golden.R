library(stablr)

# ── Golden-path rendering tests for the diagnostic plot functions ────────────
#
# vdiffr is deliberately NOT used.  Instead each test asserts that the plot
# (a) builds without error and returns a ggplot object, and (b) carries the
# correct data: layer data / built plot data are checked against the values
# stored on the fitted "stabl_fit" object (or the raw input matrices).

# ── Shared fixtures ───────────────────────────────────────────────────────────

# Small real STABL fit (n = 60, p = 8, n_bootstraps = 15, explicit seed).
# f1 and f2 carry linear signal so the hard-threshold support is non-trivial.
# suppressWarnings mirrors the convention in test-phase7.R / test-stabl-fit.R.
.make_golden_fit <- function(seed = 11L,
                             artificial_type = NULL,
                             hard_threshold = 0.5) {
  set.seed(seed)
  n <- 60L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y <- setNames(2 * x[, "f1"] - 1.5 * x[, "f2"] + rnorm(n, sd = 0.5),
                rownames(x))
  lam_grid <- data.frame(lambda = seq(0.05, 0.3, length.out = 4L))
  suppressWarnings(stabl_fit(
    x               = x,
    y               = y,
    lambda_grid     = lam_grid,
    family          = "gaussian",
    n_bootstraps    = 15L,
    artificial_type = artificial_type,
    hard_threshold  = hard_threshold,
    random_state    = seed
  ))
}

# Predictor matrix plus binary factor and continuous outcomes for the
# feature-distribution plots (n = 60, p = 8, explicit seed).
.make_golden_xy <- function(seed = 12L) {
  set.seed(seed)
  n <- 60L; p <- 8L
  x <- matrix(rnorm(n * p), n, p,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("f", seq_len(p))))
  y_factor   <- factor(rep(c("ctrl", "case"), length.out = n))
  y_numeric  <- setNames(2 * x[, "f1"] - 1.5 * x[, "f2"] + rnorm(n, sd = 0.5),
                         rownames(x))
  list(x = x, y_factor = y_factor, y_numeric = y_numeric)
}

# Collect the pre-build data of all line layers of a ggplot object.
.line_layer_data <- function(p) {
  line_layers <- Filter(function(l) inherits(l$geom, "GeomLine"), p$layers)
  lapply(line_layers, `[[`, "data")
}

# Extract the built data of the first layer whose built columns include
# `column` (same ggplot_build-based approach as test-phase7.R).
.built_layer_with <- function(p, column) {
  built <- ggplot2::ggplot_build(p)
  hits <- vapply(built$data, function(d) column %in% names(d), TRUE)
  built$data[[which(hits)[1]]]
}

# ── plot_stabl_path ───────────────────────────────────────────────────────────

test_that("plot_stabl_path renders the hard-threshold stability path", {
  skip_if_not_installed("ggplot2")
  fit <- .make_golden_fit(seed = 11L, artificial_type = NULL,
                          hard_threshold = 0.5)
  p   <- plot_stabl_path(fit)
  expect_s3_class(p, "ggplot")

  # One row per (feature, lambda) across the line layers, with no duplicates.
  line_data <- do.call(rbind, .line_layer_data(p))
  n_lambda  <- ncol(fit$stabl_scores_)
  expect_equal(nrow(line_data), length(fit$feature_names) * n_lambda)
  expect_equal(anyDuplicated(line_data[, c("feature", "lambda_idx")]), 0L)
  expect_setequal(unique(line_data$feature), fit$feature_names)

  # Fitted stabl_scores_ values, per feature, ordered along the canonical
  # decreasing-lambda path (the plot direction must not depend on the row
  # order of the supplied grid).
  ord <- order(-fit$fitted_lambda_grid$lambda)
  for (i in seq_along(fit$feature_names)) {
    rows <- line_data[line_data$feature == fit$feature_names[i], ]
    rows <- rows[order(rows$lambda_idx), ]
    expect_equal(rows$lambda_idx, seq_len(n_lambda))
    expect_equal(rows$score, unname(fit$stabl_scores_[i, ord]))
  }
  expect_setequal(unique(line_data$lambda), fit$fitted_lambda_grid$lambda)
  # The rendered lambda values decrease along the path.
  lam_path <- unique(line_data[order(line_data$lambda_idx), "lambda"])
  expect_true(all(diff(lam_path) < 0))

  # The selected flag on each row matches get_support() for that feature.
  support_map <- setNames(get_support(fit), fit$feature_names)
  expect_identical(unname(support_map[line_data$feature]), line_data$selected)

  # Threshold geometry: dashed hline at the fitted hard threshold plus its
  # text annotation.
  hline <- .built_layer_with(p, "yintercept")
  expect_equal(hline$yintercept, fit$hard_threshold)
  label <- .built_layer_with(p, "label")
  expect_equal(label$label, sprintf("Hard threshold = %.2f", fit$hard_threshold))

  # Axis labels and title.
  expect_equal(p$labels$title, "STABL Stability Path")
  expect_equal(p$labels$y, "Frequency of selection")
})

test_that("plot_stabl_path honours new_hard_threshold override", {
  skip_if_not_installed("ggplot2")
  fit <- .make_golden_fit(seed = 11L, artificial_type = NULL,
                          hard_threshold = 0.5)
  p   <- plot_stabl_path(fit, new_hard_threshold = 0.9)
  expect_s3_class(p, "ggplot")

  hline <- .built_layer_with(p, "yintercept")
  expect_equal(hline$yintercept, 0.9)
  label <- .built_layer_with(p, "label")
  expect_equal(label$label, "Threshold = 0.90")

  # The selected flag follows the overridden support, not the stored one.
  line_data   <- do.call(rbind, .line_layer_data(p))
  support_map <- setNames(get_support(fit, new_hard_threshold = 0.9),
                          fit$feature_names)
  expect_identical(unname(support_map[line_data$feature]), line_data$selected)
})

test_that("plot_stabl_path overlays artificial features and FDP+ threshold", {
  skip_if_not_installed("ggplot2")
  fit <- .make_golden_fit(seed = 13L, artificial_type = "random_permutation",
                          hard_threshold = NULL)
  p   <- plot_stabl_path(fit)
  expect_s3_class(p, "ggplot")

  line_data <- .line_layer_data(p)
  real_data <- do.call(rbind, Filter(function(d) "selected" %in% names(d),
                                     line_data))
  art_data  <- do.call(rbind, Filter(function(d) !"selected" %in% names(d),
                                     line_data))

  # Real features: one row per (feature, lambda), selected flag from support.
  n_lambda <- ncol(fit$stabl_scores_)
  ord <- order(-fit$fitted_lambda_grid$lambda)
  expect_equal(nrow(real_data), length(fit$feature_names) * n_lambda)
  support_map <- setNames(get_support(fit), fit$feature_names)
  expect_identical(unname(support_map[real_data$feature]), real_data$selected)

  # Artificial features are overlaid as their own dotted lines, one row per
  # (artificial feature, lambda), named "artificial.<i>".
  n_art <- nrow(fit$stabl_scores_artificial_)
  expect_equal(nrow(art_data), n_art * n_lambda)
  expect_setequal(unique(art_data$feature),
                  paste0("artificial.", seq_len(n_art)))
  for (i in seq_len(n_art)) {
    rows <- art_data[art_data$feature == paste0("artificial.", i), ]
    rows <- rows[order(rows$lambda_idx), ]
    expect_equal(rows$score, unname(fit$stabl_scores_artificial_[i, ord]))
  }

  # Threshold geometry: dashed hline at the FDP+-optimal threshold.
  hline <- .built_layer_with(p, "yintercept")
  expect_equal(hline$yintercept, fit$fdr_min_threshold_)
  label <- .built_layer_with(p, "label")
  expect_equal(label$label,
               sprintf("FDP+ threshold = %.2f", fit$fdr_min_threshold_))
})

# ── plot_fdr_graph ────────────────────────────────────────────────────────────

test_that("plot_stabl_path direction is independent of grid row order", {
  skip_if_not_installed("ggplot2")
  # The same lambdas supplied in ascending and descending row order must give
  # the same rendered path (canonical decreasing-lambda direction).
  make_fit <- function(lam_grid) {
    set.seed(11L)
    n <- 60L; p <- 8L
    x <- matrix(rnorm(n * p), n, p,
                 dimnames = list(paste0("s", seq_len(n)),
                                 paste0("f", seq_len(p))))
    y <- setNames(2 * x[, "f1"] - 1.5 * x[, "f2"] + rnorm(n, sd = 0.5),
                  rownames(x))
    suppressWarnings(stabl_fit(x, y, lambda_grid = lam_grid, family = "gaussian",
                               n_bootstraps = 15L, random_state = 11L))
  }
  lam_asc  <- data.frame(lambda = seq(0.05, 0.3, length.out = 4L))
  p_asc  <- plot_stabl_path(make_fit(lam_asc))
  p_desc <- plot_stabl_path(make_fit(lam_asc[rev(seq_len(nrow(lam_asc))), ,
                                         drop = FALSE]))
  b_asc  <- ggplot2::ggplot_build(p_asc)$data
  b_desc <- ggplot2::ggplot_build(p_desc)$data
  expect_equal(length(b_asc), length(b_desc))
  for (i in seq_along(b_asc)) {
    expect_equal(b_asc[[i]], b_desc[[i]],
                 info = paste("plot layer", i, "differs across grid order"))
  }
})

test_that("plot_fdr_graph plots the FDP+ curve stored in the fit", {
  skip_if_not_installed("ggplot2")
  fit <- .make_golden_fit(seed = 13L, artificial_type = "random_permutation",
                          hard_threshold = NULL)
  p   <- plot_fdr_graph(fit)
  expect_s3_class(p, "ggplot")

  # The plot data is exactly the (threshold, FDR) curve from the fit.
  expect_equal(nrow(p$data), length(fit$fdr_threshold_range))
  expect_equal(p$data$threshold, fit$fdr_threshold_range)
  expect_equal(p$data$FDR, fit$FDRs_)

  # Vertical dashed line at the FDP+-optimal threshold (first minimum of the
  # curve, matching the fit's stored fdr_min_threshold_).
  vline <- .built_layer_with(p, "xintercept")
  expect_equal(vline$xintercept, fit$fdr_min_threshold_)
  expect_equal(vline$xintercept,
               fit$fdr_threshold_range[which.min(fit$FDRs_)])

  # Horizontal dashed reference line at the default fdr_target.
  hline <- .built_layer_with(p, "yintercept")
  expect_equal(hline$yintercept, 0.05)

  # Axis labels and title.
  expect_equal(p$labels$title, "FDR Estimate")
  expect_equal(p$labels$x, "Stability threshold")
  expect_equal(p$labels$y, "FDR estimate")
})

test_that("plot_fdr_graph omits the target line when fdr_target = NULL", {
  skip_if_not_installed("ggplot2")
  fit <- .make_golden_fit(seed = 13L, artificial_type = "random_permutation",
                          hard_threshold = NULL)
  p   <- plot_fdr_graph(fit, fdr_target = NULL)
  expect_s3_class(p, "ggplot")
  n_hline <- sum(vapply(p$layers, function(l) inherits(l$geom, "GeomHline"),
                        TRUE))
  expect_equal(n_hline, 0L)
})

# ── boxplot_features ──────────────────────────────────────────────────────────

test_that("boxplot_features maps requested feature values and facets by feature", {
  skip_if_not_installed("ggplot2")
  d <- .make_golden_xy(seed = 12L)
  p <- boxplot_features(c("f1", "f2"), d$x, d$y_factor)
  expect_s3_class(p, "ggplot")

  # Long data: one row per (feature, sample) with the raw matrix values.
  expect_equal(nrow(p$data), 2L * nrow(d$x))
  expect_setequal(unique(p$data$feature), c("f1", "f2"))
  expect_equal(p$data$value[p$data$feature == "f1"], as.numeric(d$x[, "f1"]))
  expect_equal(p$data$value[p$data$feature == "f2"], as.numeric(d$x[, "f2"]))
  expect_equal(as.character(p$data$outcome[p$data$feature == "f1"]),
               as.character(d$y_factor))
  expect_equal(as.character(p$data$outcome[p$data$feature == "f2"]),
               as.character(d$y_factor))

  # Aesthetic mapping and facet variable.
  expect_equal(rlang::as_label(p$mapping$x), "outcome")
  expect_equal(rlang::as_label(p$mapping$y), "value")
  expect_equal(rlang::as_label(p$mapping$fill), "outcome")
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(rlang::as_label(p$facet$params$facets[[1]]), "feature")

  # Axis labels.
  expect_equal(p$labels$y, "Value")
})

test_that("boxplot_features silently skips features missing from x", {
  skip_if_not_installed("ggplot2")
  d <- .make_golden_xy(seed = 12L)
  p <- boxplot_features(c("f1", "f99"), d$x, d$y_factor)
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$feature), "f1")
  expect_equal(nrow(p$data), nrow(d$x))
})

# ── scatterplot_features ──────────────────────────────────────────────────────

test_that("scatterplot_features maps the two requested features", {
  skip_if_not_installed("ggplot2")
  d <- .make_golden_xy(seed = 12L)
  p <- scatterplot_features(c("f1", "f2"), d$x, d$y_numeric)
  expect_s3_class(p, "ggplot")

  # Long data: one row per (feature, sample) with the raw matrix values and
  # the continuous outcome.
  expect_equal(nrow(p$data), 2L * nrow(d$x))
  expect_setequal(unique(p$data$feature), c("f1", "f2"))
  expect_equal(p$data$value[p$data$feature == "f1"], as.numeric(d$x[, "f1"]))
  expect_equal(p$data$value[p$data$feature == "f2"], as.numeric(d$x[, "f2"]))
  expect_equal(p$data$outcome[p$data$feature == "f1"],
               as.numeric(d$y_numeric))
  expect_equal(p$data$outcome[p$data$feature == "f2"],
               as.numeric(d$y_numeric))

  # Aesthetic mapping: feature value on x, outcome on y, faceted by feature.
  expect_equal(rlang::as_label(p$mapping$x), "value")
  expect_equal(rlang::as_label(p$mapping$y), "outcome")
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(rlang::as_label(p$facet$params$facets[[1]]), "feature")

  # Axis labels.
  expect_equal(p$labels$x, "Feature value")
  expect_equal(p$labels$y, "Outcome")
})
