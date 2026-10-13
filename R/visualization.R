# visualization.R — ggplot2 diagnostic and result plots for stablr
# R port of the plot functions in stabl/stabl.py and stabl/visualization.py

.require_ggplot2 <- function() .require_pkg("ggplot2", "for plotting")

# Build a long data frame from a scores matrix (one row per feature × lambda).
# feat_names: character vector of length n.
# scores:     n × n_lambda numeric matrix.
# If `selected` is non-NULL it must be a logical vector of length n; a column
# is appended with the per-feature selection indicator (recycled across lambdas).
.scores_to_long_df <- function(feat_names, scores, lambda_vals, lam_grid, has_alpha,
                                selected = NULL) {
  n     <- length(feat_names)
  n_lam <- ncol(scores)
  lst   <- vector("list", n)
  for (i in seq_len(n)) {
    lst[[i]] <- data.frame(
      feature    = feat_names[i],
      lambda_idx = seq_len(n_lam),
      lambda     = lambda_vals,
      score      = unname(scores[i, ]),
      alpha      = if (has_alpha) unname(lam_grid$alpha) else NA_real_,
      stringsAsFactors = FALSE
    )
    if (!is.null(selected)) lst[[i]]$selected <- unname(selected[i])
  }
  do.call(rbind, lst)
}

# Color palette matching the Python STABL implementation
.stablr_colors <- list(
  selected   = "#C41E3A",  # cardinal red for stable/selected features
  noise      = "#4D4F53",  # dark grey for noise features
  artificial = "#9E9E9E",  # medium grey for artificial features
  threshold  = "#000000",  # black for threshold lines
  chance     = "#4D4F53",  # dark grey for chance diagonal
  ci         = "#e3819d"   # pink for CI bands
)

# ── Stability path ────────────────────────────────────────────────────────────

#' Plot the STABL Stability Path
#'
#' Produces a ggplot2 line chart showing the stability score (frequency of
#' selection) of each feature across the regularisation path.  Selected
#' features are highlighted in red; noise features are shown in dark grey.
#' When artificial features were used during fitting, their stability path is
#' overlaid as thin dotted grey lines.  The FDP+-optimal (or hard) threshold
#' is shown as a dashed horizontal line.
#'
#' When the lambda grid contains an `alpha` column (elastic-net mixed-alpha
#' path), the plot is automatically faceted by `alpha`.
#'
#' The path is always drawn in decreasing-lambda order (high penalty on the
#' left, low penalty on the right) within each `alpha` level, independently of
#' the row order of the supplied lambda grid.  The x-axis is the position along
#' that ordered path, so auto grids (log-spaced lambda) are evenly spaced.
#'
#' @param object A fitted `"stabl_fit"` object from [stabl_fit()].
#' @param new_hard_threshold Numeric in `(0, 1]` or `NULL`.  When supplied,
#'   overrides the threshold stored in `object`.
#' @param title Character; plot title (default `"STABL Stability Path"`).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(40 * 6), 40, 6,
#'                dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
#'   y <- setNames(rnorm(40), rownames(x))
#'   fit <- stabl_fit(x, y,
#'                    lambda_grid = data.frame(lambda = c(0.3, 0.1, 0.05)),
#'                    n_bootstraps = 6L, hard_threshold = 0.3, random_state = 1L)
#'   plot_stabl_path(fit)
#' }
#' @export
plot_stabl_path <- function(object, new_hard_threshold = NULL,
                            title = "STABL Stability Path") {
  .check_fitted_stabl(object)
  .require_ggplot2()
  title <- .validate_nonempty_scalar_character(title, "title")

  scores     <- object$stabl_scores_
  feat_names <- object$feature_names
  lam_grid   <- object$fitted_lambda_grid
  n_lambda   <- ncol(scores)

  threshold <- .resolve_stabl_threshold(object, new_hard_threshold, on_missing = "null")
  support   <- get_support(object, new_hard_threshold = new_hard_threshold)

  # Determine whether to facet by alpha
  has_alpha <- "alpha" %in% names(lam_grid)

  # Canonical left-to-right order: decreasing lambda (sparse -> dense) within
  # each alpha level.  The path direction must not depend on the row order of
  # the supplied grid; `auto_lambda_grid()` already returns decreasing lambda,
  # so auto grids are unchanged.
  ord <- if (has_alpha) {
    order(lam_grid$alpha, -lam_grid$lambda)
  } else {
    order(-lam_grid$lambda)
  }
  scores   <- scores[, ord, drop = FALSE]
  lam_grid <- lam_grid[ord, , drop = FALSE]

  lambda_vals <- unname(lam_grid$lambda)
  df_real <- .scores_to_long_df(feat_names, scores, lambda_vals, lam_grid,
                                 has_alpha, selected = support)

  df_art <- NULL
  if (!is.null(object$stabl_scores_artificial_) &&
      !is.null(object$artificial_type)) {
    art_scores <- object$stabl_scores_artificial_[, ord, drop = FALSE]
    art_names  <- paste0("artificial.", seq_len(nrow(art_scores)))
    df_art <- .scores_to_long_df(art_names, art_scores, lambda_vals, lam_grid,
                                  has_alpha)
  }

  p <- ggplot2::ggplot()

  # Artificial feature lines (drawn first, below real features)
  if (!is.null(df_art)) {
    p <- p + ggplot2::geom_line(
      data    = df_art,
      mapping = ggplot2::aes(x = lambda_idx, y = score, group = feature),
      colour  = .stablr_colors$artificial,
      alpha   = 0.4,
      linewidth = 0.4,
      linetype = "dotted"
    )
  }

  # Noise features
  df_noise <- df_real[!df_real$selected, , drop = FALSE]
  if (nrow(df_noise) > 0L) {
    p <- p + ggplot2::geom_line(
      data    = df_noise,
      mapping = ggplot2::aes(x = lambda_idx, y = score, group = feature),
      colour  = .stablr_colors$noise,
      alpha   = 0.8,
      linewidth = 0.5
    )
  }

  # Selected features (on top)
  df_sel <- df_real[df_real$selected, , drop = FALSE]
  if (nrow(df_sel) > 0L) {
    p <- p + ggplot2::geom_line(
      data    = df_sel,
      mapping = ggplot2::aes(x = lambda_idx, y = score, group = feature),
      colour  = .stablr_colors$selected,
      alpha   = 1,
      linewidth = 1
    )
  }

  # Threshold line
  if (!is.null(threshold)) {
    thresh_label <- if (!is.null(object$hard_threshold) &&
                        is.null(new_hard_threshold)) {
      sprintf("Hard threshold = %.2f", threshold)
    } else if (!is.null(object$fdr_min_threshold_) &&
               is.null(new_hard_threshold)) {
      sprintf("FDP+ threshold = %.2f", threshold)
    } else {
      sprintf("Threshold = %.2f", threshold)
    }
    p <- p + ggplot2::geom_hline(
      yintercept = threshold,
      linetype   = "dashed",
      colour     = .stablr_colors$threshold,
      linewidth  = 0.7
    ) + ggplot2::annotate(
      "text",
      x     = Inf, y = threshold,
      label = thresh_label,
      hjust = 1.05, vjust = -0.4,
      size  = 3, colour = .stablr_colors$threshold
    )
  }

  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(add = 0.02)) +
    ggplot2::labs(
      title = title,
      x     = expression(lambda ~ "index"),
      y     = "Frequency of selection"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )

  # Facet by alpha if multiple alpha values present
  if (has_alpha && length(unique(lam_grid$alpha)) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ~ alpha,
      labeller = ggplot2::label_both,
      scales   = "free_x"
    )
  }

  p
}

# ── FDR diagnostic graph ──────────────────────────────────────────────────────

#' Plot the FDP+ FDR Estimate Curve
#'
#' Displays how the estimated False Discovery Proportion (FDP+) changes across
#' the full range of candidate stability thresholds, and marks the threshold
#' that achieves the minimum FDR estimate.
#'
#' This diagnostic shows why a particular stability threshold was chosen during
#' fitting: the current rule selects the first minimum of the displayed FDP+
#' curve. Curve shape is descriptive and does not by itself establish FDP
#' control, calibration, or signal/noise separation. Interpretation depends on
#' the artificial-feature assumptions and the fitted data.
#'
#' Requires that `object` was fitted with `artificial_type` set (not `NULL`).
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param title Character scalar; plot title.  Default `"FDR Estimate"`.
#' @param fdr_target Numeric scalar or `NULL`; a visual horizontal reference
#'   line only. It is not used for fitting or threshold selection. Use `NULL`
#'   to omit the line. Default `0.05`.
#'
#' @return A `ggplot` object.  The curve shows the FDP+ estimate at each
#'   candidate threshold; a vertical dashed line marks the optimal threshold
#'   stored in `object$fdr_min_threshold_`; a horizontal dashed line marks
#'   `fdr_target` when supplied.
#'
#' @seealso [stabl_fit()], [plot_stabl_path()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(60 * 8), 60, 8,
#'                dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#'   y <- setNames(rnorm(60), rownames(x))
#'   fit <- stabl_fit(x, y,
#'                    lambda_grid        = data.frame(lambda = c(0.3, 0.1, 0.05)),
#'                    n_bootstraps       = 6L,
#'                    hard_threshold     = NULL,
#'                    artificial_type    = "random_permutation",
#'                    random_state       = 1L)
#'   plot_fdr_graph(fit)
#' }
#' }
#' @export
plot_fdr_graph <- function(object, title = "FDR Estimate", fdr_target = 0.05) {
  .check_fitted_stabl(object)
  .require_ggplot2()
  title <- .validate_nonempty_scalar_character(title, "title")
  fdr_target <- .validate_scalar_numeric(
    fdr_target,
    "fdr_target",
    min = 0,
    allow_null = TRUE
  )
  if (is.null(object$FDRs_) || is.null(object$fdr_threshold_range)) {
    stop(
      "`object` was fitted without artificial features (artificial_type = NULL). ",
      "FDR graph requires FDRs_ and fdr_threshold_range to be present.",
      call. = FALSE
    )
  }

  thresh_grid <- object$fdr_threshold_range
  fdrs        <- object$FDRs_
  df          <- data.frame(threshold = thresh_grid, FDR = fdrs)

  # Optimal threshold
  if (!is.null(object$min_fdr_) && object$min_fdr_ > 1) {
    optimal_thresh <- 1.0
    opt_label      <- "No optimal threshold (min FDR > 1)"
  } else {
    optimal_thresh <- thresh_grid[which.min(fdrs)]  # ties: first-index wins (matches Python argmin)
    opt_label      <- sprintf("Optimal threshold = %.2f", optimal_thresh)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = FDR)) +
    ggplot2::geom_line(colour = .stablr_colors$noise, linewidth = 1) +
    ggplot2::geom_vline(
      xintercept = optimal_thresh,
      linetype   = "dashed",
      colour     = .stablr_colors$selected,
      linewidth  = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x     = optimal_thresh,
      y     = max(fdrs) * 0.95,
      label = opt_label,
      hjust = -0.05,
      size  = 3,
      colour = .stablr_colors$selected
    ) +
    ggplot2::labs(
      title = title,
      x     = "Stability threshold",
      y     = "FDR estimate"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank()
    )

  if (!is.null(fdr_target)) {
    p <- p + ggplot2::geom_hline(
      yintercept = fdr_target,
      linetype   = "dashed",
      colour     = .stablr_colors$threshold,
      linewidth  = 0.6
    )
  }

  p
}

# ── ROC curve ─────────────────────────────────────────────────────────────────

#' Plot a ROC Curve
#'
#' Generates a ggplot2 ROC (Receiver Operating Characteristic) curve from
#' binary outcome labels and predicted probabilities, together with the
#' chance-level diagonal for reference.
#'
#' ROC curves visualise the trade-off between sensitivity (TPR) and
#' specificity (1 - FPR) across all possible classification thresholds.  Use
#' this plot to assess overall discriminative ability after applying STABL
#' feature selection and fitting a downstream classifier on the selected
#' features.  Tied prediction scores are aggregated at a single threshold
#' before AUC computation, making the result invariant to row order.  The Area
#' Under the Curve (AUC) is shown in the caption.
#'
#' @param y_true Numeric 0/1 or logical vector of binary outcomes (1/`TRUE`
#'   for the positive class, 0/`FALSE` for the negative class).  Must contain
#'   both classes and be the same length as `y_preds`.  Character or factor
#'   labels should be recoded explicitly so the positive class is unambiguous.
#' @param y_preds Numeric vector of predicted probabilities for the positive
#'   class.  Values must be finite and in \eqn{[0, 1]}.
#' @param title Character scalar; plot title.  Default `"ROC Curve"`.
#'
#' @return A `ggplot` object.  The AUC is shown as a caption.
#'
#' @seealso [plot_prc()] for precision-recall curves (preferred when classes
#'   are severely imbalanced).
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   y_true  <- sample(0:1, 50, replace = TRUE)
#'   y_preds <- pmin(pmax(y_true * 0.6 + rnorm(50, 0, 0.3), 0), 1)
#'   plot_roc(y_true, y_preds)
#' }
#' @export
plot_roc <- function(y_true, y_preds, title = "ROC Curve") {
  .require_ggplot2()
  title <- .validate_nonempty_scalar_character(title, "title")
  inputs  <- .coerce_binary_inputs(y_true, y_preds)
  y_true  <- inputs$y_true
  y_preds <- inputs$y_preds

  roc_df  <- .roc_curve(y_true, y_preds)
  auc_val <- .trapz(roc_df$fpr, roc_df$tpr)

  ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(
      slope     = 1, intercept = 0,
      linetype  = "dashed",
      colour    = .stablr_colors$chance,
      linewidth = 0.8,
      alpha     = 0.7
    ) +
    ggplot2::geom_line(colour = .stablr_colors$selected, linewidth = 1.2) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(0)) +
    ggplot2::labs(
      title   = title,
      x       = "False positive rate",
      y       = "True positive rate",
      caption = sprintf("AUC = %.3f", auc_val)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# ── PRC curve ─────────────────────────────────────────────────────────────────

#' Plot a Precision-Recall Curve
#'
#' Generates a ggplot2 Precision-Recall Curve (PRC) from binary outcome labels
#' and predicted probabilities, optionally with iso-F1 contour lines.
#'
#' Precision-recall curves are preferred over ROC curves when the positive
#' class is rare (class imbalance), because they focus on the model's
#' performance on positive predictions without being diluted by the large
#' number of true negatives.  Use this plot after applying STABL feature
#' selection and fitting a downstream binary classifier.  Tied prediction
#' scores are aggregated at a single threshold before stepwise average
#' precision computation, making the result invariant to row order and avoiding
#' optimistic linear interpolation across precision drops.  Average precision
#' is shown in the caption.
#'
#' @param y_true Numeric 0/1 or logical vector of binary outcomes (1/`TRUE`
#'   for the positive class, 0/`FALSE` for the negative class).  Must contain
#'   both classes and be the same length as `y_preds`.  Character or factor
#'   labels should be recoded explicitly so the positive class is unambiguous.
#' @param y_preds Numeric vector of predicted probabilities for the positive
#'   class.  Values must be finite and in \eqn{[0, 1]}.
#' @param show_iso Logical; if `TRUE` (default), four iso-F1 contour lines
#'   at F1 = 0.2, 0.4, 0.6, 0.8 are drawn as reference guides.
#' @param title Character scalar; plot title.  Default
#'   `"Precision-Recall Curve"`.
#'
#' @return A `ggplot` object.  Average precision is shown as a caption.
#'
#' @seealso [plot_roc()] for the ROC curve alternative.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   y_true  <- sample(0:1, 50, replace = TRUE)
#'   y_preds <- pmin(pmax(y_true * 0.6 + rnorm(50, 0, 0.3), 0), 1)
#'   plot_prc(y_true, y_preds)
#' }
#' @export
plot_prc <- function(y_true, y_preds, show_iso = TRUE, title = "Precision-Recall Curve") {
  .require_ggplot2()
  if (!is.logical(show_iso) || length(show_iso) != 1L || is.na(show_iso)) {
    stop("`show_iso` must be a non-missing logical scalar.", call. = FALSE)
  }
  title <- .validate_nonempty_scalar_character(title, "title")
  inputs  <- .coerce_binary_inputs(y_true, y_preds)
  y_true  <- inputs$y_true
  y_preds <- inputs$y_preds

  prc_df  <- .prc_curve(y_true, y_preds)
  ap_val <- .average_precision(prc_df)

  p <- ggplot2::ggplot()

  # Iso-F1 lines
  if (show_iso) {
    recall_seq <- seq(0.01, 1, length.out = 100)
    for (f1_level in c(0.2, 0.4, 0.6, 0.8)) {
      iso_prec <- f1_level * recall_seq / (2 * recall_seq - f1_level)
      iso_prec[iso_prec < 0 | iso_prec > 1] <- NA_real_
      iso_df <- data.frame(recall = recall_seq, precision = iso_prec)
      p <- p + ggplot2::geom_line(
        data    = iso_df[!is.na(iso_df$precision), ],
        mapping = ggplot2::aes(x = recall, y = precision),
        colour  = .stablr_colors$ci,
        alpha   = 0.5,
        linetype = "dotted",
        linewidth = 0.5
      )
    }
  }

  p <- p +
    ggplot2::geom_line(
      data    = prc_df,
      mapping = ggplot2::aes(x = recall, y = precision),
      colour  = .stablr_colors$selected,
      linewidth = 1.2
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(0)) +
    ggplot2::labs(
      title   = title,
      x       = "Recall",
      y       = "Precision",
      caption = sprintf("Average precision = %.3f", ap_val)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p
}

# ── Feature distribution plots ────────────────────────────────────────────────

#' Boxplots of Selected Features Grouped by Outcome
#'
#' Produces a ggplot2 figure with one facet per selected feature, showing the
#' distribution of each feature's values stratified by the outcome class.
#' Intended for binary and multiclass classification tasks after STABL feature
#' selection; complements [scatterplot_features()] for regression tasks.
#'
#' @param features Character vector of feature names to plot.  Features not
#'   present as columns in `x` are silently skipped.
#' @param x Numeric matrix or data frame of predictors.  Column names must
#'   include all elements of `features`.
#' @param y Factor, character, or integer vector of class labels.  Must have
#'   the same number of elements as `nrow(x)`.
#' @param title Character scalar; plot title.  Default `"Selected Features"`.
#' @param ncol Positive integer; number of facet columns in the grid.
#'   Default 3.
#'
#' @return A `ggplot` object with one facet panel per feature.
#'
#' @seealso [scatterplot_features()] for regression tasks,
#'   [save_stabl_results()] which calls this automatically.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(40 * 4), 40, 4,
#'                dimnames = list(paste0("s", 1:40), c("A", "B", "C", "D")))
#'   y <- factor(rep(c("ctrl", "case"), 20))
#'   boxplot_features(c("A", "B"), x, y)
#' }
#' @export
boxplot_features <- function(features, x, y, title = "Selected Features", ncol = 3L) {
  .require_ggplot2()
  title <- .validate_nonempty_scalar_character(title, "title")
  ncol <- .validate_scalar_integer_like(ncol, "ncol", min = 1L)
  features <- .filter_features(features, x)

  df <- .features_long(features, x, y)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = outcome, y = value, fill = outcome)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.8, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
    ggplot2::facet_wrap(~ feature, scales = "free_y", ncol = ncol) +
    ggplot2::labs(title = title, x = NULL, y = "Value") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position  = "none",
      strip.text       = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Scatterplots of Selected Features Against a Continuous Outcome
#'
#' Produces a ggplot2 figure with one facet per selected feature, showing the
#' raw data points together with a LOESS smooth and 95% confidence band.  The
#' smooth helps reveal non-linear relationships that might be missed by
#' reporting correlation coefficients alone.  Use after STABL feature
#' selection for regression tasks.
#'
#' @param features Character vector of feature names to plot.  Features not
#'   present as columns in `x` are silently skipped.
#' @param x Numeric matrix or data frame of predictors.  Column names must
#'   include all elements of `features`.
#' @param y Numeric vector of continuous outcome values.  Must have the same
#'   number of elements as `nrow(x)`.
#' @param title Character scalar; plot title.  Default `"Selected Features"`.
#' @param ncol Positive integer; number of facet columns in the grid.
#'   Default 3.
#'
#' @return A `ggplot` object with one facet panel per feature.
#'
#' @seealso [boxplot_features()] for classification tasks,
#'   [save_stabl_results()] which calls this automatically.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(40 * 4), 40, 4,
#'                dimnames = list(paste0("s", 1:40), c("A", "B", "C", "D")))
#'   y <- x[, "A"] * 0.8 + rnorm(40, 0, 0.5)
#'   scatterplot_features(c("A", "B"), x, y)
#' }
#' @export
scatterplot_features <- function(features, x, y, title = "Selected Features", ncol = 3L) {
  .require_ggplot2()
  title <- .validate_nonempty_scalar_character(title, "title")
  ncol <- .validate_scalar_integer_like(ncol, "ncol", min = 1L)
  features <- .filter_features(features, x)

  y_num <- as.numeric(y)
  df_list <- lapply(features, function(f) {
    data.frame(
      feature = f,
      value   = as.numeric(x[, f]),
      outcome = y_num,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df_list)

  ggplot2::ggplot(df, ggplot2::aes(x = value, y = outcome)) +
    ggplot2::geom_point(
      colour = .stablr_colors$noise, alpha = 0.6, size = 1
    ) +
    ggplot2::geom_smooth(
      method  = "loess",
      formula = y ~ x,
      colour  = .stablr_colors$selected,
      se      = TRUE,
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(~ feature, scales = "free_x", ncol = ncol) +
    ggplot2::labs(title = title, x = "Feature value", y = "Outcome") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ── Internal computation helpers ──────────────────────────────────────────────

# Coerce and validate binary classifier inputs; returns list(y_true, y_preds).
.coerce_binary_inputs <- function(y_true, y_preds) {
  if (length(y_true) != length(y_preds)) {
    stop("`y_true` and `y_preds` must have the same length.", call. = FALSE)
  }

  if (length(y_true) == 0L) {
    stop("`y_true` and `y_preds` must be non-empty.", call. = FALSE)
  }

  if (!(is.logical(y_true) || is.numeric(y_true))) {
    stop(
      "`y_true` must be a binary 0/1 numeric or logical vector. ",
      "Recode factor or character labels explicitly so the positive class is unambiguous.",
      call. = FALSE
    )
  }

  if (anyNA(y_true) || any(!is.finite(as.numeric(y_true)))) {
    stop("`y_true` must not contain missing or non-finite values.", call. = FALSE)
  }

  y_true_numeric <- as.numeric(y_true)
  if (!all(y_true_numeric %in% c(0, 1))) {
    stop("`y_true` must be binary, containing only 0/1 or FALSE/TRUE values.", call. = FALSE)
  }

  y_true <- as.integer(y_true_numeric)
  if (length(unique(y_true)) < 2L) {
    stop("`y_true` must contain both classes for ROC/PRC plotting.", call. = FALSE)
  }

  if (!is.numeric(y_preds)) {
    stop("`y_preds` must be a numeric vector.", call. = FALSE)
  }

  y_preds <- as.numeric(y_preds)
  if (anyNA(y_preds) || any(!is.finite(y_preds))) {
    stop("`y_preds` must contain only finite values.", call. = FALSE)
  }

  if (any(y_preds < 0 | y_preds > 1)) {
    stop("`y_preds` values must be between 0 and 1.", call. = FALSE)
  }

  list(y_true = y_true, y_preds = y_preds)
}

# Intersect features with x columns; stop if none remain.
.filter_features <- function(features, x) {
  features <- intersect(features, colnames(x))
  if (length(features) == 0L) {
    stop("None of the requested `features` are columns of `x`.", call. = FALSE)
  }
  features
}

# Compute ROC curve: returns data.frame(fpr, tpr)
.roc_curve <- function(y_true, y_preds) {
  ord    <- order(y_preds, decreasing = TRUE)
  y_s    <- y_true[ord]
  pred_s <- y_preds[ord]
  n_pos  <- sum(y_true)
  n_neg  <- length(y_true) - n_pos
  threshold_ends <- cumsum(rle(pred_s)$lengths)
  tp_cum <- cumsum(y_s)[threshold_ends]
  fp_cum <- cumsum(1L - y_s)[threshold_ends]
  data.frame(
    fpr = c(0, fp_cum / n_neg),
    tpr = c(0, tp_cum / n_pos)
  )
}

# Compute PRC curve: returns data.frame(precision, recall)
.prc_curve <- function(y_true, y_preds) {
  ord    <- order(y_preds, decreasing = TRUE)
  y_s    <- y_true[ord]
  pred_s <- y_preds[ord]
  n_pos  <- sum(y_true)
  threshold_ends <- cumsum(rle(pred_s)$lengths)
  tp_cum <- cumsum(y_s)[threshold_ends]
  prec   <- tp_cum / threshold_ends
  rec    <- tp_cum / max(n_pos, 1L)
  data.frame(
    precision = c(1, prec),
    recall    = c(0, rec)
  )
}

# Trapezoidal numerical integration
.trapz <- function(x, y) {
  ord <- order(x)
  x   <- x[ord]
  y   <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}

# Stepwise area under the precision-recall curve, equivalent to average
# precision for grouped prediction thresholds.
.average_precision <- function(prc_df) {
  sum(diff(prc_df$recall) * prc_df$precision[-1])
}

# Build long data frame for boxplot_features
.features_long <- function(features, x, y) {
  out_levels <- if (is.factor(y)) levels(y) else sort(unique(as.character(y)))
  df_list <- lapply(features, function(f) {
    data.frame(
      feature = f,
      value   = as.numeric(x[, f]),
      outcome = factor(as.character(y), levels = out_levels),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, df_list)
}

.validate_nonempty_scalar_character <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || identical(x, "")) {
    stop(sprintf("`%s` must be a non-empty character scalar.", arg), call. = FALSE)
  }
  x
}
