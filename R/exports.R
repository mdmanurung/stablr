# exports.R — CSV/disk export helpers for fitted STABL objects
# R port of export_stabl_to_csv() and save_stabl_results() from stabl/stabl.py

#' Export STABL Stability Scores to CSV
#'
#' Writes the full stability-score matrix, the per-feature maximum stability
#' scores, and (when artificial features were used) the corresponding
#' artificial-feature scores to CSV files in the directory specified by
#' `path`.
#'
#' Files written:
#' - `STABL scores.csv` — full p x L stability score matrix.
#' - `Max STABL scores.csv` — per-feature maximum stability score, sorted
#'   descending.
#' - `STABL artificial scores.csv` — artificial feature scores (only when
#'   `object$artificial_type` is not `NULL`).
#' - `Max STABL artificial scores.csv` — artificial max scores, sorted
#'   descending (only when artificial features were used).
#'
#' @param object A fitted `"stabl_fit"` object from [stabl_fit()].
#' @param path Character string; path to an existing (or to-be-created)
#'   directory.
#'
#' @return Invisibly returns `path` (as a normalized string).
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' x <- matrix(rnorm(40 * 6), 40, 6,
#'              dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
#' y <- setNames(rnorm(40), rownames(x))
#' fit <- stabl_fit(x, y,
#'                  lambda_grid  = data.frame(lambda = c(0.3, 0.1, 0.05)),
#'                  n_bootstraps = 6L, hard_threshold = 0.3, random_state = 1L)
#' out_dir <- file.path(tempdir(), "stabl_csv_export")
#' export_stabl_to_csv(fit, out_dir)
#' list.files(out_dir)
#' }
#' @export
export_stabl_to_csv <- function(object, path) {
  .check_fitted_stabl(object)
  path <- .validate_nonempty_scalar_character(path, "path")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  feat_names <- object$feature_names

  # Build column labels from the lambda grid rows
  col_labels <- .lambda_grid_row_labels(object$fitted_lambda_grid)

  # Full stability score matrix
  scores_df <- as.data.frame(object$stabl_scores_)
  rownames(scores_df) <- feat_names
  colnames(scores_df) <- col_labels
  utils::write.csv(scores_df, file = file.path(path, "STABL scores.csv"))

  # Max stability scores per feature (sorted descending)
  .write_max_scores_csv(rowMaxs(object$stabl_scores_), feat_names, path,
                        "Max STABL scores.csv")

  # Artificial feature scores (only when present)
  if (!is.null(object$stabl_scores_artificial_)) {
    art_scores <- object$stabl_scores_artificial_
    n_art      <- nrow(art_scores)
    art_names  <- paste0("artificial.", seq_len(n_art))

    art_df <- as.data.frame(art_scores)
    rownames(art_df) <- art_names
    colnames(art_df) <- col_labels
    utils::write.csv(art_df,
                     file = file.path(path, "STABL artificial scores.csv"))

    .write_max_scores_csv(rowMaxs(art_scores), art_names, path,
                          "Max STABL artificial scores.csv")
  }

  invisible(normalizePath(path))
}

#' Save All STABL Results to Disk
#'
#' Writes a complete set of STABL result artefacts — stability scores, the
#' list of selected features, the FDR diagnostic graph, the stability path
#' plot, and per-feature distribution plots — into a single output directory.
#' This is the recommended way to persist and share results from a single
#' [stabl_fit()] call.
#'
#' The function orchestrates several lower-level helpers:
#' [export_stabl_to_csv()] for tabular data, [plot_fdr_graph()] and
#' [plot_stabl_path()] for diagnostic plots, and [boxplot_features()] or
#' [scatterplot_features()] for per-feature distribution visualisations.
#' All graphics are saved using `ggplot2::ggsave()` in the format specified
#' by `figure_fmt`.
#'
#' @param object A fitted `"stabl_fit"` object returned by [stabl_fit()].
#' @param path Character string; path to the output directory.  Created
#'   recursively if it does not already exist.  Raises an error if the
#'   directory exists and `override = FALSE`.
#' @param x Numeric matrix of predictors used to fit `object`.  Required
#'   for per-feature distribution plots; must contain at least the selected
#'   feature columns.
#' @param y Outcome vector, factor, or `survival::Surv` object used to fit
#'   `object`.  Required for per-feature distribution plots.
#' @param figure_fmt Character; graphics device extension used when saving
#'   plots.  Common choices: `"pdf"` (default, publication quality),
#'   `"png"`, `"svg"`.
#' @param new_hard_threshold Numeric in `(0, 1]` or `NULL`.  When supplied,
#'   overrides the stored threshold for support extraction and stability-path
#'   plotting.
#' @param task_type Character; one of `"binary"` (default), `"multiclass"`,
#'   or `"regression"`.  Determines which per-feature plot is generated:
#'   grouped boxplots for classification, LOESS scatter for regression.
#' @param override Logical; if `TRUE`, an existing directory with the same
#'   path is accepted and its contents may be overwritten.  Default `FALSE`.
#'
#' @return Invisibly returns `path` as a normalised absolute string.
#'
#' @seealso [export_stabl_to_csv()] for the CSV-only variant,
#'   [plot_stabl_path()], [plot_fdr_graph()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(40 * 6), 40, 6,
#'                dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
#'   y <- setNames(rnorm(40), rownames(x))
#'   fit <- stabl_fit(x, y,
#'                    lambda_grid  = data.frame(lambda = c(0.3, 0.1, 0.05)),
#'                    n_bootstraps = 6L, hard_threshold = 0.3, random_state = 1L)
#'   out_dir <- file.path(tempdir(), "stabl_results")
#'   save_stabl_results(fit, out_dir, x = x, y = y,
#'                      task_type = "regression", figure_fmt = "png")
#'   list.files(out_dir)
#' }
#' }
#' @export
save_stabl_results <- function(
    object,
    path,
    x,
    y,
    figure_fmt          = "pdf",
    new_hard_threshold  = NULL,
    task_type           = "binary",
    override            = FALSE
) {
  .check_fitted_stabl(object)
  path <- .validate_nonempty_scalar_character(path, "path")
  figure_fmt <- .validate_nonempty_scalar_character(figure_fmt, "figure_fmt")
  if (!is.logical(override) || length(override) != 1L || is.na(override)) {
    stop("`override` must be a non-missing logical scalar.", call. = FALSE)
  }
  task_type <- match.arg(task_type, c("binary", "multiclass", "regression"))

  # Create output directory (fail if exists and override = FALSE)
  if (dir.exists(path) && !override) {
    stop(
      "Directory already exists: ", path,
      ". Set `override = TRUE` to overwrite.",
      call. = FALSE
    )
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  .require_ggplot2()

  # 1. CSV exports
  export_stabl_to_csv(object, path)

  # 2. FDR graph (only when artificial features were used)
  if (!is.null(object$artificial_type) && !is.null(object$FDRs_)) {
    p_fdr <- plot_fdr_graph(object)
    ggplot2::ggsave(
      filename = file.path(path, paste0("FDR Graph.", figure_fmt)),
      plot     = p_fdr,
      width    = 8, height = 4
    )
  }

  # 3. Stability path
  p_path <- plot_stabl_path(object, new_hard_threshold = new_hard_threshold)
  ggplot2::ggsave(
    filename = file.path(path, paste0("Stability Path.", figure_fmt)),
    plot     = p_path,
    width    = 6, height = 8
  )

  # 4. Selected features list → CSV
  sel_features <- get_feature_names(object,
                                        new_hard_threshold = new_hard_threshold)
  sel_dir <- file.path(path, "Selected Features")
  dir.create(sel_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(sel_features) > 0L) {
    sel_df <- data.frame(
      "Feature Name" = sel_features,
      row.names      = paste0("Feature n\u00b0", seq_along(sel_features)),
      check.names    = FALSE
    )
  } else {
    sel_df <- data.frame("Feature Name" = character(0L), check.names = FALSE)
  }
  utils::write.csv(sel_df, file = file.path(sel_dir, "Selected features.csv"))

  # 5. Per-feature distribution plots
  if (length(sel_features) > 0L) {
    p_feats <- if (task_type %in% c("binary", "multiclass")) {
      boxplot_features(features = sel_features, x = x, y = y)
    } else if (task_type == "regression") {
      scatterplot_features(features = sel_features, x = x, y = y)
    }
    if (!is.null(p_feats)) {
      ggplot2::ggsave(
        filename = file.path(sel_dir,
                             paste0("Feature distributions.", figure_fmt)),
        plot   = p_feats,
        width  = max(4, 3 * min(length(sel_features), 4)),
        height = ceiling(length(sel_features) / 4) * 3
      )
    }
  }

  invisible(normalizePath(path))
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Build, sort-descending, and write a max-scores CSV.
.write_max_scores_csv <- function(scores, feat_names, path, filename) {
  max_df <- data.frame("Max Proba" = scores, row.names = feat_names,
                       check.names = FALSE)
  max_df <- max_df[order(max_df[["Max Proba"]], decreasing = TRUE), , drop = FALSE]
  utils::write.csv(max_df, file = file.path(path, filename))
}

# Row-wise max wrapper.
# Special-cases nrow=0 and ncol=1 to preserve historical edge behavior.
rowMaxs <- function(m) {
  if (nrow(m) == 0L) return(numeric(0L))
  if (ncol(m) == 1L) return(as.numeric(m[, 1L]))
  matrixStats::rowMaxs(as.matrix(m))
}

# Build readable column labels from the lambda grid (one label per row)
.lambda_grid_row_labels <- function(lambda_grid) {
  if (nrow(lambda_grid) == 0L) return(character(0L))
  apply(lambda_grid, 1L, function(row) {
    pairs <- paste0(names(row), "=", signif(as.numeric(row), 4L))
    paste(pairs, collapse = ", ")
  })
}
