#!/usr/bin/env Rscript

.find_package_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    desc <- file.path(path, "DESCRIPTION")
    if (file.exists(desc)) {
      fields <- read.dcf(desc)
      if (identical(unname(fields[1L, "Package"]), "stablr")) return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) return(NA_character_)
    path <- parent
  }
}

.load_stablr_source <- function() {
  root <- .find_package_root()
  if (!is.na(root) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
out <- NULL
if (length(args) >= 2L) {
  hit <- which(args == "--out")
  if (length(hit) == 1L && hit < length(args)) {
    out <- args[[hit + 1L]]
  }
}
if (is.null(out) || !nzchar(out)) {
  stop("Usage: stacked_multi_omic_benchmark.R --out /tmp/output-dir", call. = FALSE)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
.load_stablr_source()

.slow_stacked_multi_omic_scalar_reference <- function(predictions, y, task_type,
                                                       n_iter,
                                                       random_state = NULL) {
  predictions <- as.matrix(predictions)
  n_omics <- ncol(predictions)
  n_samples <- nrow(predictions)
  is_obs_preds <- !is.na(predictions)
  y_na_mask <- !is.na(y)
  score_fn <- if (identical(task_type, "binary")) stablr:::.r_auc else stablr:::.r_squared

  stablr:::.with_local_seed(
    if (!is.null(random_state)) as.integer(random_state),
    {
      best_score <- -Inf
      best_weights <- rep(1 / n_omics, n_omics)
      best_probs <- rep(NA_real_, n_samples)

      for (i in seq_len(n_iter)) {
        weights <- stats::runif(n_omics, 0, 10)
        weighted_probs <- stablr:::.weighted_masked_mean(
          predictions,
          weights,
          is_obs = is_obs_preds
        )
        complete_idx <- !is.na(weighted_probs) & y_na_mask
        if (sum(complete_idx) < 2L) next

        score <- tryCatch(
          score_fn(y[complete_idx], weighted_probs[complete_idx]),
          error = function(e) NA_real_
        )
        if (!is.na(score) && score > best_score) {
          best_score <- score
          best_weights <- weights
          best_probs <- weighted_probs
        }
      }

      list(weights = best_weights, predictions = best_probs, score = best_score)
    }
  )
}

set.seed(270627L)
n <- 500L
y <- rbinom(n, 1L, 0.5)
predictions <- data.frame(
  rna = y + rnorm(n, sd = 0.4),
  protein = y + rnorm(n, sd = 0.7),
  methylation = y + rnorm(n, sd = 1.0),
  metabolomics = y + rnorm(n, sd = 1.3),
  atac = y + rnorm(n, sd = 1.6),
  mirna = y + rnorm(n, sd = 1.9)
)
n_iter <- 20000L

old_elapsed <- system.time({
  old <- .slow_stacked_multi_omic_scalar_reference(
    predictions = predictions,
    y = y,
    task_type = "binary",
    n_iter = n_iter,
    random_state = 270627L
  )
})

new_elapsed <- system.time({
  new <- stablr::stacked_multi_omic(
    predictions = predictions,
    y = y,
    task_type = "binary",
    n_iter = n_iter,
    random_state = 270627L
  )
})

score_abs_error <- abs(new$score - old$score)
weight_max_abs_error <- max(abs(new$weights$Associated_weight - old$weights))
prediction_max_abs_error <- max(abs(
  new$predictions[["Stacked Gen. Predictions"]] - old$predictions
), na.rm = TRUE)

summary <- data.frame(
  n_samples = n,
  n_omics = ncol(predictions),
  n_iter = n_iter,
  score = new$score,
  scalar_elapsed_sec = unname(old_elapsed[["elapsed"]]),
  batched_elapsed_sec = unname(new_elapsed[["elapsed"]]),
  speedup = unname(old_elapsed[["elapsed"]] / new_elapsed[["elapsed"]]),
  score_abs_error = score_abs_error,
  weight_max_abs_error = weight_max_abs_error,
  prediction_max_abs_error = prediction_max_abs_error,
  parity_ok = score_abs_error <= 1e-12 &&
    weight_max_abs_error <= 1e-12 &&
    prediction_max_abs_error <= 1e-12
)
utils::write.csv(summary, file.path(out, "stacked_multi_omic_benchmark.csv"), row.names = FALSE)
