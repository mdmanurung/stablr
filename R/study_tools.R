# Study-level reliability tools (0.3.0 cross-study round).
#
# B2 (inferred scope): per-study batch diagnostics — outcome composition and
# feature location shifts per study, used to decide whether study-disjoint
# CV (see stabl_multiomic_nested_cv) is warranted and to describe the shift
# magnitude in benchmark outputs.
#
# B3 (inferred scope): a study-disjoint train/valid split helper feeding
# stabl_multiomic_train_validate(). stabl_multiomic_cv() already supports
# group-disjoint assessment folds through its `groups` argument; this helper
# completes the train_validate side of the workflow.

#' Per-Study Batch Diagnostics for Multi-Omic STABL Inputs
#'
#' Summarises each study's outcome composition and the location shift of its
#' features relative to the pooled rest of the data. The diagnostics quantify
#' how much study-level structure a sample-i.i.d. cross-validation scheme
#' would leak across folds, and feed the study-disjoint benchmark
#' (`inst/benchmarks/b1_study_disjoint_cv.R`).
#'
#' Feature shifts are reported as absolute standardized mean differences
#' (SMD) per feature: `(mean_study - mean_rest) / pooled_sd` with
#' `pooled_sd = sqrt((sd_study^2 + sd_rest^2) / 2)`. Features with zero
#' pooled SD are excluded from the summaries and counted in `n_zero_sd`.
#'
#' @param x_list Named list of omic matrices or data frames with identical
#'   sample row names.
#' @param y Named outcome vector aligned to `x_list` rows. Factors or
#'   characters are summarised as per-level counts and prevalences; numeric
#'   outcomes as per-study mean and SD.
#' @param study_labels Named vector assigning each sample to a study
#'   (cohort, site, or batch); names must cover exactly the sample IDs. Each
#'   study must contain at least 2 samples.
#'
#' @return A list with two data frames: `study_summary` (one row per study
#'   with sample counts and outcome composition) and `feature_shift` (one row
#'   per study and omic block with `n_features`, `smd_mean_abs`,
#'   `smd_p95_abs`, `smd_max_abs`, and `n_zero_sd`).
#'
#' @seealso [stabl_multiomic_nested_cv()] for study-disjoint
#'   cross-validation and [stabl_multiomic_study_split()] for a
#'   study-disjoint train/validation split.
#'
#' @examples
#' set.seed(5L)
#' ids <- paste0("s", 1:12)
#' y <- setNames(factor(rep(c("A", "B"), times = 6L)), ids)
#' x <- matrix(rnorm(12 * 4), 12, 4, dimnames = list(ids, paste0("g", 1:4)))
#' study_labels <- setNames(rep(paste0("study", 1:3), each = 4L), ids)
#' # study1 carries a location shift on the first two features
#' x[study_labels == "study1", 1:2] <- x[study_labels == "study1", 1:2] + 2
#' diag <- stabl_multiomic_study_diagnostics(
#'   x_list = list(rna = x), y = y, study_labels = study_labels
#' )
#' diag$study_summary
#' diag$feature_shift
#' @export
stabl_multiomic_study_diagnostics <- function(x_list, y, study_labels) {
  validate_multiomic_inputs(x_list = x_list, y = y)
  sample_ids <- rownames(x_list[[1L]])
  y <- .subset_outcome_by_ids(y, sample_ids)
  study_labels <- .validate_study_labels(study_labels, sample_ids)

  studies <- sort(unique(as.character(study_labels)))
  study_of_sample <- as.character(study_labels)

  study_summary <- .study_outcome_summary(y, study_of_sample, studies)
  feature_shift <- .study_feature_shift(x_list, study_of_sample, studies)

  list(study_summary = study_summary, feature_shift = feature_shift)
}

.study_outcome_summary <- function(y, study_of_sample, studies) {
  if (is.numeric(y)) {
    rows <- lapply(studies, function(s) {
      values <- y[study_of_sample == s]
      data.frame(
        study = s,
        n_samples = length(values),
        outcome_mean = mean(values),
        outcome_sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    return(do.call(rbind, rows))
  }

  y <- factor(y)
  levels_y <- levels(y)
  rows <- lapply(studies, function(s) {
    values <- y[study_of_sample == s]
    counts <- as.integer(table(factor(as.character(values), levels = levels_y)))
    row <- data.frame(
      study = s,
      n_samples = length(values),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(levels_y)) {
      row[[paste0("n_", levels_y[[j]])]] <- counts[[j]]
      row[[paste0("prop_", levels_y[[j]])]] <- counts[[j]] / length(values)
    }
    row
  })
  do.call(rbind, rows)
}

.study_feature_shift <- function(x_list, study_of_sample, studies) {
  rows <- list()
  k <- 0L
  for (block in names(x_list)) {
    x <- x_list[[block]]
    x <- as.matrix(x)
    if (!is.numeric(x)) {
      stop(sprintf(
        "Block `%s` must be numeric after conversion to a matrix.", block
      ), call. = FALSE)
    }
    for (s in studies) {
      idx_study <- which(study_of_sample == s)
      idx_rest <- which(study_of_sample != s)
      smd_abs <- .study_block_smd_abs(x[idx_study, , drop = FALSE],
                                      x[idx_rest, , drop = FALSE])
      zero_sd <- length(smd_abs) - length(stats::na.omit(smd_abs))
      smd_abs <- stats::na.omit(smd_abs)
      k <- k + 1L
      rows[[k]] <- data.frame(
        study = s,
        block = block,
        n_features = length(smd_abs),
        smd_mean_abs = if (length(smd_abs)) mean(smd_abs) else NA_real_,
        smd_p95_abs = if (length(smd_abs)) stats::quantile(smd_abs, 0.95, names = FALSE) else NA_real_,
        smd_max_abs = if (length(smd_abs)) max(smd_abs) else NA_real_,
        n_zero_sd = zero_sd,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.study_block_smd_abs <- function(x_study, x_rest) {
  mu_study <- colMeans(x_study)
  mu_rest <- colMeans(x_rest)
  sd_study <- apply(x_study, 2L, stats::sd)
  sd_rest <- apply(x_rest, 2L, stats::sd)
  pooled_sd <- sqrt((sd_study^2 + sd_rest^2) / 2)
  smd <- (mu_study - mu_rest) / pooled_sd
  smd[pooled_sd == 0] <- NA_real_
  abs(smd)
}

#' Study-Disjoint Train/Validation Split
#'
#' Assigns whole studies to a training and a validation split so that
#' [stabl_multiomic_train_validate()] evaluates a cross-study deployment
#' setting rather than a sample-i.i.d. one. Studies are assigned largest-first
#' by a deterministic greedy pass targeting `prop_train` of the samples while
#' keeping outcome class proportions as balanced as possible across the two
#' sides; a seeded shuffle fixes the order of size ties.
#'
#' `prop_train` is a target, not a guarantee: study granularity governs the
#' achievable split, and extreme values can fail the fail-closed checks below.
#'
#' @param x_list Named list of omic matrices or data frames with identical
#'   sample row names.
#' @param y Named outcome vector aligned to `x_list` rows.
#' @param study_labels Named vector assigning each sample to a study; names
#'   must cover exactly the sample IDs. Each study must contain at least 2
#'   samples, and at least 2 studies are required.
#' @param prop_train Target proportion of samples in the training split
#'   (strictly between 0 and 1).
#' @param stratified Logical. When `TRUE` (default), class-proportion
#'   imbalance is the secondary assignment criterion after the size target.
#' @param random_state Optional integer seed for deterministic splits.
#'
#' @return A list with `train_ids`, `valid_ids` (sample IDs in the original
#'   row order), `study_assignment` (named vector mapping each study to
#'   `"train"` or `"valid"`), `class_counts` (per-side held-out outcome
#'   counts), and `fold_scheme = "study_disjoint"`.
#'
#' @examples
#' set.seed(6L)
#' ids <- paste0("s", 1:12)
#' y <- setNames(factor(rep(c("A", "B"), times = 6L)), ids)
#' x <- matrix(rnorm(12 * 4), 12, 4, dimnames = list(ids, paste0("g", 1:4)))
#' study_labels <- setNames(rep(paste0("study", 1:3), each = 4L), ids)
#' split <- stabl_multiomic_study_split(
#'   x_list = list(rna = x), y = y, study_labels = study_labels,
#'   prop_train = 2 / 3, random_state = 6L
#' )
#' split$study_assignment
#' @export
stabl_multiomic_study_split <- function(x_list, y, study_labels,
                                        prop_train = 0.75,
                                        stratified = TRUE,
                                        random_state = NULL) {
  if (!is.numeric(prop_train) || length(prop_train) != 1L ||
      is.na(prop_train) || prop_train <= 0 || prop_train >= 1) {
    stop("`prop_train` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }
  validate_multiomic_inputs(x_list = x_list, y = y)
  sample_ids <- rownames(x_list[[1L]])
  y <- .subset_outcome_by_ids(y, sample_ids)
  study_labels <- .validate_study_labels(study_labels, sample_ids)

  studies <- sort(unique(as.character(study_labels)))
  if (length(studies) < 2L) {
    stop("Study-disjoint splitting requires at least 2 studies.", call. = FALSE)
  }

  study_of_sample <- as.character(study_labels)
  class_tab <- table(
    study = factor(study_of_sample, levels = studies),
    class = factor(as.character(y))
  )

  assignment <- .greedy_split_assignment(
    class_tab = class_tab,
    prop_train = prop_train,
    stratified = stratified,
    random_state = random_state
  )

  train_studies <- names(assignment)[assignment == "train"]
  valid_studies <- names(assignment)[assignment == "valid"]
  if (length(train_studies) == 0L || length(valid_studies) == 0L) {
    stop(sprintf(
      paste0(
        "No study-disjoint split realises prop_train = %s for these study ",
        "sizes; every study would land on one side. Move `prop_train` ",
        "closer to 0.5 or supply more studies."
      ),
      format(prop_train)
    ), call. = FALSE)
  }

  train_ids <- names(y)[study_of_sample %in% train_studies]
  valid_ids <- names(y)[study_of_sample %in% valid_studies]

  if (!is.numeric(y)) {
    for (side in list(list("training", train_ids), list("validation", valid_ids))) {
      side_counts <- table(factor(as.character(y[side[[2L]]])))
      if (sum(side_counts > 0L) < 2L) {
        stop(sprintf(
          paste0(
            "The %s split contains a single outcome class (counts: %s; ",
            "studies: %s); adjust `prop_train` or supply more studies."
          ),
          side[[1L]],
          paste(sprintf("%s=%d", names(side_counts), as.integer(side_counts)),
                collapse = ", "),
          paste(if (side[[1L]] == "training") train_studies else valid_studies,
                collapse = ", ")
        ), call. = FALSE)
      }
    }
  }

  list(
    train_ids = train_ids,
    valid_ids = valid_ids,
    study_assignment = assignment,
    class_counts = list(
      train = table(y[train_ids]),
      valid = table(y[valid_ids])
    ),
    fold_scheme = "study_disjoint"
  )
}

.greedy_split_assignment <- function(class_tab, prop_train, stratified,
                                     random_state) {
  studies <- rownames(class_tab)
  sizes <- as.integer(rowSums(class_tab))
  n_total <- sum(sizes)
  target_train <- prop_train * n_total
  target_valid <- (1 - prop_train) * n_total
  overall_prop <- colSums(class_tab) / n_total

  ordered_studies <- .with_local_seed(
    if (!is.null(random_state)) as.integer(random_state),
    {
      shuffled <- sample(studies)
      shuffled_sizes <- sizes[shuffled]
      shuffled[order(-shuffled_sizes, method = "radix")]
    }
  )

  train_counts <- rep(0L, ncol(class_tab))
  valid_counts <- rep(0L, ncol(class_tab))
  train_size <- 0L
  valid_size <- 0L
  train_n_studies <- 0L
  valid_n_studies <- 0L
  assignment <- setNames(character(length(ordered_studies)), ordered_studies)

  for (study in ordered_studies) {
    study_counts <- as.integer(class_tab[study, ])
    study_size <- sum(study_counts)
    d_train <- abs(train_size + study_size - target_train)
    d_valid <- abs(valid_size + study_size - target_valid)

    if (d_train < d_valid) {
      side <- "train"
    } else if (d_valid < d_train) {
      side <- "valid"
    } else if (!isTRUE(stratified)) {
      side <- if (train_n_studies <= valid_n_studies) "train" else "valid"
    } else {
      imb_train <- sum(abs((train_counts + study_counts) /
                             (train_size + study_size) - overall_prop))
      imb_valid <- sum(abs((valid_counts + study_counts) /
                             (valid_size + study_size) - overall_prop))
      if (imb_train < imb_valid) {
        side <- "train"
      } else if (imb_valid < imb_train) {
        side <- "valid"
      } else {
        side <- if (train_n_studies <= valid_n_studies) "train" else "valid"
      }
    }

    if (side == "train") {
      train_counts <- train_counts + study_counts
      train_size <- train_size + study_size
      train_n_studies <- train_n_studies + 1L
    } else {
      valid_counts <- valid_counts + study_counts
      valid_size <- valid_size + study_size
      valid_n_studies <- valid_n_studies + 1L
    }
    assignment[study] <- side
  }
  assignment
}
