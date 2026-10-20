# Cohort replication and heterogeneity profiles (1.0.0 spec, A3).
#
# stabl_replication_profile() characterizes whether selection replicates
# across cohorts (studies, sites, batches). It distinguishes four
# estimands explicitly -- they answer different questions and agreement
# on selection alone says nothing causal:
#
#   (i)   pooled selection stability   -- the pooled stabl_fit() result;
#   (ii)  training-cohort replication  -- fraction of leave-one-study-out
#         (LOSO) training fits selecting the feature, conditional on the
#         feature being eligible (measured and non-constant) in that
#         training partition;
#   (iii) direction/effect replication -- sign concordance and an
#         inverse-variance pooled estimate with Cochran's Q
#         heterogeneity, computed ONLY from per-study effect estimates
#         the user supplies (documented estimation design is the user's
#         responsibility; the package never invents effect estimates);
#   (iv)  held-out cohort performance  -- per-study predictive
#         contribution via a documented downstream-refit design: for
#         each LOSO fold, a stats::glm on the training partition's
#         selected features is evaluated on the held-out study.
#         stabl_fit() stores no reusable model, so this is explicitly a
#         downstream-refit assessment, not a property of the STABL fit
#         itself.
#
# Implementation rules honored (spec A3): every LOSO/per-study fit is a
# fresh stabl_fit() call on the relevant training partition, so
# preprocessing, filtering, and tuning are recomputed inside training
# data by construction; no held-out cohort ever chooses a threshold.
# Uncertainty is reported as cohort-level observations (counts and
# proportions over studies); folds sharing studies are dependent, so no
# binomial CIs over folds are computed.

#' Cohort Replication and Heterogeneity Profile
#'
#' Characterizes cross-cohort replication of feature selection: pooled
#' selection stability, training-cohort replication across
#' leave-one-study-out fits, leave-one-study-out sensitivity of the
#' pooled panel, optional direction/effect replication from
#' user-supplied per-study effect estimates, and optional held-out
#' cohort performance via a documented downstream-refit design.
#'
#' @param x Numeric predictor matrix with non-empty row names (sample
#'   ids) and column names (feature ids), as for [stabl_fit()].
#' @param y Outcome vector, named to match `rownames(x)`.
#' @param study_labels Named vector (factor-coercible) assigning each
#'   sample to a study/cohort; names must match `rownames(x)` exactly.
#'   At least three studies are required.
#' @param fit_spec Named list of further arguments for [stabl_fit()],
#'   e.g. `list(n_bootstraps = 50, artificial_type = "knockoff_equi")`.
#'   The per-fit `random_state` is always set by
#'   `stabl_replication_profile`; `lambda_grid` defaults to `"auto"`
#'   when not supplied.
#' @param pooled_fit Optional pre-computed pooled [stabl_fit()] result
#'   on `(x, y)`; validated for feature-name agreement and refitted if
#'   absent.
#' @param effects Optional data frame with columns `feature`, `study`,
#'   `estimate`, `se` giving per-study effect estimates for direction
#'   replication. Estimates must come from a documented per-study
#'   estimation design; they are used as supplied. Duplicated
#'   (feature, study) rows are an error.
#' @param heldout_performance Logical; if `TRUE`, evaluate each LOSO
#'   fold's selected panel on the held-out study via a downstream
#'   [stats::glm()] refit (gaussian: MSE; binomial: balanced error
#'   rate). Empty panels contribute `NA` with reason
#'   `"empty_support"`.
#' @param seed Optional base seed for the per-fit derived seeds
#'   (recorded; fully replayable). If `NULL`, drawn from the current
#'   RNG state and recorded.
#' @param ... Reserved; unused arguments are an error.
#'
#' @return A `stabl_replication_profile` object: a list with
#'   \describe{
#'     \item{features}{Per-feature data frame: `n_measured_cohorts`,
#'       `n_eligible_cohorts`, `within_study_frequency`,
#'       `training_replication`, `fold_frequency`,
#'       `in_pooled_support`, `loso_retention`,
#'       `dropped_on_removal_of`, and (when `effects` is supplied)
#'       `n_estimates`, `sign_concordance`, `pooled_estimate`,
#'       `pooled_se`, `pooled_ci_lo`, `pooled_ci_hi`, `het_Q`,
#'       `het_df`, `het_p`.}
#'     \item{estimands}{Named list describing the four estimands and
#'       where each lives in the output.}
#'     \item{pooled_support}{Character vector of the pooled selected
#'       panel.}
#'     \item{heldout}{When `heldout_performance = TRUE`, a per-study
#'       data frame of held-out metrics with `reason` for `NA`s.}
#'     \item{eligibility}{Per feature-by-cohort eligibility with
#'       `reason` for ineligibility (`"all_missing"`,
#'       `"constant_in_cohort"`).}
#'     \item{study_labels, fit_spec, seed, fit_seeds, call}{Provenance.}
#'   }
#'
#' @section Interpretation caveat:
#'   The four estimands answer different questions. Selection agreement
#'   across cohorts is not evidence of causal or directional agreement;
#'   direction/effect replication is only as valid as the supplied
#'   per-study estimation design. Proportions are cohort-level
#'   observations over dependent folds -- no binomial confidence
#'   intervals are computed.
#'
#' @examples
#' set.seed(1)
#' ids <- paste0("s", 1:60)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(ids, paste0("f", 1:8)))
#' y <- setNames(x[, 1] * 3 + rnorm(60, sd = 0.8), ids)
#' studies <- setNames(rep(paste0("study", 1:3), each = 20L), ids)
#' rp <- stabl_replication_profile(
#'   x, y, studies,
#'   fit_spec = list(n_bootstraps = 50,
#'                   artificial_type = "random_permutation"),
#'   seed = 5
#' )
#' head(rp$features)
#' @export
stabl_replication_profile <- function(x, y, study_labels,
                                      fit_spec = list(),
                                      pooled_fit = NULL,
                                      effects = NULL,
                                      heldout_performance = FALSE,
                                      seed = NULL,
                                      ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    stop("Unused arguments: ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }
  if (!is.matrix(x) || is.null(rownames(x)) || is.null(colnames(x))) {
    stop("`x` must be a matrix with sample-id row names and feature ",
         "column names.", call. = FALSE)
  }
  if (length(y) != nrow(x) || is.null(names(y)) ||
      !identical(names(y), rownames(x))) {
    stop("`y` must be named to match `rownames(x)`.", call. = FALSE)
  }
  if (length(study_labels) != nrow(x)) {
    stop("`study_labels` must have length nrow(x).", call. = FALSE)
  }
  if (is.null(names(study_labels)) ||
      !identical(names(study_labels), rownames(x))) {
    stop("`study_labels` must be named to match `rownames(x)`.",
         call. = FALSE)
  }
  study_labels <- as.factor(study_labels)
  studies <- levels(study_labels)
  if (length(studies) < 3L) {
    stop("Need at least three studies for a replication profile.",
         call. = FALSE)
  }
  if (!is.list(fit_spec)) {
    stop("`fit_spec` must be a named list of stabl_fit() arguments.",
         call. = FALSE)
  }
  if (length(fit_spec) > 0L &&
      (is.null(names(fit_spec)) || any(!nzchar(names(fit_spec))))) {
    stop("`fit_spec` must be a fully named list.", call. = FALSE)
  }
  forbidden <- intersect(names(fit_spec), c("x", "y", "random_state"))
  if (length(forbidden) > 0L) {
    stop("`fit_spec` must not set ", paste(forbidden, collapse = ", "),
         " -- data and per-fit seeds are managed by ",
         "stabl_replication_profile().", call. = FALSE)
  }
  if (!"lambda_grid" %in% names(fit_spec)) {
    fit_spec$lambda_grid <- "auto"
  }
  if (!is.null(pooled_fit)) {
    .check_fitted_stabl(pooled_fit)
    if (!identical(pooled_fit$feature_names, colnames(x))) {
      stop("`pooled_fit` feature names do not match `colnames(x)`.",
           call. = FALSE)
    }
  }
  if (!is.null(effects)) {
    effects <- .validate_effects(effects, studies, colnames(x))
  }
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  seed <- .validate_scalar_integer_like(
    seed, "seed", min = 1L, max = .Machine$integer.max - 1L
  )

  # Per-fit seeds: one stream for within-study fits, one for LOSO fits,
  # one for the pooled fit (only consumed when pooled_fit is NULL).
  n_fits <- 2L * length(studies) + 1L
  fit_seeds <- .derive_stream_seeds(seed, n_fits)
  names(fit_seeds) <- c(paste0("within_", studies),
                        paste0("loso_", studies), "pooled")

  # ---- eligibility per feature x cohort ----------------------------------
  eligibility <- .replication_eligibility(x, study_labels, studies)

  # ---- fits ---------------------------------------------------------------
  # stabl_fit() is hermetic given random_state (PR-2 contract), so the
  # derived per-fit seed is sufficient; no ambient RNG is consumed.
  run_fit <- function(idx, sub_x, sub_y) {
    do.call(
      stabl_fit,
      c(list(x = sub_x, y = sub_y, random_state = fit_seeds[[idx]]),
        fit_spec)
    )
  }

  # (i) pooled
  if (is.null(pooled_fit)) {
    pooled_fit <- run_fit("pooled", x, y)
  }
  pooled_support <- names(which(get_support(pooled_fit)))
  family <- pooled_fit$family %||% "gaussian"
  y_levels <- if (identical(family, "binomial")) {
    levels(as.factor(y))
  } else {
    NULL
  }

  # within-study fits (estimand ii input: within-study selection
  # frequency over eligible cohorts)
  within_supports <- list()
  for (s in studies) {
    idx <- study_labels == s
    within_supports[[s]] <- names(which(get_support(run_fit(
      paste0("within_", s), x[idx, , drop = FALSE], y[idx]
    ))))
  }

  # LOSO training fits (estimands ii and iv)
  loso_supports <- list()
  heldout_rows <- list()
  for (s in studies) {
    train <- study_labels != s
    fit_s <- run_fit(
      paste0("loso_", s), x[train, , drop = FALSE], y[train]
    )
    loso_supports[[s]] <- names(which(get_support(fit_s)))
    if (isTRUE(heldout_performance)) {
      heldout_rows[[s]] <- .replication_heldout_row(
        s, fit_s, x, y, train, family, y_levels
      )
    }
  }

  features <- .assemble_replication_features(
    feature_names = colnames(x),
    studies = studies,
    eligibility = eligibility,
    within_supports = within_supports,
    loso_supports = loso_supports,
    pooled_support = pooled_support,
    effects = effects
  )

  structure(
    list(
      features = features,
      estimands = list(
        pooled_selection_stability =
          "pooled_support / in_pooled_support column",
        training_cohort_replication =
          paste0("training_replication (LOSO training fits, ",
                 "eligibility-conditioned)"),
        direction_effect_replication =
          paste0("sign_concordance / pooled_estimate columns ",
                 "(effects supplied only)"),
        heldout_cohort_performance =
          "heldout data frame (downstream glm refit; NA with reason)"
      ),
      pooled_support = pooled_support,
      heldout = if (length(heldout_rows) > 0L) {
        do.call(rbind, heldout_rows)
      } else {
        NULL
      },
      eligibility = eligibility,
      study_labels = study_labels,
      fit_spec = fit_spec,
      seed = seed,
      fit_seeds = fit_seeds,
      call = match.call()
    ),
    class = "stabl_replication_profile"
  )
}

# ---- helpers --------------------------------------------------------------

.validate_effects <- function(effects, studies, feature_names) {
  need <- c("feature", "study", "estimate", "se")
  if (!is.data.frame(effects) || !all(need %in% names(effects))) {
    stop("`effects` must be a data frame with columns: ",
         paste(need, collapse = ", "), call. = FALSE)
  }
  effects$feature <- as.character(effects$feature)
  effects$study <- as.character(effects$study)
  if (any(!is.finite(effects$estimate)) || any(!is.finite(effects$se)) ||
      any(effects$se <= 0)) {
    stop("`effects` needs finite estimates and strictly positive SEs.",
         call. = FALSE)
  }
  bad_feature <- setdiff(unique(effects$feature), feature_names)
  bad_study <- setdiff(unique(effects$study), studies)
  if (length(bad_feature) > 0L || length(bad_study) > 0L) {
    stop("`effects` references unknown features or studies: ",
         paste(c(bad_feature, bad_study), collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(paste(effects$feature, effects$study)) > 0L) {
    stop("`effects` must not contain duplicated (feature, study) rows.",
         call. = FALSE)
  }
  effects
}

.replication_eligibility <- function(x, study_labels, studies) {
  rows <- list()
  for (f in colnames(x)) {
    for (s in studies) {
      v <- x[study_labels == s, f]
      reason <- if (all(is.na(v))) {
        "all_missing"
      } else if (length(unique(v[!is.na(v)])) < 2L) {
        "constant_in_cohort"
      } else {
        NA_character_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        feature = f, study = s, eligible = is.na(reason),
        reason = reason, stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.replication_heldout_row <- function(s, fit_s, x, y, train, family,
                                     y_levels) {
  support <- names(which(get_support(fit_s)))
  if (length(support) == 0L) {
    return(data.frame(
      study = s, metric = NA_real_, metric_name = NA_character_,
      n_panel = 0L, reason = "empty_support", stringsAsFactors = FALSE
    ))
  }
  # Documented downstream-refit design: stabl_fit stores no reusable
  # model, so a stats::glm is fit on the training partition's selected
  # features and evaluated on the held-out study. This is a property of
  # the refit downstream model, not of the STABL selection itself. The
  # held-out study never influences the panel or any threshold.
  df_train <- data.frame(y = y[train], x[train, support, drop = FALSE])
  df_test <- data.frame(y = y[!train], x[!train, support, drop = FALSE])
  if (identical(family, "binomial")) {
    # Pin both partitions to the full-data level order so the response
    # coding cannot flip when a level is absent from one partition.
    df_train$y <- factor(df_train$y, levels = y_levels)
    fam <- stats::binomial()
  } else {
    fam <- stats::gaussian()
  }
  refit <- tryCatch(
    stats::glm(y ~ ., data = df_train, family = fam),
    error = function(e) e
  )
  if (inherits(refit, "error")) {
    return(data.frame(
      study = s, metric = NA_real_, metric_name = NA_character_,
      n_panel = length(support), reason = "refit_failed",
      stringsAsFactors = FALSE
    ))
  }
  # Rows with NA in the outcome or any selected predictor cannot be
  # scored: predict() yields NA (or drops them, depending on na.action),
  # which would silently misalign or NA-out the held-out metric
  # (review E4).  Restrict truth and predictions to complete rows; if no
  # complete rows remain, record the reason rather than a bogus metric.
  complete <- stats::complete.cases(df_test)
  if (!any(complete)) {
    return(data.frame(
      study = s, metric = NA_real_, metric_name = NA_character_,
      n_panel = length(support), reason = "incomplete_predictors",
      stringsAsFactors = FALSE
    ))
  }
  df_test <- df_test[complete, , drop = FALSE]
  pred <- stats::predict(refit, newdata = df_test, type = "response")
  if (identical(family, "binomial")) {
    cls <- as.integer(pred >= 0.5)
    truth <- as.integer(factor(df_test$y, levels = y_levels)) - 1L
    tpr <- if (any(truth == 1L)) {
      mean(cls[truth == 1L] == 1L)
    } else {
      NA_real_
    }
    tnr <- if (any(truth == 0L)) {
      mean(cls[truth == 0L] == 0L)
    } else {
      NA_real_
    }
    metric <- 1 - mean(c(tpr, tnr), na.rm = TRUE)
    metric_name <- "balanced_error_rate"
  } else {
    metric <- mean((df_test$y - pred)^2)
    metric_name <- "mse"
  }
  data.frame(
    study = s, metric = metric, metric_name = metric_name,
    n_panel = length(support), reason = NA_character_,
    stringsAsFactors = FALSE
  )
}

.assemble_replication_features <- function(feature_names, studies,
                                           eligibility, within_supports,
                                           loso_supports, pooled_support,
                                           effects) {
  rows <- list()
  for (f in feature_names) {
    elig <- eligibility[eligibility$feature == f, , drop = FALSE]
    n_measured <- sum(elig$reason != "all_missing" | is.na(elig$reason))
    n_eligible <- sum(elig$eligible)

    # within-study selection frequency over eligible cohorts
    elig_studies <- elig$study[elig$eligible]
    within_hits <- sum(vapply(
      elig_studies, function(s) f %in% within_supports[[s]], logical(1L)
    ))
    within_freq <- if (length(elig_studies) > 0L) {
      within_hits / length(elig_studies)
    } else {
      NA_real_
    }

    # LOSO: feature eligible in a training partition iff eligible in at
    # least one retained cohort
    loso_eligible <- vapply(studies, function(s) {
      any(elig$eligible & elig$study != s)
    }, logical(1L))
    loso_hits <- vapply(
      studies, function(s) f %in% loso_supports[[s]], logical(1L)
    )
    training_rep <- if (any(loso_eligible)) {
      mean(loso_hits[loso_eligible])
    } else {
      NA_real_
    }
    fold_freq <- mean(loso_hits)

    in_pooled <- f %in% pooled_support
    dropped <- if (in_pooled) {
      drop_studies <- studies[!loso_hits]
      if (length(drop_studies) > 0L) {
        paste(drop_studies, collapse = ",")
      } else {
        NA_character_
      }
    } else {
      NA_character_
    }

    row <- data.frame(
      feature = f,
      n_measured_cohorts = n_measured,
      n_eligible_cohorts = n_eligible,
      within_study_frequency = within_freq,
      training_replication = training_rep,
      fold_frequency = fold_freq,
      in_pooled_support = in_pooled,
      loso_retention = if (in_pooled) mean(loso_hits) else NA_real_,
      dropped_on_removal_of = dropped,
      stringsAsFactors = FALSE
    )

    if (!is.null(effects)) {
      row <- cbind(row, .replication_direction_row(f, effects))
    }
    rows[[length(rows) + 1L]] <- row
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.replication_direction_row <- function(f, effects) {
  e <- effects[effects$feature == f, , drop = FALSE]
  if (nrow(e) < 2L) {
    return(data.frame(
      n_estimates = nrow(e), sign_concordance = NA_real_,
      pooled_estimate = NA_real_, pooled_se = NA_real_,
      pooled_ci_lo = NA_real_, pooled_ci_hi = NA_real_,
      het_Q = NA_real_, het_df = NA_integer_, het_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  w <- 1 / e$se^2
  pooled <- sum(w * e$estimate) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  q <- sum(w * (e$estimate - pooled)^2)
  df <- nrow(e) - 1L
  data.frame(
    n_estimates = nrow(e),
    sign_concordance = max(mean(e$estimate > 0), mean(e$estimate < 0)),
    pooled_estimate = pooled,
    pooled_se = pooled_se,
    pooled_ci_lo = pooled - 1.96 * pooled_se,
    pooled_ci_hi = pooled + 1.96 * pooled_se,
    het_Q = q,
    het_df = df,
    het_p = stats::pchisq(q, df = df, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Print a Cohort Replication Profile
#'
#' @param x A `stabl_replication_profile` object.
#' @param ... Ignored; present for S3 compatibility.
#'
#' @return The input `x`, invisibly.
#' @examples
#' set.seed(1)
#' ids <- paste0("s", 1:60)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(ids, paste0("f", 1:8)))
#' y <- setNames(x[, 1] * 3 + rnorm(60, sd = 0.8), ids)
#' studies <- setNames(rep(paste0("study", 1:3), each = 20L), ids)
#' rp <- stabl_replication_profile(
#'   x, y, studies,
#'   fit_spec = list(n_bootstraps = 50,
#'                   artificial_type = "random_permutation"),
#'   seed = 5
#' )
#' print(rp)
#' @export
print.stabl_replication_profile <- function(x, ...) {
  cat("STABL cohort replication profile\n")
  cat(sprintf(
    "  studies: %d | features: %d | pooled support: %d\n",
    length(levels(x$study_labels)), nrow(x$features),
    length(x$pooled_support)
  ))
  f <- x$features
  tr <- f$training_replication
  if (all(is.na(tr))) {
    cat("  training replication: no eligible feature-fold pairs\n")
  } else {
    cat(sprintf(
      paste0("  training replication: median %.2f (range %.2f-%.2f)",
             " over features\n"),
      stats::median(tr, na.rm = TRUE),
      min(tr, na.rm = TRUE),
      max(tr, na.rm = TRUE)
    ))
  }
  if (length(x$pooled_support) > 0L) {
    keep <- f[f$in_pooled_support, , drop = FALSE]
    cat("  pooled panel LOSO retention:\n")
    for (i in seq_len(nrow(keep))) {
      cat(sprintf(
        "    %s: %.2f%s\n", keep$feature[[i]], keep$loso_retention[[i]],
        if (!is.na(keep$dropped_on_removal_of[[i]])) {
          paste0(" (dropped on removal of: ",
                 keep$dropped_on_removal_of[[i]], ")")
        } else {
          ""
        }
      ))
    }
  }
  if (!is.null(x$heldout)) {
    mn <- x$heldout$metric_name[!is.na(x$heldout$metric_name)]
    mn <- if (length(mn) > 0L) mn[[1L]] else "n/a"
    cat(sprintf(
      "  held-out performance (%s): median %.3g over %d studies\n",
      mn, stats::median(x$heldout$metric, na.rm = TRUE), nrow(x$heldout)
    ))
  }
  cat(paste0(
    "  note: four distinct estimands (pooled stability, training-cohort\n",
    "  replication, direction replication, held-out performance);\n",
    "  selection agreement is not causal agreement. Cohort-level\n",
    "  observations only -- no binomial CIs over dependent folds.\n"
  ))
  invisible(x)
}
