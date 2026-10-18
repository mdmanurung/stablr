# Outcome-null calibration workbench (1.0.0 spec, A2).
#
# stabl_null_check() repeatedly reruns the COMPLETE stabl_fit() pipeline on
# outcome-null versions of the data and summarizes what the selection
# machinery does when nothing is signal. Null generation is prescribed by an
# explicit exchangeability unit:
#
#   * "independent"    -- global permutation of y; legal only when samples are
#                        independent. Refused when block/cluster structure is
#                        supplied (that is exactly the case where global
#                        permutation silently breaks the design).
#   * "within_block"   -- permutation of y within blocks (e.g., sites,
#                        batches); requires `blocks`.
#   * "whole_cluster"  -- permutation of whole cluster outcome vectors across
#                        clusters (within-cluster outcome dependence is kept
#                        intact); requires `clusters`.
#   * "simulated"      -- an explicitly simulated null via a user-supplied
#                        `null_simulator`; the only design under which the
#                        truth is synthetic by construction.
#
# Unsupported or structurally contradicted designs are refused, never
# silently replaced by a global permutation.
#
# Truth-status distinction (spec A2): under a globally null DGP every
# selected feature is a false positive by construction. Permutation designs
# INDUCE such a global null and therefore quantify empirical selection
# behavior under it; they cannot validate FDP in arbitrary partly non-null
# observational data, and feature-by-feature label permutation is not a
# replacement for correctly modeled confounders. The result carries
# `truth_status` so downstream reporting cannot blur this.

#' Outcome-Null Calibration Check for STABL
#'
#' Repeatedly reruns the complete [stabl_fit()] pipeline on outcome-null
#' versions of the data (permuted or explicitly simulated) and summarizes
#' null selection behavior: the selected-panel-size distribution, the
#' probability of selecting anything at all, per-feature selection
#' frequencies, reported-FDP+ versus by-construction realized FDP, and
#' null-saturation proportions.
#'
#' @param x Numeric predictor matrix with non-empty row names (sample ids)
#'   and column names (feature ids), as for [stabl_fit()].
#' @param y Outcome vector, named to match `rownames(x)`.
#' @param design Exchangeability unit prescribing the null: one of
#'   `"independent"` (global permutation; refused when `blocks` or `clusters`
#'   are supplied), `"within_block"` (permute within `blocks`),
#'   `"whole_cluster"` (permute whole cluster outcome vectors across
#'   `clusters`), or `"simulated"` (use `null_simulator`).
#' @param fit_spec Named list of further arguments for [stabl_fit()], e.g.
#'   `list(n_bootstraps = 50, artificial_type = "knockoff_equi")`. The
#'   per-replicate `random_state` is always set by `stabl_null_check`;
#'   `lambda_grid` defaults to `"auto"` when not supplied.
#' @param n_reps Number of null replicates (permutations or simulated
#'   nulls). Total fits = `n_reps * decoy_reps`.
#' @param blocks Optional factor/vector of length `nrow(x)` giving block
#'   membership; required for `design = "within_block"`.
#' @param clusters Optional factor/vector of length `nrow(x)` giving cluster
#'   membership; required for `design = "whole_cluster"`.
#' @param null_simulator Required for `design = "simulated"`: a function
#'   `function(y, x, rep_id, seed)` returning a null outcome vector (named
#'   like `y`). Under this design the truth is synthetic: every selected
#'   feature is a false positive by construction.
#' @param decoy_reps Number of fit replicates per null replicate with
#'   distinct decoy-generator seeds; values > 1 separate sensitivity to the
#'   random decoy realization from sensitivity to the null permutation.
#' @param seed Optional base seed. Per-replicate permutation and fit seeds
#'   are derived from it by the package's hashed multi-stream derivation and
#'   stored in the result, so any replicate is exactly replayable. If NULL,
#'   a seed is drawn from the current RNG state and recorded.
#' @param max_failures Bound on tolerated fit failures. Every failure is
#'   recorded with its replicate id, seeds, and error message; if failures
#'   exceed `max_failures` the run aborts (fail-closed) rather than
#'   summarizing a mostly-failed design.
#' @param ... Reserved; unused arguments are an error.
#'
#' @return A `stabl_null_validation` object with elements:
#'   \describe{
#'     \item{design, truth_status}{The exchangeability unit and whether the
#'       null truth is `"synthetic_global_null"` (design `"simulated"`) or
#'       `"induced_global_null_permutation"` (permutation designs).}
#'     \item{n_reps, decoy_reps, n_failures, failures}{Replication counts
#'       and a data frame of bounded failure records.}
#'     \item{replicates}{Per-replicate results: seeds, panel size, selected
#'       features, reported FDP+ minimum and threshold, at-selection
#'       real/artificial proportions, and decoy-generator provenance.}
#'     \item{panel_sizes}{Integer vector of selected-panel sizes.}
#'     \item{p_any_selection}{Probability the null pipeline selects at
#'       least one feature.}
#'     \item{feature_frequency}{Named, decreasing vector of per-feature null
#'       selection frequencies.}
#'     \item{reported_fdp_plus}{Reported FDP+ minima across replicates;
#'       under the global null the realized FDP is 1 whenever the panel is
#'       non-empty, by construction.}
#'     \item{fit_spec, call}{The executed fit specification and call.}
#'   }
#'
#' @section Interpretation caveat:
#'   Null resampling quantifies empirical selection under a globally null
#'   DGP. It does not validate FDP control in partly non-null observational
#'   data, and label permutation is not a substitute for modeled
#'   confounders.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#' y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
#' nv <- stabl_null_check(
#'   x, y,
#'   design = "independent",
#'   fit_spec = list(n_bootstraps = 10,
#'                   artificial_type = "random_permutation"),
#'   n_reps = 3, seed = 7
#' )
#' nv$p_any_selection
#' @export
stabl_null_check <- function(x, y,
                             design = c("independent", "within_block",
                                        "whole_cluster", "simulated"),
                             fit_spec = list(),
                             n_reps = 20L,
                             blocks = NULL,
                             clusters = NULL,
                             null_simulator = NULL,
                             decoy_reps = 1L,
                             seed = NULL,
                             max_failures = NULL,
                             ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    stop("Unused arguments: ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }
  design <- match.arg(design)
  n_reps <- .validate_scalar_integer_like(n_reps, "n_reps", min = 1L)
  decoy_reps <- .validate_scalar_integer_like(
    decoy_reps, "decoy_reps", min = 1L
  )
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
         " -- data and per-replicate seeds are managed by stabl_null_check().",
         call. = FALSE)
  }
  if (!"lambda_grid" %in% names(fit_spec)) {
    fit_spec$lambda_grid <- "auto"
  }

  # ---- design legality: refuse, never silently substitute ------------------
  if (!is.null(blocks)) {
    blocks <- .validate_group_vector(blocks, x, "blocks")
  }
  if (!is.null(clusters)) {
    clusters <- .validate_group_vector(clusters, x, "clusters")
  }
  if (design == "independent" && (!is.null(blocks) || !is.null(clusters))) {
    stop(
      paste0(
        "design = \"independent\" applies a GLOBAL permutation of the ",
        "outcome, which is not exchangeable under block/cluster structure. ",
        "Use design = \"within_block\" or design = \"whole_cluster\", or ",
        "drop the structure only if you can defend sample independence."
      ),
      call. = FALSE
    )
  }
  if (design == "within_block") {
    if (is.null(blocks)) {
      stop("design = \"within_block\" requires `blocks`.", call. = FALSE)
    }
    if (length(unique(blocks)) < 2L) {
      stop("design = \"within_block\" needs at least two blocks; a single ",
           "block is a global permutation in disguise.", call. = FALSE)
    }
    if (any(table(blocks) < 2L)) {
      stop("Every block needs at least two samples for within-block ",
           "permutation.", call. = FALSE)
    }
  }
  if (design == "whole_cluster") {
    if (is.null(clusters)) {
      stop("design = \"whole_cluster\" requires `clusters`.", call. = FALSE)
    }
    if (length(unique(clusters)) < 3L) {
      stop("design = \"whole_cluster\" needs at least three clusters to ",
           "permute whole outcome vectors across clusters.", call. = FALSE)
    }
    if (length(unique(table(clusters))) != 1L) {
      stop("design = \"whole_cluster\" permutes whole cluster outcome ",
           "vectors across clusters and therefore requires equal-sized ",
           "clusters; use design = \"within_block\" or a simulated null ",
           "for unequal cluster sizes.", call. = FALSE)
    }
  }
  if (design == "simulated" && !is.function(null_simulator)) {
    stop("design = \"simulated\" requires `null_simulator`, a function ",
         "(y, x, rep_id, seed) returning a null outcome vector. An ",
         "explicit simulator is the point of this design: the package ",
         "will not invent a null DGP for you.", call. = FALSE)
  }

  if (is.null(max_failures)) {
    max_failures <- max(1L, floor(n_reps * decoy_reps / 2))
  }
  max_failures <- .validate_scalar_integer_like(
    max_failures, "max_failures", min = 0L
  )

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  seed <- .validate_scalar_integer_like(
    seed, "seed", min = 1L, max = .Machine$integer.max - 1L
  )
  perm_seeds <- .derive_stream_seeds(seed, n_reps)
  fit_seeds <- matrix(
    .derive_stream_seeds(seed + 1L, n_reps * decoy_reps),
    nrow = n_reps, ncol = decoy_reps
  )

  truth_status <- if (design == "simulated") {
    "synthetic_global_null"
  } else {
    "induced_global_null_permutation"
  }

  replicates <- list()
  failures <- list()
  n_done <- 0L
  for (rep_id in seq_len(n_reps)) {
    y_null <- tryCatch(
      .generate_null_outcome(
        y = y, x = x, design = design, blocks = blocks,
        clusters = clusters, null_simulator = null_simulator,
        rep_id = rep_id, seed = perm_seeds[[rep_id]]
      ),
      error = function(e) e
    )
    if (inherits(y_null, "error")) {
      failures[[length(failures) + 1L]] <- .null_failure_record(
        rep_id, NA_integer_, perm_seeds[[rep_id]], NA_integer_,
        "null_generation", conditionMessage(y_null)
      )
      next
    }
    for (decoy_id in seq_len(decoy_reps)) {
      n_done <- n_done + 1L
      fit_seed <- fit_seeds[rep_id, decoy_id]
      fit <- tryCatch(
        .with_local_seed(fit_seed, {
          do.call(
            stabl_fit,
            c(list(x = x, y = y_null, random_state = fit_seed), fit_spec)
          )
        }),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        failures[[length(failures) + 1L]] <- .null_failure_record(
          rep_id, decoy_id, perm_seeds[[rep_id]], fit_seed,
          "stabl_fit", conditionMessage(fit)
        )
      } else {
        replicates[[length(replicates) + 1L]] <-
          .null_replicate_record(fit, rep_id, decoy_id,
                                 perm_seeds[[rep_id]], fit_seed)
      }
      if (length(failures) > max_failures) {
        stop(
          sprintf(
            paste0(
              "stabl_null_check aborted: %d failures exceed max_failures ",
              "= %d after %d attempt(s). First failure: %s"
            ),
            length(failures), max_failures, n_done,
            failures[[1L]]$message
          ),
          call. = FALSE
        )
      }
    }
  }

  .assemble_null_validation(
    replicates = replicates,
    failures = failures,
    design = design,
    truth_status = truth_status,
    n_reps = n_reps,
    decoy_reps = decoy_reps,
    seed = seed,
    fit_spec = fit_spec,
    feature_names = colnames(x),
    call = match.call()
  )
}

# ---- null generators --------------------------------------------------------

.validate_group_vector <- function(g, x, arg) {
  if (length(g) != nrow(x)) {
    stop("`", arg, "` must have length nrow(x).", call. = FALSE)
  }
  g <- as.factor(g)
  if (anyNA(g)) {
    stop("`", arg, "` must not contain NA.", call. = FALSE)
  }
  g
}

.generate_null_outcome <- function(y, x, design, blocks, clusters,
                                   null_simulator, rep_id, seed) {
  if (design == "simulated") {
    out <- .with_local_seed(seed, null_simulator(y, x, rep_id, seed))
    if (length(out) != length(y)) {
      stop("`null_simulator` must return a vector of length length(y).",
           call. = FALSE)
    }
    if (is.null(names(out))) names(out) <- names(y)
    return(out)
  }
  .with_local_seed(seed, {
    if (design == "independent") {
      out <- sample(y)
    } else if (design == "within_block") {
      out <- y
      for (b in levels(blocks)) {
        idx <- which(blocks == b)
        out[idx] <- sample(y[idx])
      }
    } else {
      # whole_cluster: reassign whole cluster outcome vectors across clusters
      cl <- levels(clusters)
      perm <- sample(cl)
      out <- y
      for (i in seq_along(cl)) {
        out[clusters == perm[[i]]] <- y[clusters == cl[[i]]]
      }
    }
    out
  })
}

# ---- per-replicate records --------------------------------------------------

.null_failure_record <- function(rep_id, decoy_id, perm_seed, fit_seed,
                                 stage, message) {
  data.frame(
    rep_id = rep_id, decoy_id = decoy_id, perm_seed = perm_seed,
    fit_seed = fit_seed, stage = stage, message = message,
    stringsAsFactors = FALSE
  )
}

.null_replicate_record <- function(fit, rep_id, decoy_id, perm_seed,
                                   fit_seed) {
  support_mask <- get_support(fit)
  support <- names(support_mask)[support_mask]
  selected_threshold <- .resolve_stabl_threshold(fit)
  max_real <- rowMaxs(fit$stabl_scores_)
  prop_real <- sum(max_real > selected_threshold) / length(max_real)
  prop_art <- NA_real_
  if (!is.null(fit$stabl_scores_artificial_)) {
    max_art <- rowMaxs(fit$stabl_scores_artificial_)
    prop_art <- sum(max_art > selected_threshold) / length(max_art)
  }
  prov <- fit$artificial_provenance
  list(
    rep_id = rep_id,
    decoy_id = decoy_id,
    perm_seed = perm_seed,
    fit_seed = fit_seed,
    panel_size = length(support),
    selected = support,
    reported_fdp_plus = fit$min_fdr_,
    selected_threshold = selected_threshold,
    fdp_thresholds = fit$fdr_threshold_range,
    fdp_curve = fit$FDRs_,
    prop_real_at_selection = prop_real,
    prop_artificial_at_selection = prop_art,
    decoy_requested = prov$requested_type %||% NA_character_,
    decoy_actual = prov$actual_type %||% NA_character_,
    decoy_fallbacks = length(prov$fallback_history %||% list())
  )
}

.assemble_null_validation <- function(replicates, failures, design,
                                      truth_status, n_reps, decoy_reps,
                                      seed, fit_spec, feature_names, call) {
  panel_sizes <- vapply(
    replicates, function(r) r$panel_size, integer(1L)
  )
  freq <- numeric(0L)
  if (length(replicates) > 0L) {
    all_sel <- unlist(lapply(replicates, function(r) r$selected))
    if (length(all_sel) > 0L) {
      freq <- sort(table(all_sel) / length(replicates), decreasing = TRUE)
    }
  }
  structure(
    list(
      design = design,
      truth_status = truth_status,
      n_reps = n_reps,
      decoy_reps = decoy_reps,
      n_successful = length(replicates),
      n_failures = length(failures),
      failures = if (length(failures) > 0L) {
        do.call(rbind, failures)
      } else {
        data.frame(
          rep_id = integer(0), decoy_id = integer(0),
          perm_seed = integer(0), fit_seed = integer(0),
          stage = character(0), message = character(0),
          stringsAsFactors = FALSE
        )
      },
      replicates = replicates,
      panel_sizes = panel_sizes,
      p_any_selection = if (length(panel_sizes) > 0L) {
        mean(panel_sizes > 0L)
      } else {
        NA_real_
      },
      feature_frequency = freq,
      reported_fdp_plus = vapply(
        replicates, function(r) r$reported_fdp_plus %||% NA_real_,
        numeric(1L)
      ),
      seed = seed,
      fit_spec = fit_spec,
      feature_names = feature_names,
      call = call
    ),
    class = "stabl_null_validation"
  )
}

#' Print an Outcome-Null Validation Summary
#'
#' @param x A `stabl_null_validation` object from [stabl_null_check()].
#' @param ... Ignored; present for S3 compatibility.
#'
#' @return The input `x`, invisibly.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#' y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
#' nv <- stabl_null_check(
#'   x, y,
#'   fit_spec = list(n_bootstraps = 10,
#'                   artificial_type = "random_permutation"),
#'   n_reps = 3, seed = 7
#' )
#' print(nv)
#' @export
print.stabl_null_validation <- function(x, ...) {
  cat("STABL outcome-null validation\n")
  cat(sprintf(
    "  design: %s | truth: %s\n", x$design, x$truth_status
  ))
  cat(sprintf(
    "  replicates: %d successful / %d failed (of %d attempted)\n",
    x$n_successful, x$n_failures, x$n_reps * x$decoy_reps
  ))
  cat(sprintf(
    "  P(any selection on null): %.3f | median panel size: %s\n",
    x$p_any_selection,
    if (length(x$panel_sizes) > 0L) {
      format(stats::median(x$panel_sizes))
    } else {
      "n/a"
    }
  ))
  rep_fdp <- x$reported_fdp_plus[!is.na(x$reported_fdp_plus)]
  if (length(rep_fdp) > 0L) {
    cat(sprintf(
      paste0(
        "  reported FDP+ minima: median %.3g (range %.3g-%.3g); realized ",
        "FDP on this null is 1 whenever the panel is non-empty\n"
      ),
      stats::median(rep_fdp), min(rep_fdp), max(rep_fdp)
    ))
  }
  if (length(x$feature_frequency) > 0L) {
    top <- head(x$feature_frequency, 5L)
    cat("  most-selected null features:\n")
    for (nm in names(top)) {
      cat(sprintf("    %s: %.2f\n", nm, unname(top[[nm]])))
    }
  } else {
    cat("  no feature was ever selected on the null replicates\n")
  }
  if (x$truth_status == "induced_global_null_permutation") {
    cat(paste0(
      "  note: permutation induces a global null; this quantifies null ",
      "selection behavior,\n  not FDP validity for partly non-null ",
      "observational data.\n"
    ))
  }
  invisible(x)
}

#' Plot Outcome-Null Validation Results
#'
#' @param x A `stabl_null_validation` object from [stabl_null_check()].
#' @param which Character vector of panels: `"panel_sizes"` (null
#'   panel-size distribution), `"feature_frequency"` (per-feature null
#'   selection frequency), `"reported_vs_realized"` (reported FDP+ minima
#'   against the by-construction realized FDP of 1 on non-empty null
#'   panels), and `"saturation"` (real vs artificial at-selection
#'   proportions per replicate).
#' @param ... Ignored.
#'
#' @return Invisibly, a named list of the drawn `ggplot` objects.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- matrix(rnorm(60 * 8), 60, 8,
#'               dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#'   y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
#'   nv <- stabl_null_check(
#'     x, y,
#'     fit_spec = list(n_bootstraps = 10,
#'                     artificial_type = "random_permutation"),
#'     n_reps = 3, seed = 7
#'   )
#'   plot(nv, which = "panel_sizes")
#' }
#' @export
plot.stabl_null_validation <- function(
    x,
    which = c("panel_sizes", "feature_frequency", "reported_vs_realized",
              "saturation"),
    ...) {
  .require_ggplot2()
  which <- match.arg(which, several.ok = TRUE)
  plots <- list()

  if ("panel_sizes" %in% which) {
    df <- data.frame(panel_size = x$panel_sizes)
    plots$panel_sizes <- ggplot2::ggplot(
      df, ggplot2::aes(x = panel_size)
    ) +
      ggplot2::geom_histogram(
        binwidth = 1, boundary = -0.5, color = "grey30", fill = "grey70"
      ) +
      ggplot2::labs(
        title = "Null selected-panel sizes",
        x = "panel size", y = "replicates"
      )
  }

  if ("feature_frequency" %in% which && length(x$feature_frequency) > 0L) {
    df <- data.frame(
      feature = factor(
        names(x$feature_frequency),
        levels = rev(names(x$feature_frequency))
      ),
      frequency = as.numeric(x$feature_frequency)
    )
    plots$feature_frequency <- ggplot2::ggplot(
      df, ggplot2::aes(x = feature, y = frequency)
    ) +
      ggplot2::geom_col(fill = "grey50") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Per-feature null selection frequency",
        x = NULL, y = "frequency"
      )
  }

  if ("reported_vs_realized" %in% which && length(x$panel_sizes) > 0L) {
    df <- data.frame(
      rep = seq_along(x$panel_sizes),
      reported = x$reported_fdp_plus,
      realized = ifelse(x$panel_sizes > 0L, 1, NA_real_)
    )
    plots$reported_vs_realized <- ggplot2::ggplot(
      df, ggplot2::aes(x = rep)
    ) +
      ggplot2::geom_point(ggplot2::aes(y = reported), color = "grey30") +
      ggplot2::geom_point(
        data = df[!is.na(df$realized), , drop = FALSE],
        ggplot2::aes(y = realized), shape = 4, color = "#B00"
      ) +
      ggplot2::labs(
        title = "Reported FDP+ minima vs by-construction realized FDP",
        x = "replicate", y = "FDP"
      )
  }

  if ("saturation" %in% which && length(x$replicates) > 0L) {
    df <- data.frame(
      rep = seq_along(x$replicates),
      real = vapply(
        x$replicates, function(r) r$prop_real_at_selection, numeric(1L)
      ),
      artificial = vapply(
        x$replicates,
        function(r) r$prop_artificial_at_selection %||% NA_real_,
        numeric(1L)
      )
    )
    df <- df[stats::complete.cases(df), , drop = FALSE]
    if (nrow(df) > 0L) {
      plots$saturation <- ggplot2::ggplot(
        df, ggplot2::aes(x = real, y = artificial)
      ) +
        ggplot2::geom_point() +
        ggplot2::labs(
          title = "Null-saturation at the selected threshold",
          x = "proportion of real features selected",
          y = "proportion of artificial features selected"
        )
    }
  }

  for (p in plots) print(p)
  invisible(plots)
}
