#' Make Random-Permutation Artificial Features
#'
#' Randomly selects `n_injected` columns from `x`, copies them, and shuffles
#' each copy independently.  Mirrors the `"random_permutation"` branch of
#' `Stabl._make_artificial_features()` in the Python STABL library.
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of artificial columns to generate.
#'
#' @return Named list:
#'   \describe{
#'     \item{x_augmented}{Original matrix with artificial columns appended.}
#'     \item{noise_col_indices}{Integer vector (1-based) of original column
#'       indices selected as sources for the artificial block.}
#'     \item{artificial_provenance}{List describing the requested generator,
#'       `actual_type` used by the selected artificial columns, all
#'       `actual_types`, and the complete ordered `fallback_history`.}
#'   }
#' @keywords internal
make_rp_features <- function(x, n_injected) {
  n_features <- ncol(x)
  indices    <- sample.int(n = n_features, size = n_injected, replace = FALSE)
  x_art      <- x[, indices, drop = FALSE]
  for (i in seq_len(ncol(x_art))) {
    x_art[, i] <- sample(x_art[, i])
  }
  list(
    x_augmented      = cbind(x, x_art),
    noise_col_indices = indices,
    artificial_provenance = .make_artificial_provenance(
      requested_type = "random_permutation",
      selected_types = rep.int("random_permutation", n_injected),
      chunks = .artificial_chunk_record(
        chunk = 1L,
        requested_type = "random_permutation",
        actual_type = "random_permutation",
        n_columns = n_injected,
        fallback_reason = NA_character_
      ),
      fallback_histories = list(.empty_artificial_fallback_history())
    )
  )
}

.artificial_type_levels <- c(
  "random_permutation",
  "knockoff",
  "knockoff_equi",
  "knockoff_mvr"
)

.count_artificial_types <- function(types) {
  counts <- table(factor(types, levels = .artificial_type_levels))
  stats::setNames(as.integer(counts), names(counts))
}

.artificial_chunk_record <- function(chunk, requested_type, actual_type,
                                     n_columns, fallback_reason) {
  data.frame(
    chunk = as.integer(chunk),
    requested_type = requested_type,
    actual_type = actual_type,
    n_columns = as.integer(n_columns),
    fallback = !identical(requested_type, actual_type),
    fallback_reason = fallback_reason,
    stringsAsFactors = FALSE
  )
}

.empty_artificial_fallback_history <- function() {
  data.frame(
    step = integer(),
    from_type = character(),
    to_type = character(),
    reason = character(),
    condition_class = character(),
    stringsAsFactors = FALSE
  )
}

.artificial_fallback_event <- function(from_type, to_type, condition) {
  data.frame(
    step = 1L,
    from_type = from_type,
    to_type = to_type,
    reason = conditionMessage(condition),
    condition_class = paste(class(condition), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

.combine_artificial_fallback_histories <- function(histories, chunks) {
  if (length(histories) != nrow(chunks)) {
    stop("Artificial fallback history must have one entry per chunk.",
         call. = FALSE)
  }
  rows <- lapply(seq_along(histories), function(i) {
    history <- histories[[i]]
    if (!is.data.frame(history) ||
        !identical(names(history), names(.empty_artificial_fallback_history()))) {
      stop("Artificial fallback history has an invalid schema.", call. = FALSE)
    }
    if (!nrow(history)) return(NULL)
    history$step <- seq_len(nrow(history))
    cbind(
      data.frame(chunk = chunks$chunk[[i]], stringsAsFactors = FALSE),
      history
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      event = integer(),
      chunk = integer(),
      .empty_artificial_fallback_history(),
      stringsAsFactors = FALSE
    ))
  }
  history <- do.call(rbind, rows)
  rownames(history) <- NULL
  history <- cbind(
    data.frame(event = seq_len(nrow(history)), stringsAsFactors = FALSE),
    history
  )
  history
}

.make_artificial_provenance <- function(requested_type, selected_types, chunks,
                                        fallback_histories) {
  actual_types <- .artificial_type_levels[
    .count_artificial_types(selected_types) > 0L
  ]
  actual_type <- if (length(actual_types) == 1L) {
    actual_types[[1L]]
  } else {
    "mixed"
  }
  fallback_history <- .combine_artificial_fallback_histories(
    fallback_histories,
    chunks
  )
  list(
    requested_type = requested_type,
    actual_type = actual_type,
    actual_types = actual_types,
    n_generated = as.integer(length(selected_types)),
    n_chunks = as.integer(nrow(chunks)),
    chunk_type_counts = .count_artificial_types(chunks$actual_type),
    selected_type_counts = .count_artificial_types(selected_types),
    fallback_counts = .count_artificial_types(chunks$actual_type[chunks$fallback]),
    fallback_history = fallback_history,
    chunks = chunks
  )
}

.maximum_supported_knockoff_features <- 3000L

.assert_supported_knockoff_dimension <- function(x, requested_type) {
  n_features <- ncol(x)
  if (n_features > .maximum_supported_knockoff_features) {
    .abort_stablr(
      "stablr_high_dimensional_knockoff_unsupported",
      paste0(
        "`", requested_type, "` is scientifically supported for at most ",
        .maximum_supported_knockoff_features,
        " predictors; chunked knockoff construction is not supported."
      ),
      requested_type = requested_type,
      n_features = n_features,
      maximum_supported_features = .maximum_supported_knockoff_features
    )
  }
  invisible(x)
}

# Standardize columns for fixed-X knockoff construction.
#
# `knockoff::create.fixed()` documents its `x` argument as a "normalized
# n-by-p matrix of original variables", so centering/scaling is the caller's
#' Make Knockoff Artificial Features
#'
#' Generates **fixed-X** knockoff features via `knockoff::create.fixed()`, with
#' up to 3 000 predictors. Larger inputs fail closed because the former
#' chunked construction did not establish a valid global knockoff design.
#' Falls back to random-permutation features when the supported-dimensional
#' knockoff constructor fails (e.g., rank-deficient input).
#'
#' `create.fixed()` column-normalizes its input internally (scale only, no
#' centering) and returns knockoffs on that internal unit-norm scale.  The
#' knockoffs are therefore rescaled by the input column norms so that the
#' shipped pair `(x, xk)` satisfies the fixed-X Gram condition on the input's
#' own scale: `t(xk) xk = t(x) x` and `t(x) xk = t(x) x - diag(s)` with
#' `s_j >= 0`.  The normalization actually applied is recorded in
#' `artificial_provenance$standardization`.
#'
#' @section Caveat on bootstrap reuse:
#' Fixed-X knockoffs are generated once on the full design, so their
#' exchangeability guarantee does not survive the row subsampling of the
#' bootstrap loop: a knockoff row is only exchangeable with the matching
#' original row of that same realized design.  When correlation structure
#' matters (in particular with `sample_fraction < 1`), prefer the model-X
#' modes `"knockoff_equi"` or `"knockoff_mvr"`, whose knockoff rows are
#' i.i.d. draws from a fitted Gaussian model and remain valid under row
#' subsampling.
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details.
#'   `artificial_provenance$standardization` records the fixed-X scale
#'   contract (`knockoff_rescaled_to_input_scale = TRUE` when the fixed-X
#'   construction ran) and carries `standardized = FALSE` with NULL
#'   `means`/`sds` when the random-permutation fallback was taken.
#' @keywords internal
make_knockoff_features <- function(x, n_injected, random_state = NULL,
                                   knockoff_fallback = c("stop", "allow")) {
  knockoff_fallback <- match.arg(knockoff_fallback)
  .assert_supported_knockoff_dimension(x, "knockoff")
  .require_pkg("knockoff", "for artificial_type = \"knockoff\"")

  # NOTE: Seeding is the dispatcher's responsibility (see
  # `make_artificial_features`).  Re-seeding here would mask any RNG
  # consumed by the caller and is intentionally omitted (audit M-5).
  n_features <- ncol(x)

  .make_ko_chunk <- function(x_chunk) {
    tryCatch(
      {
        if (nrow(x_chunk) <= ncol(x_chunk)) {
          .abort_numerical_infeasibility(
            "stablr_knockoff_infeasible",
            "Input X must have dimensions n > p",
            requested_type = "knockoff",
            n_samples = nrow(x_chunk),
            n_features = ncol(x_chunk)
          )
        }
        augmented_rows <- FALSE
        knockoff_result <- .with_linear_algebra_numerical_infeasibility(
          withCallingHandlers(
            knockoff::create.fixed(x_chunk, sigma = 1),
            warning = function(w) {
              if (grepl("Augmenting the model with extra rows",
                        conditionMessage(w), fixed = TRUE)) {
                augmented_rows <<- TRUE
                invokeRestart("muffleWarning")
              }
            }
          ),
          source = "knockoff::create.fixed",
          subclass = "stablr_knockoff_infeasible"
        )
        # Scale contract (0.2.2): create.fixed() column-normalizes its input
        # (normc, center = FALSE) and returns BOTH $X (the normalized design)
        # and $Xk (knockoffs on that unit-norm scale, entries ~ 1/sqrt(n)).
        # The fixed-X Gram condition is preserved under diagonal rescaling of
        # BOTH blocks (s_j -> s_j * norm_j^2), so rescaling $Xk by the INPUT
        # column norms back-transforms the knockoffs to the input scale
        # exactly.  The shipped pair (x, xk) then satisfies, on x's own
        # scale, t(xk) xk = t(x) x and t(x) xk = t(x) x - diag(s * norms^2).
        # (0.2.0-0.2.1 centered/scaled the input first and shipped the
        # knockoffs on the standardized scale next to the raw design -- a
        # mismatched pair whenever input columns were heterogeneous.)
        col_norms <- sqrt(colSums(x_chunk^2))
        xk <- knockoff_result$Xk %*% diag(col_norms)
        colnames(xk) <- colnames(x_chunk)
        if (isTRUE(augmented_rows) || nrow(xk) != nrow(x_chunk)) {
          .abort_numerical_infeasibility(
            "stablr_knockoff_infeasible",
            paste0(
              "knockoff::create.fixed augmented the design from ",
              nrow(x_chunk), " to ", nrow(xk),
              " rows; stablr cannot align augmented knockoffs with the ",
              "unaugmented outcome."
            ),
            requested_type = "knockoff",
            n_samples = nrow(x_chunk),
            returned_rows = nrow(xk)
          )
        }
        list(
          x_art = xk,
          actual_type = "knockoff",
          fallback_reason = NA_character_,
          fallback_history = .empty_artificial_fallback_history(),
          standardization = list(
            standardized = FALSE,
            means = NULL,
            sds = NULL,
            normalization = "scale_only_internal_to_create.fixed",
            knockoff_rescaled_to_input_scale = TRUE
          )
        )
      },
      stablr_numerical_infeasibility = function(e) {
        reason <- conditionMessage(e)
        if (identical(knockoff_fallback, "stop")) {
          .abort_strict_knockoff_fallback("knockoff", reason)
        }
        warning(
          "knockoff::create.fixed failed; falling back to random ",
          "permutation for this chunk. Reason: ", reason,
          call. = FALSE
        )
        list(
          x_art = make_rp_features(x_chunk, ncol(x_chunk))$x_augmented[
            , ncol(x_chunk) + seq_len(ncol(x_chunk)), drop = FALSE
          ],
          actual_type = "random_permutation",
          fallback_reason = reason,
          fallback_history = .artificial_fallback_event(
            "knockoff",
            "random_permutation",
            e
          ),
          standardization = list(
            standardized = FALSE,
            means = NULL,
            sds = NULL
          )
        )
      }
    )
  }

  chunk_result <- .make_ko_chunk(x)
  x_art_full <- chunk_result$x_art
  orig_map <- seq_len(n_features)
  type_map <- rep.int(chunk_result$actual_type, n_features)
  chunks <- .artificial_chunk_record(
    chunk = 1L,
    requested_type = "knockoff",
    actual_type = chunk_result$actual_type,
    n_columns = ncol(chunk_result$x_art),
    fallback_reason = chunk_result$fallback_reason
  )
  fallback_histories <- list(chunk_result$fallback_history)

  # Aligned selection: injected column k must be the artificial counterpart
  # of feature sel_idx[k] -- the bootstrap accumulates scores positionally
  # (real_rows vs art_rows) and the fixed-X / model-X knockoff exchangeability
  # requires each knockoff to stay paired with its own feature.  A random
  # FEATURE subset (sorted) keeps the null pool valid; a random column
  # permutation of the knockoff matrix does not.
  sel_idx <- sort(sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE))
  x_art   <- x_art_full[, sel_idx, drop = FALSE]
  selected_types <- type_map[sel_idx]

  provenance <- .make_artificial_provenance(
    requested_type = "knockoff",
    selected_types = selected_types,
    chunks = chunks,
    fallback_histories = fallback_histories
  )
  # Record the fixed-X normalization/rescale contract actually applied
  # (NULL means/sds and no normalization tag when the random-permutation
  # fallback was taken, since that path uses no knockoff construction).
  provenance$standardization <- chunk_result$standardization

  list(
    x_augmented      = cbind(x, x_art),
    # Return original-feature indices (not x_art_full indices) so that
    # .append_noise_groups in stabl_fit.R can look up SGL groups correctly.
    noise_col_indices = orig_map[sel_idx],
    artificial_provenance = provenance
  )
}

# ── Shared helpers for model-X generators ─────────────────────────────────────

# A1 (0.2.2) strict fallback policy.  A knockoff -> random-permutation
# fallback silently changes the statistical guarantees: permutation features
# do not support FDP+ calibration, so a fit that requested a knockoff type
# but silently received permuted nulls reports FDP+ estimates whose target
# estimand changed.  The fallback is therefore an error unless the caller
# explicitly opts in with knockoff_fallback = "allow".  (The benign
# knockoff_mvr -> knockoff_equi fallback is NOT affected: both are model-X
# constructions with valid FDP+ calibration, so it stays a warning.)
.abort_strict_knockoff_fallback <- function(requested_type, reason) {
  stop(
    "knockoff construction failed for artificial_type = \"", requested_type,
    "\" and the random-permutation fallback is disabled ",
    "(knockoff_fallback = \"stop\", the strict default). Reason: ", reason,
    " Random-permutation features do not support FDP+ calibration; repair ",
    "the input (e.g., reduce p, remove constant or duplicate columns) or ",
    "explicitly opt in with knockoff_fallback = \"allow\".",
    call. = FALSE
  )
}

# Estimate a positive-definite covariance matrix from x.
# Falls back to corpcor shrinkage when cov(x) is not PD.
.estimate_pd_sigma <- function(x) {
  Sigma <- cov(x)
  min_eig <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig > 0) return(Sigma)
  # Rank-deficient (p > n) or near-singular: use Ledoit-Wolf shrinkage if available
  if (requireNamespace("corpcor", quietly = TRUE)) {
    Sigma <- suppressMessages(
      corpcor::make.positive.definite(Sigma, tol = 1e-3)
    )
  } else {
    diag(Sigma) <- diag(Sigma) + (1e-5 - min_eig)
  }
  Sigma
}

#' Make Model-X Equicorrelated Knockoff Artificial Features
#'
#' Generates **model-X equicorrelated** knockoff features via
#' `knockoff::create.gaussian(..., method = "equi")`.  This matches the
#' `GaussianSampler(X, method='equicorrelated')` call used by the Python STABL
#' library, making it the parity-correct knockoff type for cross-language
#' comparisons. Inputs above 3 000 predictors fail closed because a globally
#' valid chunked construction has not been established. Supported-dimensional
#' constructor failures fall back to random-permutation features.
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details.
#' @keywords internal
make_knockoff_equi_features <- function(x, n_injected, random_state = NULL,
                                        knockoff_fallback = c("stop", "allow")) {
  knockoff_fallback <- match.arg(knockoff_fallback)
  .assert_supported_knockoff_dimension(x, "knockoff_equi")
  .require_pkg("knockoff", "for artificial_type = \"knockoff_equi\"")

  n_features <- ncol(x)

  .make_equi_chunk <- function(x_chunk) {
    tryCatch(
      {
        mu    <- colMeans(x_chunk)
        Sigma <- .with_linear_algebra_numerical_infeasibility(
          .estimate_pd_sigma(x_chunk),
          source = "equicorrelated covariance estimation",
          subclass = "stablr_knockoff_infeasible"
        )
        list(
          x_art = .with_linear_algebra_numerical_infeasibility(
            knockoff::create.gaussian(x_chunk, mu, Sigma, method = "equi"),
            source = "knockoff::create.gaussian(equi)",
            subclass = "stablr_knockoff_infeasible"
          ),
          actual_type = "knockoff_equi",
          fallback_reason = NA_character_,
          fallback_history = .empty_artificial_fallback_history()
        )
      },
      stablr_numerical_infeasibility = function(e) {
        reason <- conditionMessage(e)
        if (identical(knockoff_fallback, "stop")) {
          .abort_strict_knockoff_fallback("knockoff_equi", reason)
        }
        warning(
          "knockoff_equi: create.gaussian failed; falling back to random ",
          "permutation for this chunk. Reason: ", reason,
          call. = FALSE
        )
        list(
          x_art = make_rp_features(x_chunk, ncol(x_chunk))$x_augmented[
            , ncol(x_chunk) + seq_len(ncol(x_chunk)), drop = FALSE
          ],
          actual_type = "random_permutation",
          fallback_reason = reason,
          fallback_history = .artificial_fallback_event(
            "knockoff_equi",
            "random_permutation",
            e
          )
        )
      }
    )
  }

  chunk_result <- .make_equi_chunk(x)
  x_art_full <- chunk_result$x_art
  orig_map <- seq_len(n_features)
  type_map <- rep.int(chunk_result$actual_type, n_features)
  chunks <- .artificial_chunk_record(
    chunk = 1L,
    requested_type = "knockoff_equi",
    actual_type = chunk_result$actual_type,
    n_columns = ncol(chunk_result$x_art),
    fallback_reason = chunk_result$fallback_reason
  )
  fallback_histories <- list(chunk_result$fallback_history)

  # Aligned selection: injected column k must be the artificial counterpart
  # of feature sel_idx[k] -- the bootstrap accumulates scores positionally
  # (real_rows vs art_rows) and the fixed-X / model-X knockoff exchangeability
  # requires each knockoff to stay paired with its own feature.  A random
  # FEATURE subset (sorted) keeps the null pool valid; a random column
  # permutation of the knockoff matrix does not.
  sel_idx <- sort(sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE))
  x_art   <- x_art_full[, sel_idx, drop = FALSE]
  selected_types <- type_map[sel_idx]

  list(
    x_augmented      = cbind(x, x_art),
    noise_col_indices = orig_map[sel_idx],
    artificial_provenance = .make_artificial_provenance(
      requested_type = "knockoff_equi",
      selected_types = selected_types,
      chunks = chunks,
      fallback_histories = fallback_histories
    )
  )
}

#' Make Model-X MVR Knockoff Artificial Features
#'
#' Generates **model-X MVR (minimum-variance-reconstructability)** knockoff
#' features.  The S-matrix is solved via `solve_mvr()` (an R-native coordinate-
#' descent port of `knockpy.mrc._solve_mvr_ungrouped`), then the knockoff
#' sample is drawn with `knockoff::create.gaussian(..., diag_s = S)`.
#' This is a novel feature exclusive to `stablr` — the Python STABL library
#' does not implement MVR knockoffs. Inputs above 3 000 predictors fail closed
#' because global exchangeability for a chunked construction has not been
#' established. For supported dimensions, MVR-solver failure falls back to
#' equi; `create.gaussian` failure falls back to random permutation.
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed; passed to `solve_mvr()` for the
#'   coordinate-shuffle RNG.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details. The retained
#'   `artificial_provenance$mvr_chunking` compatibility field always records
#'   that chunking was not applied.
#' @keywords internal
make_knockoff_mvr_features <- function(x, n_injected, random_state = NULL,
                                       knockoff_fallback = c("stop", "allow")) {
  knockoff_fallback <- match.arg(knockoff_fallback)
  .assert_supported_knockoff_dimension(x, "knockoff_mvr")
  .require_pkg("knockoff", "for artificial_type = \"knockoff_mvr\"")

  n_features <- ncol(x)

  .make_mvr_chunk <- function(x_chunk) {
    fallback_history <- .empty_artificial_fallback_history()
    tryCatch(
      {
        mu    <- colMeans(x_chunk)
        Sigma <- .with_linear_algebra_numerical_infeasibility(
          .estimate_pd_sigma(x_chunk),
          source = "MVR covariance estimation",
          subclass = "stablr_mvr_infeasible"
        )
        equi_fallback_reason <- NA_character_

        # Attempt MVR S-solve; fall back to equi on solver failure
        S_diag <- tryCatch(
          .with_linear_algebra_numerical_infeasibility(
            solve_mvr(Sigma, random_state = random_state),
            source = "MVR S-matrix solver",
            subclass = "stablr_mvr_infeasible"
          ),
          stablr_numerical_infeasibility = function(e) {
            equi_fallback_reason <<- conditionMessage(e)
            fallback_history <<- rbind(
              fallback_history,
              .artificial_fallback_event(
                "knockoff_mvr",
                "knockoff_equi",
                e
              )
            )
            warning("solve_mvr failed; using equi S for this chunk. Reason: ",
                    equi_fallback_reason, call. = FALSE)
            NULL
          }
        )

        if (is.null(S_diag)) {
          list(
            x_art = .with_linear_algebra_numerical_infeasibility(
              knockoff::create.gaussian(
                x_chunk,
                mu,
                Sigma,
                method = "equi"
              ),
              source = "knockoff::create.gaussian(equi fallback)",
              subclass = "stablr_knockoff_infeasible"
            ),
            actual_type = "knockoff_equi",
            fallback_reason = equi_fallback_reason,
            fallback_history = fallback_history
          )
        } else {
          list(
            x_art = .with_linear_algebra_numerical_infeasibility(
              knockoff::create.gaussian(
                x_chunk,
                mu,
                Sigma,
                diag_s = S_diag
              ),
              source = "knockoff::create.gaussian(MVR)",
              subclass = "stablr_knockoff_infeasible"
            ),
            actual_type = "knockoff_mvr",
            fallback_reason = NA_character_,
            fallback_history = fallback_history
          )
        }
      },
      stablr_numerical_infeasibility = function(e) {
        reason <- conditionMessage(e)
        if (identical(knockoff_fallback, "stop")) {
          .abort_strict_knockoff_fallback("knockoff_mvr", reason)
        }
        from_type <- if (nrow(fallback_history)) {
          fallback_history$to_type[[nrow(fallback_history)]]
        } else {
          "knockoff_mvr"
        }
        fallback_history <- rbind(
          fallback_history,
          .artificial_fallback_event(
            from_type,
            "random_permutation",
            e
          )
        )
        warning(
          "knockoff_mvr: create.gaussian failed; falling back to random ",
          "permutation for this chunk. Reason: ", reason,
          call. = FALSE
        )
        list(
          x_art = make_rp_features(x_chunk, ncol(x_chunk))$x_augmented[
            , ncol(x_chunk) + seq_len(ncol(x_chunk)), drop = FALSE
          ],
          actual_type = "random_permutation",
          fallback_reason = reason,
          fallback_history = fallback_history
        )
      }
    )
  }

  chunk_result <- .make_mvr_chunk(x)
  x_art_full <- chunk_result$x_art
  orig_map <- seq_len(n_features)
  type_map <- rep.int(chunk_result$actual_type, n_features)
  chunks <- .artificial_chunk_record(
    chunk = 1L,
    requested_type = "knockoff_mvr",
    actual_type = chunk_result$actual_type,
    n_columns = ncol(chunk_result$x_art),
    fallback_reason = chunk_result$fallback_reason
  )
  fallback_histories <- list(chunk_result$fallback_history)

  # Aligned selection: injected column k must be the artificial counterpart
  # of feature sel_idx[k] -- the bootstrap accumulates scores positionally
  # (real_rows vs art_rows) and the fixed-X / model-X knockoff exchangeability
  # requires each knockoff to stay paired with its own feature.  A random
  # FEATURE subset (sorted) keeps the null pool valid; a random column
  # permutation of the knockoff matrix does not.
  sel_idx <- sort(sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE))
  x_art   <- x_art_full[, sel_idx, drop = FALSE]
  selected_types <- type_map[sel_idx]

  provenance <- .make_artificial_provenance(
    requested_type = "knockoff_mvr",
    selected_types = selected_types,
    chunks = chunks,
    fallback_histories = fallback_histories
  )
  provenance$mvr_chunking <- list(
    applied = FALSE,
    chunk_size = .maximum_supported_knockoff_features,
    approximate = FALSE,
    global_exchangeability_established = NA
  )

  list(
    x_augmented      = cbind(x, x_art),
    noise_col_indices = orig_map[sel_idx],
    artificial_provenance = provenance
  )
}

#' Dispatcher for Artificial Feature Generation
#'
#' Selects and calls the appropriate artificial-feature generator based on
#' `type`, returning the augmented predictor matrix together with the column
#' indices of the injected noise block.
#'
#' Injecting artificial features supports STABL's FDP+ diagnostic calibration:
#' by mixing known-noise columns into the predictor matrix alongside real
#' features, STABL can empirically estimate how often a variable of pure noise
#' is selected at a given stability threshold.  This observed noise-selection
#' rate drives the FDP+ estimate computed in [compute_fdp_plus()]. The current
#' rule chooses its minimum and does not imply universal error-rate control.
#'
#' Four noise strategies are supported:
#' \describe{
#'   \item{`"random_permutation"`}{Copies `n_injected` randomly chosen real
#'     columns and shuffles each copy independently, breaking all signal while
#'     preserving marginal distributions.  Fast and broadly applicable.}
#'   \item{`"knockoff"`}{Generates **fixed-X** knockoffs via
#'     [knockoff::create.fixed()], which preserve the second-moment (Gram)
#'     structure of the original design on the input's own scale:
#'     `t(xk) xk = t(x) x` and `t(x) xk = t(x) x - diag(s)` with `s >= 0`
#'     (0.2.2 scale contract).  Kept for backward compatibility; note that
#'     fixed-X knockoffs are generated once on the full design, so their
#'     exchangeability does not survive the row subsampling of the bootstrap
#'     loop — prefer `"knockoff_equi"` or `"knockoff_mvr"` when correlation
#'     structure matters.}
#'   \item{`"knockoff_equi"`}{Generates **model-X equicorrelated** knockoffs
#'     via `knockoff::create.gaussian(..., method = "equi")`.  This matches the
#'     `GaussianSampler(method='equicorrelated')` call in Python STABL and is
#'     the parity-correct knockoff type for cross-language comparisons.}
#'   \item{`"knockoff_mvr"`}{Generates **model-X MVR** (minimum-variance-
#'     reconstructability) knockoffs.  The S-matrix is solved by `solve_mvr()`
#'     (an R-native port of `knockpy.mrc`); sampling uses
#'     `knockoff::create.gaussian(..., diag_s = S)`.  This is a novel feature
#'     exclusive to `stablr`.}
#' }
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#'   Must have more columns than `n_injected` for random permutation; for
#'   knockoffs, constructor failure (e.g., rank-deficient input) is governed
#'   by `knockoff_fallback`: by default (`"stop"`) it is an error, and with
#'   `knockoff_fallback = "allow"` a fallback to random permutation is
#'   attempted. Fixed-X knockoffs that require row augmentation also fail
#'   (or fall back under `"allow"`) because augmented knockoff rows cannot
#'   be aligned with the original outcome. All knockoff modes fail closed
#'   above 3 000 predictors; random permutation has no such dimensional
#'   limit.
#' @param n_injected Positive integer; number of artificial columns to append.
#'   Typically `round(ncol(x) * artificial_proportion)` as computed in
#'   [stabl_fit()].
#' @param type Character string; one of `"random_permutation"`, `"knockoff"`,
#'   `"knockoff_equi"`, or `"knockoff_mvr"`.
#' @param random_state Optional integer; passed to [set.seed()] before any
#'   random operations for reproducibility.  `NULL` leaves the RNG unchanged.
#'   Seeding happens exactly once in this dispatcher; downstream generators
#'   inherit the seeded RNG state and do not re-seed (audit M-5).
#' @param knockoff_fallback Character; the fallback policy when a knockoff
#'   constructor fails on a supported-dimensional input (audit A1, strict
#'   fallback policy).  `"stop"` (the default) raises an error: a
#'   knockoff-to-random-permutation fallback silently changes the FDP+
#'   calibration target, because permutation features do not support FDP+
#'   calibration.  `"allow"` restores the pre-0.2.2 behavior (warn and fall
#'   back to random-permutation features for the affected chunk; the
#'   `fallback_history` provenance records the event).  The benign
#'   `knockoff_mvr` -> `knockoff_equi` fallback is unaffected: both are
#'   model-X constructions with valid FDP+ calibration, so it remains a
#'   warning under either policy.
#'
#' @return Named list with three elements:
#'   \describe{
#'     \item{`x_augmented`}{Numeric matrix of size
#'       (nrow(x)) \eqn{\times} (ncol(x) + n_injected) with the artificial
#'       columns appended after the original features.}
#'     \item{`noise_col_indices`}{Integer vector of length `n_injected`
#'       containing the 1-based indices into the **original** `x` columns
#'       (not into the artificial block) that identify which source features
#'       were used to build each artificial column.  Used by [stabl_fit()]
#'       to look up sparse-group-lasso group memberships for the artificial
#'       block via `.append_noise_groups`.}
#'     \item{`artificial_provenance`}{Additive metadata describing the
#'       requested artificial-feature type, `actual_type` used for scientific
#'       classification, all generator modes used by the selected artificial
#'       columns, and the ordered structured `fallback_history`. Per-chunk
#'       counts and reasons remain available for 0.1.x compatibility.}
#'   }
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(
#'   rnorm(20 * 4), 20, 4,
#'   dimnames = list(paste0("s", 1:20), paste0("f", 1:4))
#' )
#' generated <- make_artificial_features(
#'   x, n_injected = 2L, type = "random_permutation", random_state = 1L
#' )
#' dim(generated$x_augmented)
#' generated$artificial_provenance$actual_type
#' try(make_artificial_features(x, 0L, "random_permutation"))
#'
#' @seealso [compute_fdp_plus()] which consumes the artificial-feature scores,
#'   [stabl_fit()] which calls this function automatically.
#' @export
make_artificial_features <- function(x, n_injected, type, random_state = NULL,
                                     knockoff_fallback = c("stop", "allow")) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite, non-missing values.", call. = FALSE)
  }
  n_injected <- .validate_scalar_integer_like(
    n_injected,
    "n_injected",
    min = 1L,
    max = ncol(x)
  )
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("`type` must be a single character string.", call. = FALSE)
  }
  if (!type %in% .artificial_type_levels) {
    stop(
      "`type` must be one of \"random_permutation\", \"knockoff\", ",
      "\"knockoff_equi\", or \"knockoff_mvr\", got: ", type,
      call. = FALSE
    )
  }
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }
  if (startsWith(type, "knockoff")) {
    .assert_supported_knockoff_dimension(x, type)
  }
  knockoff_fallback <- match.arg(knockoff_fallback)
  if (!is.null(random_state)) set.seed(random_state)
  switch(
    type,
    random_permutation = make_rp_features(x, n_injected),
    knockoff           = make_knockoff_features(x, n_injected,
                                                random_state = random_state,
                                                knockoff_fallback),
    knockoff_equi      = make_knockoff_equi_features(x, n_injected,
                                                     random_state = random_state,
                                                     knockoff_fallback),
    knockoff_mvr       = make_knockoff_mvr_features(x, n_injected,
                                                    random_state = random_state,
                                                    knockoff_fallback)
  )
}
