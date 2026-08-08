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

#' Make Knockoff Artificial Features
#'
#' Generates **fixed-X** knockoff features via `knockoff::create.fixed()`, with
#' column-chunking for datasets that exceed 3 000 features (mirroring the
#' Python STABL implementation that chunks calls to `GaussianSampler`).
#' Falls back to random-permutation features when the knockoff constructor
#' fails (e.g., rank-deficient input).
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details.
#' @keywords internal
make_knockoff_features <- function(x, n_injected, random_state = NULL) {
  .require_pkg("knockoff", "for artificial_type = \"knockoff\"")

  # NOTE: Seeding is the dispatcher's responsibility (see
  # `make_artificial_features`).  Re-seeding here would mask any RNG
  # consumed by the caller and is intentionally omitted (audit M-5).
  n_features <- ncol(x)
  chunk_size <- 3000L

  .make_ko_chunk <- function(x_chunk) {
    tryCatch(
      {
        augmented_rows <- FALSE
        xk <- withCallingHandlers(
          knockoff::create.fixed(x_chunk, sigma = 1)$Xk,
          warning = function(w) {
            if (grepl("Augmenting the model with extra rows",
                      conditionMessage(w), fixed = TRUE)) {
              augmented_rows <<- TRUE
              invokeRestart("muffleWarning")
            }
          }
        )
        if (isTRUE(augmented_rows) || nrow(xk) != nrow(x_chunk)) {
          stop(
            "knockoff::create.fixed augmented the design from ",
            nrow(x_chunk), " to ", nrow(xk), " rows; stablr cannot align ",
            "augmented knockoffs with the unaugmented outcome.",
            call. = FALSE
          )
        }
        list(
          x_art = xk,
          actual_type = "knockoff",
          fallback_reason = NA_character_,
          fallback_history = .empty_artificial_fallback_history()
        )
      },
      error = function(e) {
        reason <- conditionMessage(e)
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
          )
        )
      }
    )
  }

  if (n_features > chunk_size) {
    n_chunks          <- ceiling(n_features / chunk_size)
    ko_blocks         <- vector("list", n_chunks)
    orig_map_blocks   <- vector("list", n_chunks)  # track source original-feature indices
    type_blocks       <- vector("list", n_chunks)
    chunk_records     <- vector("list", n_chunks)
    fallback_histories <- vector("list", n_chunks)
    for (i in seq_len(n_chunks)) {
      col_idx              <- sample.int(n_features, size = min(chunk_size, n_features),
                                        replace = FALSE)
      chunk_result         <- .make_ko_chunk(x[, col_idx, drop = FALSE])
      ko_blocks[[i]]       <- chunk_result$x_art
      orig_map_blocks[[i]] <- col_idx  # j-th column of ko_blocks[[i]] is knockoff of col_idx[j]
      type_blocks[[i]]     <- rep.int(chunk_result$actual_type, ncol(chunk_result$x_art))
      fallback_histories[[i]] <- chunk_result$fallback_history
      chunk_records[[i]]   <- .artificial_chunk_record(
        chunk = i,
        requested_type = "knockoff",
        actual_type = chunk_result$actual_type,
        n_columns = ncol(chunk_result$x_art),
        fallback_reason = chunk_result$fallback_reason
      )
    }
    x_art_full <- do.call(cbind, ko_blocks)   # n_samples × (n_chunks * chunk_size)
    orig_map   <- unlist(orig_map_blocks)      # maps each x_art_full col -> original feature idx
    type_map   <- unlist(type_blocks, use.names = FALSE)
    chunks     <- do.call(rbind, chunk_records)
    keep_idx   <- sample.int(ncol(x_art_full), size = n_features, replace = FALSE)
    x_art_full <- x_art_full[, keep_idx, drop = FALSE]
    orig_map   <- orig_map[keep_idx]           # keep map in sync after trim
    type_map   <- type_map[keep_idx]
  } else {
    chunk_result <- .make_ko_chunk(x)
    x_art_full <- chunk_result$x_art  # n_samples × n_features, same column order as x
    orig_map   <- seq_len(n_features)  # identity mapping: col j is knockoff of feature j
    type_map   <- rep.int(chunk_result$actual_type, n_features)
    chunks     <- .artificial_chunk_record(
      chunk = 1L,
      requested_type = "knockoff",
      actual_type = chunk_result$actual_type,
      n_columns = ncol(chunk_result$x_art),
      fallback_reason = chunk_result$fallback_reason
    )
    fallback_histories <- list(chunk_result$fallback_history)
  }

  sel_idx <- sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE)
  x_art   <- x_art_full[, sel_idx, drop = FALSE]
  selected_types <- type_map[sel_idx]

  list(
    x_augmented      = cbind(x, x_art),
    # Return original-feature indices (not x_art_full indices) so that
    # .append_noise_groups in stabl_fit.R can look up SGL groups correctly.
    noise_col_indices = orig_map[sel_idx],
    artificial_provenance = .make_artificial_provenance(
      requested_type = "knockoff",
      selected_types = selected_types,
      chunks = chunks,
      fallback_histories = fallback_histories
    )
  )
}

# ── Shared helpers for model-X generators ─────────────────────────────────────

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
#' comparisons.  Column-chunking for datasets that exceed 3 000 features is
#' applied (same as [make_knockoff_features()]).  Falls back to random-permutation
#' features when the knockoff constructor fails.
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details.
#' @keywords internal
make_knockoff_equi_features <- function(x, n_injected, random_state = NULL) {
  .require_pkg("knockoff", "for artificial_type = \"knockoff_equi\"")

  n_features <- ncol(x)
  chunk_size <- 3000L

  .make_equi_chunk <- function(x_chunk) {
    tryCatch(
      {
        mu    <- colMeans(x_chunk)
        Sigma <- .estimate_pd_sigma(x_chunk)
        list(
          x_art = knockoff::create.gaussian(x_chunk, mu, Sigma, method = "equi"),
          actual_type = "knockoff_equi",
          fallback_reason = NA_character_,
          fallback_history = .empty_artificial_fallback_history()
        )
      },
      error = function(e) {
        reason <- conditionMessage(e)
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

  if (n_features > chunk_size) {
    n_chunks        <- ceiling(n_features / chunk_size)
    ko_blocks       <- vector("list", n_chunks)
    orig_map_blocks <- vector("list", n_chunks)
    type_blocks     <- vector("list", n_chunks)
    chunk_records   <- vector("list", n_chunks)
    fallback_histories <- vector("list", n_chunks)
    for (i in seq_len(n_chunks)) {
      col_idx              <- sample.int(n_features,
                                         size = min(chunk_size, n_features),
                                         replace = FALSE)
      chunk_result         <- .make_equi_chunk(x[, col_idx, drop = FALSE])
      ko_blocks[[i]]       <- chunk_result$x_art
      orig_map_blocks[[i]] <- col_idx
      type_blocks[[i]]     <- rep.int(chunk_result$actual_type, ncol(chunk_result$x_art))
      fallback_histories[[i]] <- chunk_result$fallback_history
      chunk_records[[i]]   <- .artificial_chunk_record(
        chunk = i,
        requested_type = "knockoff_equi",
        actual_type = chunk_result$actual_type,
        n_columns = ncol(chunk_result$x_art),
        fallback_reason = chunk_result$fallback_reason
      )
    }
    x_art_full <- do.call(cbind, ko_blocks)
    orig_map   <- unlist(orig_map_blocks)
    type_map   <- unlist(type_blocks, use.names = FALSE)
    chunks     <- do.call(rbind, chunk_records)
    keep_idx   <- sample.int(ncol(x_art_full), size = n_features, replace = FALSE)
    x_art_full <- x_art_full[, keep_idx, drop = FALSE]
    orig_map   <- orig_map[keep_idx]
    type_map   <- type_map[keep_idx]
  } else {
    chunk_result <- .make_equi_chunk(x)
    x_art_full <- chunk_result$x_art
    orig_map   <- seq_len(n_features)
    type_map   <- rep.int(chunk_result$actual_type, n_features)
    chunks     <- .artificial_chunk_record(
      chunk = 1L,
      requested_type = "knockoff_equi",
      actual_type = chunk_result$actual_type,
      n_columns = ncol(chunk_result$x_art),
      fallback_reason = chunk_result$fallback_reason
    )
    fallback_histories <- list(chunk_result$fallback_history)
  }

  sel_idx <- sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE)
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
#' does not implement MVR knockoffs. High-dimensional chunking is approximate:
#' global exchangeability across chunks has not been established. Chunking and
#' fallback behaviour mirror [make_knockoff_equi_features()]: MVR-solver failure falls back to equi;
#' `create.gaussian` failure falls back to random permutation.
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed; passed to `solve_mvr()` for the
#'   coordinate-shuffle RNG.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details. For MVR,
#'   `artificial_provenance$mvr_chunking` records whether the approximate
#'   high-dimensional chunking path was applied and explicitly reports that
#'   global exchangeability across chunks has not been established.
#' @keywords internal
make_knockoff_mvr_features <- function(x, n_injected, random_state = NULL) {
  .require_pkg("knockoff", "for artificial_type = \"knockoff_mvr\"")

  n_features <- ncol(x)
  chunk_size <- 3000L

  .make_mvr_chunk <- function(x_chunk) {
    fallback_history <- .empty_artificial_fallback_history()
    tryCatch(
      {
        mu    <- colMeans(x_chunk)
        Sigma <- .estimate_pd_sigma(x_chunk)
        equi_fallback_reason <- NA_character_

        # Attempt MVR S-solve; fall back to equi on solver failure
        S_diag <- tryCatch(
          solve_mvr(Sigma, random_state = random_state),
          error   = function(e) {
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
            x_art = knockoff::create.gaussian(x_chunk, mu, Sigma, method = "equi"),
            actual_type = "knockoff_equi",
            fallback_reason = equi_fallback_reason,
            fallback_history = fallback_history
          )
        } else {
          list(
            x_art = knockoff::create.gaussian(x_chunk, mu, Sigma, diag_s = S_diag),
            actual_type = "knockoff_mvr",
            fallback_reason = NA_character_,
            fallback_history = fallback_history
          )
        }
      },
      error = function(e) {
        reason <- conditionMessage(e)
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

  if (n_features > chunk_size) {
    n_chunks        <- ceiling(n_features / chunk_size)
    ko_blocks       <- vector("list", n_chunks)
    orig_map_blocks <- vector("list", n_chunks)
    type_blocks     <- vector("list", n_chunks)
    chunk_records   <- vector("list", n_chunks)
    fallback_histories <- vector("list", n_chunks)
    for (i in seq_len(n_chunks)) {
      col_idx              <- sample.int(n_features,
                                         size = min(chunk_size, n_features),
                                         replace = FALSE)
      chunk_result         <- .make_mvr_chunk(x[, col_idx, drop = FALSE])
      ko_blocks[[i]]       <- chunk_result$x_art
      orig_map_blocks[[i]] <- col_idx
      type_blocks[[i]]     <- rep.int(chunk_result$actual_type, ncol(chunk_result$x_art))
      fallback_histories[[i]] <- chunk_result$fallback_history
      chunk_records[[i]]   <- .artificial_chunk_record(
        chunk = i,
        requested_type = "knockoff_mvr",
        actual_type = chunk_result$actual_type,
        n_columns = ncol(chunk_result$x_art),
        fallback_reason = chunk_result$fallback_reason
      )
    }
    x_art_full <- do.call(cbind, ko_blocks)
    orig_map   <- unlist(orig_map_blocks)
    type_map   <- unlist(type_blocks, use.names = FALSE)
    chunks     <- do.call(rbind, chunk_records)
    keep_idx   <- sample.int(ncol(x_art_full), size = n_features, replace = FALSE)
    x_art_full <- x_art_full[, keep_idx, drop = FALSE]
    orig_map   <- orig_map[keep_idx]
    type_map   <- type_map[keep_idx]
  } else {
    chunk_result <- .make_mvr_chunk(x)
    x_art_full <- chunk_result$x_art
    orig_map   <- seq_len(n_features)
    type_map   <- rep.int(chunk_result$actual_type, n_features)
    chunks     <- .artificial_chunk_record(
      chunk = 1L,
      requested_type = "knockoff_mvr",
      actual_type = chunk_result$actual_type,
      n_columns = ncol(chunk_result$x_art),
      fallback_reason = chunk_result$fallback_reason
    )
    fallback_histories <- list(chunk_result$fallback_history)
  }

  sel_idx <- sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE)
  x_art   <- x_art_full[, sel_idx, drop = FALSE]
  selected_types <- type_map[sel_idx]

  provenance <- .make_artificial_provenance(
    requested_type = "knockoff_mvr",
    selected_types = selected_types,
    chunks = chunks,
    fallback_histories = fallback_histories
  )
  is_chunked <- n_features > chunk_size
  provenance$mvr_chunking <- list(
    applied = is_chunked,
    chunk_size = chunk_size,
    approximate = is_chunked,
    global_exchangeability_established = if (is_chunked) FALSE else NA
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
#'     [knockoff::create.fixed()], which preserve the covariance structure of
#'     the original features under the fixed-design assumption.  Kept for
#'     backward compatibility.}
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
#'   knockoffs, a fallback to random permutation is attempted when the
#'   knockoff constructor fails (e.g., rank-deficient input). Fixed-X
#'   knockoffs that require row augmentation also fall back because augmented
#'   knockoff rows cannot be aligned with the original outcome.
#' @param n_injected Positive integer; number of artificial columns to append.
#'   Typically `round(ncol(x) * artificial_proportion)` as computed in
#'   [stabl_fit()].
#' @param type Character string; one of `"random_permutation"`, `"knockoff"`,
#'   `"knockoff_equi"`, or `"knockoff_mvr"`.
#' @param random_state Optional integer; passed to [set.seed()] before any
#'   random operations for reproducibility.  `NULL` leaves the RNG unchanged.
#'   Seeding happens exactly once in this dispatcher; downstream generators
#'   inherit the seeded RNG state and do not re-seed (audit M-5).
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
#' @seealso [compute_fdp_plus()] which consumes the artificial-feature scores,
#'   [stabl_fit()] which calls this function automatically.
#' @export
make_artificial_features <- function(x, n_injected, type, random_state = NULL) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix.", call. = FALSE)
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
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }
  if (!is.null(random_state)) set.seed(random_state)
  switch(
    type,
    random_permutation = make_rp_features(x, n_injected),
    knockoff           = make_knockoff_features(x, n_injected,
                                                random_state = random_state),
    knockoff_equi      = make_knockoff_equi_features(x, n_injected,
                                                     random_state = random_state),
    knockoff_mvr       = make_knockoff_mvr_features(x, n_injected,
                                                    random_state = random_state),
    stop(
      "`type` must be one of \"random_permutation\", \"knockoff\", ",
      "\"knockoff_equi\", or \"knockoff_mvr\", got: ", type,
      call. = FALSE
    )
  )
}
