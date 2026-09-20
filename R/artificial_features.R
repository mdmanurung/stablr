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

#' Make Knockoff Artificial Features
#'
#' Generates **fixed-X** knockoff features via `knockoff::create.fixed()`, with
#' up to 3 000 predictors. Larger inputs fail closed because the former
#' chunked construction did not establish a valid global knockoff design.
#' Falls back to random-permutation features when the supported-dimensional
#' knockoff constructor fails (e.g., rank-deficient input).
#'
#' @param x Numeric matrix of predictors (samples \eqn{\times} features).
#' @param n_injected Integer; number of knockoff columns to select.
#' @param random_state Optional integer seed.
#'
#' @return Named list with elements `x_augmented`, `noise_col_indices`, and
#'   `artificial_provenance`; see [make_rp_features()] for details.
#' @keywords internal
make_knockoff_features <- function(x, n_injected, random_state = NULL) {
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
        xk <- .with_linear_algebra_numerical_infeasibility(
          withCallingHandlers(
            knockoff::create.fixed(x_chunk, sigma = 1)$Xk,
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
          fallback_history = .empty_artificial_fallback_history()
        )
      },
      stablr_numerical_infeasibility = function(e) {
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

# Estimate a positive-definite covariance matrix for the MVR generator.
# This historical helper is intentionally retained unchanged: v0.1.1 changes
# only the covariance semantics of knockoff_equi.
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

# Reproduce sklearn.covariance.LedoitWolf 1.7.1 with assume_centered = FALSE.
# In particular, sklearn uses the maximum-likelihood (1 / n) empirical
# covariance convention rather than stats::cov()'s (1 / (n - 1)) convention.
.estimate_ledoit_wolf_sigma <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L ||
      anyNA(x) || any(!is.finite(x))) {
    stop("`x` must be a non-empty finite numeric matrix.", call. = FALSE)
  }

  n_samples <- nrow(x)
  n_features <- ncol(x)
  centered <- sweep(x, 2L, colMeans(x), "-")
  empirical_covariance <- crossprod(centered) / n_samples

  if (n_features == 1L) {
    shrinkage <- 0
  } else {
    squared <- centered^2
    empirical_trace <- colSums(squared) / n_samples
    mu <- sum(empirical_trace) / n_features
    delta_raw <- sum(crossprod(centered)^2) / n_samples^2
    beta_raw <- sum(crossprod(squared))
    beta <- (beta_raw / n_samples - delta_raw) /
      (n_features * n_samples)
    delta <- (
      delta_raw - 2 * mu * sum(empirical_trace) + n_features * mu^2
    ) / n_features
    beta <- min(beta, delta)
    shrinkage <- if (beta == 0) 0 else beta / delta
  }

  mu <- sum(diag(empirical_covariance)) / n_features
  Sigma <- (1 - shrinkage) * empirical_covariance
  diag(Sigma) <- diag(Sigma) + shrinkage * mu
  Sigma <- (Sigma + t(Sigma)) / 2

  eigenvalues <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(Sigma)) || any(!is.finite(eigenvalues)) ||
      min(eigenvalues) <= 0) {
    .abort_numerical_infeasibility(
      "stablr_knockoff_infeasible",
      paste0(
        "Ledoit-Wolf covariance estimation did not produce a finite ",
        "positive-definite matrix."
      ),
      source = "sklearn-compatible Ledoit-Wolf covariance",
      n_samples = n_samples,
      n_features = n_features,
      minimum_eigenvalue = min(eigenvalues)
    )
  }

  list(
    covariance = Sigma,
    shrinkage = unname(shrinkage),
    location = unname(colMeans(x)),
    minimum_eigenvalue = unname(min(eigenvalues))
  )
}

# Match knockpy 1.3.5's ungrouped equicorrelated S-matrix construction. It
# solves on the correlation scale, applies the package's conservative
# ten-iteration PSD search at tolerance 1e-5, then rescales to covariance.
.knockpy_equicorrelated_s <- function(Sigma, tol = 1e-5, num_iter = 10L) {
  if (!is.matrix(Sigma) || !is.numeric(Sigma) || nrow(Sigma) != ncol(Sigma) ||
      !nrow(Sigma) || anyNA(Sigma) || any(!is.finite(Sigma))) {
    stop("`Sigma` must be a finite numeric square matrix.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be one finite positive number.", call. = FALSE)
  }
  num_iter <- .validate_scalar_integer_like(num_iter, "num_iter", min = 1L)

  variances <- diag(Sigma)
  if (any(variances <= 0)) {
    .abort_numerical_infeasibility(
      "stablr_knockoff_infeasible",
      "Equicorrelated knockoffs require strictly positive variances.",
      source = "knockpy-compatible equicorrelated S-matrix"
    )
  }
  scales <- sqrt(variances)
  correlation <- Sigma / outer(scales, scales)
  correlation <- (correlation + t(correlation)) / 2
  minimum_eigenvalue <- min(
    eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  )
  initial_gamma <- min(2 * minimum_eigenvalue, 1)
  if (!is.finite(initial_gamma) || initial_gamma <= 0) {
    .abort_numerical_infeasibility(
      "stablr_knockoff_infeasible",
      "The equicorrelated S-matrix is not positive definite.",
      source = "knockpy-compatible equicorrelated S-matrix",
      minimum_correlation_eigenvalue = minimum_eigenvalue
    )
  }

  initial_s <- diag(initial_gamma, nrow(Sigma))
  lower <- 0
  upper <- 1
  identity <- diag(nrow(Sigma))
  for (iteration in seq_len(num_iter)) {
    scale <- (lower + upper) / 2
    feasible <- tryCatch(
      {
        chol(2 * correlation - scale * initial_s - tol * identity)
        TRUE
      },
      error = function(...) FALSE
    )
    if (feasible) lower <- scale else upper <- scale
  }

  diag_s <- lower * initial_gamma * variances
  conditional_minimum_eigenvalue <- min(eigen(
    2 * Sigma - diag(diag_s, nrow(Sigma)),
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  if (any(!is.finite(diag_s)) || any(diag_s <= 0) ||
      !is.finite(conditional_minimum_eigenvalue) ||
      conditional_minimum_eigenvalue < tol) {
    .abort_numerical_infeasibility(
      "stablr_knockoff_infeasible",
      "The equicorrelated S-matrix failed numerical validation.",
      source = "knockpy-compatible equicorrelated S-matrix",
      minimum_eigenvalue = conditional_minimum_eigenvalue
    )
  }

  list(
    diag_s = unname(diag_s),
    initial_gamma = unname(initial_gamma),
    psd_scale = unname(lower),
    minimum_correlation_eigenvalue = unname(minimum_eigenvalue),
    conditional_minimum_eigenvalue = unname(conditional_minimum_eigenvalue),
    tolerance = unname(tol),
    iterations = as.integer(num_iter)
  )
}

#' Make Model-X Equicorrelated Knockoff Artificial Features
#'
#' Generates **model-X equicorrelated** knockoff features using a
#' scikit-learn 1.7.1-compatible Ledoit-Wolf covariance estimate and a knockpy
#' 1.3.5-compatible equicorrelated S-matrix. Sampling uses
#' `knockoff::create.gaussian(..., diag_s = S)`. This matches the
#' `GaussianSampler(X, method='equicorrelated')` call used by Python STABL.
#' Inputs above 3 000 predictors fail closed because a globally
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
make_knockoff_equi_features <- function(x, n_injected, random_state = NULL) {
  .assert_supported_knockoff_dimension(x, "knockoff_equi")
  .require_pkg("knockoff", "for artificial_type = \"knockoff_equi\"")

  n_features <- ncol(x)

  .make_equi_chunk <- function(x_chunk) {
    covariance_details <- list(
      name = "sklearn.covariance.LedoitWolf",
      compatibility_version = "1.7.1",
      assume_centered = FALSE,
      empirical_covariance_divisor = "n_samples",
      shrinkage = NA_real_
    )
    s_matrix_details <- list(
      name = "knockpy.smatrix.equicorrelated",
      compatibility_version = "1.3.5",
      tolerance = 1e-5,
      iterations = 10L
    )
    tryCatch(
      {
        mu    <- colMeans(x_chunk)
        covariance_fit <- .with_linear_algebra_numerical_infeasibility(
          .estimate_ledoit_wolf_sigma(x_chunk),
          source = "sklearn-compatible Ledoit-Wolf covariance estimation",
          subclass = "stablr_knockoff_infeasible"
        )
        Sigma <- covariance_fit$covariance
        covariance_details$shrinkage <- covariance_fit$shrinkage
        s_fit <- .with_linear_algebra_numerical_infeasibility(
          .knockpy_equicorrelated_s(Sigma),
          source = "knockpy-compatible equicorrelated S-matrix",
          subclass = "stablr_knockoff_infeasible"
        )
        s_matrix_details$initial_gamma <- s_fit$initial_gamma
        s_matrix_details$psd_scale <- s_fit$psd_scale
        s_matrix_details$minimum_correlation_eigenvalue <-
          s_fit$minimum_correlation_eigenvalue
        s_matrix_details$conditional_minimum_eigenvalue <-
          s_fit$conditional_minimum_eigenvalue
        list(
          x_art = .with_linear_algebra_numerical_infeasibility(
            knockoff::create.gaussian(
              x_chunk,
              mu,
              Sigma,
              diag_s = s_fit$diag_s
            ),
            source = "knockoff::create.gaussian(equi)",
            subclass = "stablr_knockoff_infeasible"
          ),
          actual_type = "knockoff_equi",
          fallback_reason = NA_character_,
          fallback_history = .empty_artificial_fallback_history(),
          covariance_estimator = covariance_details,
          equicorrelated_s_matrix = s_matrix_details
        )
      },
      stablr_numerical_infeasibility = function(e) {
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
          ),
          covariance_estimator = covariance_details,
          equicorrelated_s_matrix = s_matrix_details
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

  sel_idx <- sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE)
  x_art   <- x_art_full[, sel_idx, drop = FALSE]
  selected_types <- type_map[sel_idx]

  provenance <- .make_artificial_provenance(
    requested_type = "knockoff_equi",
    selected_types = selected_types,
    chunks = chunks,
    fallback_histories = fallback_histories
  )
  provenance$covariance_estimator <- chunk_result$covariance_estimator
  provenance$equicorrelated_s_matrix <- chunk_result$equicorrelated_s_matrix

  list(
    x_augmented      = cbind(x, x_art),
    noise_col_indices = orig_map[sel_idx],
    artificial_provenance = provenance
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
make_knockoff_mvr_features <- function(x, n_injected, random_state = NULL) {
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

  sel_idx <- sample.int(n = ncol(x_art_full), size = n_injected, replace = FALSE)
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
#'     [knockoff::create.fixed()], which preserve the covariance structure of
#'     the original features under the fixed-design assumption.  Kept for
#'     backward compatibility.}
#'   \item{`"knockoff_equi"`}{Generates **model-X equicorrelated** knockoffs
#'     using a scikit-learn 1.7.1-compatible Ledoit-Wolf covariance estimate
#'     and knockpy 1.3.5-compatible equicorrelated S-matrix. Sampling uses
#'     `knockoff::create.gaussian(..., diag_s = S)`. This matches the
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
#'   knockoff rows cannot be aligned with the original outcome. All knockoff
#'   modes fail closed above 3 000 predictors; random permutation has no such
#'   dimensional limit.
#' @param n_injected Positive integer; number of artificial columns to append.
#'   Typically `round(ncol(x) * artificial_proportion)` as computed in
#'   [stabl_fit()].
#' @param type Character string; one of `"random_permutation"`, `"knockoff"`,
#'   `"knockoff_equi"`, or `"knockoff_mvr"`.
#' @param random_state Optional integer used as a locally scoped seed for
#'   reproducibility. The caller's RNG state is restored on exit. `NULL`
#'   leaves generation on the caller's current RNG stream. Seeding happens
#'   exactly once in this dispatcher; downstream generators do not re-seed.
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
make_artificial_features <- function(x, n_injected, type, random_state = NULL) {
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
  generate <- function() {
    switch(
      type,
      random_permutation = make_rp_features(x, n_injected),
      knockoff           = make_knockoff_features(x, n_injected,
                                                  random_state = random_state),
      knockoff_equi      = make_knockoff_equi_features(
        x,
        n_injected,
        random_state = random_state
      ),
      knockoff_mvr       = make_knockoff_mvr_features(
        x,
        n_injected,
        random_state = random_state
      )
    )
  }
  if (is.null(random_state)) return(generate())

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(random_state)
  generate()
}
