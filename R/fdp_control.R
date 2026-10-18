#' Compute FDP+ (False Discovery Proportion Upper Bound)
#'
#' Sweeps a grid of stability-score thresholds and computes the FDP+ estimate
#' at each one, then identifies the threshold that minimises it.  This is a
#' direct port of `Stabl._compute_FDPplus()` from the Python STABL library.
#'
#' The FDP+ at threshold \eqn{\tau} is estimated as:
#' \deqn{\widehat{\text{FDP}}(\tau) =
#'   \frac{(1/\pi) \cdot |\{j : \hat{q}_j^{\text{art}} > \tau\}| + 1}
#'        {\max(1,\, |\{j : \hat{q}_j > \tau\}|)}}
#' where \eqn{\pi} is `artificial_proportion`, \eqn{\hat{q}_j} are the
#' maximum-over-lambda stability scores for real feature \eqn{j}, and
#' \eqn{\hat{q}_j^{\text{art}}} are the scores for artificial feature \eqn{j}.
#'
#' @section Comparator and tie policy (contract):
#' Exceedance counts use a \strong{strict} `>` comparator: a feature whose
#' score equals the threshold exactly is \emph{not} counted.  This matches
#' the Python STABL reference implementation and is part of the parity
#' contract; do not "fix" it to `>=`.  When several thresholds tie for the
#' minimum FDP+, the \strong{first} (lowest) threshold on the grid wins,
#' matching Python's `argmin` first-index semantics.  `min_fdr` is a
#' decoy-based diagnostic (a stochastic upper bound under exchangeability of
#' the artificial features), not a controlled FDR guarantee; when it exceeds
#' 1, `fdr_min_threshold` is capped at 1 and downstream selection yields an
#' empty support.
#'
#' @param stabl_scores Numeric matrix (features \eqn{\times} lambdas) of
#'   stability scores for the real features.
#' @param stabl_scores_artificial Numeric matrix (artificial features
#'   \eqn{\times} lambdas) of stability scores for the injected noise features.
#' @param artificial_proportion Positive numeric scalar; the fraction of
#'   artificial features relative to real features (used as \eqn{\pi}).
#' @param fdr_threshold_range Numeric vector of threshold values to sweep.
#'   Default: `seq(0, 0.99, by = 0.01)`, matching Python STABL's
#'   `np.arange(0., 1., .01)`.
#'
#' @return Named list:
#'   \describe{
#'     \item{FDRs}{Numeric vector of FDP+ estimates, one per threshold.}
#'     \item{min_fdr}{Minimum FDP+ achieved across the threshold range.}
#'     \item{fdr_min_threshold}{Threshold value achieving `min_fdr`, capped
#'       at 1.}
#'     \item{fdrs_table}{Numeric matrix (lambdas \eqn{\times} thresholds) of
#'       per-lambda FDP+ values.}
#'   }
#'
#' @examples
#' # Small synthetic stability-score matrices (3 real + 3 artificial features,
#' # 2 lambda values)
#' scores_real <- matrix(c(0.8, 0.6, 0.1, 0.9, 0.5, 0.05), nrow = 3)
#' scores_art  <- matrix(c(0.1, 0.2, 0.3, 0.15, 0.25, 0.1), nrow = 3)
#' result <- compute_fdp_plus(scores_real, scores_art,
#'                            artificial_proportion = 1.0)
#' result$min_fdr           # minimum FDP+ estimate
#' result$fdr_min_threshold # threshold that achieves it
#'
#' @seealso [stabl_fit()] which calls this internally when artificial features
#'   are used, [plot_fdr_graph()] to visualise the FDP+ curve.
#' @export
compute_fdp_plus <- function(
    stabl_scores,
    stabl_scores_artificial,
    artificial_proportion,
    fdr_threshold_range = seq(0, 0.99, by = 0.01)
) {
  score_inputs <- list(
    stabl_scores = stabl_scores,
    stabl_scores_artificial = stabl_scores_artificial
  )
  for (name in names(score_inputs)) {
    scores <- score_inputs[[name]]
    if (!is.matrix(scores) || !is.numeric(scores)) {
      stop("`", name, "` must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(scores) == 0L || ncol(scores) == 0L) {
      stop("`", name, "` must have at least one row and one column.", call. = FALSE)
    }
    if (any(!is.finite(scores))) {
      stop("`", name, "` must contain only finite values.", call. = FALSE)
    }
    if (any(scores < 0 | scores > 1)) {
      stop("`", name, "` values must lie in [0, 1].", call. = FALSE)
    }
  }
  if (ncol(stabl_scores) != ncol(stabl_scores_artificial)) {
    stop(
      "`stabl_scores` and `stabl_scores_artificial` must have the same number of columns.",
      call. = FALSE
    )
  }
  artificial_proportion <- .validate_proportion(
    artificial_proportion,
    "artificial_proportion"
  )
  if (!is.numeric(fdr_threshold_range) || !is.null(dim(fdr_threshold_range)) ||
      length(fdr_threshold_range) == 0L ||
      any(!is.finite(fdr_threshold_range)) ||
      any(fdr_threshold_range < 0 | fdr_threshold_range > 1)) {
    stop(
      "`fdr_threshold_range` must be a non-empty finite numeric vector in [0, 1].",
      call. = FALSE
    )
  }

  inv_prop    <- 1.0 / artificial_proportion
  max_scores  <- rowMaxs(stabl_scores)
  max_art     <- rowMaxs(stabl_scores_artificial)

  # Sort once, then count values strictly above each threshold using binary
  # interval lookup. This preserves strict `>` ties while avoiding the large
  # feature-by-threshold logical matrices produced by outer().
  n_real <- .count_strict_exceedances(max_scores, fdr_threshold_range)
  n_art  <- .count_strict_exceedances(max_art, fdr_threshold_range)
  FDRs   <- (inv_prop * n_art + 1.0) / pmax(1.0, n_real)

  # Per-lambda table uses the same strict counting implementation column-wise.
  nr_table <- .count_matrix_strict_exceedances(
    stabl_scores,
    fdr_threshold_range
  )
  na_table <- .count_matrix_strict_exceedances(
    stabl_scores_artificial,
    fdr_threshold_range
  )
  fdrs_table <- (inv_prop * na_table + 1.0) / pmax(1.0, nr_table)

  min_fdr <- min(FDRs)
  fdr_min_threshold <- if (min_fdr > 1.0) {
    1.0
  } else {
    min(fdr_threshold_range[which.min(FDRs)], 1.0)  # ties: first-index wins (matches Python argmin)
  }

  list(
    FDRs              = FDRs,
    min_fdr           = min_fdr,
    fdr_min_threshold = fdr_min_threshold,
    fdrs_table        = fdrs_table
  )
}

.count_strict_exceedances <- function(values, thresholds) {
  sorted <- sort(unname(values), method = "quick")
  counts <- as.numeric(length(sorted) - findInterval(thresholds, sorted))
  names(counts) <- names(thresholds)
  counts
}

.count_matrix_strict_exceedances <- function(scores, thresholds) {
  counts <- vapply(
    seq_len(ncol(scores)),
    function(j) .count_strict_exceedances(scores[, j], thresholds),
    numeric(length(thresholds))
  )
  counts <- matrix(
    counts,
    nrow = length(thresholds),
    ncol = ncol(scores)
  )
  unname(t(counts))
}
