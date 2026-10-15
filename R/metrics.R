# metrics.R — Selection-stability similarity measures
# R port of stabl/metrics.py
#
# These functions compare lists of selected features across repeat runs of
# STABL (or any other feature-selection method) to assess reproducibility.
# Three families of similarity are provided:
#   - Jaccard: set-overlap fraction, no correction for chance.
#   - Adjusted: chance-corrected overlap (analogous to Cohen's kappa for sets).
#   - Pearson: correlation-inspired chance correction, normalized by set size.
# For simulation benchmarks where a ground-truth feature set is known, the
# file also provides FDR, TPR, and F-score functions.

# ── Jaccard ──────────────────────────────────────────────────────────────────

#' Jaccard Similarity Between Two Feature Sets
#'
#' Measures the overlap between two sets of selected features as a simple
#' fraction of shared features out of all features that appeared in either
#' set.  This is the most interpretable pairwise similarity when no prior
#' expectation of set sizes exists.
#'
#' Use Jaccard when you want a quick, easy-to-explain measure.  Prefer
#' [adjusted_similarity()] or [pearson_similarity()] when comparing across
#' runs that select very different numbers of features, because Jaccard does
#' not correct for the chance overlap expected between small or large sets.
#'
#' When both sets are empty the function returns 0 (convention matching the
#' Python implementation).
#'
#' @param list1 Character or integer vector of selected feature identifiers.
#' @param list2 Character or integer vector of selected feature identifiers.
#'
#' @return Numeric scalar in \eqn{[0, 1]}.  A value of 1 means the two sets
#'   are identical; 0 means they share no features.
#'
#' @seealso [jaccard_similarity_matrix()] for computing all pairwise Jaccard similarities
#'   at once.
#'
#' @examples
#' jaccard_similarity(c("A", "B", "C"), c("B", "C", "D"))  # 0.5
#' jaccard_similarity(c("A", "B"),      c("A", "B"))        # 1
#' jaccard_similarity(character(0),     character(0))       # 0 (both empty)
#' @export
jaccard_similarity <- function(list1, list2) {
  s1 <- unique(list1)
  s2 <- unique(list2)
  inter <- length(intersect(s1, s2))
  union_size <- length(s1) + length(s2) - inter
  if (inter == 0L && union_size == 0L) return(0)
  inter / union_size
}

#' Pairwise Jaccard Matrix from a List of Feature Sets
#'
#' Computes all N\eqn{\times}N pairwise Jaccard similarities in one call.
#' Useful for visualising the reproducibility of STABL across multiple
#' cross-validation folds or repeated runs: a boxplot of the off-diagonal
#' values summarises how consistently the same features are selected.
#'
#' By default the self-similarity diagonal (always 1) is removed, so the
#' returned matrix has \eqn{N \times (N - 1)} columns and each row contains
#' the \eqn{N - 1} similarities of run \eqn{i} with every other run.
#'
#' @param list_of_lists A list of character/integer vectors, one per STABL run
#'   (e.g. one per cross-validation fold).
#' @param remove_diag Logical; if `TRUE` (default) the diagonal column
#'   (self-similarity = 1) is removed from the output.
#'
#' @return Numeric matrix of dimension N\eqn{\times}N (or N\eqn{\times}(N-1)
#'   when `remove_diag = TRUE`).  Row/column order matches `list_of_lists`.
#'
#' @seealso [jaccard_similarity()], [adjusted_similarity_pairwise()],
#'   [pearson_similarity_pairwise()]
#'
#' @examples
#' sets <- list(c("A","B","C"), c("B","C","D"), c("A","C","E"))
#' jaccard_similarity_matrix(sets)           # 3 x 2 (diagonal removed)
#' jaccard_similarity_matrix(sets, remove_diag = FALSE)  # 3 x 3
#' @export
jaccard_similarity_matrix <- function(list_of_lists, remove_diag = TRUE) {
  n   <- length(list_of_lists)
  mat <- matrix(0, nrow = n, ncol = n)
  diag(mat) <- 1
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      v         <- jaccard_similarity(list_of_lists[[i]], list_of_lists[[j]])
      mat[i, j] <- v
      mat[j, i] <- v
    }
  }
  if (remove_diag) {
    out <- matrix(0, nrow = n, ncol = n - 1L)
    for (i in seq_len(n)) {
      out[i, ] <- mat[i, setdiff(seq_len(n), i)]
    }
    mat <- out
  }
  mat
}

# ── Adjusted similarity ───────────────────────────────────────────────────────

#' Adjusted Similarity Between Two Feature Sets
#'
#' Computes a chance-corrected similarity between two feature sets.  Unlike
#' [jaccard_similarity()], the adjusted measure accounts for the expected
#' random overlap given the sizes of both sets and the total feature universe,
#' so it does not systematically penalise methods that select many features.
#'
#' The formula is analogous to Cohen's kappa for sets:
#' \deqn{S_{\text{adj}}(A, B) =
#'   \frac{r - \mathbb{E}[r]}{\min(k_1, k_2) - \max(0, k_1 + k_2 - d)}}
#' where \eqn{r = |A \cap B|}, \eqn{k_i = |A_i|}, \eqn{d =} `n_universe`,
#' and \eqn{\mathbb{E}[r] = k_1 k_2 / d}.
#'
#' Returns 0 when either set is empty or equals the full universe (edge cases
#' where the correction denominator is zero).
#'
#' @param set1 Character or integer vector of selected feature identifiers.
#' @param set2 Character or integer vector of selected feature identifiers.
#' @param n_universe Integer; total number of candidate features in the
#'   universe (i.e. the number of columns in the original predictor matrix).
#' @param list1 Deprecated; use `set1` instead.  If supplied by name, a
#'   deprecation warning is emitted and the value is used as `set1` when the
#'   latter is `NULL`.  It will be removed in a future release.
#' @param list2 Deprecated; use `set2` instead.  If supplied by name, a
#'   deprecation warning is emitted and the value is used as `set2` when the
#'   latter is `NULL`.  It will be removed in a future release.
#' @param nb_total_elements Deprecated; use `n_universe` instead.  If supplied
#'   by name, a deprecation warning is emitted and the value is used as
#'   `n_universe` when the latter is `NULL`.  It will be removed in a future
#'   release.
#'
#' @return Numeric scalar in \eqn{(-1, 1]}.  Values above 0 indicate more
#'   overlap than expected by chance; 1 means perfect agreement; negative
#'   values indicate less overlap than chance.
#'
#' @seealso [adjusted_similarity_pairwise()] for computing all pairwise values at
#'   once, [adjusted_similarity_summary()] for a summary statistic,
#'   [jaccard_similarity()] for the non-chance-corrected alternative.
#'
#' @examples
#' # 10-feature universe; sets A and B share 2 out of 3 features each
#' adjusted_similarity(c("f1","f2","f3"), c("f2","f3","f4"), n_universe = 10L)
#' adjusted_similarity(c("f1","f2"),      c("f1","f2"),      n_universe = 10L) # 0.8
#' adjusted_similarity(character(0),      c("f1"),           n_universe = 10L) # 0
#' @export
adjusted_similarity <- function(set1 = NULL, set2 = NULL, n_universe = NULL,
                                list1 = NULL, list2 = NULL,
                                nb_total_elements = NULL) {
  # Deprecation shim: `list1`/`list2`/`nb_total_elements` were renamed to
  # `set1`/`set2`/`n_universe`.  missing() (not is.null()) so that explicit
  # NULLs also warn; purely positional calls are unaffected.
  if (!missing(list1) || !missing(list2) || !missing(nb_total_elements)) {
    warning(
      "`list1`, `list2` and `nb_total_elements` are deprecated; use ",
      "`set1`, `set2` and `n_universe` instead.",
      call. = FALSE
    )
    if (is.null(set1)) set1 <- list1
    if (is.null(set2)) set2 <- list2
    if (is.null(n_universe)) n_universe <- nb_total_elements
  }
  if (is.null(n_universe)) {
    stop("argument \"n_universe\" is missing, with no default", call. = FALSE)
  }

  s1 <- unique(set1)
  s2 <- unique(set2)
  r  <- length(intersect(s1, s2))
  k1 <- length(s1)
  k2 <- length(s2)
  u  <- k1 + k2 - r

  d <- as.integer(n_universe)
  if (u > d) {
    stop(
      "Union cardinal (", u, ") exceeds n_universe (", d, ").",
      call. = FALSE
    )
  }
  if (k1 == d || k2 == d || k1 == 0L || k2 == 0L) return(0)

  num   <- r - (k1 * k2) / d
  denom <- min(k1, k2) - max(0L, k1 + k2 - d)
  num / denom
}

#' Upper-Triangle Adjusted Similarity Values
#'
#' Computes the full pairwise adjusted-similarity matrix for N feature sets and
#' returns only the N\eqn{\times}(N-1)/2 off-diagonal upper-triangle values as
#' a flat vector.  This format is convenient for computing summary statistics
#' (see [adjusted_similarity_summary()]) or for Wilcoxon/permutation tests
#' comparing reproducibility across methods.
#'
#' @param list_of_lists A list of character/integer vectors, one per STABL run.
#' @param nb_total_elements Integer; total number of candidate features.
#'
#' @return Numeric vector of length N*(N-1)/2 containing the pairwise
#'   adjusted-similarity values for all unique pairs (row-major upper-triangle
#'   order, matching the Python convention).
#'
#' @seealso [adjusted_similarity()] for the pairwise function,
#'   [adjusted_similarity_summary()] for a one-number summary,
#'   [pearson_similarity_pairwise()] for an alternative chance-correction.
#'
#' @examples
#' sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
#' adjusted_similarity_pairwise(sets, nb_total_elements = 10L)
#' @export
adjusted_similarity_pairwise <- function(list_of_lists, nb_total_elements) {
  .pairwise_upper_triangle(
    list_of_lists,
    function(a, b) adjusted_similarity(a, b, n_universe = nb_total_elements)
  )
}

#' Summary Statistic of Adjusted Similarity Values
#'
#' Convenience wrapper that computes all pairwise adjusted similarities
#' (via [adjusted_similarity_pairwise()]) and reduces them to a single
#' location statistic with an associated spread measure.  Useful for
#' reporting a single reproducibility number per STABL configuration in
#' benchmarking tables.
#'
#' @param list_of_lists A list of character/integer vectors, one per run.
#' @param nb_total_elements Integer; total number of candidate features.
#' @param stat Character; `"median"` (default, robust to outliers) or
#'   `"mean"`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`statistic`}{The median (or mean) of adjusted-similarity values.}
#'     \item{`err`}{For `"median"`: the 25th and 75th percentile vector
#'       (IQR bounds).  For `"mean"`: the root-mean-squared deviation
#'       (RMSD / population SD).}
#'   }
#'
#' @seealso [adjusted_similarity_pairwise()], [pearson_similarity_summary()]
#'
#' @examples
#' sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
#' adjusted_similarity_summary(sets, nb_total_elements = 10L)
#' adjusted_similarity_summary(sets, nb_total_elements = 10L, stat = "mean")
#' @export
adjusted_similarity_summary <- function(list_of_lists, nb_total_elements,
                                        stat = "median") {
  vals <- adjusted_similarity_pairwise(list_of_lists, nb_total_elements)
  .similarity_summary(vals, stat)
}

# ── Pearson similarity ────────────────────────────────────────────────────────

#' Pearson-Corrected Similarity Between Two Feature Sets
#'
#' Computes a Pearson-correlation-inspired similarity that corrects for
#' expected random intersection.  This is a second chance-correction approach
#' (alongside [adjusted_similarity()]) that normalises by the geometric mean
#' of the within-set variances under independent Bernoulli sampling.
#'
#' The formula is:
#' \deqn{S_{\text{Pearson}}(A, B) =
#'   \frac{r - k_i k_j / d}{d \cdot \upsilon_i \upsilon_j}}
#' where \eqn{r = |A \cap B|}, \eqn{k_i = |A_i|}, \eqn{d} is the universe
#' size, and \eqn{\upsilon_i = \sqrt{\pi_i (1 - \pi_i)}} with
#' \eqn{\pi_i = k_i / d}.
#'
#' Edge cases: returns 1 when both sets are empty or both equal the universe;
#' returns 0 when one set is empty or equals the universe.
#'
#' @param set1 Character or integer vector of selected feature identifiers.
#' @param set2 Character or integer vector of selected feature identifiers.
#' @param n_universe Integer; total number of candidate features in the universe.
#' @param list_i Deprecated; use `set1` instead.  If supplied by name, a
#'   deprecation warning is emitted and the value is used as `set1` when the
#'   latter is `NULL`.  It will be removed in a future release.
#' @param list_j Deprecated; use `set2` instead.  If supplied by name, a
#'   deprecation warning is emitted and the value is used as `set2` when the
#'   latter is `NULL`.  It will be removed in a future release.
#' @param d Deprecated; use `n_universe` instead.  If supplied by name, a
#'   deprecation warning is emitted and the value is used as `n_universe` when
#'   the latter is `NULL`.  It will be removed in a future release.
#'
#' @return Numeric scalar.  Positive values indicate more overlap than chance;
#'   the maximum is typically close to 1 for perfectly matching sets.
#'
#' @seealso [pearson_similarity_pairwise()] for computing all pairwise values at
#'   once, [pearson_similarity_summary()] for a summary statistic,
#'   [adjusted_similarity()] for an alternative chance-correction.
#'
#' @examples
#' pearson_similarity(c("f1","f2","f3"), c("f2","f3","f4"), n_universe = 10L)
#' pearson_similarity(c("f1","f2"),      c("f1","f2"),      n_universe = 10L) # near 1
#' pearson_similarity(character(0),      c("f1","f2"),      n_universe = 10L) # 0
#' @export
pearson_similarity <- function(set1 = NULL, set2 = NULL, n_universe = NULL,
                               list_i = NULL, list_j = NULL, d = NULL) {
  # Deprecation shim: `list_i`/`list_j`/`d` were renamed to
  # `set1`/`set2`/`n_universe`.  missing() (not is.null()) so that explicit
  # NULLs also warn; purely positional calls are unaffected.
  if (!missing(list_i) || !missing(list_j) || !missing(d)) {
    warning(
      "`list_i`, `list_j` and `d` are deprecated; use ",
      "`set1`, `set2` and `n_universe` instead.",
      call. = FALSE
    )
    if (is.null(set1)) set1 <- list_i
    if (is.null(set2)) set2 <- list_j
    if (is.null(n_universe)) n_universe <- d
  }
  if (is.null(n_universe)) {
    stop("argument \"n_universe\" is missing, with no default", call. = FALSE)
  }

  si <- unique(set1)
  sj <- unique(set2)
  r  <- length(intersect(si, sj))
  ki <- length(si)
  kj <- length(sj)
  d  <- as.integer(n_universe)

  if ((ki == d && kj == d) || (ki == 0L && kj == 0L)) {
    return(1)
  }
  if (ki == 0L || kj == 0L || ki == d || kj == d) {
    return(0)
  }

  exp_r      <- (ki * kj) / d
  pi_i       <- ki / d
  pi_j       <- kj / d
  upsilon_i  <- sqrt(pi_i * (1 - pi_i))
  upsilon_j  <- sqrt(pi_j * (1 - pi_j))
  (r - exp_r) / (d * upsilon_i * upsilon_j)
}

#' Upper-Triangle Pearson Similarity Values
#'
#' Computes all pairwise Pearson-corrected similarities and returns the
#' N\eqn{\times}(N-1)/2 upper-triangle values as a flat vector (same layout
#' as [adjusted_similarity_pairwise()], enabling direct comparison between the
#' two metrics on the same data).
#'
#' @param list_of_lists A list of character/integer vectors, one per run.
#' @param nb_total_elements Integer; total number of candidate features.
#' @param d Deprecated; use `nb_total_elements` instead.
#'
#' @return Numeric vector of length N*(N-1)/2.
#'
#' @seealso [pearson_similarity()] for the pairwise function,
#'   [pearson_similarity_summary()] for a one-number summary,
#'   [adjusted_similarity_pairwise()] for an alternative chance-correction.
#'
#' @examples
#' sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
#' pearson_similarity_pairwise(sets, nb_total_elements = 10L)
#' @export
pearson_similarity_pairwise <- function(list_of_lists, nb_total_elements,
                                        d = NULL) {
  # Deprecation shim: `d` was renamed to `nb_total_elements` for consistency
  # with the adjusted-similarity family.  missing() (not is.null()) so that
  # explicit NULLs also warn; purely positional calls are unaffected.
  if (!missing(d) && !is.null(d)) {
    warning(
      "`d` is deprecated; use `nb_total_elements` instead.", call. = FALSE
    )
    if (is.null(nb_total_elements)) nb_total_elements <- d
  }
  .pairwise_upper_triangle(
    list_of_lists,
    function(a, b) pearson_similarity(a, b, n_universe = nb_total_elements)
  )
}

#' Summary Statistic of Pearson Similarity Values
#'
#' Convenience wrapper analogous to [adjusted_similarity_summary()] but
#' using the Pearson-corrected similarity.  See [pearson_similarity()] for a
#' description of the underlying metric and when to prefer it over the
#' adjusted similarity.
#'
#' @param list_of_lists A list of character/integer vectors, one per run.
#' @param nb_total_elements Integer; total number of candidate features.
#' @param stat Character; `"median"` (default) or `"mean"`.
#' @param d Deprecated; use `nb_total_elements` instead.
#'
#' @return A named list with `statistic` and `err` (see
#'   [adjusted_similarity_summary()] for the exact definitions).
#'
#' @seealso [pearson_similarity_pairwise()], [adjusted_similarity_summary()]
#'
#' @examples
#' sets <- list(c("f1","f2","f3"), c("f2","f3","f4"), c("f1","f3","f5"))
#' pearson_similarity_summary(sets, nb_total_elements = 10L)
#' pearson_similarity_summary(sets, nb_total_elements = 10L, stat = "mean")
#' @export
pearson_similarity_summary <- function(list_of_lists, nb_total_elements,
                                       stat = "median", d = NULL) {
  # Deprecation shim: see pearson_similarity_pairwise().
  if (!missing(d) && !is.null(d)) {
    warning(
      "`d` is deprecated; use `nb_total_elements` instead.", call. = FALSE
    )
    if (is.null(nb_total_elements)) nb_total_elements <- d
  }
  vals <- pearson_similarity_pairwise(list_of_lists, nb_total_elements)
  .similarity_summary(vals, stat)
}

# ── Classification metrics ────────────────────────────────────────────────────

#' FDR Between Two Feature Sets
#'
#' Computes the False Discovery Rate when `list1` is treated as the predicted
#' selection and `list2` as the known ground-truth.  This metric is only
#' meaningful in simulation benchmarks where the true signal features are
#' known in advance.
#'
#' @param list1 Predicted feature set (character or integer vector).
#' @param list2 True feature set (ground truth).
#'
#' @return Numeric scalar in \eqn{[0, 1]}.  Returns 0 when the predicted set
#'   is empty (no false discoveries possible).
#'
#' @seealso [tpr_similarity()], [fscore_similarity()]
#'
#' @examples
#' true_signal <- c("f1", "f2", "f3")
#' predicted   <- c("f1", "f2", "f4", "f5")  # f4, f5 are false discoveries
#' fdr_similarity(predicted, true_signal)      # 2/4 = 0.5
#' fdr_similarity(character(0), true_signal)   # 0 (empty prediction)
#' @export
fdr_similarity <- function(list1, list2) {
  tp <- length(intersect(list1, list2))
  fp <- length(setdiff(list1,  list2))
  if (fp + tp == 0L) return(0)
  fp / (tp + fp)
}

#' TPR Between Two Feature Sets
#'
#' Computes the True Positive Rate (sensitivity / recall) when `list1` is the
#' predicted selection and `list2` is the known ground-truth.  Use together
#' with [fdr_similarity()] and [fscore_similarity()] to characterise the
#' precision-recall trade-off in simulation benchmarks.
#'
#' @param list1 Predicted feature set (character or integer vector).
#' @param list2 True feature set (ground truth).
#'
#' @return Numeric scalar in \eqn{[0, 1]}.  Returns 0 when the true set is
#'   empty.
#'
#' @seealso [fdr_similarity()], [fscore_similarity()]
#'
#' @examples
#' true_signal <- c("f1", "f2", "f3")
#' predicted   <- c("f1", "f2", "f4")  # misses f3
#' tpr_similarity(predicted, true_signal)       # 2/3
#' tpr_similarity(character(0), true_signal)    # 0
#' tpr_similarity(c("f1","f2","f3"), true_signal)  # 1
#' @export
tpr_similarity <- function(list1, list2) {
  tp <- length(intersect(list1, list2))
  fn <- length(setdiff(list2,  list1))
  if (fn + tp == 0L) return(0)
  tp / (tp + fn)
}

#' F-Score Between Two Feature Sets
#'
#' Computes the F\eqn{_\beta} score between a predicted and a ground-truth
#' feature set.  The default `beta = 1` gives the standard F1 score (harmonic
#' mean of precision and recall).  Larger `beta` values weight recall more
#' heavily; smaller values emphasise precision.
#'
#' Use [fdr_similarity()] and [tpr_similarity()] when you want to report
#' precision and recall separately; use `fscore_similarity()` when you need
#' a single number that balances both.
#'
#' @param list1 Predicted feature set (character or integer vector).
#' @param list2 True feature set (ground truth).
#' @param beta Positive numeric; controls trade-off between precision and
#'   recall.  Default 1 gives the F1 score.
#'
#' @return Numeric scalar in \eqn{[0, 1]}.  Returns 0 when both the predicted
#'   and true sets are empty.
#'
#' @seealso [fdr_similarity()], [tpr_similarity()]
#'
#' @examples
#' true_signal <- c("f1", "f2", "f3")
#' predicted   <- c("f1", "f2", "f4")  # 2 TP, 1 FP, 1 FN
#' fscore_similarity(predicted, true_signal)           # F1
#' fscore_similarity(predicted, true_signal, beta = 2) # recall-weighted
#' fscore_similarity(character(0), character(0))       # 0
#' @export
fscore_similarity <- function(list1, list2, beta = 1) {
  tp <- length(intersect(list1, list2))
  fp <- length(setdiff(list1,  list2))
  fn <- length(setdiff(list2,  list1))
  num <- (1 + beta^2) * tp
  den <- (1 + beta^2) * tp + beta^2 * fn + fp
  if (den == 0L) return(0)
  num / den
}

# ── Internal helpers ──────────────────────────────────────────────────────────

.pairwise_upper_triangle <- function(list_of_lists, fn) {
  n    <- length(list_of_lists)
  vals <- numeric(n * (n - 1L) / 2L)
  if (n < 2L) return(vals)
  k <- 1L
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      vals[k] <- fn(list_of_lists[[i]], list_of_lists[[j]])
      k <- k + 1L
    }
  }
  vals
}

.similarity_summary <- function(vals, stat) {
  if (stat == "median") {
    list(
      statistic = stats::median(vals),
      err       = as.numeric(stats::quantile(vals, probs = c(0.25, 0.75)))
    )
  } else if (stat == "mean") {
    mu <- mean(vals)
    list(
      statistic = mu,
      err       = sqrt(mean((vals - mu)^2))
    )
  } else {
    stop(
      "`stat` must be \"median\" or \"mean\", got \"", stat, "\".",
      call. = FALSE
    )
  }
}

# ── Classification metrics (used by multiomic nested-CV workflows) ────────────

# Compute per-class and aggregate classification metrics from truth/predicted
# factor vectors.  Called by both nested_cv.R and multiomic_workflows.R so it
# lives here rather than in either caller file.
.classification_metrics <- function(truth, predicted) {
  if (length(truth) != length(predicted) || length(truth) == 0L) {
    stop("`truth` and `predicted` must be non-empty and have equal length.",
         call. = FALSE)
  }
  if (anyNA(truth) || anyNA(predicted)) {
    stop("Classification metrics require complete truth and predictions.",
         call. = FALSE)
  }
  truth <- if (is.factor(truth)) {
    factor(truth, levels = levels(truth))
  } else {
    factor(truth)
  }
  predicted_values <- as.character(predicted)
  unseen <- setdiff(unique(predicted_values), levels(truth))
  if (length(unseen)) {
    stop(
      "Predicted classes must use the truth levels; unseen: ",
      paste(unseen, collapse = ", "), ".",
      call. = FALSE
    )
  }
  predicted <- factor(predicted_values, levels = levels(truth))
  tab       <- table(truth = truth, predicted = predicted)
  recalls    <- diag(tab) / pmax(rowSums(tab), 1L)
  precisions <- diag(tab) / pmax(colSums(tab), 1L)
  f1 <- ifelse(
    precisions + recalls > 0,
    2 * precisions * recalls / (precisions + recalls),
    0
  )
  accuracy <- sum(diag(tab)) / sum(tab)
  # Matthews Correlation Coefficient (Gorodkin 2004, multiclass generalisation):
  #   MCC = (n * sum_k(C[k,k]) - sum_k(row_k * col_k)) /
  #         sqrt((n^2 - sum(col_k^2)) * (n^2 - sum(row_k^2)))
  # Returns 0 when the denominator is 0 (degenerate predictions).
  n         <- sum(tab)
  row_k     <- rowSums(tab)
  col_k     <- colSums(tab)
  mcc_num   <- n * sum(diag(tab)) - sum(row_k * col_k)
  mcc_denom <- sqrt((n^2 - sum(col_k^2)) * (n^2 - sum(row_k^2)))
  mcc       <- if (mcc_denom == 0) 0 else mcc_num / mcc_denom
  list(
    accuracy           = unname(accuracy),
    balanced_error_rate = unname(1 - mean(recalls)),
    per_class_recall   = recalls,
    macro_f1           = unname(mean(f1)),
    mcc                = unname(mcc),
    confusion          = tab
  )
}

# ── Binary probability metrics (validated workflows) ────────────────────────

# Brier score for P(event): mean((y_event - p)^2).  Lower is better.
.brier_score <- function(y_event, probability) {
  mean((y_event - probability)^2)
}

# Clipped log loss for P(event); probabilities are clipped to [eps, 1 - eps]
# so a confidently wrong prediction is heavily but finitely penalised.
.log_loss_clipped <- function(y_event, probability, eps = 1e-15) {
  p <- pmin(pmax(probability, eps), 1 - eps)
  -mean(y_event * log(p) + (1 - y_event) * log(1 - p))
}

# Rank-based AUROC (Mann-Whitney U with average ranks for ties).  Returns
# NA_real_ when the truth contains a single class: the AUC is undefined
# there, and NA is the explicit diagnostic rather than a fabricated value.
.auroc <- function(y_event, probability) {
  if (length(unique(y_event)) < 2L) return(NA_real_)
  r <- rank(probability, ties.method = "average")
  n1 <- sum(y_event == 1L)
  n0 <- length(y_event) - n1
  unname((sum(r[y_event == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0))
}

# Binary probability metrics with explicit event-class mapping.  `truth` must
# have exactly two levels; `probability` is P(truth == event_class).
.binary_probability_metrics <- function(truth, probability, event_class) {
  # Preserve the declared level set: a single-class assessment fold still has
  # two levels, and AUROC must return its diagnostic NA rather than error.
  truth <- if (is.factor(truth)) {
    factor(truth, levels = levels(truth))
  } else {
    factor(truth)
  }
  if (nlevels(truth) != 2L) {
    stop("Binary probability metrics require exactly two outcome classes.",
         call. = FALSE)
  }
  event_class <- as.character(event_class)
  if (!event_class %in% levels(truth)) {
    stop("`event_class` must be one of the outcome levels.", call. = FALSE)
  }
  if (length(probability) != length(truth) || anyNA(probability)) {
    stop("`probability` must be complete and aligned to `truth`.",
         call. = FALSE)
  }
  if (any(probability < 0 | probability > 1)) {
    stop("`probability` must lie in [0, 1].", call. = FALSE)
  }
  y_event <- as.integer(truth == event_class)
  list(
    brier = .brier_score(y_event, probability),
    log_loss = .log_loss_clipped(y_event, probability),
    auroc = .auroc(y_event, probability),
    event_class = event_class,
    n = length(truth),
    n_events = sum(y_event)
  )
}

# ── Deprecation shims (0.3.0 API overhaul) ────────────────────────────────────
# The 0.3.0 rename scheme: `*_similarity()` for a pair, `*_similarity_matrix()`
# for an all-pairs matrix, `*_similarity_pairwise()` for the vector of unique
# pairwise values, `*_similarity_summary()` for a one-number summary.  The old
# names remain as shims emitting a one-time-per-call deprecation warning.

#' Deprecated: use [jaccard_similarity_matrix()]
#'
#' `jaccard_matrix()` is deprecated and will be removed in a future release;
#' use [jaccard_similarity_matrix()] instead.
#' @param list_of_lists See [jaccard_similarity_matrix()].
#' @examples
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' # warns; use jaccard_similarity_matrix()
#' same <- suppressWarnings(identical(
#'   jaccard_matrix(sets),
#'   jaccard_similarity_matrix(sets)
#' ))
#' stopifnot(same)
#' @keywords internal
#' @export
jaccard_matrix <- function(list_of_lists, remove_diag = TRUE) {
  warning(
    "`jaccard_matrix()` is deprecated; use `jaccard_similarity_matrix()` ",
    "instead.", call. = FALSE
  )
  jaccard_similarity_matrix(list_of_lists, remove_diag = remove_diag)
}

#' Deprecated: use [adjusted_similarity_pairwise()]
#'
#' `adjusted_similarity_values()` is deprecated and will be removed in a
#' future release; use [adjusted_similarity_pairwise()] instead.
#' @examples
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' # warns; use adjusted_similarity_pairwise()
#' same <- suppressWarnings(identical(
#'   adjusted_similarity_values(sets, nb_total_elements = 10L),
#'   adjusted_similarity_pairwise(sets, nb_total_elements = 10L)
#' ))
#' stopifnot(same)
#' @keywords internal
#' @export
adjusted_similarity_values <- function(list_of_lists, nb_total_elements) {
  warning(
    "`adjusted_similarity_values()` is deprecated; use ",
    "`adjusted_similarity_pairwise()` instead.", call. = FALSE
  )
  adjusted_similarity_pairwise(list_of_lists, nb_total_elements)
}

#' Deprecated: use [adjusted_similarity_summary()]
#'
#' `adjusted_similarity_measure()` is deprecated and will be removed in a
#' future release; use [adjusted_similarity_summary()] instead.
#' @param list_of_lists,nb_total_elements,stat See
#'   [adjusted_similarity_summary()].
#' @examples
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' # warns; use adjusted_similarity_summary()
#' same <- suppressWarnings(identical(
#'   adjusted_similarity_measure(sets, nb_total_elements = 10L),
#'   adjusted_similarity_summary(sets, nb_total_elements = 10L)
#' ))
#' stopifnot(same)
#' @keywords internal
#' @export
adjusted_similarity_measure <- function(list_of_lists, nb_total_elements,
                                        stat = "median") {
  warning(
    "`adjusted_similarity_measure()` is deprecated; use ",
    "`adjusted_similarity_summary()` instead.", call. = FALSE
  )
  adjusted_similarity_summary(list_of_lists, nb_total_elements, stat = stat)
}

#' Deprecated: use [pearson_similarity_pairwise()]
#'
#' `pearson_similarity_values()` is deprecated and will be removed in a
#' future release; use [pearson_similarity_pairwise()] instead.
#' @examples
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' # warns; use pearson_similarity_pairwise()
#' same <- suppressWarnings(identical(
#'   pearson_similarity_values(sets, d = 10L),
#'   pearson_similarity_pairwise(sets, nb_total_elements = 10L)
#' ))
#' stopifnot(same)
#' @keywords internal
#' @export
pearson_similarity_values <- function(list_of_lists, d) {
  warning(
    "`pearson_similarity_values()` is deprecated; use ",
    "`pearson_similarity_pairwise()` instead.", call. = FALSE
  )
  pearson_similarity_pairwise(list_of_lists, d)
}

#' Deprecated: use [pearson_similarity_summary()]
#'
#' `pearson_similarity_measure()` is deprecated and will be removed in a
#' future release; use [pearson_similarity_summary()] instead.
#' @examples
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' # warns; use pearson_similarity_summary()
#' same <- suppressWarnings(identical(
#'   pearson_similarity_measure(sets, d = 10L),
#'   pearson_similarity_summary(sets, nb_total_elements = 10L)
#' ))
#' stopifnot(same)
#' @keywords internal
#' @export
pearson_similarity_measure <- function(list_of_lists, d, stat = "median") {
  warning(
    "`pearson_similarity_measure()` is deprecated; use ",
    "`pearson_similarity_summary()` instead.", call. = FALSE
  )
  pearson_similarity_summary(list_of_lists, d, stat = stat)
}
