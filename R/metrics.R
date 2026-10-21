# metrics.R — Selection-stability similarity measures
# R port of stabl/metrics.py
#
# These functions compare lists of selected features across repeat runs of
# STABL (or any other feature-selection method) to assess reproducibility.
# The user-facing entry points are two functions, in the mixOmics style of a
# small verb-based API:
#   - stabl_similarity():         scalar or pairwise similarity under one of
#                                 three methods (Jaccard, adjusted, Pearson);
#   - stabl_similarity_summary(): one-number summary of pairwise similarity.
# The eight historical 0.3.0 function names (jaccard_similarity,
# jaccard_similarity_matrix, adjusted_similarity, adjusted_similarity_pairwise,
# adjusted_similarity_summary, pearson_similarity, pearson_similarity_pairwise,
# pearson_similarity_summary) were removed ahead of version 1.0.0; use the
# two entry points above.
# For simulation benchmarks where a ground-truth feature set is known, the
# file also provides FDR, TPR, and F-score functions.

# ── Similarity: unified interface ────────────────────────────────────────────

#' Similarity Between Feature Sets
#'
#' Computes the similarity between two sets of selected features, or between
#' all pairs of sets in a list, using one of three methods.  Use this to
#' assess the reproducibility of STABL (or any other feature-selection
#' method) across cross-validation folds or repeated runs.
#'
#' Three methods are available:
#' \describe{
#'   \item{`"jaccard"` (default)}{Set-overlap fraction
#'     \eqn{|A \cap B| / |A \cup B|}.  The most interpretable measure; it
#'     does not correct for the chance overlap expected between small or
#'     large sets.  When both sets are empty the value is 0 (convention
#'     matching the Python implementation).}
#'   \item{`"adjusted"`}{Chance-corrected overlap, analogous to Cohen's
#'     kappa for sets:
#'     \deqn{S_{\text{adj}}(A, B) =
#'       \frac{r - \mathbb{E}[r]}{\min(k_1, k_2) - \max(0, k_1 + k_2 - d)}}
#'     where \eqn{r = |A \cap B|}, \eqn{k_i = |A_i|}, \eqn{d =}
#'     `n_universe`, and \eqn{\mathbb{E}[r] = k_1 k_2 / d}.  Returns 0 when
#'     either set is empty or equals the full universe.}
#'   \item{`"pearson"`}{Pearson-correlation-inspired chance correction (the
#'     phi coefficient):
#'     \deqn{S_{\text{Pearson}}(A, B) =
#'       \frac{r - k_i k_j / d}{d \cdot \upsilon_i \upsilon_j}}
#'     with \eqn{\upsilon_i = \sqrt{\pi_i (1 - \pi_i)}} and
#'     \eqn{\pi_i = k_i / d}.  Returns 1 when both sets are empty or both
#'     equal the universe; 0 when exactly one of them is.}
#' }
#'
#' @param x Either a character/integer vector of selected feature
#'   identifiers (then `y` must also be supplied and the result is a single
#'   scalar), or a list of such vectors, one per STABL run (then `y` must be
#'   `NULL` and all pairwise similarities are computed).
#' @param y `NULL` (default) or a character/integer vector of selected
#'   feature identifiers.  Supply `y` for scalar mode, or a list as `x` for
#'   pairwise mode — not both.
#' @param method Character; one of `"jaccard"` (default), `"adjusted"`, or
#'   `"pearson"`.
#' @param n_universe Integer; total number of candidate features in the
#'   universe (i.e. the number of columns in the original predictor matrix).
#'   Required for `method = "adjusted"` and `method = "pearson"`; ignored
#'   for `method = "jaccard"`.
#' @param as_matrix Logical; pairwise mode only.  If `FALSE` (default),
#'   returns the \eqn{N (N - 1) / 2} off-diagonal upper-triangle values as a
#'   flat vector in row-major order (matching the Python convention).  If
#'   `TRUE`, returns the full symmetric \eqn{N \times N} matrix.
#' @param remove_diag Logical; only used when `as_matrix = TRUE`.  If `TRUE`
#'   the self-similarity diagonal is removed row-wise, yielding an
#'   \eqn{N \times (N - 1)} matrix in which row \eqn{i} contains the
#'   similarities of run \eqn{i} with every other run.
#'
#' @return In scalar mode, a numeric scalar.  In pairwise mode with
#'   `as_matrix = FALSE`, a numeric vector of length \eqn{N (N - 1) / 2}.
#'   In pairwise mode with `as_matrix = TRUE`, a numeric matrix of dimension
#'   \eqn{N \times N} (or \eqn{N \times (N - 1)} when `remove_diag = TRUE`).
#'
#'   The matrix diagonal holds the self-similarity of each set: exactly 1
#'   for `"jaccard"` and `"pearson"`, and the computed self-value for
#'   `"adjusted"` (\eqn{1 - k / d} for a set of size \eqn{k \le d / 2}),
#'   which is below 1 by construction because the measure is
#'   chance-corrected.
#'
#' @seealso [stabl_similarity_summary()] for reducing pairwise values to a
#'   single location statistic.
#'
#' @examples
#' # Scalar mode: two sets
#' stabl_similarity(c("A", "B", "C"), c("B", "C", "D"))                # 0.5
#' stabl_similarity(c("f1", "f2", "f3"), c("f2", "f3", "f4"),
#'                  method = "adjusted", n_universe = 10L)
#'
#' # Pairwise mode: a list of sets, one per run
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' stabl_similarity(sets)                                        # flat vector
#' stabl_similarity(sets, as_matrix = TRUE)                      # 3 x 3
#' stabl_similarity(sets, as_matrix = TRUE, remove_diag = TRUE)  # 3 x 2
#' stabl_similarity(sets, method = "adjusted", n_universe = 10L)
#' @export
stabl_similarity <- function(x, y = NULL,
                             method = c("jaccard", "adjusted", "pearson"),
                             n_universe = NULL,
                             as_matrix = FALSE,
                             remove_diag = FALSE) {
  method <- match.arg(method)
  if (method %in% c("adjusted", "pearson") && is.null(n_universe)) {
    stop("argument \"n_universe\" is missing, with no default", call. = FALSE)
  }
  if (is.list(x) && is.null(y)) {
    if (!as_matrix) {
      return(.pairwise_upper_triangle(
        x, .similarity_scalar_fn(method, n_universe)
      ))
    }
    return(.similarity_matrix(x, method, n_universe, remove_diag = remove_diag))
  }
  if (!is.list(x) && !is.null(y)) {
    return(.similarity_scalar_fn(method, n_universe)(x, y))
  }
  stop(
    "Supply either two feature sets (`x` and `y` as vectors) for a single ",
    "similarity value, or a list of feature sets as `x` (with `y = NULL`) ",
    "for all pairwise values.",
    call. = FALSE
  )
}

#' Summary Statistic of Pairwise Similarity Values
#'
#' Computes all pairwise similarities between the feature sets in `x` (via
#' [stabl_similarity()]) and reduces them to a single location statistic
#' with an associated spread measure.  Useful for reporting a single
#' reproducibility number per STABL configuration in benchmarking tables.
#'
#' @param x A list of character/integer vectors, one per STABL run.
#' @param method Character; one of `"jaccard"` (default), `"adjusted"`, or
#'   `"pearson"`.  See [stabl_similarity()] for method descriptions.
#' @param n_universe Integer; total number of candidate features in the
#'   universe.  Required for `method = "adjusted"` and
#'   `method = "pearson"`; ignored for `method = "jaccard"`.
#' @param stat Character; `"median"` (default, robust to outliers) or
#'   `"mean"`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`statistic`}{The median (or mean) of the pairwise similarity
#'       values.}
#'     \item{`err`}{For `"median"`: the 25th and 75th percentile vector
#'       (IQR bounds).  For `"mean"`: the root-mean-squared deviation
#'       (RMSD / population SD).}
#'   }
#'
#' @seealso [stabl_similarity()]
#'
#' @examples
#' sets <- list(c("f1", "f2", "f3"), c("f2", "f3", "f4"), c("f1", "f3", "f5"))
#' stabl_similarity_summary(sets)
#' stabl_similarity_summary(sets, method = "adjusted", n_universe = 10L)
#' stabl_similarity_summary(sets, method = "pearson", n_universe = 10L,
#'                          stat = "mean")
#' @export
stabl_similarity_summary <- function(x,
                                     method = c("jaccard", "adjusted", "pearson"),
                                     n_universe = NULL,
                                     stat = c("median", "mean")) {
  method <- match.arg(method)
  stat   <- match.arg(stat)
  if (method %in% c("adjusted", "pearson") && is.null(n_universe)) {
    stop("argument \"n_universe\" is missing, with no default", call. = FALSE)
  }
  if (length(x) < 2L) {
    stop("`x` must contain at least two sets to summarize pairwise ",
         "similarity.", call. = FALSE)
  }
  vals <- .pairwise_upper_triangle(x, .similarity_scalar_fn(method, n_universe))
  .similarity_summary(vals, stat)
}

# ── Internal similarity engines ───────────────────────────────────────────────

# Scalar Jaccard overlap between two sets (0 when both are empty).
.jaccard_scalar <- function(set1, set2) {
  s1 <- unique(set1)
  s2 <- unique(set2)
  inter <- length(intersect(s1, s2))
  union_size <- length(s1) + length(s2) - inter
  if (inter == 0L && union_size == 0L) return(0)
  inter / union_size
}

# Scalar chance-corrected (Cohen's-kappa-like) similarity between two sets.
.adjusted_scalar <- function(set1, set2, n_universe) {
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

# Scalar Pearson-corrected (phi-coefficient) similarity between two sets.
.pearson_scalar <- function(set1, set2, n_universe) {
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

# Binary scalar function for the chosen method, with n_universe bound in.
.similarity_scalar_fn <- function(method, n_universe = NULL) {
  switch(
    method,
    jaccard  = .jaccard_scalar,
    adjusted = function(a, b) .adjusted_scalar(a, b, n_universe),
    pearson  = function(a, b) .pearson_scalar(a, b, n_universe)
  )
}

# Full symmetric N x N pairwise similarity matrix.  The diagonal is 1 for
# the jaccard and pearson methods (self-similarity is always 1 under their
# edge rules) and the honest computed self-value for the adjusted method.
.similarity_matrix <- function(list_of_lists, method, n_universe,
                               remove_diag = FALSE) {
  n <- length(list_of_lists)
  if (n < 1L) {
    stop("`list_of_lists` must contain at least one set.", call. = FALSE)
  }
  fn  <- .similarity_scalar_fn(method, n_universe)
  mat <- matrix(0, nrow = n, ncol = n)
  if (method == "adjusted") {
    for (i in seq_len(n)) {
      mat[i, i] <- fn(list_of_lists[[i]], list_of_lists[[i]])
    }
  } else {
    diag(mat) <- 1
  }
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        v         <- fn(list_of_lists[[i]], list_of_lists[[j]])
        mat[i, j] <- v
        mat[j, i] <- v
      }
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
#' stablr:::fdr_similarity(predicted, true_signal)      # 2/4 = 0.5
#' stablr:::fdr_similarity(character(0), true_signal)   # 0 (empty prediction)
#' @keywords internal
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
#' stablr:::tpr_similarity(predicted, true_signal)       # 2/3
#' stablr:::tpr_similarity(character(0), true_signal)    # 0
#' stablr:::tpr_similarity(c("f1","f2","f3"), true_signal)  # 1
#' @keywords internal
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
#' stablr:::fscore_similarity(predicted, true_signal)           # F1
#' stablr:::fscore_similarity(predicted, true_signal, beta = 2) # recall-weighted
#' stablr:::fscore_similarity(character(0), character(0))       # 0
#' @keywords internal
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
  # Consistent with .classification_metrics: NA truth is a caller error, not
  # a value to propagate silently into brier/log_loss/auroc.
  if (anyNA(truth)) {
    stop("`truth` must not contain NA values.", call. = FALSE)
  }
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
