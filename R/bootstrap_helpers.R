#' Classic Bootstrap Sampler Indices
#'
#' Draws a random set of row indices from a training set for a single STABL
#' bootstrap iteration.  This is the standard (non-grouped) sampler used when
#' all samples are independent.
#'
#' STABL accumulates how often each feature is selected across many bootstrap
#' subsamples; these subsamples must be drawn consistently to avoid bias.
#' This function centralises that sampling so that class-imbalance handling
#' and the "degenerate bootstrap" guard (see below) are applied uniformly.
#'
#' @details
#' **Class weighting:** When `class_weights` is supplied each sample's
#' inclusion probability is proportional to the weight of its class.  This is
#' useful for imbalanced binary or multi-class outcomes where unweighted
#' subsampling would frequently produce single-class bootstraps.
#'
#' **Degenerate bootstrap guard:** When the drawn subsample contains only one
#' unique class but the full outcome contains at least two, the function
#' retries (without a seed) until a class-diverse subsample is found.  This
#' prevents downstream model failures in classifiers that require at least two
#' classes.
#'
#' @param y Outcome vector whose length equals the number of training samples.
#'   Values are only examined to check class diversity and to apply
#'   `class_weights`; the type is arbitrary (numeric, factor, character).
#' @param n_subsamples Positive integer; number of row indices to draw.  Must
#'   not exceed `length(y)` when `replace = FALSE`.
#' @param replace Logical; sample with replacement?  Default `TRUE`.
#' @param class_weights Optional named numeric vector keyed by class labels
#'   (as produced by `as.character(y)`).  When `NULL` all samples are equally
#'   likely.
#' @param stratify Logical; if `TRUE`, draw approximately the same class
#'   proportions as `y` by sampling within each outcome class.  This is useful
#'   for small or imbalanced classification tasks where unstratified
#'   subsampling can produce very sparse minority-class bootstraps.  Default
#'   `FALSE` preserves the original STABL sampling behavior.
#' @param strata Optional categorical stratification design.  Provide a vector
#'   for one stratification factor, or a `data.frame`/matrix/list for a joint
#'   design across multiple factors.  Bootstrap samples are drawn within the
#'   interaction of all supplied columns.  When `NULL`, no stratification is
#'   used unless `stratify = TRUE`, in which case `y` is used as the strata.
#' @param seed Optional integer; passed to [set.seed()] before sampling for
#'   reproducibility.  `NULL` leaves the RNG state unchanged.
#'
#' @return Integer vector of length `n_subsamples` containing 1-based row
#'   indices into the training set.
#'
#' @seealso [group_bootstrap_indices()] for repeated-measures/grouped data,
#'   [stabl_fit()] which calls this sampler automatically.
#'
#' @examples
#' set.seed(42L)
#' y <- c(rep(0, 15), rep(1, 15))  # balanced binary outcome
#' idx <- classic_bootstrap_indices(y, n_subsamples = 20L, seed = 1L)
#' table(y[idx])  # class distribution in the bootstrap
#'
#' # Stratified sampling preserves class proportions
#' idx_str <- classic_bootstrap_indices(y, n_subsamples = 20L,
#'                                      stratify = TRUE, seed = 1L)
#' table(y[idx_str])
#' @export
classic_bootstrap_indices <- function(
  y,
  n_subsamples,
  replace = TRUE,
  class_weights = NULL,
  stratify = FALSE,
  strata = NULL,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- length(y)
  if (!replace && n_subsamples > n) {
    stop("`n_subsamples` cannot exceed sample size when replace = FALSE.", call. = FALSE)
  }

  if (!is.logical(stratify) || length(stratify) != 1L || is.na(stratify)) {
    stop("`stratify` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(stratify) && !is.null(class_weights)) {
    stop("`class_weights` cannot be combined with `stratify = TRUE`.", call. = FALSE)
  }
  if (!is.null(strata) && !is.null(class_weights)) {
    stop("`class_weights` cannot be combined with `strata`.", call. = FALSE)
  }

  probs <- NULL
  if (!is.null(class_weights)) {
    class_ids <- as.character(y)
    if (is.null(names(class_weights))) {
      stop("`class_weights` must be a named numeric vector.", call. = FALSE)
    }
    if (!all(unique(class_ids) %in% names(class_weights))) {
      stop("All observed classes must have an entry in `class_weights`.", call. = FALSE)
    }
    probs <- unname(class_weights[class_ids])
    probs <- probs / sum(probs)
  }

  strata_ids <- .bootstrap_strata_ids(
    strata = if (isTRUE(stratify) && is.null(strata)) y else strata,
    n = n,
    arg = "strata"
  )

  draw_once <- if (!is.null(strata_ids)) {
    function() .stratified_bootstrap_indices(strata_ids, n_subsamples, replace)
  } else {
    function() sample.int(n = n, size = n_subsamples, replace = replace, prob = probs)
  }

  idx <- draw_once()

  idx <- .guard_diverse_bootstrap(
    y, idx, draw_once,
    caller = "classic_bootstrap_indices",
    hint   = "Check class balance or increase `sample_fraction`."
  )

  idx
}

#' Group-Aware Bootstrap Sampler Indices
#'
#' Draws row indices for a single STABL bootstrap iteration by sampling
#' **groups** rather than individual samples.  Use this instead of
#' [classic_bootstrap_indices()] whenever the training data contains repeated
#' measurements or known clusters (e.g. multiple time-points per subject,
#' technical replicates, or family members).
#'
#' @details
#' **Why group-level sampling prevents leakage:** If individual rows from the
#' same subject appear in both the bootstrap subsample and its complement, the
#' stability score becomes inflated because the learner can partially memorise
#' subject-level patterns.  By sampling complete groups, the subsample and its
#' "holdout" are subject-disjoint, preserving the independence assumption
#' underlying STABL's FDP+ diagnostic assumptions.
#'
#' **How groups are sampled:** Groups are drawn one at a time (without
#' replacement by default) until the running tally of included rows reaches
#' or exceeds `n_subsamples`.  Whole groups are always kept intact: if the
#' last added group causes the count to exceed `n_subsamples`, the surplus
#' rows are *not* trimmed.  This guarantees the strongest leakage prevention
#' (the subsample and its complement are always group-disjoint), at the cost
#' of a slightly variable subsample size.  When `replace = TRUE` the same
#' group may be drawn multiple times (useful when there are very few groups).
#'
#' **Degenerate bootstrap guard:** As with [classic_bootstrap_indices()], if
#' the final subsample contains only one unique class the function retries
#' until a class-diverse sample is obtained.
#'
#' @param y Outcome vector with the same length as `groups`.
#' @param groups Vector of group identifiers (same length as `y`).  Values may
#'   be any type that supports `unique()` and equality comparison.
#' @param n_subsamples Positive integer; target number of rows in the
#'   subsample.  The actual count may differ slightly due to whole-group
#'   granularity (whole groups are always preserved; the realised count is
#'   the smallest group-aligned value `>= n_subsamples`).
#' @param replace Logical; may the same group be sampled more than once?
#'   Default `FALSE`.
#' @param stratify Logical; if `TRUE`, sample whole groups within outcome
#'   strata so the realised subsample retains class representation where the
#'   group structure permits it.  Each group must map to exactly one outcome
#'   class.  Whole groups are still preserved, so realised counts can exceed
#'   per-class targets.  Default `FALSE`.
#' @param strata Optional categorical stratification design.  Provide a vector
#'   for one stratification factor, or a `data.frame`/matrix/list for a joint
#'   design across multiple factors.  Each group must map to exactly one
#'   realised stratum.
#' @param seed Optional integer; passed to [set.seed()] before sampling.
#'   `NULL` leaves the RNG state unchanged.
#'
#' @return Integer vector of 1-based row indices into the training set.
#'
#' @seealso [classic_bootstrap_indices()] for independent-sample data,
#'   [stabl_fit()] which calls this sampler automatically when `groups` is
#'   supplied.
#'
#' @examples
#' # 3 groups with 4 samples each; binary outcome
#' y      <- c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
#' groups <- rep(c("G1", "G2", "G3"), each = 4L)
#' idx    <- group_bootstrap_indices(y, groups, n_subsamples = 8L, seed = 1L)
#' table(groups[idx])  # whole groups only; count >= 8
#' @export
group_bootstrap_indices <- function(y, groups, n_subsamples, replace = FALSE,
                                    stratify = FALSE, strata = NULL,
                                    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(groups) != length(y)) {
    stop("`groups` must be the same length as `y`.", call. = FALSE)
  }

  n <- length(y)
  if (!replace && n_subsamples > n) {
    stop("`n_subsamples` cannot exceed sample size when replace = FALSE.", call. = FALSE)
  }
  if (!is.logical(stratify) || length(stratify) != 1L || is.na(stratify)) {
    stop("`stratify` must be TRUE or FALSE.", call. = FALSE)
  }
  strata_ids <- .bootstrap_strata_ids(
    strata = if (isTRUE(stratify) && is.null(strata)) y else strata,
    n = n,
    arg = "strata"
  )

  group_levels <- unique(groups)
  draw_once <- if (!is.null(strata_ids)) {
    function() .stratified_group_bootstrap_indices(
      strata_ids = strata_ids,
      groups = groups,
      group_levels = group_levels,
      n_subsamples = n_subsamples,
      replace = replace
    )
  } else {
    function() .unstratified_group_bootstrap_indices(
      groups = groups,
      group_levels = group_levels,
      n_subsamples = n_subsamples,
      replace = replace
    )
  }

  sampled_idx <- draw_once()
  # Whole-group invariant (D3): never trim partial groups, even if the
  # final tally exceeds n_subsamples.

  sampled_idx <- .guard_diverse_bootstrap(
    y, sampled_idx, draw_once,
    caller = "group_bootstrap_indices",
    hint   = "Check class balance or group structure."
  )

  sampled_idx
}

.guard_diverse_bootstrap <- function(y, idx, draw_once, caller, hint,
                                     max_retries = 1000L) {
  if (length(unique(y)) < 2L) return(idx)
  attempt <- 0L
  while (length(unique(y[idx])) < 2L) {
    attempt <- attempt + 1L
    if (attempt > max_retries) {
      stop(
        caller, ": could not draw a class-diverse subsample after ",
        max_retries, " attempts. ", hint,
        call. = FALSE
      )
    }
    idx <- draw_once()
  }
  idx
}

.bootstrap_strata_ids <- function(strata, n, arg = "strata") {
  if (is.null(strata)) {
    return(NULL)
  }

  if (is.data.frame(strata)) {
    design <- strata
  } else if (is.matrix(strata)) {
    design <- as.data.frame(strata, stringsAsFactors = FALSE)
  } else if (is.list(strata) && !is.null(strata)) {
    design <- as.data.frame(strata, stringsAsFactors = FALSE)
  } else {
    if (length(strata) != n) {
      stop("`", arg, "` must have length ", n, ".", call. = FALSE)
    }
    design <- data.frame(stratum = strata, check.names = FALSE)
  }

  if (nrow(design) != n) {
    stop("`", arg, "` must have ", n, " rows.", call. = FALSE)
  }
  if (ncol(design) == 0L) {
    stop("`", arg, "` must contain at least one column.", call. = FALSE)
  }
  if (anyNA(design)) {
    stop("`", arg, "` cannot contain missing values.", call. = FALSE)
  }

  pieces <- lapply(design, function(x) as.character(x))
  do.call(interaction, c(pieces, list(drop = TRUE, sep = "\r", lex.order = TRUE)))
}

.subset_bootstrap_strata_by_ids <- function(strata, sample_ids,
                                            arg = "bootstrap_strata") {
  if (is.null(strata)) {
    return(NULL)
  }

  n <- length(sample_ids)
  if (is.matrix(strata)) {
    strata <- as.data.frame(strata, stringsAsFactors = FALSE)
  }
  if (is.data.frame(strata)) {
    out <- strata
    rn <- rownames(out)
    has_ids <- .bootstrap_strata_has_row_ids(rn, sample_ids, nrow(out))
    if (has_ids) {
      if (!all(sample_ids %in% rn)) {
        stop("Sample mismatch between `x` row names and `", arg, "` row names.",
             call. = FALSE)
      }
      return(out[sample_ids, , drop = FALSE])
    }
    if (nrow(out) != n) {
      stop("`", arg, "` must have one row per sample or row names matching `x`.",
           call. = FALSE)
    }
    rownames(out) <- sample_ids
    return(out)
  }

  if (is.list(strata) && !is.null(strata) && !is.atomic(strata)) {
    return(.subset_bootstrap_strata_by_ids(
      as.data.frame(strata, stringsAsFactors = FALSE),
      sample_ids = sample_ids,
      arg = arg
    ))
  }

  if (is.null(names(strata))) {
    if (length(strata) != n) {
      stop("`", arg, "` must have one value per sample or names matching `x`.",
           call. = FALSE)
    }
    return(stats::setNames(strata, sample_ids))
  }
  if (!all(sample_ids %in% names(strata))) {
    stop("Sample mismatch between `x` row names and `", arg, "` names.",
         call. = FALSE)
  }
  strata[sample_ids]
}

.bootstrap_strata_has_row_ids <- function(row_ids, sample_ids, n_rows) {
  if (is.null(row_ids) || anyNA(row_ids) || any(!nzchar(row_ids))) {
    return(FALSE)
  }
  default_row_ids <- as.character(seq_len(n_rows))
  !identical(row_ids, default_row_ids) || setequal(row_ids, sample_ids)
}

.stratified_counts <- function(strata_ids, n_subsamples, replace,
                               min_per_stratum = 1L) {
  strata_tab <- table(strata_ids)
  n_strata <- length(strata_tab)

  min_per_stratum <- .validate_scalar_integer_like(
    min_per_stratum,
    "min_per_stratum",
    min = 1L
  )

  if (n_strata == 0L) {
    stop("`strata` must contain at least one sample.", call. = FALSE)
  }
  minimum_total <- n_strata * min_per_stratum
  if (n_subsamples < minimum_total) {
    if (min_per_stratum == 1L) {
      stop(
        "`n_subsamples` must be at least the number of realised strata.",
        call. = FALSE
      )
    }
    stop(
      "`n_subsamples` must be at least ", minimum_total,
      " to allocate ", min_per_stratum,
      " observation(s) to each realised stratum.",
      call. = FALSE
    )
  }
  if (!replace && any(as.integer(strata_tab) < min_per_stratum)) {
    stop(
      "Each realised stratum must contain at least ", min_per_stratum,
      " observation(s) when `replace = FALSE`.",
      call. = FALSE
    )
  }
  if (n_strata == 1L) {
    out <- n_subsamples
    names(out) <- names(strata_tab)
    out <- as.integer(out)
    names(out) <- names(strata_tab)
    return(out)
  }
  exact <- as.numeric(strata_tab) / sum(strata_tab) * n_subsamples
  counts <- floor(exact)
  names(counts) <- names(strata_tab)

  below_minimum <- counts < min_per_stratum
  if (any(below_minimum)) {
    counts[below_minimum] <- min_per_stratum
  }

  remainders <- exact - floor(exact)
  names(remainders) <- names(strata_tab)
  needed <- n_subsamples - sum(counts)
  if (needed > 0L) {
    order_names <- names(sort(remainders, decreasing = TRUE))
    for (i in seq_len(needed)) {
      pick <- order_names[[((i - 1L) %% length(order_names)) + 1L]]
      counts[[pick]] <- counts[[pick]] + 1L
    }
  }
  while (sum(counts) > n_subsamples) {
    candidates <- names(counts[counts > min_per_stratum])
    if (length(candidates) == 0L) {
      stop("Could not allocate stratified bootstrap counts.", call. = FALSE)
    }
    pick <- candidates[[which.max(counts[candidates])]]
    counts[[pick]] <- counts[[pick]] - 1L
  }

  if (!replace && any(counts > as.integer(strata_tab[names(counts)]))) {
    stop(
      "Stratified bootstrap count exceeds an observed stratum size with ",
      "`replace = FALSE`. Increase `sample_fraction` or use replacement.",
      call. = FALSE
    )
  }

  out <- as.integer(counts)
  names(out) <- names(counts)
  out
}

.stratified_bootstrap_indices <- function(strata_ids, n_subsamples, replace,
                                          min_per_stratum = 1L) {
  strata_ids <- as.character(strata_ids)
  counts <- .stratified_counts(
    strata_ids,
    n_subsamples,
    replace,
    min_per_stratum = min_per_stratum
  )
  idx_by_stratum <- split(seq_along(strata_ids), strata_ids)

  idx <- unlist(
    lapply(names(counts), function(stratum) {
      sample(idx_by_stratum[[stratum]], size = counts[[stratum]], replace = replace)
    }),
    use.names = FALSE
  )
  sample(idx, length(idx), replace = FALSE)
}

.unstratified_group_bootstrap_indices <- function(groups, group_levels,
                                                  n_subsamples, replace) {
  .sample_groups_until(group_levels, groups, n_subsamples, replace)
}

# Shared group-sampling while-loop: draw groups until target_n indices are
# accumulated; use index-then-subset to avoid R's sample(x,1) length-1 pitfall.
.sample_groups_until <- function(group_levels, groups, target_n, replace) {
  remaining  <- group_levels
  idx_chunks <- list()
  total_len  <- 0L

  while (total_len < target_n && length(remaining) > 0L) {
    pick_pos  <- sample.int(length(remaining), size = 1L)
    g         <- remaining[[pick_pos]]
    remaining <- if (replace) remaining else remaining[-pick_pos]
    chunk     <- which(groups == g)
    idx_chunks[[length(idx_chunks) + 1L]] <- chunk
    total_len <- total_len + length(chunk)
  }

  idx <- unlist(idx_chunks, use.names = FALSE)
  if (!replace) unique(idx) else idx
}

.stratified_group_bootstrap_indices <- function(strata_ids, groups, group_levels,
                                                n_subsamples, replace,
                                                min_per_stratum = 1L) {
  group_strata <- vapply(group_levels, function(g) {
    stratum <- unique(as.character(strata_ids[groups == g]))
    if (length(stratum) != 1L) {
      stop(
        "Grouped stratified bootstrap requires each group to map to exactly ",
        "one realised stratum.",
        call. = FALSE
      )
    }
    stratum
  }, character(1L))

  target_counts <- .stratified_counts(
    strata_ids,
    n_subsamples,
    replace,
    min_per_stratum = min_per_stratum
  )
  # List-collect across strata to avoid quadratic c()-in-loop growth;
  # unlist() once at the end for both the inner (per-stratum) and outer
  # (cross-stratum) accumulations.
  outer_chunks <- list()

  for (stratum in names(target_counts)) {
    stratum_idx <- .sample_groups_until(
      group_levels[group_strata == stratum], groups, target_counts[[stratum]], replace
    )

    if (length(stratum_idx) < target_counts[[stratum]]) {
      stop(
        "Could not satisfy grouped stratified bootstrap target for stratum `",
        stratum, "`.",
        call. = FALSE
      )
    }
    outer_chunks[[length(outer_chunks) + 1L]] <- stratum_idx
  }

  sampled_idx <- unlist(outer_chunks, use.names = FALSE)
  if (!replace) sampled_idx <- unique(sampled_idx)

  sample(sampled_idx, length(sampled_idx), replace = FALSE)
}
