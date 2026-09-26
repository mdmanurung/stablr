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
#'   likely.  Not available in [group_bootstrap_indices()]: weighting there would
#'   have to be defined over groups rather than samples, and group-level
#'   imbalance is already addressed exactly by `stratify = TRUE`, which samples
#'   whole groups within outcome strata.
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
#' @param seed Optional integer; seeds this call for reproducibility.  The
#'   caller's RNG state is saved and restored, so a seeded call does not
#'   advance the global stream.  `NULL` draws from the ambient stream.
#'
#' @section Seeding in a loop:
#' Because `seed` scopes to a single call and the prior RNG state is restored
#' afterwards, calling this function repeatedly with the same `seed` returns
#' the same indices every time. To draw many distinct replicates, either leave
#' `seed = NULL` and seed once outside the loop, or vary the seed per iteration.
#' [stabl_fit()] does the latter: it never passes `seed` here and instead
#' derives one seed for the whole bootstrap loop.
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

  idx <- .with_local_seed(seed, {
    idx <- draw_once()

    .guard_diverse_bootstrap(
      y, idx, draw_once,
      caller = "classic_bootstrap_indices",
      hint   = "Check class balance or increase `sample_fraction`."
    )
  })

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
#' **Why group-level sampling:** STABL has no holdout.  The bootstrap subsample
#' is the only data any replicate is fitted on, so a subject that appears on
#' both sides of the draw does not leak into a model -- there is no model
#' fitted on the complement.  What whole-group sampling buys is an
#' exchangeable draw: the rows of one subject are not independent replicates, so
#' a row-level subsample contains correlated rows and the resulting stability
#' scores concentrate on subject-specific patterns instead of on features that
#' generalise.  Drawing whole subjects keeps the subsample a set of independent
#' units.  This is orthogonal to the FDP+ calibration, which assumes the
#' artificial features are exchangeable with the real ones and says nothing
#' about subject disjointness.
#'
#' **How groups are sampled:** Groups are drawn one at a time (without
#' replacement by default) until the running tally of included rows reaches
#' or exceeds `n_subsamples`.  Whole groups are always kept intact: if the
#' last added group causes the count to exceed `n_subsamples`, the surplus
#' rows are *not* trimmed.  The realised count is therefore the smallest
#' group-aligned value `>= n_subsamples`, which can be materially larger than
#' `n_subsamples` when groups are large; a draw that cannot reach the target at
#' all is an error, never a short subsample.  When `replace = TRUE` the same
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
#'   subsample.  The actual count may be larger due to whole-group
#'   granularity (whole groups are always preserved; the realised count is
#'   the smallest group-aligned value `>= n_subsamples`).  A target that
#'   cannot be reached from the available groups is an error.
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
#' @param seed Optional integer; seeds this call for reproducibility.  The
#'   caller's RNG state is saved and restored, so a seeded call does not
#'   advance the global stream.  `NULL` draws from the ambient stream.  See the
#'   *Seeding in a loop* section of [classic_bootstrap_indices()] for the
#'   consequence of reusing one seed across many calls.
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
  if (length(groups) != length(y)) {
    stop("`groups` must be the same length as `y`.", call. = FALSE)
  }
  .validate_group_ids(groups)

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

  # Whole-group invariant (D3): never trim partial groups, even if the
  # final tally exceeds n_subsamples.

  sampled_idx <- .with_local_seed(seed, {
    idx <- draw_once()

    .guard_diverse_bootstrap(
      y, idx, draw_once,
      caller = "group_bootstrap_indices",
      hint   = "Check class balance or group structure."
    )
  })

  sampled_idx
}

.guard_diverse_bootstrap <- function(y, idx, draw_once, caller, hint,
                                     max_retries = 1000L) {
  if (length(unique(y)) < 2L) return(idx)
  attempt <- 0L
  while (length(unique(y[idx])) < 2L) {
    attempt <- attempt + 1L
    if (attempt > max_retries) {
      .abort_stablr(
        "stablr_bootstrap_infeasible",
        paste0(
          caller, ": could not draw a class-diverse subsample after ",
          max_retries, " attempts. ", hint
        ),
        caller = caller,
        max_retries = max_retries
      )
    }
    idx <- draw_once()
  }
  idx
}

# A group label that matches no row (`NA`/`NaN`) contributes an empty chunk to
# `.sample_groups_until()`, so the group pool can drain without ever reaching the
# row target.  Rejecting it here names the cause instead of letting the draw
# silently come up short.
.validate_group_ids <- function(groups) {
  if (anyNA(groups)) {
    stop(
      "`groups` cannot contain missing values: a group label with no matching ",
      "row contributes no samples to the bootstrap subsample.",
      call. = FALSE
    )
  }
  invisible(TRUE)
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
  .sample_groups_until(group_levels, groups, n_subsamples, replace,
                       label = "Grouped bootstrap")
}

# Shared group-sampling while-loop: draw groups until target_n indices are
# accumulated; use index-then-subset to avoid R's sample(x,1) length-1 pitfall.
#
# `target_n` is a lower bound, not an exact count: whole groups are never
# trimmed, so the loop legitimately overshoots.  Coming up *short* is a
# different matter -- the group pool ran dry -- and is an error rather than a
# silently smaller subsample, because nothing downstream records the realised
# length and the caller would report the requested target as if it were met.
.sample_groups_until <- function(group_levels, groups, target_n, replace,
                                 label = "Grouped bootstrap") {
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
  realised <- if (!replace) length(unique(idx)) else length(idx)
  if (realised < target_n) {
    .abort_stablr(
      "stablr_bootstrap_infeasible",
      paste0(
        label, " could only reach ", realised,
        " of the requested ", target_n, " rows from ",
        length(group_levels), " group(s). Lower `n_subsamples` (or raise ",
        "`sample_fraction` granularity) or allow `replace = TRUE`."
      ),
      target_n = target_n,
      realised = realised,
      n_groups = length(group_levels)
    )
  }
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
      group_levels[group_strata == stratum],
      groups,
      target_counts[[stratum]],
      replace,
      label = paste0("Grouped stratified bootstrap for stratum `", stratum, "`")
    )
    outer_chunks[[length(outer_chunks) + 1L]] <- stratum_idx
  }

  sampled_idx <- unlist(outer_chunks, use.names = FALSE)
  if (!replace) sampled_idx <- unique(sampled_idx)

  sample(sampled_idx, length(sampled_idx), replace = FALSE)
}
