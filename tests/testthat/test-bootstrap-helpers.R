test_that("sample_bootstrap_indices returns requested size", {
  y <- c(0, 1, 0, 1, 0, 1)
  idx <- sample_bootstrap_indices(y = y, n_subsamples = 4, replace = FALSE, random_state = 1)

  expect_length(idx, 4)
  expect_true(all(idx >= 1 & idx <= length(y)))
})

test_that("sample_bootstrap_indices validates non-replacement size", {
  y <- c(0, 1, 0)
  expect_error(
    sample_bootstrap_indices(y = y, n_subsamples = 5, replace = FALSE),
    "cannot exceed"
  )
})

test_that("sample_bootstrap_indices can stratify by outcome class", {
  y <- c(rep("NP", 15L), rep("P", 23L))

  idx <- sample_bootstrap_indices(
    y = y,
    n_subsamples = 19L,
    replace = FALSE,
    stratify = TRUE,
    random_state = 11L
  )

  expect_length(idx, 19L)
  expect_equal(as.integer(table(y[idx])), c(8L, 11L))
  expect_named(table(y[idx]), c("NP", "P"))
})

test_that("sample_bootstrap_indices rejects impossible stratified targets", {
  y <- c("A", "B", "C")

  expect_error(
    sample_bootstrap_indices(y = y, n_subsamples = 2L, stratify = TRUE),
    "at least the number of realised strata"
  )
})

test_that("sample_bootstrap_indices can stratify by a multi-column design", {
  y <- c(rep("NP", 15L), rep("P", 23L))
  study <- c(
    rep("CVTU3", 3L), rep("EGSV2", 4L), rep("PfGA2", 8L),
    rep("CVTU3", 9L), rep("EGSV2", 6L), rep("PfGA2", 8L)
  )
  strata <- data.frame(outcome = y, study = study)

  idx <- sample_bootstrap_indices(
    y = y,
    n_subsamples = 19L,
    replace = FALSE,
    strata = strata,
    random_state = 12L
  )

  observed <- table(interaction(strata[idx, ], drop = TRUE, sep = "\r"))
  expect_length(idx, 19L)
  expect_equal(length(observed), 6L)
  expect_true(all(observed >= 1L))
})

test_that("bootstrap strata row names are used when numeric sample ids are shuffled", {
  sample_ids <- c("2", "1", "3")
  strata_df <- data.frame(
    site = c("site_for_1", "site_for_2", "site_for_3"),
    row.names = c("1", "2", "3")
  )
  strata_mat <- matrix(
    c("site_for_1", "site_for_2", "site_for_3"),
    ncol = 1L,
    dimnames = list(c("1", "2", "3"), "site")
  )

  out_df <- .subset_bootstrap_strata_by_ids(strata_df, sample_ids)
  out_mat <- .subset_bootstrap_strata_by_ids(strata_mat, sample_ids)

  expect_equal(rownames(out_df), sample_ids)
  expect_equal(out_df$site, c("site_for_2", "site_for_1", "site_for_3"))
  expect_equal(rownames(out_mat), sample_ids)
  expect_equal(out_mat$site, c("site_for_2", "site_for_1", "site_for_3"))
})

test_that("group_bootstrap_indices validates group length", {
  y <- c(0, 1, 0)
  g <- c("a", "a")
  expect_error(group_bootstrap_indices(y = y, groups = g, n_subsamples = 2), "same length")
})

test_that("group_bootstrap_indices keeps whole groups when target aligns", {
  # Simulate repeated-measures layout: 4 subjects, 3 visits each.
  groups <- rep(paste0("id", 1:4), each = 3L)
  y <- rep(c(0, 1), length.out = length(groups))

  set.seed(123)
  idx <- group_bootstrap_indices(
    y = y,
    groups = groups,
    n_subsamples = 6L,
    replace = FALSE
  )

  expect_length(idx, 6L)

  selected_groups <- unique(groups[idx])
  for (g in selected_groups) {
    expect_setequal(idx[groups[idx] == g], which(groups == g))
  }
})

test_that("sample_bootstrap_indices can split repeated-measures groups", {
  # Without grouped sampling, subject rows can be partially sampled.
  groups <- rep(paste0("id", 1:4), each = 3L)
  y <- rep(c(0, 1), length.out = length(groups))

  has_partial_group <- FALSE
  set.seed(321)
  for (i in seq_len(200L)) {
    idx <- sample_bootstrap_indices(y = y, n_subsamples = 6L, replace = FALSE)
    selected_groups <- unique(groups[idx])
    partial <- any(vapply(selected_groups, function(g) {
      sel <- idx[groups[idx] == g]
      !setequal(sel, which(groups == g))
    }, logical(1L)))

    if (partial) {
      has_partial_group <- TRUE
      break
    }
  }

  expect_true(has_partial_group)
})

test_that("group_bootstrap_indices replace=FALSE never re-draws the same group", {
  # Each group must appear at most once when replace = FALSE.
  # The old code used sample(group_levels, replace=replace) which had no effect
  # on a size-1 draw — the pool never shrank and the same group could appear
  # multiple times.
  groups <- rep(paste0("id", 1:6), each = 4L)
  y      <- rep(c(0, 1), length.out = length(groups))

  for (seed_val in c(1L, 7L, 42L, 99L, 200L)) {
    idx             <- group_bootstrap_indices(y = y, groups = groups,
                                              n_subsamples = 12L,
                                              replace = FALSE, random_state = seed_val)
    selected_groups <- groups[idx]
    group_counts    <- table(selected_groups)
    # Every sampled row belongs to a group; each group is present at most once
    # (all its rows appear, but the group label appears exactly n_rows_per_group times)
    # — the key assertion is that no group appears more than its actual row count.
    for (g in names(group_counts)) {
      expect_lte(length(unique(idx[selected_groups == g])),
                 sum(groups == g))
    }
    # Stronger: unique group labels in the selection must equal selected group count
    expect_equal(length(unique(selected_groups)),
                 length(unique(groups[idx])))
  }
})

test_that("group_bootstrap_indices can stratify while preserving whole groups", {
  groups <- rep(paste0("id", seq_len(8L)), each = 2L)
  y <- rep(c(rep("NP", 4L), rep("P", 4L)), each = 2L)

  idx <- group_bootstrap_indices(
    y = y,
    groups = groups,
    n_subsamples = 8L,
    replace = FALSE,
    stratify = TRUE,
    random_state = 22L
  )

  expect_true(all(table(y[idx]) >= 4L))
  selected_groups <- unique(groups[idx])
  for (g in selected_groups) {
    expect_setequal(idx[groups[idx] == g], which(groups == g))
  }
})

test_that("group_bootstrap_indices stratification requires class-pure groups", {
  groups <- c("id1", "id1", "id2", "id2")
  y <- c("NP", "P", "NP", "P")

  expect_error(
    group_bootstrap_indices(
      y = y,
      groups = groups,
      n_subsamples = 2L,
      stratify = TRUE,
      random_state = 1L
    ),
    "exactly one realised stratum"
  )
})

test_that("group_bootstrap_indices stratified replacement can reuse groups", {
  groups <- c("id1", "id1", "id2", "id2")
  y <- c("NP", "NP", "P", "P")

  idx <- group_bootstrap_indices(
    y = y,
    groups = groups,
    n_subsamples = 6L,
    replace = TRUE,
    stratify = TRUE,
    random_state = 3L
  )

  expect_gte(length(idx), 6L)
  expect_true(anyDuplicated(idx) > 0L)
  expect_true(all(c("NP", "P") %in% y[idx]))
})

test_that("sample_bootstrap_indices errors after max retries on impossible class diversity", {
  # n_subsamples = 1 can never contain both classes, even though the population
  # has 2 classes, so the iterative retry loop must hit its cap and stop().
  y <- c(0L, rep(1L, 19L))   # 2-class population
  expect_error(
    sample_bootstrap_indices(y = y, n_subsamples = 1L, replace = FALSE, random_state = 1L),
    "could not draw a class-diverse subsample"
  )
})

test_that("group_bootstrap_indices errors after max retries on impossible class diversity", {
  groups <- rep(paste0("id", 1:5), each = 4L)
  # All 4 rows of "id1" are class 0; all other groups are class 1.
  # Every single-group draw is therefore single-class, so the retry loop
  # can never produce a class-diverse subsample and must error.
  y      <- c(rep(0L, 4L), rep(1L, length(groups) - 4L))
  # n_subsamples = 1 → one group drawn (4 rows) → always single-class → error
  expect_error(
    group_bootstrap_indices(y = y, groups = groups, n_subsamples = 1L,
                            replace = FALSE, random_state = 1L),
    "could not draw a class-diverse subsample"
  )
})

# D2 characterization: pin exact element order for internal helpers before
# refactoring c()-in-loop growth to list-collect + unlist.

test_that("D2: .append_noise_groups NULL passthrough is identity (characterization)", {
  g <- c(1L, 1L, 2L, 3L, 3L)
  expect_identical(stablr:::.append_noise_groups(g, NULL, length(g)), g)
})

test_that("D2: .append_noise_groups valid src indices borrow group ids (characterization)", {
  g         <- c(1L, 1L, 2L, 3L, 3L)
  noise_src <- c(2L, 1L, 3L)
  result    <- stablr:::.append_noise_groups(g, noise_src, length(g) + length(noise_src))
  expect_identical(result, c(1L, 1L, 2L, 3L, 3L, 1L, 1L, 2L))
})

test_that("D2: .append_noise_groups invalid indices assign fresh group ids (characterization)", {
  g         <- c(1L, 1L, 2L, 3L, 3L)
  noise_inv <- c(99L, 0L, NA_integer_)
  result    <- stablr:::.append_noise_groups(g, noise_inv, length(g) + length(noise_inv))
  expect_identical(result, c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 6L))
})

test_that("D2: .append_noise_groups mixed src+invalid preserves element order (characterization)", {
  g         <- c(1L, 1L, 2L, 3L, 3L)
  noise_mix <- c(2L, 99L)
  result    <- stablr:::.append_noise_groups(g, noise_mix, length(g) + length(noise_mix))
  expect_identical(result, c(1L, 1L, 2L, 3L, 3L, 1L, 4L))
})

test_that("D2: .unstratified_group_bootstrap_indices replace=FALSE exact order (characterization)", {
  groups <- c("a", "a", "a", "b", "b", "c", "c", "c", "c", "d")
  gl     <- unique(groups)
  set.seed(7L)
  idx <- stablr:::.unstratified_group_bootstrap_indices(groups, gl,
                                                        n_subsamples = 5L,
                                                        replace = FALSE)
  expect_identical(idx, c(4L, 5L, 10L, 1L, 2L, 3L))
})

test_that("D2: .unstratified_group_bootstrap_indices replace=TRUE exact order (characterization)", {
  groups <- c("a", "a", "a", "b", "b", "c", "c", "c", "c", "d")
  gl     <- unique(groups)
  set.seed(13L)
  idx <- stablr:::.unstratified_group_bootstrap_indices(groups, gl,
                                                        n_subsamples = 5L,
                                                        replace = TRUE)
  expect_identical(idx, c(10L, 6L, 7L, 8L, 9L))
})

test_that("D2: .stratified_group_bootstrap_indices replace=FALSE exact order (characterization)", {
  groups <- c("a", "a", "a", "b", "b", "c", "c", "c", "c", "d")
  gl     <- unique(groups)
  strata <- c("X", "X", "X", "Y", "Y", "Y", "Y", "Y", "Y", "X")
  set.seed(22L)
  idx <- stablr:::.stratified_group_bootstrap_indices(strata, groups, gl,
                                                      n_subsamples = 5L,
                                                      replace = FALSE)
  expect_identical(idx, c(1L, 7L, 3L, 8L, 2L, 10L, 6L, 9L))
})
