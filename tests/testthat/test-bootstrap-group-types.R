# Audit findings H-1 (V-1, V-7) and Q5 (whole-group invariant, D3).
# These tests pin two invariants of group_bootstrap_indices:
#   1. The function never fabricates a group label, regardless of the storage
#      type of `groups` (integer, numeric, character, factor).  This guards
#      against the R `sample(x, 1L)` length-1 numeric pitfall: when
#      length(remaining) == 1 and remaining is numeric, sample(remaining, 1L)
#      silently behaves as sample.int(remaining, 1L).
#   2. Whole groups are never split between the bootstrap subsample and its
#      complement (D3 decision: drop the trim policy).  Per STABL.md §Step 2
#      this is the strongest leakage-prevention contract.

.assert_no_fabricated_labels <- function(idx, groups, label) {
  observed <- unique(groups[idx])
  expected_universe <- unique(groups)
  expect_true(
    all(observed %in% expected_universe),
    label = paste0("[", label, "] sampled rows reference only existing group labels")
  )
}

.assert_whole_groups_only <- function(idx, groups, label) {
  for (g in unique(groups[idx])) {
    rows_in_g          <- which(groups == g)
    rows_in_g_sampled  <- intersect(idx, rows_in_g)
    expect_setequal(rows_in_g_sampled, rows_in_g)
  }
}

# ----- Parametrised by group-label storage type --------------------------------
group_specs <- list(
  list(name = "integer",   make = function() rep(1:6, each = 5L)),
  list(name = "numeric",   make = function() rep(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), each = 5L)),
  list(name = "character", make = function() rep(paste0("g", 1:6), each = 5L)),
  list(name = "factor",    make = function() factor(rep(paste0("g", 1:6), each = 5L)))
)

for (spec in group_specs) {
  local({
    nm     <- spec$name
    maker  <- spec$make

    test_that(paste0("group_bootstrap_indices never fabricates labels (", nm, " groups)"), {
      groups <- maker()
      y      <- rep(c(0L, 1L), length.out = length(groups))

      # Hammer many seeds; pre-fix `sample(remaining, 1L)` triggers the pitfall
      # on numeric/integer groups whenever the pool shrinks to length 1.
      for (s in seq_len(50L)) {
        idx <- group_bootstrap_indices(y = y, groups = groups,
                                       n_subsamples = 18L,
                                       replace = FALSE, random_state = s)
        .assert_no_fabricated_labels(idx, groups, paste(nm, "seed", s))
      }
    })

    test_that(paste0("group_bootstrap_indices preserves whole groups (", nm, " groups)"), {
      groups <- maker()
      y      <- rep(c(0L, 1L), length.out = length(groups))

      for (s in seq_len(20L)) {
        idx <- group_bootstrap_indices(y = y, groups = groups,
                                       n_subsamples = 18L,
                                       replace = FALSE, random_state = s)
        .assert_whole_groups_only(idx, groups, paste(nm, "seed", s))
      }
    })
  })
}

# ----- Regression guard: explicit numeric-pitfall reproduction ----------------
test_that("group_bootstrap_indices does not fabricate labels when the eligible pool shrinks to 1 (numeric)", {
  # Three large numeric-labeled groups.  With n_subsamples == sum-of-two-groups
  # the third draw is forced; right before that draw `remaining` has length 1
  # and is numeric — exactly the configuration where the pre-fix code returns
  # sample.int(remaining, 1L).
  groups <- rep(c(10, 20, 30), each = 4L)
  y      <- rep(c(0L, 1L), length.out = length(groups))

  for (s in seq_len(200L)) {
    idx <- group_bootstrap_indices(y = y, groups = groups,
                                   n_subsamples = 12L,  # = total: must hit all 3
                                   replace = FALSE, random_state = s)
    .assert_no_fabricated_labels(idx, groups, paste("forced-singleton seed", s))
  }
})

# ----- Regression guard: whole-group invariant when n_subsamples is mid-group -
test_that("group_bootstrap_indices does not split groups when target falls inside a group", {
  # 5 groups of 4 rows each; target 6 = group-of-4 + half-of-next.  Pre-fix
  # code would trim the second group to 2 rows; post-fix keeps it whole.
  groups <- rep(paste0("g", 1:5), each = 4L)
  y      <- rep(c(0L, 1L), length.out = length(groups))

  for (s in seq_len(50L)) {
    idx <- group_bootstrap_indices(y = y, groups = groups,
                                   n_subsamples = 6L,
                                   replace = FALSE, random_state = s)
    .assert_whole_groups_only(idx, groups, paste("mid-group seed", s))
    # Realised count must be a multiple of group size (4), and >= target
    expect_true(length(idx) %% 4L == 0L)
    expect_gte(length(idx), 6L)
  }
})
