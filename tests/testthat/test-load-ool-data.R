# E1: Tests for load_ool_data() — bundled extdata fixtures must always be present.

test_that("E1: load_ool_data train returns three-element list with correct names", {
  ool <- load_ool_data(split = "train")

  expect_type(ool, "list")
  expect_named(ool, c("x_list", "y", "ids"))
  expect_named(ool$x_list, c("cytof", "proteomics"))
})

test_that("E1: load_ool_data train matrices have 150 rows and 100 columns", {
  ool <- load_ool_data(split = "train")

  expect_equal(nrow(ool$x_list$cytof),        150L)
  expect_equal(ncol(ool$x_list$cytof),        100L)
  expect_equal(nrow(ool$x_list$proteomics),   150L)
  expect_equal(ncol(ool$x_list$proteomics),   100L)
  expect_length(ool$y,                        150L)
  expect_length(ool$ids,                      150L)
})

test_that("E1: load_ool_data train rows are aligned across omics and outcome", {
  ool <- load_ool_data(split = "train")

  # Row names of both matrices and names of y must all equal ids
  expect_equal(rownames(ool$x_list$cytof),      ool$ids)
  expect_equal(rownames(ool$x_list$proteomics), ool$ids)
  expect_equal(names(ool$y),                    ool$ids)
})

test_that("E1: load_ool_data train matrices are numeric with no NAs", {
  ool <- load_ool_data(split = "train")

  expect_true(is.numeric(ool$x_list$cytof))
  expect_true(is.numeric(ool$x_list$proteomics))
  expect_false(anyNA(ool$x_list$cytof))
  expect_false(anyNA(ool$x_list$proteomics))
  expect_false(anyNA(ool$y))
})

test_that("E1: load_ool_data valid split loads without error and is non-empty", {
  ool_v <- load_ool_data(split = "valid")

  expect_named(ool_v, c("x_list", "y", "ids"))
  expect_gt(length(ool_v$ids), 0L)
  # Rows are aligned
  expect_equal(rownames(ool_v$x_list$cytof),      ool_v$ids)
  expect_equal(rownames(ool_v$x_list$proteomics), ool_v$ids)
  expect_equal(names(ool_v$y),                    ool_v$ids)
})

test_that("E1: load_ool_data rejects invalid split argument", {
  expect_error(load_ool_data(split = "test"), regexp = "arg")
})

test_that("load_ool_data fails closed on duplicate sample IDs (A16 regression)", {
  dup_df <- data.frame(id = c("a", "a", "b"), f1 = 1:3)
  testthat::with_mocked_bindings(
    read.csv = function(...) dup_df,
    expect_error(load_ool_data("train"), "Duplicate sample IDs")
  )
})

test_that("load_ool_data fails closed on non-numeric feature columns (A16 regression)", {
  bad_df <- data.frame(id = c("a", "b"), f1 = c("x", "y"),
                       stringsAsFactors = FALSE)
  testthat::with_mocked_bindings(
    read.csv = function(...) bad_df,
    expect_error(load_ool_data("train"), "Non-numeric")
  )
})
