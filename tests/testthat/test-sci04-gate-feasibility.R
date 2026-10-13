# Audit E1 (0.2.2): deterministic reference tests for the SCI-04 v2 draft
# feasibility arithmetic.  The draft proposes the distribution-free one-sided
# Hoeffding upper bound U = min(1, mean(X) + sqrt(log(1/alpha) / (2*n))) with
# alpha = 0.05 and release criterion U <= 0.10.  The criterion is achievable
# only if the margin alone is below 0.10; at the previously considered
# n = 100 it is not (margin 0.1224 > 0.10), making the gate vacuous.  These
# tests pin the arithmetic from the formula and predeclared settings only --
# no locked-run data is involved.

.hoeffding_margin <- function(n, alpha = 0.05) {
  sqrt(log(1 / alpha) / (2 * n))
}

test_that("SCI-04 draft: Hoeffding margin formula matches the contract", {
  expect_equal(.hoeffding_margin(100), 0.1224, tolerance = 1e-3)
  expect_equal(.hoeffding_margin(150), 0.0999, tolerance = 1e-3)
  expect_equal(.hoeffding_margin(200), 0.0865, tolerance = 1e-3)
  expect_equal(.hoeffding_margin(250), 0.0774, tolerance = 1e-3)
})

test_that("SCI-04 draft: n = 100 makes the U <= 0.10 criterion vacuous", {
  # Even a perfect run (mean(X) = 0) exceeds the criterion at n = 100.
  u_perfect <- min(1, 0 + .hoeffding_margin(100))
  expect_gt(u_perfect, 0.10)
})

test_that("SCI-04 draft: the predeclared n must be >= 150 for feasibility", {
  expect_lt(.hoeffding_margin(150), 0.10)
  expect_lt(.hoeffding_margin(200), 0.10)
  # The feasibility boundary: n = 149 is still infeasible, n = 150 is not.
  expect_gt(.hoeffding_margin(149), 0.10)
})

test_that("SCI-04 draft: contract file carries the feasibility clause", {
  contract_path <- system.file(
    "analysis", "contracts", "SCI-04-null-selected-fraction-v2.md",
    package = "stablr"
  )
  skip_if(identical(contract_path, ""), "contract file not installed")
  text <- paste(readLines(contract_path), collapse = "\n")
  expect_match(text, "Feasibility of the release criterion", fixed = TRUE)
  expect_match(text, "n >= 150", fixed = TRUE)
  expect_match(text, "0.1224", fixed = TRUE)
})
