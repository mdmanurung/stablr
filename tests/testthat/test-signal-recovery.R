# Audit V-5: with a clean signal (3 strong features in p=30, n=200) STABL
# must recover the planted set with a compact FDP+-thresholded support.

test_that("stabl_fit recovers a 3-feature gaussian signal with compact support", {
  skip_on_cran()
  withr::local_seed(7)

  n <- 200L; p <- 30L
  x <- matrix(rnorm(n * p), nrow = n,
               dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  beta <- numeric(p); beta[1:3] <- c(2.0, -1.5, 1.0)
  y <- setNames(as.numeric(x %*% beta + rnorm(n, 0, 0.3)), rownames(x))

  fit <- stabl_fit(x, y, lambda_grid = "auto", n_lambda = 20L,
                   n_bootstraps = 100L,
                   artificial_type = "random_permutation",
                   random_state = 7L)

  selected <- get_feature_names(fit)
  # All three planted features must be selected; allow a small number of
  # extras to accommodate bootstrap/path variability.
  expect_true(all(c("f1", "f2", "f3") %in% selected))
  expect_lte(length(selected), 8L)
})
