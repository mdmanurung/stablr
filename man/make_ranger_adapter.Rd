% Generated manually for the stablr improvement plan (item D1)
\name{make_ranger_adapter}
\alias{make_ranger_adapter}
\title{Build a Random-Forest Learner Adapter (Calibrated, Scoped)}
\usage{
make_ranger_adapter(num_trees = 500L)
}
\arguments{
\item{num_trees}{Positive integer; number of trees per ranger forest.}
}
\value{
A function with signature
\code{function(x, y, lambda_val) -> logical vector of length ncol(x)},
mirroring \code{\link[=make_glmnet_adapter]{make_glmnet_adapter()}}: the
mask selects features whose permutation importance meets the threshold
carried by \code{lambda_val}.
}
\description{
Learner adapter extending STABL to non-linear signals with random forests
(\code{ranger}).  One forest is fit per bootstrap subsample;
the "penalty path" is an importance-threshold path built from quantiles of
the pooled permutation-importance distribution (see
\code{\link[=auto_rf_threshold_grid]{auto_rf_threshold_grid()}}), so the
selection set grows monotonically as the threshold decreases — the same
sparse-to-dense nesting the lambda path provides for glmnet learners.
}
\details{
\strong{Calibrated within a measured scope (WS-2, 60 simulation runs):}
with \code{knockoff_equi} decoys the mean empirical FDP is 0.025 on a
linear design and 0.065 on a low-dimensional nonlinear design (p = 10),
meeting the predeclared criterion (FDP <= 0.15, TPR >= 0.8). The criterion
is NOT met on a high-dimensional correlated design at p = 50 (mean FDP
0.44) or for binary outcomes (0.26); avoid this pairing there. This is
simulation evidence, not a guarantee for your dataset.

\strong{Statistical guardrails:} permutation importance is used (not Gini).
Random forests must be paired with model-X knockoff artificial features
(\code{knockoff_equi} / \code{knockoff_mvr}): permuted nulls do not receive
the correlation leakage that real null features receive in a multivariate
learner, so \code{random_permutation} nulls make FDP+ anti-conservative —
stablr warns when this combination is requested.
}
\examples{
if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(x[, 1] * x[, 2] * 0.5 + rnorm(60, sd = 0.3), rownames(x))
  adapter <- make_ranger_adapter(num_trees = 100L)
  grid <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  mask <- adapter(x, y, grid[1L, ])  # sparsest end of the path
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
}
\seealso{
\code{\link[=auto_rf_threshold_grid]{auto_rf_threshold_grid()}},
\code{\link[=make_glmnet_adapter]{make_glmnet_adapter()}},
\code{\link[=stabl_fit]{stabl_fit()}} (\code{base_learner = "rf"}).
}
