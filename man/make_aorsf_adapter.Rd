% Generated manually for the stablr improvement plan (item D2)
\name{make_aorsf_adapter}
\alias{make_aorsf_adapter}
\title{Build an Oblique Random-Forest Learner Adapter (Fenced)}
\usage{
make_aorsf_adapter(num_trees = 500L)
}
\arguments{
\item{num_trees}{Positive integer; number of trees per aorsf oblique forest.}
}
\value{
A function with signature
\code{function(x, y, lambda_val) -> logical vector of length ncol(x)},
mirroring \code{\link[=make_glmnet_adapter]{make_glmnet_adapter()}}: the
mask selects features whose negation importance meets the threshold carried
by \code{lambda_val}.
}
\description{
Learner adapter extending STABL to non-linear signals with
oblique random forests (\code{aorsf}).  One forest is fit per bootstrap
subsample; the "penalty path" is an importance-threshold path built from
quantiles of the pooled negation-importance distribution (see
\code{\link[=auto_oblique_threshold_grid]{auto_oblique_threshold_grid()}}),
so the selection set grows monotonically as the threshold decreases — the
same sparse-to-dense nesting the lambda path provides for glmnet learners.
}
\details{
\strong{Fenced (WS-3, 40 simulation runs):} no predeclared remedy reached
the calibration criterion (mean FDP <= 0.15 on both designs; best remedy
0.18 / 0.39). Knockoff decoy importance is inflated for oblique forests, so
FDP+ is anti-conservative. Use \code{base_learner = "rf"} with
\code{knockoff_equi} for knockoff-calibrated selection instead.

\strong{Statistical guardrails:} negation importance is used (not permutation
importance) because it discriminates relevant from irrelevant predictors
better in oblique forests (Jaeger et al. 2023).  Oblique random forests must
be paired with model-X knockoff artificial features (\code{knockoff_equi} /
\code{knockoff_mvr}): permuted nulls do not receive the correlation leakage
that real null features receive in a multivariate learner, so
\code{random_permutation} nulls make FDP+ anti-conservative — stablr warns
when this combination is requested.

\strong{D2 caveat:} importance is spread across the linear combinations used
in each split, so per-feature selection frequencies are noisier than for
axis-aligned random forests (\code{\link[=make_ranger_adapter]{make_ranger_adapter()}}); small
stability-score differences between features are less reliable and more
bootstraps are recommended.
}
\examples{
if (requireNamespace("aorsf", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(0.8 * x[, 1] + 0.8 * x[, 2] + 0.5 * x[, 1] * x[, 2] +
                rnorm(60, sd = 0.3), rownames(x))
  adapter <- make_aorsf_adapter(num_trees = 100L)
  grid <- auto_oblique_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  mask <- adapter(x, y, grid[1L, ])  # sparsest end of the path
  cat("selected:", sum(mask), "of", ncol(x), "features\n")
}
}
\seealso{
\code{\link[=auto_oblique_threshold_grid]{auto_oblique_threshold_grid()}},
\code{\link[=make_ranger_adapter]{make_ranger_adapter()}},
\code{\link[=make_glmnet_adapter]{make_glmnet_adapter()}},
\code{\link[=stabl_fit]{stabl_fit()}} (\code{base_learner = "oblique_rf"}).
}
