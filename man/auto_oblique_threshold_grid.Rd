% Generated manually for the stablr improvement plan (item D2)
\name{auto_oblique_threshold_grid}
\alias{auto_oblique_threshold_grid}
\title{Auto-Generate an Importance-Threshold Grid for Oblique-Random-Forest STABL (Experimental)}
\usage{
auto_oblique_threshold_grid(x, y, n_grid = 30L, num_trees = 500L, seed = NULL)
}
\arguments{
\item{x}{A \code{matrix} or \code{data.frame} of features (samples in rows).}
\item{y}{Outcome vector.}
\item{n_grid}{Positive integer; number of grid points.}
\item{num_trees}{Positive integer; number of trees per aorsf oblique forest.}
\item{seed}{Optional integer seed for the forest fit.}
}
\value{
A \code{data.frame} with columns \code{lambda} (the importance threshold)
and \code{threshold} (identical), sorted in \emph{decreasing} order so the
path runs sparse-to-dense, matching the lambda-path convention used by
\code{\link[=stabl_fit]{stabl_fit()}}.
}
\description{
Builds the oblique-random-forest analog of
\code{\link[=auto_rf_threshold_grid]{auto_rf_threshold_grid()}}: one aorsf
oblique forest with negation importance is fit on the full data, and the
grid is taken from quantiles of the pooled importance distribution
(\code{probs = seq(0.5, 0.999, length.out = n_grid)}).  Experimental
extension beyond Python STABL.
}
\details{
Negation importance is used instead of permutation importance because it
discriminates relevant from irrelevant predictors better in oblique forests
(Jaeger et al. 2023).  Because importance is spread across the linear
combinations used in each split, the resulting per-feature selection
frequencies are noisier than for axis-aligned random forests — prefer more
bootstraps when using \code{base_learner = "oblique_rf"} in
\code{\link[=stabl_fit]{stabl_fit()}}.
}
\examples{
if (requireNamespace("aorsf", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(0.8 * x[, 1] + 0.8 * x[, 2] + 0.5 * x[, 1] * x[, 2] +
                rnorm(60, sd = 0.3), rownames(x))
  grid <- auto_oblique_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  grid
}
}
\seealso{
\code{\link[=make_aorsf_adapter]{make_aorsf_adapter()}},
\code{\link[=auto_rf_threshold_grid]{auto_rf_threshold_grid()}},
\code{\link[=auto_lambda_grid]{auto_lambda_grid()}}.
}
