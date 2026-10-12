% Generated manually for the stablr improvement plan (item D1)
\name{auto_rf_threshold_grid}
\alias{auto_rf_threshold_grid}
\title{Auto-Generate an Importance-Threshold Grid for Random-Forest STABL (Experimental)}
\usage{
auto_rf_threshold_grid(x, y, n_grid = 30L, num_trees = 500L, seed = NULL)
}
\arguments{
\item{x}{A \code{matrix} or \code{data.frame} of features (samples in rows).}
\item{y}{Outcome vector.}
\item{n_grid}{Positive integer; number of grid points.}
\item{num_trees}{Positive integer; number of trees per ranger forest.}
\item{seed}{Optional integer seed for the forest fit.}
}
\value{
A \code{data.frame} with columns \code{lambda} (the importance threshold)
and \code{threshold} (identical), sorted in \emph{decreasing} order so the
path runs sparse-to-dense, matching the lambda-path convention used by
\code{\link[=stabl_fit]{stabl_fit()}}.
}
\description{
Builds the random-forest analog of \code{\link[=auto_lambda_grid]{auto_lambda_grid()}}:
one ranger forest with permutation importance is fit on the full data, and
the grid is taken from quantiles of the pooled importance distribution
(\code{probs = seq(0.5, 0.999, length.out = n_grid)}).  Experimental
extension beyond Python STABL.
}
\examples{
if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(60 * 6), 60, 6,
              dimnames = list(paste0("s", 1:60), paste0("f", 1:6)))
  y <- setNames(x[, 1] * x[, 2] * 0.5 + rnorm(60, sd = 0.3), rownames(x))
  grid <- auto_rf_threshold_grid(x, y, n_grid = 5L, num_trees = 100L, seed = 1L)
  grid
}
}
\seealso{
\code{\link[=make_ranger_adapter]{make_ranger_adapter()}},
\code{\link[=auto_lambda_grid]{auto_lambda_grid()}}.
}
