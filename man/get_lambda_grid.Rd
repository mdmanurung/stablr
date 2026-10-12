% Generated manually for the stablr improvement plan (item C3)
\name{get_lambda_grid}
\alias{get_lambda_grid}
\alias{get_lambda_grid.stabl_fit}
\title{Get the Fitted Lambda Grid}
\usage{
get_lambda_grid(object)
}
\arguments{
\item{object}{A fitted \code{"stabl_fit"} object returned by \code{\link[=stabl_fit]{stabl_fit()}}.}
}
\value{
The \code{data.frame} penalty grid actually used by the fit (one row
per lambda, with columns \code{lambda} and, for elastic net, \code{alpha}).
}
\description{
Returns the penalty grid that the fit consumed, after any
\code{"auto"} resolution.  Useful for reproducing a fit exactly or for
inspecting auto-generated grids.
}
\examples{
set.seed(1L)
x <- matrix(rnorm(30 * 5), 30, 5,
            dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
y <- setNames(rnorm(30), rownames(x))
fit <- stabl_fit(x, y,
                 lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
                 n_bootstraps   = 4L,
                 hard_threshold = 0.3,
                 random_state   = 1L)
get_lambda_grid(fit)
}
\seealso{
\code{\link[=get_fdrs]{get_fdrs()}} for the FDP+ diagnostics across the swept thresholds.
}
