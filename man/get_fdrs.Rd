% Generated manually for the stablr improvement plan (item C3)
\name{get_fdrs}
\alias{get_fdrs}
\alias{get_fdrs.stabl_fit}
\title{Get FDP+ Diagnostic Estimates Across the Swept Thresholds}
\usage{
get_fdrs(object)
}
\arguments{
\item{object}{A fitted \code{"stabl_fit"} object returned by \code{\link[=stabl_fit]{stabl_fit()}}.}
}
\value{
Named numeric vector with one FDP+ estimate per swept stability
threshold (names taken from the fit's \code{fdr_threshold_range}).
}
\description{
Returns the vector of FDP+ diagnostic estimates computed by
\code{\link[=stabl_fit]{stabl_fit()}} across the swept stability thresholds.
}
\details{
These are \emph{diagnostic} estimates used to choose the stability
threshold, \bold{not} a guarantee that false discovery proportion is
controlled at any level.  Calibration depends on the artificial features
acting as fair null competitors.  Fails with an error when the fit was
produced without artificial features (\code{artificial_type = NULL}).
}
\examples{
set.seed(1L)
x <- matrix(rnorm(40 * 6), 40, 6,
            dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
y <- setNames(rnorm(40), rownames(x))
fit <- stabl_fit(x, y,
                 lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
                 n_bootstraps   = 4L,
                 artificial_type = "random_permutation",
                 random_state   = 1L)
fdrs <- get_fdrs(fit)  # named numeric, one FDP+ estimate per threshold
fdrs
}
\seealso{
\code{\link[=get_stabl_scores]{get_stabl_scores()}} for the full selection-frequency path,
\code{\link[=get_lambda_grid]{get_lambda_grid()}} for the penalty grid used.
}
