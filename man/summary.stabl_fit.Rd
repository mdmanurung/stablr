% Generated manually for the stablr improvement plan (item C4)
\name{summary.stabl_fit}
\alias{summary.stabl_fit}
\alias{print.summary.stabl_fit}
\title{Summarize a Fitted stabl_fit Object}
\usage{
\method{summary}{stabl_fit}(object, ...)

\method{print}{summary.stabl_fit}(x, ...)
}
\arguments{
\item{object}{A fitted \code{"stabl_fit"} object returned by \code{\link[=stabl_fit]{stabl_fit()}}.}
\item{x}{A \code{"summary.stabl_fit"} object as returned by \code{summary()}.}
\item{...}{Unused; present for S3 consistency.}
}
\value{
\code{summary()} returns a list of class \code{"summary.stabl_fit"} with
elements \code{selected_features}, \code{n_selected}, \code{stabl_threshold},
\code{min_fdr}, \code{lambda_grid_used}, \code{n_features_in} and
\code{n_bootstraps}.  \code{print()} returns the object invisibly.
}
\description{
Compact human-readable summary of a STABL fit: the selected feature set,
the stability threshold, the minimal FDP+ diagnostic, and the dimensions
of the fit.  The \code{min_fdr} field is an FDP+ \emph{diagnostic}, not an
FDR guarantee.
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
summary(fit)
}
\seealso{
\code{\link[=get_support]{get_support()}} for the binary selection mask,
\code{\link[=get_fdrs]{get_fdrs()}} for the full FDP+ diagnostic vector.
}
