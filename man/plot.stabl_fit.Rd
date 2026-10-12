% Generated manually for the stablr improvement plan (item C4)
\name{plot.stabl_fit}
\alias{plot.stabl_fit}
\title{Plot a Fitted stabl_fit Object}
\usage{
\method{plot}{stabl_fit}(x, ...)
}
\arguments{
\item{x}{A fitted \code{"stabl_fit"} object returned by \code{\link[=stabl_fit]{stabl_fit()}}.}
\item{...}{
Arguments forwarded to \code{\link[=plot_stabl_path]{plot_stabl_path()}}.
}
}
\value{
A \code{ggplot} object (invisibly), as returned by
\code{\link[=plot_stabl_path]{plot_stabl_path()}}.
}
\description{
Convenience S3 method so a fitted STABL model can be visualised with
\code{plot(fit)}; delegates to \code{\link[=plot_stabl_path]{plot_stabl_path()}}.
}
\examples{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(30 * 5), 30, 5,
              dimnames = list(paste0("s", 1:30), paste0("f", 1:5)))
  y <- setNames(rnorm(30), rownames(x))
  fit <- stabl_fit(x, y,
                   lambda_grid    = data.frame(lambda = c(0.2, 0.1, 0.05)),
                   n_bootstraps   = 4L,
                   hard_threshold = 0.3,
                   random_state   = 1L)
  plot(fit)
}
}
\seealso{
\code{\link[=plot_stabl_path]{plot_stabl_path()}} for the full signature and styling options.
}
