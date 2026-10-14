\name{auto_ordinal_lambda_grid}
\alias{auto_ordinal_lambda_grid}
\title{Build a Data-Driven Lambda Grid for the Ordinal Learner (Experimental)}
\usage{
auto_ordinal_lambda_grid(x, y, n_lambda = 30L, alpha = 1)
}
\arguments{
\item{x}{A \code{matrix} or \code{data.frame} of features (samples in rows).}
\item{y}{Factor outcome; the level order is the band order.  An ordered
factor is recommended.}
\item{n_lambda}{Positive integer; number of lambda values.  Default \code{30L}.}
\item{alpha}{Numeric in \code{[0, 1]}; the elastic-net mixing weight.
Default \code{1} (pure lasso).}
}
\value{
A \code{data.frame} with a single \code{lambda} column, sorted in decreasing
order so the path runs sparse-to-dense, matching the convention used by
\code{\link[=stabl_fit]{stabl_fit()}}.
}
\description{
Builds the ordinal analog of
\code{\link[=auto_lambda_grid]{auto_lambda_grid()}}: one
\code{ordinalNet} cumulative-logit model with parallel terms is fit on the
full data, and the grid is taken from its decreasing lambda path.  The grid
is computed once on the full data; each bootstrap subsample is subsequently
fit at these fixed penalties with its own model.
}
\examples{
if (requireNamespace("ordinalNet", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(400L), 40L, 10L,
              dimnames = list(NULL, paste0("f", 1:10)))
  eta <- x[, 1] - 0.5 * x[, 2]
  y   <- cut(eta + rnorm(40L, sd = 0.5), breaks = 4,
             labels = paste0("L", 1:4), ordered_result = TRUE)
  grid <- auto_ordinal_lambda_grid(x, y, n_lambda = 5L)
  head(grid)
}
}
\seealso{
\code{\link[=make_ordinal_adapter]{make_ordinal_adapter()}},
\code{\link[=auto_lambda_grid]{auto_lambda_grid()}},
\code{\link[=stabl_fit]{stabl_fit()}}
}
