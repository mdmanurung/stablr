% Generated manually for the stablr improvement plan (item D3)
\name{print.stabl_env_stability}
\alias{print.stabl_env_stability}
\title{Print an Environment-Stability Analysis}
\usage{
\method{print}{stabl_env_stability}(x, ...)
}
\arguments{
\item{x}{A \code{"stabl_env_stability"} object returned by
\code{\link[=stabl_environment_stability]{stabl_environment_stability()}}.}
\item{...}{Unused; included for S3 consistency.}
}
\value{
The input object, invisibly.
}
\description{
Compact summary of an environment-stability analysis: the environments
analyzed, the number of features, and the top features by
minimum-across-environment stability score.
}
\examples{
set.seed(1L)
n <- 90L
env <- factor(rep(c("batchA", "batchB", "batchC"), each = n / 3))
x <- matrix(rnorm(n * 6), n, 6,
            dimnames = list(paste0("s", seq_len(n)), paste0("f", 1:6)))
y <- setNames(1.2 * x[, "f1"] + 0.8 * x[, "f2"] + rnorm(n, sd = 0.5),
              rownames(x))
res <- stabl_environment_stability(
  x, y, env,
  lambda_grid   = data.frame(lambda = c(0.6, 0.4, 0.25)),
  n_bootstraps  = 10L,
  hard_threshold = 0.5,
  random_state  = 1L
)
print(res)
}
\seealso{
\code{\link[=stabl_environment_stability]{stabl_environment_stability()}}.
}
