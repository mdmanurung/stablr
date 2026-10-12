% Generated manually for the stablr improvement plan (item D3)
\name{stabl_environment_stability}
\alias{stabl_environment_stability}
\title{Environment-Stability Workflow: Reliable and Transportable Feature Sets}
\usage{
stabl_environment_stability(x, y, environments, lambda_grid = "auto",
                            random_state = NULL, ...)
}
\arguments{
\item{x}{A \code{matrix} or \code{data.frame} of features (samples in rows).}
\item{y}{Named outcome vector; sample IDs must match \code{rownames(x)}.}
\item{environments}{Factor or character vector with one label per sample
(batch, cohort, site, ...).  Must contain at least two non-empty levels
and no \code{NA}.}
\item{lambda_grid}{Passed to \code{\link[=stabl_fit]{stabl_fit()}}; either
\code{"auto"} (per-environment grids) or an explicit grid.  Note that
\code{"auto"} grids reach the saturated small-lambda regime, which makes
max-over-lambda scores lenient; explicit top-heavy grids discriminate
better.}
\item{random_state}{Optional integer seed; per-environment seeds are
derived deterministically as \code{random_state + level_index}.}
\item{...}{Additional arguments forwarded to \code{\link[=stabl_fit]{stabl_fit()}}
in every environment.}
}
\value{
A list of class \code{"stabl_env_stability"} with elements
\code{env_fits} (named list of per-environment \code{"stabl_fit"} objects),
\code{env_scores} (features x environments matrix of max-over-lambda
selection frequencies), \code{min_scores} (minimum across environments),
\code{mean_scores}, \code{environments} (levels used),
\code{pairwise_adjusted_similarity} (environment x environment similarity
of selected sets), and \code{call}.
}
\description{
Fits STABL separately in each environment and aggregates per-feature
selection frequencies by their \emph{minimum} across environments: a
feature is only as reliable as its weakest environment.  This delivers
ICP-inspired stability-across-environments (Peters et al.) as a
\emph{predictive} selection heuristic — without causal assumptions, without
an environment-variable requirement beyond labels, and without new
dependencies.
}
\details{
Environment labels may be batch, cohort, or site identifiers.  Group
integrity within each per-environment fit is preserved by passing
environment labels through the usual grouped-bootstrap machinery.  The
minimum-aggregation rule is deliberately conservative: features whose
signal is environment-specific receive a low minimum score.  This is a
stability heuristic, \bold{not} a causal-effect estimate and not an FDR
guarantee.
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
res$min_scores                      # transportability summary
get_min_scores(res)                 # same, via accessor
}
\seealso{
\code{\link[=stabl_fit]{stabl_fit()}},
\code{\link[=get_importances]{get_importances()}},
\code{\link[=adjusted_similarity]{adjusted_similarity()}}.
}
