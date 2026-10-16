\name{make_ordinal_adapter}
\alias{make_ordinal_adapter}
\title{Build an Ordinal-Regression Learner Adapter (Experimental)}
\usage{
make_ordinal_adapter(alpha = 1)
}
\arguments{
\item{alpha}{Numeric in \code{(0, 1]}; the elastic-net mixing weight inside
\code{ordinalNet}.  Default \code{1} (pure lasso).}
}
\value{
A function with signature
\code{function(x, y, lambda_val) -> logical vector of length ncol(x)}
for use inside \code{\link[=stabl_fit]{stabl_fit()}}, where \code{lambda_val}
is a 1-row \code{data.frame} from the lambda grid (see
\code{\link[=auto_ordinal_lambda_grid]{auto_ordinal_lambda_grid()}}).
}
\description{
Simulated calibration evidence (10 seeds x 2 designs, proportional-odds
ground truth): knockoff_equi mean FDP+ was 0.00 (nominal 0.1) but mean TPR
was only 0.37 on both designs -- the learner is strongly conservative and
did NOT meet the predeclared graduation bar (FDP+ <= 0.15 AND TPR >= 0.8).
It stays experimental; see benchmarks/ws9_summary.md.

Returns a closure that fits ONE elastic-net cumulative-logit model
(\code{ordinalNet}) on a single bootstrap subsample and returns a logical
selection mask over the features, mirroring the role of
\code{\link[=make_glmnet_adapter]{make_glmnet_adapter()}} for ordered
outcomes.
}
\details{
With \code{parallelTerms = TRUE} the proportional-odds model gives one slope
per feature.  A feature is selected at a grid row when the absolute value of
its slope exceeds the same small coefficient cutoff the glmnet adapters use.
The grid is a lambda path (see
\code{\link[=auto_ordinal_lambda_grid]{auto_ordinal_lambda_grid()}}); lower
penalties give supersets of the selections at higher penalties, so the
stability-accumulation and FDP+ logic in
\code{\link[=stabl_fit]{stabl_fit()}} work unchanged.

Guardrails: \code{y} must be a factor whose level order is the band order
(an ordered factor is recommended).  This learner is an extension beyond
Python STABL and has no Python-parity counterpart; it ships EXPERIMENTAL
until the WS-9 benchmark passes.  Artificial features are outcome-agnostic,
so both \code{random_permutation} and model-X knockoff decoys apply
unchanged.
}
\examples{
if (requireNamespace("ordinalNet", quietly = TRUE)) {
  set.seed(1L)
  x <- matrix(rnorm(400L), 40L, 10L,
              dimnames = list(NULL, paste0("f", 1:10)))
  eta <- x[, 1] - 0.5 * x[, 2]
  y   <- cut(eta + rnorm(40L, sd = 0.5), breaks = 4,
             labels = paste0("L", 1:4), ordered_result = TRUE)
  grid   <- auto_ordinal_lambda_grid(x, y, n_lambda = 5L)
  adapter <- make_ordinal_adapter()
  mask <- adapter(x, y, grid[1L, , drop = FALSE])
  cat("selected at sparsest lambda:", sum(mask), "of", ncol(x), "\\n")
}
}
\seealso{
\code{\link[=auto_ordinal_lambda_grid]{auto_ordinal_lambda_grid()}},
\code{\link[=make_glmnet_adapter]{make_glmnet_adapter()}},
\code{\link[=stabl_fit]{stabl_fit()}}
}
