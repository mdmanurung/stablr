# Ordinal-regression learner adapter (ordinalNet) — experimental extension
# (workstream WS-9 of the 0.3.0 feature-maturity roadmap).
#
# Statistical design:
# * The base learner is the elastic-net cumulative-logit (proportional-odds)
#   model of ordinalNet with parallelTerms = TRUE: one slope per feature,
#   shared across the outcome bands.  Selection follows the same
#   absolute-coefficient rule as the glmnet adapters.
# * The outcome is a factor whose LEVEL ORDER is the band order.  An ordered
#   factor is recommended; a plain factor is accepted with its level order.
# * The regularisation path is a lambda path (ordinalNet's own decreasing
#   sequence, or a user-supplied grid), so max-over-grid stability
#   frequencies and the FDP+ machinery work unchanged.
# * Decoys and FDP+ are outcome-agnostic: artificial features are generated
#   from `x` alone, so knockoff and permutation decoys apply unchanged.
# * No Python-parity counterpart exists (Python STABL has no ordinal
#   learner); this adapter is R-only and ships EXPERIMENTAL until the WS-9
#   benchmark (proportional-odds simulations, both decoy types, graduation
#   criterion mean FDP <= 0.15 AND TPR >= 0.8) passes.

#' Build an Ordinal-Regression Learner Adapter (Experimental)
#'
#' Returns a closure that fits ONE elastic-net cumulative-logit model
#' (`ordinalNet`) on a single bootstrap subsample and returns a logical
#' selection mask over the features, mirroring the role of
#' [make_glmnet_adapter()] for ordered outcomes.
#'
#' @section Selection rule:
#' With `parallelTerms = TRUE` the proportional-odds model gives one slope per
#' feature.  A feature is selected at a grid row when the absolute value of
#' its slope exceeds the same small coefficient cutoff the glmnet adapters
#' use.  The grid is a lambda path (see [auto_ordinal_lambda_grid()]); lower
#' penalties give supersets of the selections at higher penalties, so the
#' stability-accumulation and FDP+ logic in [stabl_fit()] work unchanged.
#'
#' @section Statistical guardrails (read before using):
#' * **Outcome type.**  `y` must be a factor; its level order is the band
#'   order.  An ordered factor is recommended.  Numeric outcomes must be
#'   banded by the caller (for example with `cut()`).
#' * **Experimental.**  This learner is an extension beyond Python STABL and
#'   has no Python-parity counterpart.  FDP+ calibration for the ordinal
#'   learner is validated by the WS-9 simulation benchmark; consult the
#'   benchmark summary before treating FDP+ as calibrated.
#' * **Decoys.**  Artificial features are outcome-agnostic; both
#'   `random_permutation` and model-X knockoffs apply unchanged.
#'
#' @param alpha Numeric in `(0, 1]`; the elastic-net mixing weight inside
#'   `ordinalNet`.  Default `1` (pure lasso).
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`
#'   for use inside [stabl_fit()], where `lambda_val` is a 1-row `data.frame`
#'   from the lambda grid (see [auto_ordinal_lambda_grid()]).
#'
#' @seealso [auto_ordinal_lambda_grid()], [make_glmnet_adapter()],
#'   [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("ordinalNet", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(400L), 40L, 10L,
#'               dimnames = list(NULL, paste0("f", 1:10)))
#'   eta <- x[, 1] - 0.5 * x[, 2]
#'   y   <- cut(eta + rnorm(40L, sd = 0.5), breaks = 4,
#'              labels = paste0("L", 1:4), ordered_result = TRUE)
#'   grid   <- auto_ordinal_lambda_grid(x, y, n_lambda = 5L)
#'   adapter <- make_ordinal_adapter()
#'   mask <- adapter(x, y, grid[1L, , drop = FALSE])
#'   cat("selected at sparsest lambda:", sum(mask), "of", ncol(x), "\n")
#' }
#' @section Calibration status (experimental): Simulated calibration evidence
#' (10 seeds x 2 proportional-odds designs) shows knockoff_equi mean FDP+ of
#' 0.00 (nominal 0.1) but mean TPR of only 0.37 on both designs -- the learner
#' is strongly conservative and did not meet the predeclared graduation bar.
#' It stays experimental; see benchmarks/ws9_summary.md.
#'
#' @export
make_ordinal_adapter <- function(alpha = 1) {
  .require_pkg("ordinalNet", "for base_learner = \"ordinal\"")
  alpha <- .validate_scalar_numeric(alpha, "alpha", min = 0, min_open = TRUE,
                                    max = 1)
  batch_adapter <- .make_ordinal_batch_adapter(alpha = alpha)

  function(x, y, lambda_val) {
    batch_adapter(x, y, lambda_val)[, 1L]
  }
}

#' Build a Data-Driven Lambda Grid for the Ordinal Learner (Experimental)
#'
#' Fits one `ordinalNet` cumulative-logit model on the full training data and
#' returns its decreasing lambda path as a pre-expanded `data.frame` suitable
#' for the `lambda_grid` argument of [stabl_fit()] with
#' `base_learner = "ordinal"`.  This is the ordinal analog of
#' [auto_lambda_grid()].
#'
#' The grid is computed once on the full data (exactly like
#' [auto_lambda_grid()]); each bootstrap subsample is subsequently fit at
#' these fixed penalties with its own model.
#'
#' @param x Numeric matrix of predictors (samples by features).  Column names
#'   are recommended; unnamed matrices are named internally by `ordinalNet`.
#' @param y Factor outcome; the level order is the band order.  An ordered
#'   factor is recommended.
#' @param n_lambda Positive integer; number of lambda values.  Default `30L`
#'   (matches [auto_lambda_grid()]'s default).
#' @param alpha Numeric in `[0, 1]`; the elastic-net mixing weight.  Default
#'   `1` (pure lasso).
#'
#' @return A `data.frame` with a single `lambda` column, sorted in decreasing
#'   order so the path runs sparse-to-dense, matching the convention used by
#'   [stabl_fit()].
#'
#' @seealso [make_ordinal_adapter()], [auto_lambda_grid()], [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("ordinalNet", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(400L), 40L, 10L,
#'               dimnames = list(NULL, paste0("f", 1:10)))
#'   eta <- x[, 1] - 0.5 * x[, 2]
#'   y   <- cut(eta + rnorm(40L, sd = 0.5), breaks = 4,
#'              labels = paste0("L", 1:4), ordered_result = TRUE)
#'   grid <- auto_ordinal_lambda_grid(x, y, n_lambda = 5L)
#'   head(grid)
#' }
#' @export
auto_ordinal_lambda_grid <- function(x, y, n_lambda = 30L, alpha = 1) {
  .require_pkg("ordinalNet", "for base_learner = \"ordinal\"")
  n_lambda <- .validate_scalar_integer_like(n_lambda, "n_lambda", min = 1L)
  alpha <- .validate_scalar_numeric(alpha, "alpha", min = 0, min_open = TRUE,
                                    max = 1)
  y <- .validate_ordinal_outcome(y)
  x <- as.matrix(x)

  fit <- .with_glmnet_numerical_infeasibility(
    ordinalNet::ordinalNet(
      x                 = x,
      y                 = y,
      alpha             = alpha,
      family            = "cumulative",
      link              = "logit",
      parallelTerms     = TRUE,
      nonparallelTerms  = FALSE,
      nLambda           = n_lambda
    ),
    source = "ordinalNet automatic lambda path"
  )
  lam_seq <- fit$lambdaVals
  if (is.null(lam_seq) || length(lam_seq) == 0L || anyNA(lam_seq)) {
    .abort_numerical_infeasibility(
      "stablr_learner_numerical_infeasibility",
      "ordinalNet automatic lambda path returned no usable lambda values.",
      source = "ordinalNet automatic lambda path"
    )
  }
  data.frame(lambda = sort(lam_seq, decreasing = TRUE))
}

# ---- Internal helpers ---------------------------------------------------------

# Validate the ordinal outcome: a factor with at least two levels and no NA.
# An unordered factor is accepted with a warning (the level order is the band
# order).  Returns the factor unchanged.
.validate_ordinal_outcome <- function(y) {
  if (!is.factor(y)) {
    stop(
      "`base_learner = \"ordinal\"` requires `y` to be a factor whose ",
      "level order is the band order. Band a numeric outcome with cut() ",
      "first.",
      call. = FALSE
    )
  }
  if (anyNA(y)) {
    stop("`y` must not contain NA values for the ordinal learner.",
         call. = FALSE)
  }
  if (length(levels(y)) < 2L) {
    stop("`y` must have at least two levels for the ordinal learner.",
         call. = FALSE)
  }
  if (!is.ordered(y)) {
    warning(
      "`y` is a factor but not an ordered factor; the level order ",
      paste0("(", paste(levels(y), collapse = " < "), ") "),
      "is used as the band order.",
      call. = FALSE
    )
  }
  y
}

# Batch adapter: ONE ordinalNet fit per bootstrap subsample; masks for the
# whole grid are read from that single coefficient path.  Signature matches
# the other batch adapters:
#   function(x, y, lambda_grid) -> logical matrix (n_features x n_lambdas)
.make_ordinal_batch_adapter <- function(alpha = 1) {
  force(alpha)
  function(x, y, lambda_grid) {
    y <- .validate_ordinal_outcome(y)
    n_lambdas  <- nrow(lambda_grid)
    n_features <- ncol(x)

    lam_seq <- sort(as.numeric(lambda_grid[["lambda"]]), decreasing = TRUE)
    if (length(lam_seq) != n_lambdas || anyNA(lam_seq) ||
        any(!is.finite(lam_seq))) {
      stop(
        "The ordinal lambda grid must contain finite numeric `lambda` values.",
        call. = FALSE
      )
    }

    fit <- .with_glmnet_numerical_infeasibility(
      ordinalNet::ordinalNet(
        x                 = x,
        y                 = y,
        alpha             = alpha,
        family            = "cumulative",
        link              = "logit",
        parallelTerms     = TRUE,
        nonparallelTerms  = FALSE,
        lambdaVals        = lam_seq
      ),
      source = "ordinalNet bootstrap path"
    )

    # Coefficient layout (parallelTerms = TRUE): fit$coefs is
    # (n_lambda x n_coef); the first (K - 1) columns are the intercepts, then
    # one slope column per feature, in the column order of `x`.
    coefs <- fit$coefs
    if (is.null(coefs) || ncol(coefs) < n_features + 1L) {
      stop(
        "ordinalNet did not return usable coefficients for this bootstrap ",
        "subsample.",
        call. = FALSE
      )
    }
    slopes <- coefs[, seq.int(ncol(coefs) - n_features + 1L, ncol(coefs)),
                    drop = FALSE]
    if (nrow(slopes) != n_lambdas) {
      # ordinalNet may trim the path; fall back to per-lambda extraction so
      # the returned matrix always has one column per grid row.
      slopes <- t(vapply(
        lam_seq,
        function(s) {
          idx <- which.min(abs(fit$lambdaVals - s))
          coefs[idx, seq.int(ncol(coefs) - n_features + 1L, ncol(coefs))]
        },
        numeric(n_features)
      ))
    }

    masks <- t(abs(slopes) > .BOOTSTRAP_COEF_THRESHOLD)
    dimnames(masks) <- NULL
    masks
  }
}
