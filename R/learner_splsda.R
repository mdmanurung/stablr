# sPLS-DA learner adapter (mixOmics) — classification-only base learner for
# the STABL bootstrap loop.  Follows the adapter contract of
# learner_adapters.R: the internal batch adapter maps a whole lambda-grid
# (here: a keepX grid) to a logical mask matrix in one call per bootstrap
# subsample; the exported factory returns the single-grid-row closure.

.SPLSDA_LOAD_TOL <- 1e-10

.stop_if_not_binomial_splsda <- function(family) {
  if (!identical(family, "binomial")) {
    stop("`base_learner = \"splsda\"` supports classification only: use ",
         "`family = \"binomial\"` (sPLS-DA has no gaussian, cox, or ",
         "ordinal variant).", call. = FALSE)
  }
}

.as_splsda_factor <- function(y) {
  y_f <- if (is.factor(y)) {
    y
  } else {
    values <- unique(as.character(y))
    if (!all(values %in% c("0", "1"))) {
      stop("`y` must be 0/1 numeric or a factor for ",
           "`base_learner = \"splsda\"`.", call. = FALSE)
    }
    factor(as.character(y), levels = c("0", "1"))
  }
  y_f <- droplevels(y_f)
  if (nlevels(y_f) < 2L) {
    stop("`y` must have at least two observed levels for ",
         "`base_learner = \"splsda\"`.", call. = FALSE)
  }
  y_f
}

# Internal: batch adapter over a keepX grid.  Signature matches
# .make_glmnet_batch_adapter: function(x, y, lambda_grid) -> logical matrix
# (ncol(x) x nrow(lambda_grid)).
.make_splsda_batch_adapter <- function(family, ncomp, bootstrap_threshold) {
  .stop_if_not_binomial_splsda(family)
  .require_pkg("mixOmics", "for base_learner = \"splsda\"")
  ncomp <- .validate_scalar_integer_like(ncomp, "ncomp", min = 1L)
  function(x, y, lambda_grid) {
    keepx <- if ("keepX" %in% names(lambda_grid)) {
      as.integer(lambda_grid[["keepX"]])
    } else {
      as.integer(lambda_grid[["lambda"]])  # `lambda` column doubles as keepX
    }
    if (anyNA(keepx) || any(keepx < 1L)) {
      stop("`lambda_grid` rows must carry keepX >= 1 (in the `keepX` column ",
           "or, falling back, in the `lambda` column).", call. = FALSE)
    }
    y_f <- .as_splsda_factor(y)
    # mixOmics requires unique sample ids and unique column identifiers; the
    # STABL augmented matrix deliberately reuses source feature names for
    # artificial columns (positions, not names, drive the accumulation).
    if (is.null(colnames(x))) {
      colnames(x) <- paste0("x.", seq_len(ncol(x)))
    } else {
      colnames(x) <- make.unique(colnames(x))
    }
    n_features <- ncol(x)
    n_lambdas  <- nrow(lambda_grid)
    result <- matrix(FALSE, nrow = n_features, ncol = n_lambdas)
    # A constant column is degenerate for PLS.  Drop such columns for this
    # bootstrap's fit and map the mask back so only the dropped columns are
    # all-FALSE — previously ONE constant column discarded the mask for ALL
    # features at EVERY grid row (the check was inside the loop but
    # lambda-invariant), biasing stability frequencies downward when rare
    # binary features are present (review E7).
    constant <- apply(x, 2L, function(col) diff(range(col)) == 0)
    x_fit <- x[, !constant, drop = FALSE]
    for (l in seq_len(n_lambdas)) {
      # An all-constant subsample is degenerate: contribute an all-FALSE
      # mask instead of failing the run.
      if (ncol(x_fit) == 0L) next
      fit <- mixOmics::splsda(x_fit, y_f,
                              keepX = rep(min(keepx[l], ncol(x_fit)), ncomp),
                              ncomp = ncomp)
      load <- fit$loadings$X
      if (is.null(dim(load))) next
      result[!constant, l] <- matrixStats::rowMaxs(abs(as.matrix(load))) >
        bootstrap_threshold
    }
    result
  }
}

#' Build a sPLS-DA Learner Adapter
#'
#' Returns a closure that fits a `mixOmics::splsda` model on a single
#' bootstrap subsample and returns a logical selection mask over the features.
#' The closure mirrors the adapter contract of [make_glmnet_adapter()]: the
#' STABL bootstrap loop supplies the subsample and one lambda-grid row, and
#' the adapter answers with the selected-feature indicator.
#'
#' The regularisation parameter swept by the STABL lambda grid is sPLS-DA's
#' `keepX` (the number of features kept per component).  Grid rows are read
#' from a `keepX` column when present and fall back to the `lambda` column
#' otherwise, so [auto_splsda_keepx_grid()] output and hand-built
#' `data.frame(lambda = ...)` grids both work.
#'
#' A feature is counted as selected in a bootstrap when the magnitude of its
#' loading on any retained component exceeds `bootstrap_threshold`.  Sparse
#' PLS loadings are exactly zero for non-retained features, so this is an
#' exact selection indicator, not a soft threshold.
#'
#' @param family Character; must be `"binomial"`.  sPLS-DA is a
#'   classification learner; there is no gaussian, cox, or ordinal variant.
#' @param ncomp Positive integer; number of latent components to retain in
#'   every fit.  Default `2L`.
#' @param bootstrap_threshold Positive numeric; loading-magnitude cutoff
#'   above which a feature is counted as selected.  Default `1e-10` (any
#'   non-zero sparse loading counts).
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`
#'   for use inside [stabl_fit()].
#'
#' @seealso [auto_splsda_keepx_grid()], [make_glmnet_adapter()],
#'   [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("mixOmics", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(120L), 30L, 20L)
#'   y <- factor(rep(c("0", "1"), each = 15L))
#'   adapter <- stablr:::make_splsda_adapter()
#'   mask <- adapter(x, y, data.frame(lambda = 5))
#'   cat("selected:", sum(mask), "of", ncol(x), "features\n")
#' }
#' @keywords internal
make_splsda_adapter <- function(family              = "binomial",
                                ncomp               = 2L,
                                bootstrap_threshold = .SPLSDA_LOAD_TOL) {
  .stop_if_not_binomial_splsda(family)
  .require_pkg("mixOmics", "for base_learner = \"splsda\"")
  batch_adapter <- .make_splsda_batch_adapter(
    family = family,
    ncomp = ncomp,
    bootstrap_threshold = bootstrap_threshold
  )

  function(x, y, lambda_val) {
    batch_adapter(x, y, lambda_val)[, 1L]
  }
}

#' Build an Automatic keepX Grid for the sPLS-DA Learner
#'
#' Returns a pre-expanded lambda grid for `base_learner = "splsda"` whose rows
#' sweep sPLS-DA's per-component `keepX` (the number of features retained per
#' latent component).  The grid is geometric in `p = ncol(x)`, spanning from
#' roughly `p / 50` up to `p / 3`, which keeps the sparse end informative for
#' FDP+ thresholding while still probing supports large enough to matter for
#' prediction.
#'
#' The returned `data.frame` carries the keepX values in both the `lambda`
#' column (the column [stabl_fit()] requires) and a dedicated `keepX` column
#' (which [make_splsda_adapter()] reads preferentially).
#'
#' @param x Numeric matrix; the real feature matrix (only `ncol(x)` is used).
#' @param n_lambda Positive integer; number of grid rows.  Default `10L`.
#'   Duplicated rounded values are dropped, so the returned grid may be
#'   shorter than `n_lambda`.
#'
#' @return A `data.frame` with columns `lambda` and `keepX` (identical
#'   positive integer values, sorted increasingly).
#'
#' @seealso [make_splsda_adapter()], [auto_lambda_grid()]
#'
#' @examples
#' grid <- auto_splsda_keepx_grid(matrix(0, 10, 150), n_lambda = 8)
#' grid
#' @export
auto_splsda_keepx_grid <- function(x, n_lambda = 10L) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`x` must be a matrix or data.frame.", call. = FALSE)
  }
  n_lambda <- .validate_scalar_integer_like(n_lambda, "n_lambda",
                                            min = 1L, max = 1000L)
  p <- ncol(x)
  lo <- max(1L, as.integer(floor(p / 50)))
  hi <- max(lo + 1L, as.integer(floor(p / 3)))
  keepx <- unique(as.integer(round(exp(seq(log(lo), log(hi),
                                          length.out = n_lambda)))))
  keepx <- sort(keepx)
  data.frame(lambda = keepx, keepX = keepx)
}
