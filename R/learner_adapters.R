#' Build a glmnet Learner Adapter
#'
#' Returns a closure that fits a `glmnet` model on a single bootstrap
#' subsample and returns a logical selection mask over the features.  The
#' closure mirrors the role of `fit_bootstrapped_sample()` in the Python STABL
#' library.
#'
#' Learner adapters decouple the modelling back-end from the STABL bootstrap
#' loop, making it easy to substitute different penalised-regression solvers
#' without changing the stability-accumulation logic.  This factory produces
#' the standard lasso / elastic-net adapter that is the default back-end for
#' [stabl_fit()].
#'
#' The returned function accepts a pre-expanded lambda-grid row (a 1-row
#' `data.frame`) and applies `glmnet` at that exact penalty value.  The
#' `alpha` elastic-net mixing parameter is resolved in priority order:
#' `alpha_fixed` argument > `alpha` column in `lambda_val` > default of 1
#' (pure lasso).
#'
#' @details
#' This adapter factory is primarily an internal backend used by [stabl_fit()].
#' Most users should configure `base_learner`, `family`, and `lambda_grid`
#' directly in [stabl_fit()] rather than calling this factory manually.
#'
#' @param family Character; the `glmnet` response family, for example
#'   `"gaussian"`, `"binomial"`, `"multinomial"`, or `"cox"`.
#' @param alpha_fixed Numeric scalar or `NULL`.  When not `NULL`, this value
#'   overrides any `alpha` column in `lambda_val`.
#' @param bootstrap_threshold Positive numeric; absolute-coefficient cutoff
#'   used to decide that a feature is selected in a given bootstrap.  Features
#'   with `|coef| > bootstrap_threshold` are counted as selected.
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`
#'   for use inside [stabl_fit()].
#'
#' @seealso [make_adaptive_lasso_adapter()], [make_sgl_adapter()],
#'   [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1L)
#'   x <- matrix(rnorm(200L), 20L, 10L)
#'   y <- rnorm(20L)
#'   adapter <- make_glmnet_adapter(family = "gaussian", bootstrap_threshold = 1e-5)
#'   mask <- adapter(x, y, data.frame(lambda = 0.1))
#'   cat("selected:", sum(mask), "of", ncol(x), "features\n")
#' }
#' @export
make_glmnet_adapter <- function(
    family              = "gaussian",
    alpha_fixed         = NULL,
    bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
) {
  .require_pkg("glmnet")
  batch_adapter <- .make_glmnet_batch_adapter(
    family = family,
    alpha_fixed = alpha_fixed,
    bootstrap_threshold = bootstrap_threshold
  )

  function(x, y, lambda_val) {
    batch_adapter(x, y, lambda_val)[, 1L]
  }
}

#' Build an Adaptive Lasso Learner Adapter
#'
#' Returns a closure that computes feature-specific penalty weights from a
#' ridge initialisation on each bootstrap subsample, then fits a lasso model
#' with those `penalty.factor` values to obtain the selected-feature mask.
#'
#' Adaptive lasso improves on standard lasso by assigning stronger penalties
#' to features with small initial coefficients (likely noise) and weaker
#' penalties to features with large initial coefficients (likely signal).
#' This asymmetric penalisation achieves the oracle property under regularity
#' conditions, selecting the true support more reliably than plain lasso when
#' signal features have moderate to large effect sizes.
#'
#' Weights are defined as:
#' \deqn{w_j = 1 / (|\hat\beta_j^{\mathrm{init}}| + \epsilon)^\gamma}
#' where \eqn{\hat\beta_j^{\mathrm{init}}} comes from a ridge regression on
#' the same bootstrap subsample.  The `epsilon` floor avoids division by zero
#' for features with near-zero ridge coefficients.
#'
#' @details
#' This adapter factory is primarily an internal backend used by [stabl_fit()].
#' For end-to-end feature selection workflows, prefer
#' `base_learner = "adaptive_lasso"` in [stabl_fit()].
#'
#' @param family Character; the `glmnet` response family, for example
#'   `"gaussian"`, `"binomial"`, `"multinomial"`, or `"cox"`.
#' @param gamma Positive numeric scalar; controls how sharply the weights
#'   down-penalise features with large ridge coefficients.  Larger values
#'   make the penalty more selective.  Default `1.0` matches the Python
#'   reference implementation.
#' @param epsilon Positive numeric scalar; added to the denominator of the
#'   weight to avoid division by zero.  Default `1e-6`.
#' @param bootstrap_threshold Positive numeric; absolute-coefficient cutoff
#'   used to decide that a feature is selected in a given bootstrap.  Features
#'   with `|coef| > bootstrap_threshold` are counted as selected.
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`.
#'
#' @seealso [make_glmnet_adapter()], [make_sgl_adapter()], [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(2L)
#'   x <- matrix(rnorm(200L), 20L, 10L)
#'   y <- rnorm(20L)
#'   adapter <- make_adaptive_lasso_adapter(family = "gaussian", gamma = 1, epsilon = 1e-6)
#'   mask <- adapter(x, y, data.frame(lambda = 0.1))
#'   cat("selected:", sum(mask), "of", ncol(x), "features\n")
#' }
#' @export
make_adaptive_lasso_adapter <- function(
    family              = "gaussian",
    gamma               = 1.0,
    epsilon             = 1e-6,
    bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
) {
  .require_pkg("glmnet")
  if (!is.numeric(gamma) || length(gamma) != 1L || gamma <= 0) {
    stop("`gamma` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive numeric scalar.", call. = FALSE)
  }
  batch_adapter <- .make_adaptive_lasso_batch_adapter(
    family = family,
    gamma = gamma,
    epsilon = epsilon,
    bootstrap_threshold = bootstrap_threshold
  )

  function(x, y, lambda_val) {
    batch_adapter(x, y, lambda_val)[, 1L]
  }
}

.stop_if_cox_sgl <- function(family) {
  if (identical(family, "cox")) {
    stop(
      "Cox family is not supported by sparse_group_lasso. ",
      "Use base_learner = \"lasso\" or \"adaptive_lasso\" with family = \"cox\".",
      call. = FALSE
    )
  }
}

#' Build a Sparse Group Lasso Learner Adapter
#'
#' Returns a closure that fits `sparsegl::sparsegl()` on a single bootstrap
#' subsample and returns a logical selection mask over the features.
#'
#' Sparse-group lasso is useful when features have a known (or inferred) block
#' structure, for example gene pathways, omic layers, or correlated feature
#' clusters.  It imposes simultaneous sparsity within and between groups: the
#' group-level penalty encourages whole groups to be zeroed out, while the
#' within-group lasso penalty allows groups to have only a sparse subset of
#' active features.  This can substantially improve stability in structured
#' high-dimensional settings.
#'
#' The `alpha` value (`asparse` in `sparsegl`) controls the balance between
#' the within-group lasso penalty and the group-level penalty: 0 = pure group
#' lasso; 1 = pure lasso (no group penalty).  The default of `0.05` used when
#' no `alpha` column is present matches the Python STABL reference.
#'
#' @details
#' This adapter factory is primarily an internal backend used by [stabl_fit()].
#' For end-to-end feature selection workflows, prefer
#' `base_learner = "sparse_group_lasso"` in [stabl_fit()] and provide
#' `feature_groups` (or `corr_group_threshold`) there.
#'
#' @param family Character; response family (`"gaussian"`, `"binomial"`, or
#'   `"multinomial"`).  Cox regression is not supported by `sparsegl`.
#' @param feature_groups Integer or factor vector of length `p` (number of
#'   features) assigning each feature to a group.  Values are coerced via
#'   `as.integer(as.factor(...))`, so any type that has a natural ordering is
#'   accepted.
#' @param alpha_fixed Numeric scalar in `[0, 1]` or `NULL`.  When not `NULL`,
#'   overrides any `alpha` column in `lambda_val`.  Controls the lasso /
#'   group-lasso mixing weight (`asparse` in `sparsegl`).
#' @param bootstrap_threshold Positive numeric; absolute-coefficient cutoff
#'   used to decide that a feature is selected in a given bootstrap.
#'
#' @return A function with signature
#'   `function(x, y, lambda_val) -> logical vector of length ncol(x)`.
#'
#' @seealso [make_glmnet_adapter()], [make_adaptive_lasso_adapter()],
#'   [stabl_fit()]
#'
#' @examples
#' if (requireNamespace("sparsegl", quietly = TRUE)) {
#'   set.seed(3L)
#'   x <- matrix(rnorm(200L), 20L, 10L)
#'   y <- rnorm(20L)
#'   groups  <- rep(1:5, each = 2L)
#'   adapter <- make_sgl_adapter(family = "gaussian", feature_groups = groups)
#'   mask <- adapter(x, y, data.frame(lambda = 0.1))
#'   cat("selected:", sum(mask), "of", ncol(x), "features\n")
#' }
#' @export
make_sgl_adapter <- function(
    family              = "gaussian",
    feature_groups,
    alpha_fixed         = NULL,
    bootstrap_threshold = .BOOTSTRAP_COEF_THRESHOLD
) {
  .require_pkg("sparsegl", "for base_learner = \"sparse_group_lasso\"")
  .stop_if_cox_sgl(family)
  if (missing(feature_groups) || is.null(feature_groups)) {
    stop("`feature_groups` must be provided for sparse group lasso.", call. = FALSE)
  }

  groups <- as.integer(as.factor(feature_groups))
  if (length(groups) == 0L) {
    stop("`feature_groups` cannot be empty.", call. = FALSE)
  }
  if (anyNA(groups)) {
    stop("`feature_groups` cannot contain NA values.", call. = FALSE)
  }
  if (!is.null(alpha_fixed) &&
      (!is.numeric(alpha_fixed) || length(alpha_fixed) != 1L ||
       alpha_fixed < 0 || alpha_fixed > 1)) {
    stop("`alpha_fixed` must be NULL or a numeric scalar in [0, 1].", call. = FALSE)
  }

  .resolve_asparse <- function(lambda_val) {
    if (!is.null(alpha_fixed)) {
      return(as.numeric(alpha_fixed))
    }
    if ("alpha" %in% names(lambda_val)) {
      alpha_val <- as.numeric(lambda_val[["alpha"]])
      if (length(alpha_val) != 1L || is.na(alpha_val) ||
          alpha_val < 0 || alpha_val > 1) {
        stop("`alpha` in `lambda_val` must be in [0, 1] for sparse_group_lasso.", call. = FALSE)
      }
      return(alpha_val)
    }
    .SGL_DEFAULT_ALPHA
  }

  function(x, y, lambda_val) {
    if (length(groups) != ncol(x)) {
      stop("`feature_groups` length must match ncol(x) for sparse_group_lasso.", call. = FALSE)
    }

    lambda_use  <- as.numeric(lambda_val[["lambda"]])
    asparse_use <- .resolve_asparse(lambda_val)

    # sparsegl requires features sorted by group (non-decreasing group index).
    sort_ord <- order(groups)
    inv_ord  <- order(sort_ord)
    x_s      <- x[, sort_ord, drop = FALSE]
    grp_s    <- groups[sort_ord]

    if (family == "multinomial") {
      if (!is.factor(y)) {
        y <- factor(y)
      }
      levs <- levels(y)
      if (length(levs) < 2L) {
        stop("Multinomial sparse-group fitting requires at least 2 classes.", call. = FALSE)
      }

      coef_mat <- matrix(0, nrow = ncol(x), ncol = length(levs))
      for (k in seq_along(levs)) {
        yk <- as.integer(y == levs[[k]])
        fit_k <- sparsegl::sparsegl(
          x = x_s,
          y = yk,
          group = grp_s,
          family = "binomial",
          lambda = lambda_use,
          asparse = asparse_use
        )
        coef_mat[, k] <- .feature_abs_coefs_sparsegl(fit_k, s = lambda_use)[inv_ord]
      }
      return(rowMaxs(coef_mat) > bootstrap_threshold)
    }

    fit <- sparsegl::sparsegl(
      x = x_s,
      y = y,
      group = grp_s,
      family = family,
      lambda = lambda_use,
      asparse = asparse_use
    )
    .feature_abs_coefs_sparsegl(fit, s = lambda_use)[inv_ord] > bootstrap_threshold
  }
}

#' Build a Data-Driven Lambda Grid
#'
#' Runs `glmnet` on the full training data to obtain a calibrated penalty
#' sequence, returning the result as a pre-expanded `data.frame` suitable for
#' use as the `lambda_grid` argument of [stabl_fit()].
#'
#' Choosing lambda values manually is error-prone: too large a lambda selects
#' nothing; too small a lambda selects everything.  This function delegates
#' the sequence computation to `glmnet`'s own warm-start path algorithm,
#' which guarantees the grid spans from near-zero sparsity to near-full
#' sparsity for the given data, family, and alpha.
#'
#' For elastic-net models, supplying a vector of `l1_ratio` values causes the
#' function to fit a separate path per alpha, row-bind the resulting grids,
#' and add an `alpha` column so that [make_glmnet_adapter()] can dispatch
#' correctly.  This mirrors `auto_mode_lambda_grid()` in the Python reference
#' implementation.
#'
#' @param x Numeric matrix of predictors (samples by features).
#' @param y Outcome vector.  For `"gaussian"` provide a numeric vector; for
#'   `"binomial"`/`"multinomial"` provide a factor or 0/1 integer vector;
#'   for `"cox"` provide a `survival::Surv` object.
#' @param family Character; `glmnet` response family (`"gaussian"`,
#'   `"binomial"`, `"multinomial"`, or `"cox"`).  Default `"gaussian"`.
#' @param n_lambda Positive integer; desired number of lambda values per alpha
#'   level.  Actual length may be slightly shorter if `glmnet` detects
#'   saturation early.  Default `30L`.
#' @param l1_ratio Numeric scalar, numeric vector, or `NULL`.  When `NULL`
#'   (default), a pure lasso path (alpha = 1) is used.  When a vector is
#'   supplied (e.g. `c(0.5, 0.75, 1.0)`), one path is fitted per value and
#'   the results are combined; an `alpha` column is added to the output.
#'
#' @return A `data.frame` with at least a `lambda` column.  An `alpha` column
#'   is included whenever `l1_ratio` is not `NULL`.
#'
#' @seealso [stabl_fit()] which calls this automatically when
#'   `lambda_grid = "auto"`.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(40 * 6), 40, 6,
#'              dimnames = list(paste0("s", 1:40), paste0("f", 1:6)))
#' y <- setNames(rnorm(40), rownames(x))
#' grid <- auto_lambda_grid(x, y, family = "gaussian", n_lambda = 10L)
#' head(grid)
#' @export
auto_lambda_grid <- function(
    x,
    y,
    family   = "gaussian",
    n_lambda = 30L,
    l1_ratio = NULL
) {
  .require_pkg("glmnet")
  cox_ties <- .COX_TIES

  alphas <- if (is.null(l1_ratio)) 1.0 else as.numeric(l1_ratio)
  add_alpha <- !is.null(l1_ratio)

  grids <- vector("list", length(alphas))
  for (i in seq_along(alphas)) {
    fit_tmp <- .with_glmnet_numerical_infeasibility(
      glmnet::glmnet(
        x = x,
        y = y,
        family = family,
        alpha = alphas[[i]],
        nlambda = n_lambda,
        cox.ties = cox_ties
      ),
      source = "glmnet automatic lambda path"
    )
    lam_seq <- fit_tmp$lambda
    grids[[i]] <- if (add_alpha) {
      data.frame(alpha = alphas[[i]], lambda = lam_seq)
    } else {
      data.frame(lambda = lam_seq)
    }
  }

  do.call(rbind, grids)
}

.feature_abs_coefs <- function(fit, s, family = "gaussian") {
  coef_obj <- glmnet::coef.glmnet(fit, s = s)
  row_sel <- if (family == "cox") TRUE else -1L

  if (inherits(coef_obj, "dgCMatrix")) {
    return(abs(as.numeric(coef_obj[row_sel, 1L])))
  }
  if (is.matrix(coef_obj)) {
    return(abs(as.numeric(coef_obj[row_sel, 1L])))
  }
  if (is.list(coef_obj)) {
    mats <- lapply(coef_obj, function(m) as.numeric(m[row_sel, 1L]))
    coef_mat <- do.call(cbind, mats)
    return(rowMaxs(abs(coef_mat)))
  }

  stop("Unsupported coefficient structure returned by glmnet.", call. = FALSE)
}

.feature_abs_coefs_sparsegl <- function(fit, s) {
  coef_obj <- stats::coef(fit, s = s)

  if (inherits(coef_obj, "dgCMatrix")) {
    return(abs(as.numeric(coef_obj[-1L, 1L])))
  }
  if (is.matrix(coef_obj)) {
    return(abs(as.numeric(coef_obj[-1L, 1L])))
  }

  stop("Unsupported coefficient structure returned by sparsegl.", call. = FALSE)
}

# ---- Batch coefficient extraction helpers ------------------------------------
# These extract coefficients for ALL lambdas in one call, returning a
# (p × n_lambda) numeric matrix of absolute coefficient values.

.near_lambda_match <- function(lambda_seq, fit_lambda,
                               tolerance = 100 * .Machine$double.eps) {
  col_idx <- match(lambda_seq, fit_lambda)
  miss <- is.na(col_idx)
  if (!any(miss) || length(fit_lambda) == 0L) {
    return(col_idx)
  }

  for (i in which(miss)) {
    if (!is.finite(lambda_seq[[i]])) {
      next
    }
    dist <- abs(fit_lambda - lambda_seq[[i]])
    nearest <- which.min(dist)
    if (length(nearest) == 0L || !is.finite(dist[[nearest]])) {
      next
    }
    scale <- max(1, abs(lambda_seq[[i]]), abs(fit_lambda[[nearest]]))
    if (dist[[nearest]] <= tolerance * scale) {
      col_idx[[i]] <- nearest
    }
  }

  col_idx
}

.feature_abs_coefs_batch <- function(fit, lambda_seq, family = "gaussian") {
  # Fast path: when every requested lambda is on fit$lambda, read fit$beta
  # directly instead of making ~n_lambda coef.glmnet() calls.  For on-grid
  # lambdas this is bit-identical to the per-lambda path (verified by
  # test-coef-batch-ongrid.R) but avoids repeated sparse-to-dense coercions.
  # Near-on-grid lambdas are accepted within a tight floating-point tolerance;
  # true off-grid values still fall back to per-lambda coef().
  col_idx <- .near_lambda_match(lambda_seq, fit$lambda)
  if (!anyNA(col_idx)) {
    if (is.list(fit$beta)) {
      # Multinomial: fit$beta is a list of p × n_lambda matrices (one per class,
      # no intercept row).  Return per-feature max absolute coef across classes.
      # Strip dimnames to match the unnamed matrix returned by the slow path.
      m <- Reduce(pmax, lapply(fit$beta, function(b) {
        abs(as.matrix(b[, col_idx, drop = FALSE]))
      }))
      dimnames(m) <- NULL
      m
    } else {
      # Gaussian / binomial / Cox: fit$beta is p × n_lambda (no intercept row).
      # Strip dimnames to match the slow path (do.call(cbind, unnamed vectors)).
      m <- abs(as.matrix(fit$beta[, col_idx, drop = FALSE]))
      dimnames(m) <- NULL
      m
    }
  } else {
    col_list <- lapply(lambda_seq, function(s) {
      .feature_abs_coefs(fit = fit, s = s, family = family)
    })
    do.call(cbind, col_list)
  }
}

.feature_abs_coefs_sparsegl_batch <- function(fit, lambda_seq) {
  col_list <- lapply(lambda_seq, function(s) .feature_abs_coefs_sparsegl(fit = fit, s = s))
  do.call(cbind, col_list)  # p x n_lambda numeric matrix
}

# ---- Internal batch adapters -------------------------------------------------
# Each factory returns a function with signature:
#   function(x, y, lambda_grid) -> logical matrix (n_features × n_lambdas)
#
# By fitting the full lambda path in one glmnet/sparsegl call per bootstrap
# (instead of one call per lambda), the bootstrap loop in stabl_fit() reduces
# from n_bootstraps × n_lambdas model fits to just n_bootstraps fits.

.make_glmnet_batch_adapter <- function(family, alpha_fixed, bootstrap_threshold) {
  cox_ties <- .COX_TIES
  function(x, y, lambda_grid) {
    n_lambdas  <- nrow(lambda_grid)
    n_features <- ncol(x)
    result     <- matrix(FALSE, nrow = n_features, ncol = n_lambdas)

    has_alpha_col <- "alpha" %in% names(lambda_grid)

    if (has_alpha_col && is.null(alpha_fixed)) {
      # Elastic net with varying alpha: one path call per unique alpha
      for (a in unique(lambda_grid[["alpha"]])) {
        row_idx           <- which(lambda_grid[["alpha"]] == a)
        lambda_seq        <- lambda_grid[["lambda"]][row_idx]
        fit <- .with_glmnet_numerical_infeasibility(
          glmnet::glmnet(
            x,
            y,
            family = family,
            alpha = a,
            lambda = sort(lambda_seq, decreasing = TRUE),
            cox.ties = cox_ties
          ),
          source = "glmnet bootstrap path"
        )
        result[, row_idx] <- .feature_abs_coefs_batch(
          fit, lambda_seq, family = family
        ) > bootstrap_threshold
      }
    } else {
      # Single alpha: use alpha_fixed if given, lasso (1.0) as default
      alpha      <- if (!is.null(alpha_fixed)) alpha_fixed else 1.0
      lambda_seq <- lambda_grid[["lambda"]]
      fit <- .with_glmnet_numerical_infeasibility(
        glmnet::glmnet(
          x,
          y,
          family = family,
          alpha = alpha,
          lambda = sort(lambda_seq, decreasing = TRUE),
          cox.ties = cox_ties
        ),
        source = "glmnet bootstrap path"
      )
      result     <- .feature_abs_coefs_batch(fit, lambda_seq, family = family) > bootstrap_threshold
    }

    result
  }
}

.make_adaptive_lasso_batch_adapter <- function(family, gamma, epsilon,
                                               bootstrap_threshold) {
  cox_ties <- .COX_TIES
  function(x, y, lambda_grid) {
    lambda_seq <- lambda_grid[["lambda"]]

    # Ridge initialization to compute adaptive penalty weights
    init_fit <- .with_glmnet_numerical_infeasibility(
      glmnet::glmnet(
        x,
        y,
        family = family,
        alpha = 0,
        nlambda = 30L,
        cox.ties = cox_ties
      ),
      source = "glmnet adaptive-lasso ridge initialization"
    )
    init_lambda    <- tail(init_fit$lambda, n = 1L)
    init_scores    <- .feature_abs_coefs(fit = init_fit, s = init_lambda,
                       family = family)
    penalty_factor <- 1.0 / ((init_scores + epsilon) ^ gamma)

    # Fit adaptive lasso across the full lambda path (one call)
    fit <- .with_glmnet_numerical_infeasibility(
      glmnet::glmnet(
        x,
        y,
        family = family,
        alpha = 1,
        lambda = sort(lambda_seq, decreasing = TRUE),
        penalty.factor = penalty_factor,
        cox.ties = cox_ties
      ),
      source = "glmnet adaptive-lasso bootstrap path"
    )
    coef_mat <- .feature_abs_coefs_batch(fit, lambda_seq, family = family)
    coef_mat > bootstrap_threshold
  }
}

.make_sgl_batch_adapter <- function(family, feature_groups, alpha_fixed,
                                    bootstrap_threshold) {
  .stop_if_cox_sgl(family)
  groups   <- as.integer(as.factor(feature_groups))
  # sparsegl requires features sorted by group (non-decreasing group index).
  # sort_ord/inv_ord/grp_s depend only on feature_groups, so hoist them here
  # and reuse across every bootstrap call.  Only x_s stays inside the closure
  # because x is the per-bootstrap sample.
  sort_ord <- order(groups)
  inv_ord  <- order(sort_ord)
  grp_s    <- groups[sort_ord]

  function(x, y, lambda_grid) {
    n_lambdas  <- nrow(lambda_grid)
    n_features <- ncol(x)
    result     <- matrix(FALSE, nrow = n_features, ncol = n_lambdas)

    has_alpha_col <- "alpha" %in% names(lambda_grid)

    # Build list of (alpha, row_idx) pairs for batching
    alpha_batches <- if (!is.null(alpha_fixed)) {
      list(list(alpha = as.numeric(alpha_fixed), row_idx = seq_len(n_lambdas)))
    } else if (has_alpha_col) {
      lapply(unique(lambda_grid[["alpha"]]), function(a) {
        list(alpha = a, row_idx = which(lambda_grid[["alpha"]] == a))
      })
    } else {
      list(list(alpha = .SGL_DEFAULT_ALPHA, row_idx = seq_len(n_lambdas)))
    }

    for (ab in alpha_batches) {
      row_idx    <- ab$row_idx
      lambda_seq <- lambda_grid[["lambda"]][row_idx]
      asparse    <- ab$alpha

      x_s <- x[, sort_ord, drop = FALSE]

      if (family == "multinomial") {
        if (!is.factor(y)) y <- factor(y)
        levs     <- levels(y)
        coef_acc <- matrix(0.0, nrow = n_features, ncol = length(row_idx))
        for (k in seq_along(levs)) {
          yk    <- as.integer(y == levs[[k]])
          fit_k <- sparsegl::sparsegl(x = x_s, y = yk, group = grp_s,
                                      family  = "binomial",
                                      lambda  = sort(lambda_seq, decreasing = TRUE),
                                      asparse = asparse)
          coef_k_sorted <- .feature_abs_coefs_sparsegl_batch(fit_k, lambda_seq)
          coef_k        <- coef_k_sorted[inv_ord, , drop = FALSE]
          coef_acc      <- pmax(coef_acc, coef_k)
        }
        result[, row_idx] <- coef_acc > bootstrap_threshold
      } else {
        fit      <- sparsegl::sparsegl(x = x_s, y = y, group = grp_s,
                                       family  = family,
                                       lambda  = sort(lambda_seq, decreasing = TRUE),
                                       asparse = asparse)
        coef_mat          <- .feature_abs_coefs_sparsegl_batch(fit, lambda_seq)
        result[, row_idx] <- coef_mat[inv_ord, , drop = FALSE] > bootstrap_threshold
      }
    }

    result
  }
}
