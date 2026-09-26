# knockoff_mvr.R — MVR (minimum-variance-reconstructability) S-matrix solver
#
# R-native port of knockpy/mrc.py (Spector & Janson 2022, arXiv 2011.14625).
# This implements coordinate descent over diag(S) to minimize the MVR loss:
#
#   L(S) = trace((2*Sigma - S)^{-1}) + trace(S^{-1})
#        = sum(1/eig(2*Sigma - S)) + sum(1/diag(S))
#
# The resulting S is then passed to knockoff::create.gaussian(X, mu, Sigma,
# diag_s=S) for the actual model-X knockoff sampling.
#
# Every infeasibility below raises stablr_mvr_infeasible (a
# stablr_numerical_infeasibility) so that make_knockoff_mvr_features() can
# fall back to equicorrelated knockoffs instead of handing a degenerate S to
# knockoff::create.gaussian().
#
# All functions are internal (not exported).  They are called by
# make_knockoff_mvr_features() in R/artificial_features.R.
#
# Reference:
#   Spector & Janson (2022). "Powerful Knockoffs via Minimizing
#   Reconstructability." arXiv:2011.14625.

# ── Rank-1 Cholesky up/down-date ─────────────────────────────────────────────
#
# Port of knockpy/mrc.py::cholupdate (the pure-numpy fallback).
# Operates on the UPPER-triangular Cholesky factor R of a matrix V:
#   add = TRUE  → returns chol(V + x x')
#   add = FALSE → returns chol(V - x x')
#
# Arguments:
#   R   : p × p upper-triangular matrix (chol(V)).
#   x   : numeric vector of length p.
#   add : logical; TRUE = rank-1 update, FALSE = rank-1 downdate.
#
# Returns the updated R (a new matrix; R is not modified in place).
#
# Algorithm: Givens rotations, transcribed column-by-column from the Python.
# x is consumed destructively within this function (local copy).
.mvr_cholupdate <- function(R, x, add) {
  p <- length(x)
  for (k in seq_len(p)) {
    Rkk <- R[k, k]
    xk  <- x[k]
    r <- .mvr_cholupdate_radicand(Rkk, xk, add = add, pivot = k)
    c <- r / Rkk
    s <- xk / Rkk
    R[k, k] <- r
    if (k < p) {
      tail_idx <- (k + 1L):p
      R_tail   <- R[k, tail_idx]   # save OLD R[k, tail] before update
      x_tail   <- x[tail_idx]
      if (add) {
        R[k, tail_idx] <- (R_tail + s * x_tail) / c
      } else {
        R[k, tail_idx] <- (R_tail - s * x_tail) / c
      }
      # x update uses the ALREADY-UPDATED R[k, tail_idx] — this matches knockpy:
      #   x[k+1:] = c * x[k+1:] - s * R[k, k+1:]   (after R was updated above)
      x[tail_idx] <- c * x_tail - s * R[k, tail_idx]
    }
  }
  R
}

# The downdate radicand Rkk^2 - xk^2 is only representable when |xk| <= |Rkk|.
# A negative value must not reach sqrt(): the resulting NaN poisons the whole
# factor and every later coordinate silently returns garbage. Values that are
# negative only at the scale of double precision sit on the PSD boundary and are
# reported separately from a genuine magnitude violation, because both leave
# no usable pivot.
.mvr_cholupdate_radicand <- function(Rkk, xk, add, pivot) {
  if (!is.finite(Rkk) || Rkk <= 0 || !is.finite(xk)) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR Cholesky factor lost a usable pivot at position ", pivot,
        "; the coordinate descent cannot continue."
      ),
      source = "MVR Cholesky rank-1 update",
      pivot = as.integer(pivot)
    )
  }
  if (add) return(sqrt(Rkk^2 + xk^2))

  radicand <- Rkk^2 - xk^2
  if (radicand < 0) {
    # The two squares each carry ~2 ulp, so near cancellation leaves a residue
    # of order 4 * eps * max(Rkk^2, xk^2) that is not a real violation.
    roundoff <- 8 * .Machine$double.eps * max(Rkk^2, xk^2)
    if (radicand < -roundoff) {
      .abort_numerical_infeasibility(
        "stablr_mvr_infeasible",
        paste0(
          "The MVR Cholesky downdate is infeasible at position ", pivot,
          " (|xk| exceeds |Rkk| by ", format(sqrt(abs(radicand))), ")."
        ),
        source = "MVR Cholesky downdate",
        pivot = as.integer(pivot),
        radicand = radicand
      )
    }
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR Cholesky downdate reached the positive-semidefinite boundary ",
        "at position ", pivot, " (radicand ", format(radicand), "); no ",
        "strictly positive pivot remains."
      ),
      source = "MVR Cholesky downdate",
      pivot = as.integer(pivot),
      radicand = radicand
    )
  }
  sqrt(radicand)
}


# ── MVR loss ─────────────────────────────────────────────────────────────────
#
# L(S) = sum(1 / (eig(2*Sigma - S) + smoothing)) + sum(1 / (diag(S) + smoothing))
# Returns Inf when S or (2*Sigma - S) has a non-positive eigenvalue.
.mvr_loss <- function(Sigma, S, smoothing = 0) {
  s_diag    <- diag(S)
  eigs_diff <- eigen(2 * Sigma - S, symmetric = TRUE, only.values = TRUE)$values
  if (min(s_diag) < 0 || min(eigs_diff) < 0) return(Inf)
  sum(1 / (eigs_diff + smoothing)) + sum(1 / (s_diag + smoothing))
}

# The convergence test computes improvement = previous_loss - loss. An infinite
# loss therefore yields improvement = -Inf, which the test reads as converged
# and the final PSD guards then shrink into a plausible-looking S. A
# non-finite loss is divergence and must abort before it can be mistaken for
# convergence.
.mvr_loss_checked <- function(Sigma, S, smoothing = 0) {
  loss <- .mvr_loss(Sigma, S, smoothing = smoothing)
  if (!is.finite(loss)) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR loss is not finite; coordinate descent diverged instead of ",
        "converging."
      ),
      source = "MVR loss evaluation",
      smoothing = smoothing
    )
  }
  loss
}


# ── Per-coordinate quadratic solve ───────────────────────────────────────────
#
# Port of knockpy/mrc.py::_solve_mvr_quadratic (scalar, no group support).
# We minimise f(delta) = 1/(sj + delta) - delta*cn / (1 - delta*cd)
# subject to delta in the open PSD-feasible interval (-sj, 1/cd).
#
# Arguments:
#   cn, cd   : scalars from the Cholesky solves (cn < 0, cd > 0).
#   sj       : current value of S[j,j].
#   smoothing: passed through from the outer solver (default 0).
#
# Returns: optimal delta (scalar).
.solve_mvr_quadratic <- function(cn, cd, sj, smoothing = 0) {
  # Quadratic coefficients for the first-order condition
  coef2 <- -cn - cd^2
  coef1 <- 2 * (-cn * (sj + smoothing) + cd)
  coef0 <- -cn * (sj + smoothing)^2 - 1

  roots <- .mvr_real_roots(coef0, coef1, coef2)

  # Keep roots inside the PSD-feasible open interval (-sj, 1/cd). The slack is
  # absolute, so when both bounds are themselves near 1e-12 the two can cross.
  lo <- -sj + .Machine$double.eps * 1e4
  hi <- 1 / cd - .Machine$double.eps * 1e4
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR per-coordinate feasible interval is empty or inverted: (",
        format(lo), ", ", format(hi), ")."
      ),
      source = "MVR per-coordinate feasible interval",
      cn = cn, cd = cd, sj = sj
    )
  }
  roots <- roots[roots > lo & roots < hi]

  if (length(roots) == 0L) {
    .mvr_abort_infeasible_quadratic(cn = cn, cd = cd, sj = sj)
  }

  # Evaluate objective and pick the minimiser
  f_obj <- function(delta) {
    1 / (sj + delta) - (delta * cn) / (1 - delta * cd)
  }
  losses <- vapply(roots, f_obj, numeric(1L))
  roots[which.min(losses)]
}

# polyroot() finds genuinely real roots with an imaginary residue proportional to
# the root itself, so the real-root test has to be relative to the root
# magnitude. An absolute floor additionally admits genuinely complex roots of
# small magnitude, which are then used as if they were stationary points.
.mvr_real_roots <- function(coef0, coef1, coef2) {
  all_roots <- polyroot(c(coef0, coef1, coef2))
  modulus <- Mod(all_roots)
  is_real <- abs(Im(all_roots)) <=
    1e-10 * pmax(modulus, .Machine$double.xmin)
  Re(all_roots[is_real])
}

# The single degradation .solve_mvr_quadratic() can report. It carries its own
# condition class plus a restart holding the substitute delta, so solve_mvr()
# recovers from exactly this condition and nothing else: an unrelated warning
# raised anywhere inside .solve_mvr_quadratic() stays an ordinary warning and
# reaches the user instead of being replaced by a non-optimal value. The
# stablr_mvr_degraded_delta restart name is the pair to this condition class.
.mvr_abort_infeasible_quadratic <- function(cn, cd, sj) {
  withRestarts(
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible_quadratic",
      paste0(
        "The MVR per-coordinate quadratic has no stationary point inside its ",
        "feasible interval."
      ),
      source = "MVR per-coordinate quadratic",
      cn = cn,
      cd = cd,
      sj = sj
    ),
    stablr_mvr_degraded_delta = function(value) value
  )
}


# ── PSD guards ───────────────────────────────────────────────────────────────
# Port of knockpy/utilities.py::shift_until_PSD and scale_until_PSD.

.shift_until_psd <- function(S, tol = 1e-5) {
  min_eig <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < tol) {
    diag(S) <- diag(S) + (tol - min_eig)
  }
  S
}

# The binary search only ever evaluates gamma in (0, 1), so `lo` keeps its
# initial value of 0 when 2 * Sigma is itself below the feasibility tolerance.
# The result is then S = 0, whose strictly-zero diagonal makes
# knockoff::create.gaussian() emit power-free knockoffs that are still selected
# and consumed by compute_fdp_plus(). The post-condition mirrors the one
# .knockpy_equicorrelated_s() applies to the equicorrelated S matrix.
.scale_until_psd <- function(Sigma, S, tol = 1e-5, num_iter = 10L) {
  # Ensure S itself is PSD first
  S <- .shift_until_psd(S, tol)

  # Binary search for the largest gamma in [0,1] such that 2*Sigma - gamma*S is PSD
  lo <- 0
  hi <- 1
  for (j in seq_len(num_iter)) {
    gamma <- (lo + hi) / 2
    V     <- 2 * Sigma - gamma * S
    min_eig_V <- min(
      .with_linear_algebra_numerical_infeasibility(
        eigen(V, symmetric = TRUE, only.values = TRUE)$values,
        source = "MVR positive-semidefinite feasibility check",
        subclass = "stablr_mvr_infeasible"
      )
    )
    if (min_eig_V >= tol) lo <- gamma else hi <- gamma
  }

  scaled <- lo * S
  diag_s <- diag(scaled)
  minimum_eigenvalue <- min(
    .with_linear_algebra_numerical_infeasibility(
      eigen(2 * Sigma - scaled, symmetric = TRUE, only.values = TRUE)$values,
      source = "MVR positive-semidefinite feasibility check",
      subclass = "stablr_mvr_infeasible"
    )
  )
  if (any(!is.finite(diag_s)) || any(diag_s <= 0) ||
      !is.finite(minimum_eigenvalue) || minimum_eigenvalue < tol) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR S-matrix failed numerical validation: the rescaled diagonal ",
        "must be strictly positive and 2 * Sigma - gamma * S must stay ",
        "positive definite at the search tolerance."
      ),
      source = "MVR positive-semidefinite rescale",
      psd_scale = lo,
      tolerance = tol,
      minimum_eigenvalue = minimum_eigenvalue
    )
  }

  list(S = scaled, gamma = lo, minimum_eigenvalue = minimum_eigenvalue)
}


# ── Solve-result guards ──────────────────────────────────────────────────────
#
# base::forwardsolve() and base::backsolve() raise only on an exactly-zero
# pivot, and a triangular matrix with a non-zero diagonal is never singular, so
# the silent path is not a singular triangular solve: it is a factor that no
# longer factors the matrix it is supposed to factor. Nothing in the arithmetic
# downstream of a rank-1 update can detect that, because cd and cn both scale
# with the factor and every derived delta stays finite and well signed. These
# guards therefore compare the factor against the matrix it claims, and check
# the recovered quantities for the sign and finiteness that their definitions
# require.
#
# The remaining algebraic relationship, (M^{-2})_{jj} >= (M^{-1})_{jj}^2, is
# deliberately not used: it is invariant to a uniform rescaling of the factor,
# so it cannot distinguish a drifted factor from a correct one, and it degrades
# towards equality exactly as A becomes near-diagonal, which would force a
# tolerance loose enough to accept anything.

.mvr_verify_cholesky_factor <- function(
  R_upper,
  M,
  source = "MVR Cholesky factorisation"
) {
  diagonal <- diag(R_upper)
  reconstructed <- crossprod(R_upper)
  scale <- max(abs(M))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  deviation <- max(abs(reconstructed - M)) / scale
  if (any(!is.finite(diagonal)) || any(diagonal <= 0) ||
      !is.finite(deviation) ||
      deviation > sqrt(.Machine$double.eps)) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR Cholesky factor no longer factors 2 * Sigma - S (relative ",
        "reconstruction deviation ", format(deviation), ")."
      ),
      source = source,
      relative_reconstruction_deviation = unname(deviation)
    )
  }
  invisible(R_upper)
}

# For M = R' R the forward solve returns vd = L^{-1} e_j, so cd = ||vd||^2 is
# exactly (M^{-1})_{jj}, and the back solve returns vn = R^{-1} vd, so
# cn = -||vn||^2 is exactly -(M^{-2})_{jj}. Both are therefore strictly signed
# for a correct factor: cd > 0 and cn < 0, and 1 / cd — the right end of the
# per-coordinate feasible interval — must stay representable.
.mvr_coordinate_quadratic <- function(R_upper, j) {
  p <- nrow(R_upper)
  diagonal <- diag(R_upper)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR Cholesky factor is singular or non-finite; the coordinate ",
        "solves would return garbage."
      ),
      source = "MVR coordinate-descent Cholesky factor"
    )
  }

  ej <- numeric(p)
  ej[j] <- 1
  vd <- .with_linear_algebra_numerical_infeasibility(
    forwardsolve(t(R_upper), ej),
    source = "MVR coordinate forward solve",
    subclass = "stablr_mvr_infeasible"
  )
  vn <- .with_linear_algebra_numerical_infeasibility(
    backsolve(R_upper, vd),
    source = "MVR coordinate back solve",
    subclass = "stablr_mvr_infeasible"
  )
  cd <- sum(vd^2)
  cn <- -sum(vn^2)
  if (any(!is.finite(vd)) || any(!is.finite(vn)) ||
      !is.finite(cd) || cd <= 0 || !is.finite(cn) || cn >= 0 ||
      !is.finite(1 / cd) || 1 / cd <= 0) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR coordinate solves did not recover (M^{-1})_{jj} > 0 and a ",
        "strictly negative cn for coordinate ", j, "."
      ),
      source = "MVR coordinate solves",
      coordinate = as.integer(j),
      cd = unname(cd),
      cn = unname(cn)
    )
  }
  list(cd = cd, cn = cn)
}


# ── Main MVR S-matrix solver ──────────────────────────────────────────────────
#
# Coordinate descent minimising the MVR loss.
# Port of knockpy/mrc.py::_solve_mvr_ungrouped + the PSD guards from solve_mvr.
#
# @param Sigma      p × p positive-definite covariance matrix (estimated
#                   from data before calling this function).
# @param num_iter   Maximum number of coordinate-descent sweeps.  Default 50.
# @param converge_tol Convergence criterion on the decayed-exponential-moving
#                   average improvement.  Default 1e-2 (matches knockpy default
#                   for _solve_mvr_ungrouped).
# @param tol        Minimum permissible eigenvalue for the final PSD guards.
# @param smoothing  Add smoothing to eigenvalues in the loss.  Default 0.
# @param random_state Optional integer seed for the coordinate-shuffle RNG;
#                   the caller's RNG state is restored after solving.
#
# @return Numeric vector of length p — the diagonal of S (covariance-scale).
#         Drop-in for the `diag_s` argument of knockoff::create.gaussian().
solve_mvr <- function(Sigma,
                      num_iter      = 50L,
                      converge_tol  = 1e-2,
                      tol           = 1e-5,
                      smoothing     = 0,
                      random_state  = NULL) {
  if (!is.null(random_state)) {
    return(.with_local_seed(
      random_state,
      solve_mvr(
        Sigma = Sigma,
        num_iter = num_iter,
        converge_tol = converge_tol,
        tol = tol,
        smoothing = smoothing,
        random_state = NULL
      )
    ))
  }

  if (!is.matrix(Sigma) || !is.numeric(Sigma) || nrow(Sigma) != ncol(Sigma) ||
      nrow(Sigma) < 1L || anyNA(Sigma) || any(!is.finite(Sigma))) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      "`Sigma` must be a finite numeric square matrix.",
      source = "MVR S-matrix solver"
    )
  }

  p <- nrow(Sigma)

  # Initialise S = min_eig(Sigma) * I  (same as knockpy)
  min_eig <- min(
    .with_linear_algebra_numerical_infeasibility(
      eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values,
      source = "MVR covariance decomposition",
      subclass = "stablr_mvr_infeasible"
    )
  )
  # MVR is defined for a positive-definite Sigma only. Shifting a singular
  # Sigma to the feasibility tolerance silently answers a different question
  # and returns an S that is unrelated to the input.
  if (min_eig <= 0) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "MVR requires a positive-definite covariance; `Sigma` has minimum ",
        "eigenvalue ", format(min_eig), "."
      ),
      source = "MVR S-matrix solver",
      minimum_eigenvalue = min_eig
    )
  }
  S <- diag(min_eig, p)

  # Initialise upper-triangular Cholesky factor of (2*Sigma - S)
  initial_matrix <- 2 * Sigma - S + smoothing * diag(p)
  R_upper <- .with_linear_algebra_numerical_infeasibility(
    chol(initial_matrix),
    source = "MVR Cholesky factorisation",
    subclass = "stablr_mvr_infeasible"
  )
  .mvr_verify_cholesky_factor(
    R_upper, initial_matrix, source = "MVR Cholesky factorisation"
  )

  loss              <- Inf
  decayed_improvement <- 10  # large initial value → don't stop early on iter 0
  degraded_coordinates <- 0L
  coordinate_updates  <- 0L

  for (i in seq_len(num_iter)) {
    # Randomise coordinate order each sweep (matches knockpy np.random.shuffle)
    perm <- sample.int(p)

    for (j_idx in perm) {
      # -- 1. Compute cd and cn via Cholesky forward/back solves ----------------
      solved <- .mvr_coordinate_quadratic(R_upper, j_idx)
      cd     <- solved$cd
      cn     <- solved$cn

      # -- 2. Solve per-coordinate quadratic ------------------------------------
      # Only stablr_mvr_infeasible_quadratic is degraded here. The substitute is
      # 1% of the interval width, which is not the minimiser, so the
      # substitution is counted and reported once the descent finishes.
      delta  <- withCallingHandlers(
        .solve_mvr_quadratic(cn = cn, cd = cd, sj = S[j_idx, j_idx],
                             smoothing = smoothing),
        stablr_mvr_infeasible_quadratic = function(cnd) {
          degraded_coordinates <<- degraded_coordinates + 1L
          invokeRestart(
            "stablr_mvr_degraded_delta",
            (1 / cd - S[j_idx, j_idx]) * 0.01
          )
        }
      )

      # -- 3. Rank-1 Cholesky update of R_upper --------------------------------
      x_vec          <- numeric(p)
      x_vec[j_idx]   <- sqrt(abs(delta))
      if (delta > 0) {
        # s_j increases → M = 2Sigma - S decreases → downdate
        R_upper <- .mvr_cholupdate(R_upper, x_vec, add = FALSE)
      } else if (delta < 0) {
        # s_j decreases → M increases → update
        R_upper <- .mvr_cholupdate(R_upper, x_vec, add = TRUE)
      }
      .mvr_verify_cholesky_factor_diagonal(R_upper)

      # -- 4. Update S ----------------------------------------------------------
      S[j_idx, j_idx] <- S[j_idx, j_idx] + delta
      coordinate_updates <- coordinate_updates + 1L
    }

    # -- Convergence check (decayed EMA of improvement) -----------------------
    prev_loss <- loss
    loss      <- .mvr_loss_checked(Sigma, S, smoothing = smoothing)

    if (i > 1L) {
      improvement         <- prev_loss - loss
      decayed_improvement <- decayed_improvement / 10 + 9 * improvement / 10
    }
    if (decayed_improvement < converge_tol) break
  }

  .mvr_verify_cholesky_factor(
    R_upper,
    2 * Sigma - diag(diag(S)) + smoothing * diag(p),
    source = "MVR Cholesky factor drift check"
  )

  if (degraded_coordinates > 0L) {
    warning(
      "knockoff_mvr: ", degraded_coordinates, " of ", coordinate_updates,
      " coordinate updates had no feasible stationary point; a 1%-of-bound ",
      "heuristic was substituted, so the returned S is not the MVR minimiser.",
      call. = FALSE
    )
  }

  # -- Final PSD guards (mirror knockpy's solve_mvr post-processing) -----------
  S <- .shift_until_psd(S, tol = tol)
  result <- .scale_until_psd(Sigma, S, tol = tol, num_iter = 10L)
  S <- result$S

  diag(S)  # return as a vector (drop-in for create.gaussian diag_s)
}

# A Cholesky factor of a positive-definite matrix has a strictly positive,
# finite diagonal. Checking it after every rank-1 update costs O(p) and catches
# a zero pivot or a NaN before the next coordinate reads it as a solve input.
.mvr_verify_cholesky_factor_diagonal <- function(R_upper) {
  diagonal <- diag(R_upper)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    .abort_numerical_infeasibility(
      "stablr_mvr_infeasible",
      paste0(
        "The MVR Cholesky factor lost a strictly positive pivot during ",
        "coordinate descent."
      ),
      source = "MVR coordinate-descent Cholesky factor"
    )
  }
  invisible(R_upper)
}
