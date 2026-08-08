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
    r <- if (add) sqrt(Rkk^2 + xk^2) else sqrt(Rkk^2 - xk^2)
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

  # Roots of the quadratic (2 roots for a degree-2 polynomial)
  roots <- Re(polyroot(c(coef0, coef1, coef2)))  # polyroot uses increasing degree
  # Keep only real roots (imaginary part exactly 0 after polyroot — not always
  # guaranteed due to floating-point, so check the complex polyroot output)
  all_roots <- polyroot(c(coef0, coef1, coef2))
  roots <- Re(all_roots[abs(Im(all_roots)) < 1e-10 * abs(Re(all_roots)) + 1e-14])

  # Keep roots inside the PSD-feasible open interval (-sj, 1/cd)
  lo <- -sj + .Machine$double.eps * 1e4
  hi <- 1 / cd - .Machine$double.eps * 1e4
  roots <- roots[roots > lo & roots < hi]

  if (length(roots) == 0L) {
    # Fallback: pick the interior point that minimises f; this should not
    # happen with well-conditioned inputs but mirrors knockpy's RuntimeError
    # by returning the centre of the feasible interval.
    warning("MVR quadratic: no feasible root; returning interval midpoint.",
            call. = FALSE)
    return((lo + hi) / 2)
  }

  # Evaluate objective and pick the minimiser
  f_obj <- function(delta) {
    1 / (sj + delta) - (delta * cn) / (1 - delta * cd)
  }
  losses <- vapply(roots, f_obj, numeric(1L))
  roots[which.min(losses)]
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

.scale_until_psd <- function(Sigma, S, tol = 1e-5, num_iter = 10L) {
  # Ensure S itself is PSD first
  S <- .shift_until_psd(S, tol)

  # Binary search for the largest gamma in [0,1] such that 2*Sigma - gamma*S is PSD
  lo <- 0
  hi <- 1
  for (j in seq_len(num_iter)) {
    gamma <- (lo + hi) / 2
    V     <- 2 * Sigma - gamma * S
    # Try Cholesky: if it succeeds, gamma is feasible
    feasible <- tryCatch(
      .with_linear_algebra_numerical_infeasibility(
        {
          min_eig_V <- min(
            eigen(V, symmetric = TRUE, only.values = TRUE)$values
          )
          min_eig_V >= tol
        },
        source = "MVR positive-semidefinite feasibility check",
        subclass = "stablr_mvr_infeasible"
      ),
      stablr_numerical_infeasibility = function(e) FALSE
    )
    if (feasible) lo <- gamma else hi <- gamma
  }
  list(S = lo * S, gamma = lo)
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

  p    <- nrow(Sigma)
  inds <- seq_len(p)

  # Initialise S = min_eig(Sigma) * I  (same as knockpy)
  min_eig <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  # Guard: if Sigma is not PD, shift it
  if (min_eig <= 0) {
    Sigma   <- Sigma + (-min_eig + tol) * diag(p)
    min_eig <- tol
  }
  S <- diag(min_eig, p)

  # Initialise upper-triangular Cholesky factor of (2*Sigma - S)
  R_upper <- chol(2 * Sigma - S + smoothing * diag(p))

  loss              <- Inf
  decayed_improvement <- 10  # large initial value → don't stop early on iter 0

  for (i in seq_len(num_iter)) {
    # Randomise coordinate order each sweep (matches knockpy np.random.shuffle)
    perm <- sample.int(p)

    for (j_idx in perm) {
      # -- 1. Compute cd and cn via Cholesky forward/back solves ----------------
      ej     <- numeric(p); ej[j_idx] <- 1
      # vd = L^{-1} e_j  (forward sub with lower-tri L = t(R_upper))
      vd     <- forwardsolve(t(R_upper), ej)
      cd     <- sum(vd^2)                      # (M^{-1})_{jj}

      # vn = L^{-T} vd = R_upper^{-1} vd  (back sub)
      vn     <- backsolve(R_upper, vd)
      cn     <- -sum(vn^2)                     # -(M^{-2})_{jj}

      # -- 2. Solve per-coordinate quadratic ------------------------------------
      delta  <- tryCatch(
        .solve_mvr_quadratic(cn = cn, cd = cd, sj = S[j_idx, j_idx],
                             smoothing = smoothing),
        warning = function(w) {
          # Degenerate case: nudge by a small positive amount toward upper bound
          (1 / cd - S[j_idx, j_idx]) * 0.01
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

      # -- 4. Update S ----------------------------------------------------------
      S[j_idx, j_idx] <- S[j_idx, j_idx] + delta
    }

    # -- Convergence check (decayed EMA of improvement) -----------------------
    prev_loss <- loss
    loss      <- .mvr_loss(Sigma, S, smoothing = smoothing)

    if (i > 1L) {
      improvement         <- prev_loss - loss
      decayed_improvement <- decayed_improvement / 10 + 9 * improvement / 10
    }
    if (decayed_improvement < converge_tol) break
  }

  # -- Final PSD guards (mirror knockpy's solve_mvr post-processing) -----------
  S <- .shift_until_psd(S, tol = tol)
  result <- .scale_until_psd(Sigma, S, tol = tol, num_iter = 10L)
  S <- result$S

  diag(S)  # return as a vector (drop-in for create.gaussian diag_s)
}
