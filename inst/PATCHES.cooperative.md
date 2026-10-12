# PATCHES.cooperative — divergence log vs CRAN multiview 1.0

The files `R/cooperative-*.R` are a vendored, near-verbatim fork of CRAN
`multiview` 1.0 (GPL-2, see `inst/COPYING.cooperative`). This file is the
machine-readable log of every deliberate divergence from upstream so that a
future upstream sync (or a CRAN reviewer) can audit the fork in minutes.

Format: one entry per divergence, with the upstream file, the vendored file,
the approximate location, the nature of the change, and whether it alters
statistical behavior under stablr's usage of the backend. stablr always calls
this backend with **unit case weights and no offsets** (the public wrappers in
`R/cooperative-backend.R` never pass `weights` or `offset`), so "No" below
means bit-identical results for every code path stablr exercises today.

---

## D1. Function renames: `.mv_` prefix

- **Upstream:** every exported function (`multiview()`, `cv.multiview()`,
  `multiview.path()`, `build_predmat()`, ...).
- **Vendored:** all files `R/cooperative-*.R`.
- **Location:** function definitions throughout (e.g. `cooperative-multiview.R:208`
  `.mv_multiview`, `cooperative-cv.multiview.R:173` `.mv_cv_multiview`).
- **Nature:** mechanical rename with `.mv_` prefix to keep the vendored
  namespace private (the functions are internal; S3 methods
  `predict.multiview`, `coef.multiview`, `plot.multiview`, `family.multiview`
  and their `cv.multiview` counterparts are registered by hand in `NAMESPACE`).
- **Statistical behavior under stablr usage:** No (pure rename).

## D2. Family restriction to gaussian/binomial with hard stops

- **Upstream:** `R/multiview.R` (supports gaussian, binomial, poisson, cox).
- **Vendored:** `R/cooperative-multiview.R:222-244`.
- **Nature:** `poisson`, `quasipoisson`, `cox` and any other family now raise
  `stop("Native cooperative learning in stablr currently supports gaussian and
  binomial only.")`. The cox code paths (`multiview.cox.path`, `cv.coxnet`,
  `build_predmat_coxnetlist`) are unreachable dead code retained for parity.
- **Statistical behavior under stablr usage:** No (stablr only requests
  gaussian/binomial; unsupported families were never reachable).

## D3. Upstream typo fix: `length(row)` -> `length(rho)`

- **Upstream:** `R/multiview.R` (argument validation of `rho`).
- **Vendored:** `R/cooperative-multiview.R:247-250`
  (`if (length(rho) != 1L || rho < 0) stop(...)`).
- **Nature:** upstream 1.0 validates `length(row)` — an undefined variable
  (typo for `rho`), which makes the check throw `object 'row' not found`
  whenever the check runs. Fixed to `length(rho)`.
- **Statistical behavior under stablr usage:** No (fixes a crash; valid scalar
  `rho` inputs behave identically).

## D4. (B1) Offset subsetting in `build_predmat`

- **Upstream:** `R/build_predmat.R:26` — same bug (unreported upstream).
- **Vendored:** `R/cooperative-build_predmat.R:26`.
- **Nature:** the fold loop had `offset_sub <- rep(offset_sub[which, drop =
  FALSE])`, which (a) self-references `offset_sub` before it exists (the
  argument is `offset` — a hard error, or silent use of a stray global), and
  (b) wraps the subset in a spurious `rep()`. Fixed to
  `offset_sub <- offset[which, drop = FALSE]`, mirroring the `x_sub_list`
  row-subsetting two lines above and the training-side idiom
  `offset[!which, drop = FALSE]` in the `cv.multiview` fold loop: `offset`
  entering `build_predmat` is a length-`nobs` vector (after the `drop()`
  normalization at the top of the function), not a per-view list, and
  `predict.multiview` forwards `newoffset` to glmnet, which expects one entry
  per prediction row.
- **Statistical behavior under stablr usage:** No (the buggy branch is only
  reached with a non-NULL offset; stablr never passes offsets, in which case
  `offset_sub` is `NULL` before and after the fix).

## D5. (B2) Halved-step warm-start residual weighting

- **Upstream:** `R/multiview.path.R:766` — same bug (unreported upstream).
- **Vendored:** `R/cooperative-multiview.path.R:769` (inside `multiview.fit`,
  the `if (halved)` rebuild of `fit$warm_fit$r`).
- **Nature:** the step-halving branch rebuilt the warm-start residual as
  `w_std * (target - g_eta) + g_offset` — the working residual weighted but
  the offset added unweighted — while the non-halved branch (line 658) and the
  initial-residual convention use `w_std * (target - g_eta + g_offset)`. The
  halved branch now uses the same expression, so the offset enters the
  weighted residual consistently.
- **Statistical behavior under stablr usage:** No. With no offset,
  `g_offset` is a zero vector (`offset <- numeric(nobs)` upstream of the
  loop), so both expressions are identical for ANY weights, including the
  unit weights stablr always uses. The divergence only matters for weighted
  fits with a non-zero offset that also trigger step-halving.

## D6. (B3) Pass `lower.limits` to `elnet.fit`

- **Upstream:** `R/multiview.path.R:680` — same asymmetry (unreported upstream).
- **Vendored:** `R/cooperative-multiview.path.R:682` (the `elnet.fit` call in
  `multiview.fit`).
- **Nature:** the call passed `upper.limits = upper.limits` but not
  `lower.limits`, so the coefficient lower bounds used by the solver on the
  first (non-warm) call were the `elnet.fit` defaults (`-Inf`, mapped to
  `-control$big` unscaled) instead of the caller's `lower.limits` (already
  transformed to `-control$big * xs` in `.mv_multiview`). `lower.limits` is a
  formal of `multiview.fit` (line 468), so it is simply forwarded now:
  `lower.limits = lower.limits`. Warm-started calls are unaffected either way
  (they reuse `cl` stored in the warm-start object).
- **Statistical behavior under stablr usage:** No. stablr never sets
  coefficient limits, so the bound is `-control$big` (9.9e35) either way —
  never binding, numerically identical solutions.

---

### Reporting upstream

D3-D6 are bugs present in CRAN multiview 1.0 (`build_predmat.R:26`,
`multiview.path.R:680`, `multiview.path.R:766`); consider reporting them to
the multiview maintainers. D4-D6 are fixed here in-place; re-apply them after
any upstream sync.
