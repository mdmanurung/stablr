# Selection reliability diagnostics for STABL fits (1.0.0 spec, A1).
#
# `stabl_diagnose()` inspects a fitted `stabl_fit` object and reports, without
# ever modifying the fit or its selected support:
#   * real/artificial stability-score distributions and decoy provenance,
#   * selected real/decoy counts along the full threshold path,
#   * the selected cut-off, the full FDP+ curve, boundary solutions, ties,
#     flat minima, and low-decoy-count instability,
#   * null-like score saturation (the RF/oblique-RF leakage signature),
#   * correlation/dependence diagnostics as *assumption probes* (never as
#     proofs of permutation exchangeability or knockoff validity),
#   * a learner x decoy x outcome x dimensionality evidence profile sourced
#     from the package's benchmark/audit artifacts,
#   * machine-readable `warning_codes`, `assumption_notes`, and
#     `calibration_scope`.
#
# Status convention for the evidence profile: `validated_in_stablr` means the
# package holds evidence that the FDP+ diagnostic is approximately calibrated
# or conservative for that regime (algorithmic parity and/or simulation
# benchmarks); `empirically_concerning` means anti-conservative or otherwise
# adverse evidence exists; `not_evaluated` means no benchmark artifact covers
# the combination; `unsupported` means the combination fails closed by design.

# Benchmark-grounded evidence table. Provenance points at artifacts in
# benchmarks/ and docs/audit/2026-10-17-adversarial/. Regime descriptions are
# part of the evidence: a row never generalizes beyond its evaluated regime.
.stabl_diagnostic_evidence <- data.frame(
  stringsAsFactors = FALSE,
  learner = c(
    "lasso", "lasso", "lasso", "lasso", "lasso",
    "elastic_net", "adaptive_lasso",
    "rf", "rf", "rf",
    "oblique_rf",
    "ordinal",
    "sparse_group_lasso",
    "splsda", "splsda", "splsda",
    "cox_learners"
  ),
  decoy = c(
    "knockoff_equi", "random_permutation", "any", "any", "knockoff_equi",
    "any", "any",
    "knockoff_equi", "knockoff_equi", "random_permutation",
    "any",
    "any",
    "any",
    "random_permutation", "knockoff_equi", "knockoff_mvr",
    "any"
  ),
  family = c(
    "gaussian", "gaussian", "gaussian", "gaussian", "binomial",
    "any", "any",
    "gaussian", "gaussian", "any",
    "any",
    "ordinal",
    "any",
    "binomial", "binomial", "binomial",
    "any"
  ),
  regime = c(
    "sparse independent, linear signal, low dimension",
    "sparse independent, linear signal, low dimension",
    "pure null (no signal)",
    "block-correlated AR(0.6), weak signal, n=100 (Phase 6 D3 conditions)",
    "high-dimensional multi-cohort (MINT: K=3 x 80, p=150)",
    "linear signal (algorithmic parity); correlated blocks degrade FDP",
    "linear signal; correlated blocks degrade FDP",
    "low-dimensional nonlinear (n=120, p=10)",
    "high-dimensional nonlinear or binomial outcome",
    "any (permutation decoys with a multivariate learner)",
    "any",
    "low-dimensional ordinal designs",
    "correlated block designs",
    "multi-cohort binomial (MINT: K=3 x 80, p=150)",
    "multi-cohort binomial (MINT: K=3 x 80, p=150)",
    "multi-cohort binomial (MINT: K=3 x 80, p=150)",
    "any"
  ),
  status = c(
    "validated_in_stablr", "empirically_concerning",
    "empirically_concerning", "validated_in_stablr",
    "empirically_concerning",
    "validated_in_stablr", "validated_in_stablr",
    "validated_in_stablr", "empirically_concerning",
    "empirically_concerning",
    "empirically_concerning",
    "empirically_concerning",
    "empirically_concerning",
    "validated_in_stablr", "validated_in_stablr", "validated_in_stablr",
    "not_evaluated"
  ),
  provenance = c(
    "Python-parity fixtures (SHA-256-manifested); audit Phase 6 D1 (realized FDP 0.458 vs diagnostic 0.473)",
    "audit Phase 6 D1 (realized FDP 0.210 vs diagnostic 0.149)",
    "audit Phase 6 D2 / F-P6-1 (realized FDP 0.92-0.93; 18-33 features selected on null data)",
    "audit Phase 6 D3 (realized FDP 0.354 vs diagnostic 0.522; conservative)",
    "benchmarks/mint_comparison partA (realized FDP 0.83 vs reported 0.10; calibration 8.0)",
    "Python-parity fixtures (SHA-256-manifested); benchmarks/learner_comparison (block DGP FDP 0.53 KE / 0.40 RP)",
    "benchmarks/learner_comparison (block DGP FDP 0.48 KE / 0.29 RP)",
    "benchmarks/f5_summary.md (full mode: FDP 0.065 at nominal 0.1)",
    "benchmarks/ws2_summary.md (FDP 0.438 high-dim, 0.265 binomial; graduation NOT MET)",
    "benchmarks/f5_summary.md + ws2_summary.md (FDP 0.32-0.44); finding F-005 mechanism",
    "benchmarks/ws3_summary.md (fenced: no predeclared remedy met the FDP criterion)",
    "benchmarks/ws9_summary.md (graduation NOT MET: FDP 0.140/0.156 vs <= 0.15 bar)",
    "benchmarks/learner_comparison/guide_table.csv (block DGP FDP 0.91 KE / 0.74 RP)",
    "benchmarks/mint_comparison partA (calibration 1.01; realized FDP 0.31 vs reported 0.30)",
    "benchmarks/mint_comparison partA (calibration 0.49; conservative)",
    "benchmarks/mint_comparison/modelx_extension (median per-seed calibration 0.38; conservative)",
    "no FDP-calibration benchmark artifact"
  ),
  note = c(
    "Algorithmic parity with Python STABL plus approximately calibrated FDP+ minimum in this regime (simulation scope).",
    "FDP+ minimum anti-conservative under permutation decoys even in the independent regime.",
    "The FDP+-minimizing rule still selects sizeable panels on null data; an empty support is the only outcome a null dataset cannot produce by this rule.",
    "FDP+ minimum conservative under block correlation at the evaluated signal level; regime-dependence itself is the finding (F-P6-2). A sibling binomial block-DGP benchmark (learner_comparison: lasso+KE calibration 2.0) does not corroborate beyond these conditions.",
    "Cohort heterogeneity breaks knockoff exchangeability; the FDP+ minimum is severely anti-conservative in this regime.",
    "Parity is algorithmic, not calibration evidence; treat correlated designs as concerning.",
    "No algorithmic-parity coverage (not in the Python reference); correlated designs degrade FDP.",
    "Empirical calibration only; simulation scope, not application-data evidence.",
    "Design-scope restriction retained: do not read the FDP+ minimum as calibrated here.",
    "Anti-conservative by construction: permuted decoys do not receive the correlation leakage real null features receive.",
    "Adapter fenced by benchmark outcome; scheduled for retirement from the core API in 1.0.0.",
    "Ordinal adapter stays experimental with documented scope.",
    "Group-structured selection on correlated blocks is far outside any calibrated regime.",
    "Best-calibrated sPLS-DA arm in the MINT benchmark; single-regime simulation evidence only.",
    "Conservative in the MINT benchmark (realized FDP below reported FDP+); single-regime simulation evidence only.",
    "Conservative in the MINT model-X extension (realized FDP below reported FDP+); single-regime simulation evidence only.",
    "Cox workflows have no FDP+ calibration evidence; a validated Cox workflow is a future conditional item."
  )
)

#' Diagnose selection reliability of a fitted STABL object
#'
#' Inspects a `stabl_fit` result and returns a structured
#' `stabl_diagnostics` object: score-distribution summaries, decoy-generator
#' provenance, the full threshold path of real/decoy selection counts, FDP+
#' curve features (selected cut-off, ties, flat minima, boundary solutions,
#' low-decoy-count instability), a null-saturation probe, a score-dependence
#' assumption probe, and a learner x decoy x outcome x dimensionality
#' evidence profile with benchmark provenance.
#'
#' The diagnostic is read-only: it never changes the fit or its selected
#' support. It distinguishes a user-supplied `hard_threshold` from the
#' minimized-FDP+ rule and reports the effective threshold comparator
#' (strict `>`) and the FDP+ denominator convention
#' `(n_artificial / artificial_proportion + 1) / pmax(1, n_real)`.
#'
#' @param fit A fitted [stabl_fit()] object (including `stabl()` results,
#'   which inherit from `stabl_fit`).
#' @param saturation_threshold Numeric in \[0, 1\]; thresholds at or below
#'   this value are considered "near-zero" for the null-saturation probe.
#'   Default `0.05`.
#' @param saturation_proportion Numeric in (0, 1\]; the probe flags
#'   null-like saturation when a near-zero threshold selects at least this
#'   proportion of real *and* of artificial features. Default `0.5`.
#' @param low_decoy_count Non-negative integer; the FDP+ minimum is flagged
#'   as unstable when the artificial-feature count at the selected threshold
#'   is below this value. Default `5L`.
#' @param ... Reserved for future use; passing arguments here is an error.
#'
#' @return An object of class `stabl_diagnostics` with components:
#'   \describe{
#'     \item{score_summary}{List of distribution summaries (quantiles, mean)
#'       for real and artificial per-feature maximum scores.}
#'     \item{decoy_provenance}{The fit's `artificial_provenance` (requested
#'       vs actual generator, fallback counts), or `NULL`.}
#'     \item{threshold_path}{`data.frame` with `threshold`, `n_real`,
#'       `n_artificial`, `prop_real`, `prop_artificial`, and `fdp_plus`
#'       along the full `fdr_threshold_range`.}
#'     \item{per_lambda}{`data.frame` of per-lambda real/artificial selected
#'       counts at the selected threshold (regularization-path view).}
#'     \item{cutoff}{List describing the selected threshold, the selection
#'       rule in effect (`"user_hard_threshold"`, `"fdp_plus_minimum"`, or
#'       `"none"`), the comparator, and the FDP+ denominator convention.}
#'     \item{fdp_curve}{List with the FDP+ curve, all minimizing thresholds,
#'       and tie/flat/boundary flags.}
#'     \item{null_saturation}{List with the null-saturation probe result.}
#'     \item{dependence_probe}{List with score-profile correlation summaries
#'       (assumption probe only).}
#'     \item{evidence_profile}{`data.frame` of matching evidence rows with
#'       `status` and `provenance`.}
#'     \item{warning_codes}{Character vector of machine-readable codes.}
#'     \item{assumption_notes}{Character vector of assumption statements.}
#'     \item{calibration_scope}{Character vector delimiting what the FDP+
#'       diagnostic does and does not establish.}
#'   }
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#' y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
#' fit <- stabl_fit(x, y, n_bootstraps = 20, lambda_grid = "auto",
#'                  artificial_type = "random_permutation")
#' d <- stabl_diagnose(fit)
#' d$warning_codes
#'
#' @seealso [compute_fdp_plus()], [plot_fdr_graph()], [get_support()]
#' @export
stabl_diagnose <- function(fit,
                           ...,
                           saturation_threshold = 0.05,
                           saturation_proportion = 0.5,
                           low_decoy_count = 5L) {
  .check_fitted_stabl(fit)
  dots <- list(...)
  if (length(dots) > 0L) {
    stop("Unused arguments: ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }
  saturation_threshold <- .validate_scalar_numeric(
    saturation_threshold, "saturation_threshold", min = 0, max = 1
  )
  saturation_proportion <- .validate_proportion(
    saturation_proportion, "saturation_proportion"
  )
  low_decoy_count <- .validate_scalar_integer_like(
    low_decoy_count, "low_decoy_count", min = 0
  )

  .diagnose_stabl_scores(
    stabl_scores            = fit$stabl_scores_,
    stabl_scores_artificial = fit$stabl_scores_artificial_,
    fdr_threshold_range     = fit$fdr_threshold_range,
    FDRs                    = fit$FDRs_,
    min_fdr                 = fit$min_fdr_,
    fdr_min_threshold       = fit$fdr_min_threshold_,
    hard_threshold          = fit$hard_threshold,
    artificial_type         = fit$artificial_type,
    artificial_provenance   = fit$artificial_provenance,
    artificial_proportion   = fit$artificial_proportion,
    effective_artificial_proportion = fit$effective_artificial_proportion,
    base_learner            = fit$base_learner,
    family                  = fit$family,
    n_bootstraps            = fit$n_bootstraps,
    feature_names           = fit$feature_names,
    saturation_threshold    = saturation_threshold,
    saturation_proportion   = saturation_proportion,
    low_decoy_count         = low_decoy_count
  )
}

# Internal workhorse operating on explicit components so that hand-worked
# score fixtures can be diagnosed without fitting a model. Not exported.
.diagnose_stabl_scores <- function(stabl_scores,
                                   stabl_scores_artificial = NULL,
                                   fdr_threshold_range = seq(0, 0.99, by = 0.01),
                                   FDRs = NULL,
                                   min_fdr = NULL,
                                   fdr_min_threshold = NULL,
                                   hard_threshold = NULL,
                                   artificial_type = NULL,
                                   artificial_provenance = NULL,
                                   artificial_proportion = NULL,
                                   effective_artificial_proportion = NULL,
                                   base_learner = NULL,
                                   family = NULL,
                                   n_bootstraps = NULL,
                                   feature_names = NULL,
                                   saturation_threshold = 0.05,
                                   saturation_proportion = 0.5,
                                   low_decoy_count = 5L) {
  has_artificial <- !is.null(stabl_scores_artificial) &&
    nrow(stabl_scores_artificial) > 0L

  max_real <- rowMaxs(stabl_scores)
  n_real_total <- length(max_real)
  max_art <- if (has_artificial) rowMaxs(stabl_scores_artificial) else NULL
  n_art_total <- length(max_art)

  warning_codes   <- character(0)
  assumption_notes <- character(0)
  calibration_scope <- character(0)

  # ---- 1. Score distributions + decoy provenance ---------------------------
  summarize_scores <- function(v) {
    list(
      n         = length(v),
      quantiles = stats::quantile(v, probs = c(0, 0.25, 0.5, 0.75, 1),
                                  names = FALSE),
      mean      = mean(v),
      prop_zero = mean(v == 0),
      prop_one  = mean(v == 1)
    )
  }
  score_summary <- list(real = summarize_scores(max_real))
  if (has_artificial) {
    score_summary$artificial <- summarize_scores(max_art)
  }

  # ---- 2. Threshold path (strict `>` counting, identical to FDP+) ----------
  n_real_path <- .count_strict_exceedances(max_real, fdr_threshold_range)
  n_art_path <- if (has_artificial) {
    .count_strict_exceedances(max_art, fdr_threshold_range)
  } else {
    rep(NA_real_, length(fdr_threshold_range))
  }
  fdp_path <- if (!is.null(FDRs)) {
    FDRs
  } else if (has_artificial && !is.null(effective_artificial_proportion)) {
    (n_art_path / effective_artificial_proportion + 1.0) /
      pmax(1.0, n_real_path)
  } else {
    rep(NA_real_, length(fdr_threshold_range))
  }
  threshold_path <- data.frame(
    threshold        = fdr_threshold_range,
    n_real           = as.numeric(n_real_path),
    n_artificial     = as.numeric(n_art_path),
    prop_real        = as.numeric(n_real_path) / n_real_total,
    prop_artificial  = if (has_artificial) {
      as.numeric(n_art_path) / n_art_total
    } else {
      rep(NA_real_, length(fdr_threshold_range))
    },
    fdp_plus         = as.numeric(fdp_path)
  )

  # ---- 3. Cutoff: rule in effect, comparator, denominator ------------------
  selection_rule <- if (!is.null(hard_threshold)) {
    "user_hard_threshold"
  } else if (!is.null(fdr_min_threshold)) {
    "fdp_plus_minimum"
  } else {
    "none"
  }
  selected_threshold <- switch(
    selection_rule,
    user_hard_threshold = hard_threshold,
    fdp_plus_minimum    = fdr_min_threshold,
    none                = NA_real_
  )
  if (selection_rule == "user_hard_threshold") {
    fdp_ref <- if (is.null(fdr_min_threshold)) "n/a" else format(fdr_min_threshold)
    assumption_notes <- c(
      assumption_notes,
      paste0(
        "A user-supplied hard_threshold (", format(hard_threshold),
        ") is in effect; the minimized-FDP+ threshold (", fdp_ref,
        ") is shown for reference only and did not drive selection."
      )
    )
    warning_codes <- c(warning_codes, "STABL_DIAG_USER_THRESHOLD")
  }
  if (selection_rule == "none") {
    warning_codes <- c(warning_codes, "STABL_DIAG_NO_THRESHOLD")
    assumption_notes <- c(
      assumption_notes,
      "No threshold is available (no artificial features and no ",
      "hard_threshold): the selected support is empty by contract."
    )
  }
  cutoff <- list(
    selected_threshold  = selected_threshold,
    selection_rule      = selection_rule,
    hard_threshold      = hard_threshold,
    fdr_min_threshold   = fdr_min_threshold,
    comparator          = "strictly greater (>)",
    denominator         = list(
      artificial_proportion = effective_artificial_proportion,
      formula = "(n_artificial / artificial_proportion + 1) / pmax(1, n_real)",
      plus_one_offset = TRUE,
      zero_real_guard   = "pmax(1, n_real)"
    )
  )

  # ---- 4. FDP+ curve features: ties, flat minima, boundary -----------------
  fdp_curve <- NULL
  if (!is.null(FDRs) && !is.null(min_fdr)) {
    at_min <- which(FDRs == min_fdr)
    minima_thresholds <- fdr_threshold_range[at_min]
    is_tied <- length(at_min) > 1L
    is_flat <- is_tied && any(diff(at_min) == 1L)
    first_min_idx <- at_min[[1L]]
    is_boundary <- first_min_idx == 1L ||
      first_min_idx == length(fdr_threshold_range)
    fdp_curve <- list(
      thresholds         = fdr_threshold_range,
      fdp_plus           = FDRs,
      min_fdr            = min_fdr,
      minima_thresholds  = minima_thresholds,
      n_minima           = length(at_min),
      is_tied            = is_tied,
      is_flat            = is_flat,
      is_boundary        = is_boundary,
      tie_policy         = "first minimizing threshold wins (matches Python argmin)"
    )
    if (is_tied) {
      warning_codes <- c(warning_codes, "STABL_DIAG_TIED_MINIMA")
    }
    if (is_flat) {
      warning_codes <- c(warning_codes, "STABL_DIAG_FLAT_MINIMUM")
    }
    if (is_boundary) {
      warning_codes <- c(warning_codes, "STABL_DIAG_BOUNDARY_MINIMUM")
      assumption_notes <- c(
        assumption_notes,
        "The first FDP+ minimum lies at the edge of fdr_threshold_range; ",
        "consider widening the range before trusting the minimum."
      )
    }
    # Low-decoy-count instability at the selected threshold.
    if (has_artificial && !is.null(fdr_min_threshold)) {
      idx <- which(fdr_threshold_range == fdr_min_threshold)
      if (length(idx) == 1L && n_art_path[[idx]] < low_decoy_count) {
        warning_codes <- c(warning_codes, "STABL_DIAG_LOW_DECOY_COUNT")
        assumption_notes <- c(
          assumption_notes,
          sprintf(
            paste0(
              "Only %d artificial feature(s) exceed the selected threshold; ",
              "the FDP+ minimum rests on very few decoys and is unstable ",
              "under decoy re-realization."
            ),
            as.integer(n_art_path[[idx]])
          )
        )
      }
    }
  }

  # ---- 5. Null-like saturation probe ---------------------------------------
  # Two readings, distinct purposes:
  #  * `near_zero_region` (data): thresholds <= `saturation_threshold` that
  #    select large proportions of BOTH reals and decoys — the literal
  #    RF/oblique-RF signature of spec A1 (their lambda path is a near-zero
  #    importance-threshold path). Reported for information; dense lasso
  #    score matrices legitimately saturate at t = 0, so this alone never
  #    warns.
  #  * `detected` (warning): the SELECTED threshold itself admits at least
  #    `saturation_proportion` of both reals and decoys — the panel the user
  #    receives is then decoy-dominated by the fit's own decoy model.
  saturation <- list(
    detected = FALSE,
    at_selected = list(prop_real = NA_real_, prop_artificial = NA_real_),
    near_zero_region = list(thresholds = numeric(0)),
    saturation_threshold = saturation_threshold,
    saturation_proportion = saturation_proportion
  )
  if (has_artificial) {
    near_zero_sat <- fdr_threshold_range <= saturation_threshold &
      threshold_path$prop_real >= saturation_proportion &
      threshold_path$prop_artificial >= saturation_proportion
    saturation$near_zero_region$thresholds <-
      fdr_threshold_range[near_zero_sat]

    if (!is.na(selected_threshold)) {
      # Computed directly at the selected threshold (strict >), so a
      # user-supplied hard_threshold off the FDP+ grid is still probed.
      prop_r <- sum(max_real > selected_threshold) / n_real_total
      prop_a <- sum(max_art > selected_threshold) / n_art_total
      saturation$at_selected$prop_real <- prop_r
      saturation$at_selected$prop_artificial <- prop_a
      if (prop_r >= saturation_proportion &&
          prop_a >= saturation_proportion) {
          saturation$detected <- TRUE
          warning_codes <- c(warning_codes, "STABL_DIAG_NULL_SATURATION")
          learner_context <- if (!is.null(base_learner)) {
            sprintf(" (learner: %s)", base_learner)
          } else {
            ""
          }
          assumption_notes <- c(
            assumption_notes,
            sprintf(
              paste0(
                "Null-like score saturation%s: the selected threshold %.3g ",
                "admits %.0f%% of real and %.0f%% of artificial features. ",
                "This is the expected signature of importance leakage ",
                "(e.g., random-forest-type learners on correlated nulls); ",
                "the selected panel is decoy-dominated by the fit's own ",
                "decoy model and the FDP+ minimum is not trustworthy."
              ),
              learner_context, selected_threshold,
              100 * prop_r, 100 * prop_a
            )
          )
      }
    }
  }

  # ---- 6. Dependence probe (assumption probe, NOT a validity test) ---------
  dependence_probe <- list(
    real = .score_profile_correlation(stabl_scores),
    artificial = if (has_artificial) {
      .score_profile_correlation(stabl_scores_artificial)
    } else {
      NULL
    },
    note = paste0(
      "Score-profile correlations are an assumption probe only: they do ",
      "not test or establish permutation exchangeability or knockoff ",
      "validity."
    )
  )
  if (!is.null(dependence_probe$real$mean_abs) &&
      dependence_probe$real$mean_abs > 0.5) {
    warning_codes <- c(warning_codes, "STABL_DIAG_HIGH_SCORE_DEPENDENCE")
    assumption_notes <- c(
      assumption_notes,
      sprintf(
        paste0(
          "Mean absolute score-profile correlation among real features is ",
          "%.2f: features co-select across the regularization path. ",
          "Interpret per-feature selection cautiously; this probe does not ",
          "by itself invalidate the decoy design."
        ),
        dependence_probe$real$mean_abs
      )
    )
  }

  # ---- 7. Evidence profile ---------------------------------------------------
  evidence_profile <- .match_diagnostic_evidence(
    base_learner = base_learner,
    artificial_type = artificial_type,
    family = family
  )
  if (nrow(evidence_profile) == 0L) {
    warning_codes <- c(warning_codes, "STABL_DIAG_EVIDENCE_NOT_EVALUATED")
    evidence_profile <- data.frame(
      learner = base_learner %||% "unknown",
      decoy = artificial_type %||% "none",
      family = family %||% "unknown",
      regime = "any",
      status = "not_evaluated",
      provenance = NA_character_,
      note = paste0(
        "No benchmark artifact covers this learner x decoy x outcome ",
        "combination; the FDP+ diagnostic is unvalidated here."
      ),
      stringsAsFactors = FALSE
    )
  } else if (all(evidence_profile$status == "not_evaluated")) {
    # Every matching regime is unevaluated: surface that absence of evidence.
    warning_codes <- c(warning_codes, "STABL_DIAG_EVIDENCE_NOT_EVALUATED")
  } else if (all(evidence_profile$status == "empirically_concerning")) {
    # Warn only when NO matching evidence row supports a calibrated or
    # conservative regime; mixed profiles are reported without a blanket
    # warning (no binary safe/unsafe score, per spec).
    warning_codes <- c(warning_codes, "STABL_DIAG_EVIDENCE_CONCERNING")
  }

  # ---- 8. Calibration scope --------------------------------------------------
  if (!has_artificial) {
    warning_codes <- c(warning_codes, "STABL_DIAG_NO_ARTIFICIALS")
    calibration_scope <- c(
      calibration_scope,
      paste0(
        "Fitted without artificial features: no FDP+ diagnostic exists for ",
        "this fit. Threshold selection is an uncalibrated hard-threshold ",
        "choice."
      )
    )
  } else {
    calibration_scope <- c(
      calibration_scope,
      paste0(
        "The FDP+ minimum is a stochastic-bound diagnostic, not a ",
        "calibrated estimate of realized FDP (F-001). On pure-null designs ",
        "the minimizing rule still selected 18-33 features with realized ",
        "FDP 0.92-0.93 while reporting 0.42-0.66 (audit Phase 6, F-P6-1): ",
        "an empty support is the only outcome a null dataset cannot ",
        "produce by this rule."
      ),
      paste0(
        "Calibration is regime-dependent in both directions (audit Phase 6, ",
        "F-P6-2): approximately calibrated for knockoff decoys on ",
        "independent signal, anti-conservative for permutation decoys, ",
        "conservative under block correlation. No single interpretation of ",
        "the FDP+ minimum is safe across regimes."
      )
    )
  }

  # ---- 9. Empty support --------------------------------------------------------
  support_empty <- if (is.na(selected_threshold)) {
    TRUE
  } else {
    !any(max_real > selected_threshold)
  }
  if (support_empty) {
    warning_codes <- c(warning_codes, "STABL_DIAG_EMPTY_SUPPORT")
  }

  structure(
    list(
      score_summary     = score_summary,
      decoy_provenance  = artificial_provenance,
      threshold_path    = threshold_path,
      per_lambda        = .per_lambda_counts(
        stabl_scores, stabl_scores_artificial, selected_threshold
      ),
      cutoff            = cutoff,
      fdp_curve         = fdp_curve,
      null_saturation   = saturation,
      dependence_probe  = dependence_probe,
      evidence_profile  = evidence_profile,
      warning_codes     = unique(warning_codes),
      assumption_notes  = assumption_notes,
      calibration_scope = calibration_scope,
      support_empty     = support_empty,
      context = list(
        base_learner = base_learner,
        family = family,
        artificial_type = artificial_type,
        n_features = n_real_total,
        n_artificial = n_art_total,
        n_lambdas = ncol(stabl_scores),
        n_bootstraps = n_bootstraps
      )
    ),
    class = "stabl_diagnostics"
  )
}

# Mean/max absolute off-diagonal correlation of per-feature score profiles
# across the regularization path. Returns NULLs when not computable.
.score_profile_correlation <- function(scores) {
  if (nrow(scores) < 2L || ncol(scores) < 3L) {
    return(list(mean_abs = NULL, max_abs = NULL, n_pairs = 0L))
  }
  # Zero-variance profiles (constant scores across the path) produce NA
  # correlations with a "standard deviation is zero" warning; those pairs are
  # filtered below, so muffle exactly that warning and no other.
  profile_cor <- withCallingHandlers(
    stats::cor(t(scores)),
    warning = function(w) {
      if (grepl("standard deviation is zero", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  off <- profile_cor[upper.tri(profile_cor)]
  off <- off[is.finite(off)]
  if (length(off) == 0L) {
    return(list(mean_abs = NULL, max_abs = NULL, n_pairs = 0L))
  }
  list(
    mean_abs = mean(abs(off)),
    max_abs  = max(abs(off)),
    n_pairs  = length(off)
  )
}

# Per-lambda selected counts at the chosen threshold (regularization path).
.per_lambda_counts <- function(stabl_scores, stabl_scores_artificial,
                               selected_threshold) {
  if (is.na(selected_threshold)) {
    return(NULL)
  }
  out <- data.frame(
    lambda_index = seq_len(ncol(stabl_scores)),
    n_real = as.numeric(colSums(stabl_scores > selected_threshold)),
    stringsAsFactors = FALSE
  )
  if (!is.null(stabl_scores_artificial) && nrow(stabl_scores_artificial) > 0L) {
    out$n_artificial <- as.numeric(
      colSums(stabl_scores_artificial > selected_threshold)
    )
  }
  out
}

# Match the fit's learner x decoy x family against the evidence table.
# "any" rows match every value; a fit with unknown metadata matches only
# "any" rows for that dimension.
.match_diagnostic_evidence <- function(base_learner, artificial_type, family) {
  tbl <- .stabl_diagnostic_evidence
  match_dim <- function(column, value) {
    column == "any" | (!is.null(value) & column == value)
  }
  keep <- match_dim(tbl$learner, base_learner) &
    match_dim(tbl$decoy, artificial_type) &
    match_dim(tbl$family, family)
  tbl[keep, , drop = FALSE]
}

#' Print STABL Selection Diagnostics
#'
#' @param x A `stabl_diagnostics` object from [stabl_diagnose()].
#' @param ... Ignored; present for S3 compatibility.
#'
#' @return The input `x`, invisibly.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#' y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
#' fit <- stabl_fit(x, y, n_bootstraps = 20, lambda_grid = "auto",
#'                  artificial_type = "random_permutation")
#' print(stabl_diagnose(fit))
#' @export
print.stabl_diagnostics <- function(x, ...) {
  ctx <- x$context
  cat("STABL selection diagnostics\n")
  cat(sprintf(
    "  learner: %s | family: %s | decoys: %s\n",
    ctx$base_learner %||% "unknown", ctx$family %||% "unknown",
    ctx$artificial_type %||% "none"
  ))
  cat(sprintf(
    "  features: %d real / %d artificial | lambdas: %d | bootstraps: %s\n",
    ctx$n_features, ctx$n_artificial %||% 0L, ctx$n_lambdas,
    if (is.null(ctx$n_bootstraps)) "unknown" else as.character(ctx$n_bootstraps)
  ))
  cat(sprintf(
    "  selection rule: %s | threshold: %s | comparator: %s\n",
    x$cutoff$selection_rule,
    if (is.na(x$cutoff$selected_threshold)) {
      "none"
    } else {
      format(x$cutoff$selected_threshold)
    },
    x$cutoff$comparator
  ))
  if (!is.null(x$fdp_curve)) {
    cat(sprintf(
      "  FDP+ minimum: %s at threshold %s (%d minimizing threshold%s)\n",
      format(signif(x$fdp_curve$min_fdr, 3L)),
      format(x$fdp_curve$minima_thresholds[[1L]]),
      x$fdp_curve$n_minima, if (x$fdp_curve$n_minima > 1L) "s" else ""
    ))
  }
  if (isTRUE(x$null_saturation$detected)) {
    cat(sprintf(
      paste0(
        "  NULL SATURATION: the selected threshold %.3g admits %.0f%% of ",
        "real and %.0f%% of artificial features\n"
      ),
      x$cutoff$selected_threshold,
      100 * x$null_saturation$at_selected$prop_real,
      100 * x$null_saturation$at_selected$prop_artificial
    ))
  }
  cat("  evidence:\n")
  for (i in seq_len(nrow(x$evidence_profile))) {
    row <- x$evidence_profile[i, ]
    cat(sprintf("    - [%s] %s x %s x %s (%s)\n",
                row$status, row$learner, row$decoy, row$family, row$regime))
  }
  if (length(x$warning_codes) > 0L) {
    cat("  warning codes:\n")
    for (code in x$warning_codes) cat("    - ", code, "\n", sep = "")
  }
  if (length(x$assumption_notes) > 0L) {
    cat("  assumption notes:\n")
    for (note in x$assumption_notes) {
      cat("    - ", paste0(note, collapse = ""), "\n", sep = "")
    }
  }
  if (length(x$calibration_scope) > 0L) {
    cat("  calibration scope:\n")
    for (scope in x$calibration_scope) {
      cat("    - ", paste0(scope, collapse = ""), "\n", sep = "")
    }
  }
  invisible(x)
}

#' Plot STABL selection diagnostics
#'
#' @param x A `stabl_diagnostics` object from [stabl_diagnose()].
#' @param which Character vector of panels to draw: `"fdp_curve"` (FDP+ along
#'   the threshold range, with the selected cut-off and all minimizing
#'   thresholds marked), `"threshold_path"` (real vs artificial selection
#'   counts along the path), and `"score_distributions"` (ECDFs of real and
#'   artificial per-feature maximum scores).
#' @param ... Passed to plotting helpers; currently ignored.
#'
#' @return Invisibly, a named list of the drawn `ggplot` objects.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- matrix(rnorm(60 * 8), 60, 8,
#'               dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#'   y <- setNames(x[, 1] * 2 + rnorm(60), rownames(x))
#'   fit <- stabl_fit(x, y, n_bootstraps = 20, lambda_grid = "auto",
#'                    artificial_type = "random_permutation")
#'   plot(stabl_diagnose(fit), which = "fdp_curve")
#' }
#' @export
plot.stabl_diagnostics <- function(
    x,
    which = c("fdp_curve", "threshold_path", "score_distributions"),
    ...) {
  .require_ggplot2()
  which <- match.arg(which, several.ok = TRUE)
  plots <- list()

  if ("fdp_curve" %in% which && !is.null(x$fdp_curve)) {
    df <- data.frame(
      threshold = x$fdp_curve$thresholds,
      fdp_plus  = x$fdp_curve$fdp_plus
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = fdp_plus)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(
        xintercept = x$fdp_curve$minima_thresholds,
        linetype = if (x$fdp_curve$n_minima > 1L) "dotted" else "dashed"
      ) +
      ggplot2::labs(
        title = "FDP+ diagnostic",
        x = "stability threshold",
        y = "FDP+"
      )
    if (!is.na(x$cutoff$selected_threshold) &&
        x$cutoff$selection_rule == "user_hard_threshold") {
      p <- p + ggplot2::geom_vline(
        xintercept = x$cutoff$selected_threshold,
        linetype = "solid", colour = "grey40"
      )
    }
    plots$fdp_curve <- p
  }

  if ("threshold_path" %in% which) {
    df <- x$threshold_path
    df_long <- data.frame(
      threshold = rep(df$threshold, 2L),
      count     = c(df$n_real, df$n_artificial),
      kind      = rep(c("real", "artificial"), each = nrow(df))
    )
    df_long <- df_long[is.finite(df_long$count), , drop = FALSE]
    plots$threshold_path <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = threshold, y = count, colour = kind)
    ) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = "Selection counts along the threshold path",
        x = "stability threshold",
        y = "features above threshold (strict >)",
        colour = NULL
      )
  }

  if ("score_distributions" %in% which) {
    # Survival view of per-feature maximum scores, real vs artificial:
    # prop_* is already P(score > t) along the threshold grid.
    path <- x$threshold_path
    df_surv <- data.frame(
      threshold = rep(path$threshold, 2L),
      proportion = c(path$prop_real, path$prop_artificial),
      kind = rep(c("real", "artificial"), each = nrow(path))
    )
    df_surv <- df_surv[is.finite(df_surv$proportion), , drop = FALSE]
    plots$score_distributions <- ggplot2::ggplot(
      df_surv,
      ggplot2::aes(x = threshold, y = proportion, colour = kind)
    ) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = "Score survival (proportion of features with score > t)",
        x = "stability score t",
        y = "proportion",
        colour = NULL
      )
  }

  for (p in plots) print(p)
  invisible(plots)
}
