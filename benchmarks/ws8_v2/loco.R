# WS-8 v2 M2 — Leakage-proof leave-one-cohort-out (LOCO) engine.
#
# Defines the 8 fixed outer folds and provides fold-safe feature transforms,
# training-only selection/sign assignment, cohort-respecting inner ridge CV,
# and immutable held-out predictions.  Every data-adaptive step is fit on the
# seven training cohorts only.  A provenance guard fails if a held-out donor ID
# appears in any upstream fit.
#
# This is an engine, not a benchmark: run.R calls it per method.  Pure
# functions plus a fold-definition writer.  Sourced by run.R / null_probe.R.
#
# Application case study — NOT calibration evidence.

suppressPackageStartupMessages(library(glmnet))

# ---- Fold definition ----------------------------------------------------------
# Eight fixed outer folds, one per cohort, in alphabetical order.  Returns a
# data.frame of donor / cohort / fold / role for the manifest.
make_loco_folds <- function(donor, cohort) {
  stopifnot(length(donor) == length(cohort))
  levels <- sort(unique(cohort))
  do.call(rbind, lapply(levels, function(h) {
    data.frame(fold = h, donor = donor[cohort == h], cohort = h,
               role = "test", stringsAsFactors = FALSE)
  }))
}

# ---- Fold-safe feature transform ----------------------------------------------
# The composition transform (add-half smoothing + logit) is per-donor and
# deterministic, so it carries no leakage.  Any LEARNED transform (imputation,
# scaling, filtering, ancestry collapsing, knockoff covariance) must be fit on
# training donors only.  This helper fits a per-fold transform recipe on the
# training matrix and applies it to both train and test.
#
# Here the only learned step we allow is dropping features that are constant
# within any TRAINING cohort (zero variance -> singular knockoff covariance).
# The drop set is computed on training data only, then applied to test.
fit_fold_transform <- function(X_train, cohort_train) {
  stopifnot(is.matrix(X_train))
  zv <- unique(unlist(lapply(sort(unique(cohort_train)), function(e) {
    Xe <- X_train[cohort_train == e, , drop = FALSE]
    v <- apply(Xe, 2, function(cc) stats::var(cc, na.rm = TRUE) == 0)
    names(v)[v]
  })))
  list(drop = zv, keep = setdiff(colnames(X_train), zv),
       n_train = nrow(X_train), train_cohorts = sort(unique(cohort_train)))
}
apply_fold_transform <- function(recipe, X) {
  X[, recipe$keep, drop = FALSE]
}

# ---- Training-only direction signs --------------------------------------------
# Sign of the per-feature association with age, estimated on TRAINING donors
# only: sign of the cohort-balanced (unweighted) mean per-cohort Spearman rho
# across training cohorts.  Indeterminate (|mean rho| < tol) -> 0.
training_signs <- function(X_train, y_train, cohort_train, tol = 1e-6) {
  feats <- colnames(X_train)
  rhos <- sapply(feats, function(f) {
    r <- sapply(sort(unique(cohort_train)), function(e) {
      idx <- cohort_train == e
      suppressWarnings(stats::cor(X_train[idx, f], y_train[idx],
                                  method = "spearman"))
    })
    mean(r, na.rm = TRUE)  # cohort-balanced mean rho
  })
  s <- ifelse(rhos > tol, 1L, ifelse(rhos < -tol, -1L, 0L))
  setNames(as.integer(s), feats)
}

# ---- Cohort-respecting inner ridge CV -----------------------------------------
# Tune ridge lambda by leave-one-TRAINING-cohort-out inner CV (each training
# cohort held out in turn).  Falls back to random donor folds only if a training
# cohort is too small to hold out (< 5 donors).  Returns the fitted cv.glmnet
# and the chosen lambda.
train_ridge_cohort_cv <- function(X_train, y_train, cohort_train) {
  cohs <- sort(unique(cohort_train))
  if (all(table(cohort_train) >= 5)) {
    # foldid assigning each training cohort its own fold -> LOCO inner CV
    foldid <- match(cohort_train, cohs)
    cvfit <- cv.glmnet(X_train, y_train, alpha = 0, foldid = foldid)
  } else {
    cvfit <- cv.glmnet(X_train, y_train, alpha = 0, nfolds = 5)
  }
  cvfit
}

# ---- Leakage guard ------------------------------------------------------------
# Fail if any test donor ID appears in the training set used for a fit.
assert_no_leakage <- function(train_ids, test_ids, context = "") {
  bad <- intersect(train_ids, test_ids)
  if (length(bad) > 0L) {
    stop(sprintf("LEAKAGE in %s: %d test donor IDs in training set (e.g. %s)",
                 context, length(bad), paste(head(bad, 3), collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

# ---- One outer fold -----------------------------------------------------------
# select_fn: function(X_train, y_train, cohort_train, covar = NULL) -> list
#   (ranked = character vector of feature names in rank order, native =
#   character vector of the native selected panel, scores = named numeric
#   training scores, optional detail = features x cohorts audit matrix).
# K: top-K budget for the matched-budget panel.
# covar: optional data.frame of training-only adjustment covariates (sex,
#   ethnicity) with rows aligned to `donor`; only selectors that declare a
#   covariate-adjusted rule (M1) use it.  Passed through untouched by RNG.
# Returns per-fold replication + prediction inputs (no evaluation here).
run_one_fold <- function(h, X, y, cohort, donor, select_fn, K = 10L,
                         covar = NULL) {
  te <- cohort == h; tr <- !te
  train_ids <- donor[tr]; test_ids <- donor[te]
  assert_no_leakage(train_ids, test_ids, context = paste0("fold ", h))

  # Fold-safe transform (training-fit feature filtering).
  recipe <- fit_fold_transform(X[tr, , drop = FALSE], cohort[tr])
  Xtr <- apply_fold_transform(recipe, X[tr, , drop = FALSE])
  Xte <- apply_fold_transform(recipe, X[te, , drop = FALSE])
  # stabl_fit requires outcomes/groups named by sample id matching rownames(x).
  ytr <- setNames(y[tr], donor[tr]); yte <- setNames(y[te], donor[te])
  ctr <- setNames(as.character(cohort[tr]), donor[tr])

  # Training-only selection/ranking (cohort vector named by sample id).
  if (!is.null(covar)) {
    stopifnot(identical(rownames(covar), donor))
    sel <- select_fn(Xtr, ytr, ctr, covar = covar[tr, , drop = FALSE])
  } else {
    sel <- select_fn(Xtr, ytr, ctr)
  }
  stopifnot(all(sel$ranked %in% colnames(Xtr)))
  panel_topk <- head(sel$ranked, K)
  panel_native <- intersect(sel$native, colnames(Xtr))

  # Training-only direction signs (for the replication endpoint).
  signs <- training_signs(Xtr, ytr, ctr)

  list(
    fold = h, test_ids = test_ids, train_ids = train_ids,
    recipe = recipe, panel_topk = panel_topk, panel_native = panel_native,
    signs = signs, scores = sel$scores, sel_full = sel,
    Xtr = Xtr, Xte = Xte, ytr = ytr, yte = yte,
    cohort_tr = ctr, cohort_te = setNames(as.character(cohort[te]), donor[te])
  )
}
