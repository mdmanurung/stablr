# interpretation.R — descriptive biomarker interpretation for fitted STABL
# objects.
#
# Evidence contract (roadmap section 3): every quantity returned here is
# *post-selection descriptive* or *apparent* unless the caller supplies
# genuinely held-out assessment data.  Nothing here re-runs STABL, adjusts
# for selection, or produces inferential p-values.  The four evidence
# dimensions — stability, conditional effect, redundancy group, and
# permutation importance — are reported separately and are never merged into
# one score.

# ── Public entry point ───────────────────────────────────────────────────────

#' Interpret a fitted STABL selection
#'
#' Builds a canonical feature-evidence table for a fitted [stabl_fit()] or
#' [stabl()] object **without re-running STABL**.  For each feature the table
#' records the stability score, the selection flag, a conditional signed
#' effect from a separately fitted downstream model, a reproducible
#' correlation/redundancy group, and a seeded permutation importance.
#'
#' @details
#' **Evidence labels.**  All effects and importances are computed on the data
#' supplied as `x`/`y`.  When that data is the fitting data, every row is
#' labelled `post_selection_descriptive` (effects) and `apparent`
#' (importance): the same samples informed both selection and the downstream
#' fit, so the quantities are optimistic and carry no inferential meaning.
#' Supplying genuinely held-out `x`/`y` relabels the importance scope, but
#' the caller is responsible for that data being untouched by fitting — this
#' function cannot verify it.
#'
#' **What each column is.**  `stability` is the max-over-lambda selection
#' frequency from [get_importances()]; it is not a coefficient magnitude or a
#' prediction importance.  `coefficient` / `standardized_effect` /
#' `effect_direction` come from a downstream [stats::glm()] fit on the
#' selected features only and are conditional, descriptive associations.
#' `correlation_group` is a reproducible high-|correlation| cluster; it is
#  not biological substitutability.  `permutation_importance` is the mean
#' increase in loss when the feature is permuted (positive = helpful);
#' negative values are retained.
#'
#' @param object A fitted `"stabl_fit"` or `"stabl"` object.
#' @param x Feature matrix (samples × features) with row and column names.
#'   Must contain every selected feature; extra columns are ignored.
#' @param y Outcome aligned to the rows of `x` (named vector or factor).
#' @param family `"gaussian"` or `"binomial"`; the downstream effect model.
#' @param assessment Logical; set `TRUE` only when `x`/`y` are genuinely
#'   held out from the fitting procedure.  Changes the `evaluation_scope`
#'   label, nothing else.
#' @param cor_method `"pearson"` or `"spearman"`.
#' @param cor_threshold Numeric in `[0, 1]`; |correlation| at which two
#'   features join one redundancy group.
#' @param n_repeats Number of permutation repeats for importance.
#' @param loss Loss for permutation importance: `"squared_error"` (gaussian),
#'   `"brier"` or `"log_loss"` (binomial).
#' @param random_state Seed for permutation importance.
#' @param event_class For binomial, the level treated as the event; defaults
#'   to the second level of the outcome factor.
#' @param ... Reserved for future use.
#'
#' @return An object of class `"stabl_interpretation"`: a list with
#'   `evidence` (the canonical data.frame), `effects_model` (the downstream
#'   fit), `correlation` (method, threshold, groups), `importance`
#'   (per-repeat matrix and settings), `event_class`, `family`, and
#'   `provenance` (call, seed, scopes).  Missing or unestimable entries are
#'   explicit `NA` with a reason column.
#'
#' @examples
#' set.seed(1L)
#' x <- matrix(rnorm(60 * 8), 60, 8,
#'             dimnames = list(paste0("s", 1:60), paste0("f", 1:8)))
#' y <- setNames(0.8 * x[, "f1"] + rnorm(60), rownames(x))
#' fit <- stabl_fit(x, y,
#'                  lambda_grid = data.frame(lambda = c(0.3, 0.1, 0.05)),
#'                  n_bootstraps = 10L, artificial_type = "random_permutation",
#'                  random_state = 1L)
#' insight <- stabl_interpret(fit, x, y, family = "gaussian", n_repeats = 3L)
#' summary(insight)
#' @export
stabl_interpret <- function(object, x, y, family = c("gaussian", "binomial"),
                            assessment = FALSE,
                            cor_method = c("pearson", "spearman"),
                            cor_threshold = 0.7,
                            n_repeats = 20L,
                            loss = NULL,
                            random_state = NULL,
                            event_class = NULL,
                            ...) {
  .check_fitted_stabl(object)
  family <- match.arg(family)
  cor_method <- match.arg(cor_method)
  cor_threshold <- .validate_proportion(cor_threshold, "cor_threshold", allow_one = TRUE)
  n_repeats <- .validate_scalar_integer_like(n_repeats, "n_repeats", min = 1L)
  if (!is.null(random_state)) {
    random_state <- .validate_scalar_integer_like(random_state, "random_state")
  }

  # ── Input validation: strict IDs, names, level order ─────────────────────
  x <- as.matrix(x)
  if (is.null(colnames(x)) || is.null(rownames(x))) {
    stop("`x` must have both row names (sample IDs) and column names ",
         "(feature names).", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("`x` must be numeric after as.matrix().", call. = FALSE)
  }
  feat_names <- object$feature_names
  missing_feats <- setdiff(feat_names, colnames(x))
  if (length(missing_feats)) {
    stop(sprintf("`x` is missing %d fitted feature(s): %s",
                 length(missing_feats),
                 paste(head(missing_feats, 5L), collapse = ", ")),
         call. = FALSE)
  }
  x <- x[, feat_names, drop = FALSE]

  y <- .interpret_align_outcome(y, rownames(x), family, event_class)
  event_class <- attr(y, "event_class")
  y_vec <- attr(y, "y_numeric")

  # ── Selection and stability (reuse accessors; semantics unchanged) ───────
  support <- get_support(object)
  importances <- get_importances(object)
  threshold <- .resolve_stabl_threshold(object, NULL, on_missing = "na")
  selected <- names(support)[support]

  scope <- if (isTRUE(assessment)) "held_out_assessment" else "fitting_data"
  evidence_type <- if (isTRUE(assessment)) "cross_validated" else
    "post_selection_descriptive"

  # ── Downstream effect model on selected features ─────────────────────────
  effects <- .interpret_effects(x, y_vec, selected, family, event_class)

  # ── Redundancy groups on the supplied data ───────────────────────────────
  corr <- .interpret_correlation_groups(x, selected, cor_method, cor_threshold)

  # ── Permutation importance ───────────────────────────────────────────────
  loss <- loss %||% switch(family, gaussian = "squared_error", binomial = "brier")
  importance <- .interpret_permutation_importance(
    effects$model, x, y_vec, selected, family, loss, n_repeats,
    random_state, event_class
  )

  # ── Canonical evidence table ─────────────────────────────────────────────
  evidence <- data.frame(
    feature = feat_names,
    selected = unname(support[feat_names]),
    stability = unname(importances[feat_names]),
    selection_threshold = threshold,
    coefficient = unname(effects$coefficient[feat_names]),
    standardized_effect = unname(effects$standardized_effect[feat_names]),
    effect_direction = unname(effects$effect_direction[feat_names]),
    effect_measure = effects$effect_measure,
    correlation_group = unname(corr$group[feat_names]),
    permutation_importance = unname(importance$mean[feat_names]),
    importance_sd = unname(importance$sd[feat_names]),
    importance_n_repeats = ifelse(feat_names %in% selected, n_repeats, 0L),
    evidence_type = evidence_type,
    evaluation_scope = scope,
    n_observations = nrow(x),
    model_id = "downstream_refit",
    stringsAsFactors = FALSE
  )
  # Unselected features have no downstream effect or importance: explicit NA.
  evidence$coefficient[!evidence$selected] <- NA_real_
  evidence$standardized_effect[!evidence$selected] <- NA_real_
  evidence$effect_direction[!evidence$selected] <- NA_character_
  evidence$permutation_importance[!evidence$selected] <- NA_real_
  evidence$importance_sd[!evidence$selected] <- NA_real_
  evidence$na_reason <- ifelse(
    evidence$selected, NA_character_, "not_selected"
  )
  # Features selected but with an unestimable effect (e.g. separation).
  unestimable <- evidence$selected & is.na(evidence$coefficient)
  evidence$na_reason[unestimable] <- "effect_not_estimable"

  out <- list(
    evidence = evidence,
    effects_model = effects$model,
    correlation = list(
      method = cor_method,
      threshold = cor_threshold,
      groups = corr$groups,
      scope = scope
    ),
    importance = list(
      loss = loss,
      repeats = importance$repeats,
      n_repeats = n_repeats,
      scope = scope
    ),
    event_class = event_class,
    family = family,
    provenance = list(
      n_observations = nrow(x),
      n_features = length(feat_names),
      n_selected = length(selected),
      random_state = random_state,
      evaluation_scope = scope,
      evidence_type = evidence_type
    )
  )
  class(out) <- "stabl_interpretation"
  out
}

# ── Outcome alignment ────────────────────────────────────────────────────────

.interpret_align_outcome <- function(y, sample_ids, family, event_class) {
  if (is.null(names(y))) {
    if (length(y) != length(sample_ids)) {
      stop("`y` must be named or have the same length as nrow(x).",
           call. = FALSE)
    }
    names(y) <- sample_ids
  }
  if (!identical(sort(names(y)), sort(sample_ids))) {
    stop("`y` names do not match the row names of `x`.", call. = FALSE)
  }
  y <- y[sample_ids]

  if (identical(family, "binomial")) {
    f <- if (is.factor(y)) y else factor(y)
    if (nlevels(f) != 2L) {
      stop("Binomial interpretation requires exactly two outcome levels.",
           call. = FALSE)
    }
    event_class <- event_class %||% levels(f)[2L]
    if (!event_class %in% levels(f)) {
      stop("`event_class` must be one of the outcome levels.", call. = FALSE)
    }
    y_numeric <- as.integer(f == event_class)
    attr(y, "event_class") <- event_class
    attr(y, "y_numeric") <- setNames(y_numeric, sample_ids)
  } else {
    if (!is.numeric(y)) {
      stop("Gaussian interpretation requires a numeric outcome.", call. = FALSE)
    }
    attr(y, "event_class") <- NA_character_
    attr(y, "y_numeric") <- setNames(as.numeric(y), sample_ids)
  }
  y
}

# ── Downstream effect model ──────────────────────────────────────────────────

.interpret_effects <- function(x, y_vec, selected, family, event_class) {
  empty <- list(
    model = NULL,
    coefficient = setNames(rep(NA_real_, ncol(x)), colnames(x)),
    standardized_effect = setNames(rep(NA_real_, ncol(x)), colnames(x)),
    effect_direction = setNames(rep(NA_character_, ncol(x)), colnames(x)),
    effect_measure = switch(family, gaussian = "mean_difference",
                            binomial = "log_odds_ratio")
  )
  if (length(selected) == 0L) return(empty)

  df <- as.data.frame(x[, selected, drop = FALSE])
  df$.y <- unname(y_vec)
  fam <- if (identical(family, "binomial")) {
    stats::binomial(link = "logit")
  } else {
    stats::gaussian()
  }
  model <- tryCatch(
    suppressWarnings(stats::glm(.y ~ ., data = df, family = fam)),
    error = function(e) e
  )
  if (inherits(model, "error")) return(empty)

  coefs <- stats::coef(model)[-1L]  # drop intercept
  # Standardized effect: coefficient * SD(feature) (per-SD change).
  sds <- apply(x[, selected, drop = FALSE], 2L, stats::sd, na.rm = TRUE)
  std <- coefs * sds[names(coefs)]
  direction <- ifelse(is.na(coefs), NA_character_,
                      ifelse(coefs > 0, "positive",
                             ifelse(coefs < 0, "negative", "zero")))

  out <- empty
  out$model <- model
  out$coefficient[names(coefs)] <- unname(coefs)
  out$standardized_effect[names(coefs)] <- unname(std)
  out$effect_direction[names(coefs)] <- unname(direction)
  out
}

# ── Redundancy groups ────────────────────────────────────────────────────────

.interpret_correlation_groups <- function(x, selected, cor_method, cor_threshold) {
  feats <- colnames(x)
  group <- setNames(rep(NA_integer_, length(feats)), feats)
  groups <- list()
  if (length(selected) < 2L) {
    if (length(selected) == 1L) {
      group[selected] <- 1L
      groups <- list(`1` = selected)
    }
    return(list(group = group, groups = groups))
  }
  xs <- x[, selected, drop = FALSE]
  cmat <- suppressWarnings(stats::cor(xs, use = "pairwise.complete.obs",
                                      method = cor_method))
  cmat[is.na(cmat)] <- 0
  # Single-linkage clustering at the |correlation| threshold: reproducible and
  # deterministic for a fixed input.
  adj <- abs(cmat) >= cor_threshold
  diag(adj) <- FALSE
  n <- length(selected)
  comp <- seq_len(n)
  repeat {
    merged <- FALSE
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        if (adj[i, j] && comp[i] != comp[j]) {
          old <- comp[j]
          comp[comp == old] <- comp[i]
          merged <- TRUE
        }
      }
    }
    if (!merged) break
  }
  uniq <- sort(unique(comp))
  gid <- match(comp, uniq)
  group[selected] <- gid
  groups <- split(selected, gid)
  names(groups) <- as.character(sort(unique(gid)))
  list(group = group, groups = groups)
}

# ── Permutation importance ───────────────────────────────────────────────────

.interpret_permutation_importance <- function(model, x, y_vec, selected, family,
                                              loss, n_repeats, random_state,
                                              event_class) {
  feats <- colnames(x)
  out <- list(
    mean = setNames(rep(NA_real_, length(feats)), feats),
    sd = setNames(rep(NA_real_, length(feats)), feats),
    repeats = NULL
  )
  if (is.null(model) || length(selected) == 0L) return(out)

  xs <- x[, selected, drop = FALSE]
  pred_fun <- function(mat) {
    nd <- as.data.frame(mat)
    colnames(nd) <- selected
    p <- stats::predict(model, newdata = nd, type = "response")
    unname(p)
  }
  loss_fun <- .interpret_loss(loss, family)

  base_pred <- pred_fun(xs)
  base_loss <- loss_fun(y_vec, base_pred)

  if (!is.null(random_state)) set.seed(random_state)
  rep_mat <- matrix(NA_real_, nrow = n_repeats, ncol = length(selected),
                    dimnames = list(NULL, selected))
  for (r in seq_len(n_repeats)) {
    for (j in selected) {
      xp <- xs
      xp[, j] <- sample(xp[, j])  # one shared permutation per column
      rep_mat[r, j] <- loss_fun(y_vec, pred_fun(xp)) - base_loss
    }
  }
  out$mean[selected] <- colMeans(rep_mat)
  out$sd[selected] <- apply(rep_mat, 2L, stats::sd)
  out$repeats <- rep_mat
  out
}

.interpret_loss <- function(loss, family) {
  if (identical(family, "gaussian") || identical(loss, "squared_error")) {
    return(function(y, p) mean((y - p)^2))
  }
  if (identical(loss, "brier")) {
    return(function(y, p) mean((y - p)^2))
  }
  if (identical(loss, "log_loss")) {
    return(function(y, p) {
      p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
      -mean(y * log(p) + (1 - y) * log(1 - p))
    })
  }
  stop("Unsupported loss: ", loss, call. = FALSE)
}

# ── S3 methods ───────────────────────────────────────────────────────────────

#' @export
print.stabl_interpretation <- function(x, ...) {
  ev <- x$evidence
  cat("<stabl_interpretation>\n")
  cat(sprintf("  family: %s | scope: %s | evidence: %s\n",
              x$family, x$provenance$evaluation_scope,
              x$provenance$evidence_type))
  cat(sprintf("  %d feature(s), %d selected, %d observations\n",
              nrow(ev), sum(ev$selected), x$provenance$n_observations))
  invisible(x)
}

#' @export
summary.stabl_interpretation <- function(object, ...) {
  ev <- object$evidence[object$evidence$selected, , drop = FALSE]
  ev <- ev[order(-ev$stability), , drop = FALSE]
  rownames(ev) <- NULL
  ev
}

#' @export
plot.stabl_interpretation <- function(x, ...) {
  .require_ggplot2()
  ev <- x$evidence[x$evidence$selected, , drop = FALSE]
  if (nrow(ev) == 0L) {
    stop("No selected features to plot.", call. = FALSE)
  }
  ev$feature <- factor(ev$feature, levels = ev$feature[order(ev$stability)])
  ggplot2::ggplot(ev, ggplot2::aes(x = ev$stability, y = ev$feature)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Stability score", y = NULL,
                  title = "Selected-feature stability") +
    ggplot2::theme_minimal(base_size = 11)
}
