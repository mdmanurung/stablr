#!/usr/bin/env Rscript
# ============================================================================
# WS-4 — Cooperative vs early vs late fusion benchmark (predeclared)
# ============================================================================
# 2 views; complementary arm (both views informative) vs redundant arm
# (view 2 pure noise); binomial + gaussian; n_train 150, n_valid 50;
# 3 fusion modes; 5 seeds; 100 bootstraps.
#
# Predeclared evaluation rule (MATURITY_PREDESIGN.md, WS-4):
#   - Held-out performance on the validation split: BER (binomial) / MSE
#     (gaussian).
#       * late fusion:  native `late_fusion$valid_predictions`
#         (probability scale for binomial, outcome scale for gaussian).
#       * cooperative:  native `cooperative_fusion$valid_predictions`
#         (response scale for binomial under the default deviance
#         type_measure; link = outcome scale for gaussian).
#       * early fusion: `stabl_fit()` stores no reusable model, so held-out
#         predictions come from a downstream `stats::glm` refit on the
#         concatenated selected features (`early_fusion$selected_train`) —
#         the same downstream convention the package's late-fusion path
#         uses per omic. Recorded as `early` in the perf_source column.
#   - Selection-set Jaccard across seeds, computed on the mode-appropriate
#     selection set (early: concatenated early-fusion selection; late:
#     per-omic stabl selections; cooperative: multiview per-view selection).
#   - View-wise selection shares: fraction of the mode's selected features
#     carrying the `v2_` prefix.
#
# Run: Rscript benchmarks/ws4_fusion_benchmark.R   (~0.5-3 h)

.libPaths(c("/workspace/Rlib-shared", .libPaths()))
if (file.exists("DESCRIPTION")) pkgload::load_all(".", compile = TRUE, quiet = TRUE) else library(stablr)

seeds <- 1:5
n_bootstraps <- 100L
out_dir <- "benchmarks"

make_data <- function(arm, family, seed) {
  set.seed(seed)
  n_tr <- 150; n_va <- 50; p <- 25
  mk_view <- function(seed_shift) {
    Sigma <- matrix(0.2, p, p); diag(Sigma) <- 1
    L <- t(chol(Sigma))
    X <- matrix(rnorm((n_tr + n_va) * p), n_tr + n_va, p) %*% L
    colnames(X) <- paste0("v", seed_shift, "_f", seq_len(p))
    rownames(X) <- paste0("s", seq_len(n_tr + n_va))
    X
  }
  x1 <- mk_view(1); x2 <- mk_view(2)
  if (arm == "complementary") {
    eta <- 0.9 * x1[, 1] + 0.9 * x2[, 1]
  } else {
    eta <- 0.9 * x1[, 1]
  }
  y <- if (family == "binomial") rbinom(n_tr + n_va, 1, 1/(1 + exp(-eta))) else eta + rnorm(n_tr + n_va, sd = 0.5)
  y <- setNames(as.numeric(y), paste0("s", seq_len(n_tr + n_va)))
  tr <- seq_len(n_tr); va <- n_tr + seq_len(n_va)
  list(x1 = x1, x2 = x2, y = y, tr = tr, va = va,
       signals_v1 = "v1_f1", signals_v2 = if (arm == "complementary") "v2_f1" else character(0))
}

jaccard <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

balanced_error <- function(truth, pred_class) {
  truth <- as.numeric(truth); pred_class <- as.numeric(pred_class)
  if (length(unique(truth)) < 2) return(NA_real_)
  0.5 * (mean(pred_class[truth == 0] != 0) + mean(pred_class[truth == 1] != 1))
}

# ---- Mode-specific extraction ----------------------------------------------
# Returns list(perf, perf_source, sel_features, sel_note)
extract_mode <- function(fit, mode, family, y_valid, d) {
  if (mode == "early") {
    ef <- fit$early_fusion
    if (is.null(ef) || is.null(ef$selected_valid) || ncol(ef$selected_valid) == 0L) {
      return(list(perf = NA_real_, perf_source = "early_downstream_glm",
                  sel = if (!is.null(ef)) ef$selected_features else character(0),
                  note = "no selected features for downstream refit"))
    }
    out <- tryCatch({
      train_df <- as.data.frame(ef$selected_train); train_df$.y <- unname(d$y[d$tr])
      if (family == "binomial") {
        m <- stats::glm(.y ~ ., data = train_df, family = stats::binomial(link = "logit"))
        s <- unname(stats::predict(m, newdata = as.data.frame(ef$selected_valid), type = "response"))
        pf <- balanced_error(y_valid, as.numeric(s >= 0.5))
      } else {
        m <- stats::glm(.y ~ ., data = train_df, family = stats::gaussian())
        s <- unname(stats::predict(m, newdata = as.data.frame(ef$selected_valid)))
        pf <- mean((s - as.numeric(y_valid))^2)
      }
      list(perf = pf, perf_source = "early_downstream_glm",
           sel = ef$selected_features, note = NA_character_)
    }, error = function(e) list(perf = NA_real_, perf_source = "early_downstream_glm",
                                sel = ef$selected_features,
                                note = paste("downstream refit failed:", conditionMessage(e))))
    return(out)
  }
  if (mode == "late") {
    lf <- fit$late_fusion
    if (is.null(lf) || is.null(lf$valid_predictions)) {
      return(list(perf = NA_real_, perf_source = "late_native",
                  sel = c(fit$selected_features$view1, fit$selected_features$view2),
                  note = "no native valid_predictions"))
    }
    s <- as.numeric(lf$valid_predictions)
    if (family == "binomial") {
      pf <- balanced_error(y_valid, as.numeric(s >= 0.5))
    } else {
      pf <- mean((s - as.numeric(y_valid))^2)
    }
    return(list(perf = pf, perf_source = "late_native",
                sel = c(fit$selected_features$view1, fit$selected_features$view2),
                note = NA_character_))
  }
  if (mode == "cooperative") {
    cf <- fit$cooperative_fusion
    if (is.null(cf) || is.null(cf$valid_predictions)) {
      sel_cf <- if (!is.null(cf)) unlist(cf$selected_features) else character(0)
      return(list(perf = NA_real_, perf_source = "cooperative_native",
                  sel = sel_cf, note = "no native valid_predictions"))
    }
    s <- as.numeric(cf$valid_predictions)
    if (family == "binomial") {
      pf <- balanced_error(y_valid, as.numeric(s >= 0.5))
    } else {
      pf <- mean((s - as.numeric(y_valid))^2)
    }
    return(list(perf = pf, perf_source = "cooperative_native",
                sel = unlist(cf$selected_features), note = NA_character_))
  }
  stop("unknown mode")
}

rows <- list()
for (arm in c("complementary", "redundant")) {
  for (family in c("binomial", "gaussian")) {
    for (mode in c("early", "late", "cooperative")) {
      selections <- list(); perf <- c(); src <- c(); notes <- c()
      for (seed in seeds) {
        d <- make_data(arm, family, seed)
        lam <- data.frame(lambda = seq(0.02, 0.2, length.out = 8))
        res <- tryCatch(
          withCallingHandlers({
            fit <- stabl_multiomic_train_validate(
              x_train_list = list(view1 = d$x1[d$tr, ], view2 = d$x2[d$tr, ]),
              y_train = d$y[d$tr],
              lambda_grid = lam,
              x_valid_list = list(view1 = d$x1[d$va, ], view2 = d$x2[d$va, ]),
              y_valid = d$y[d$va],
              base_learner = "lasso", family = family,
              n_bootstraps = n_bootstraps,
              artificial_type = "random_permutation",
              random_state = seed,
              early_fusion = mode == "early",
              late_fusion = mode == "late",
              cooperative_fusion = mode == "cooperative"
            )
            ex <- extract_mode(fit, mode, family, d$y[d$va], d)
            selections[[seed]] <- ex$sel
            perf[seed] <- ex$perf; src[seed] <- ex$perf_source; notes[seed] <- ex$note
            data.frame(seed = seed, performance = ex$perf,
                       perf_source = ex$perf_source, note = ex$note)
          },
          warning = function(w) if (grepl("anti-conservative|not validated|deprecat", conditionMessage(w)))
            invokeRestart("muffleWarning")),
          error = function(e) data.frame(seed = seed, performance = NA_real_,
                                         perf_source = "error",
                                         note = conditionMessage(e)))
        rows[[length(rows) + 1L]] <- cbind(
          data.frame(arm = arm, family = family, mode = mode,
                     mean_seed_jaccard = NA_real_, view2_share = NA_real_),
          res
        )
      }
      js <- NA_real_
      if (length(selections) == length(seeds) && length(seeds) > 1) {
        js <- mean(combn(selections, 2, FUN = function(p) jaccard(p[[1]], p[[2]])))
      }
      all_sel <- unique(unlist(selections))
      v2_share <- if (length(all_sel) == 0) NA_real_ else
        mean(grepl("^v2_", all_sel))
      rows[[length(rows) + 1L]] <- data.frame(
        arm = arm, family = family, mode = mode, seed = NA,
        performance = mean(perf, na.rm = TRUE),
        perf_source = paste(unique(src), collapse = "|"),
        note = paste(unique(notes[!is.na(notes)]), collapse = " | "),
        mean_seed_jaccard = js, view2_share = v2_share
      )
    }
  }
}
results <- do.call(rbind, rows)
write.csv(results, file.path(out_dir, "ws4_results.csv"), row.names = FALSE)

sink(file.path(out_dir, "ws4_summary.md"))
cat("# WS-4 — Fusion benchmark summary\n\n")
cat("- Predeclared: MATURITY_PREDESIGN.md; 5 seeds; 100 bootstraps; lasso; random_permutation decoys.\n")
cat("- Held-out metric: BER (binomial) / MSE (gaussian) on the 50-sample validation split.\n")
cat("- Evaluation rule: late/cooperative use native valid_predictions; early uses a downstream\n",
    "  stats::glm refit on the concatenated selected features (stabl_fit stores no reusable model).\n\n")
cat("## Mean validation performance and selection stability\n\n")
sum_cols <- intersect(c("arm", "family", "mode", "performance", "perf_source",
                        "mean_seed_jaccard", "view2_share", "note"),
                      names(results))
print(results[is.na(results$seed), sum_cols], row.names = FALSE)
cat("\n## Per-seed rows\n\n")
print(results[!is.na(results$seed), sum_cols], row.names = FALSE)
cat("\n## Caveats\n\n- Simulation-only; two-view designs; view-2 noise in the redundant arm is the deliberate null probe.\n",
    "- Cooperative fusion: built-in vendored engine, gaussian/binomial only (Cox deferred — see roadmap WS-4).\n",
    "- Early-fusion held-out performance depends on the documented downstream glm refit; it is not a\n",
    "  package-native prediction path (stabl_fit stores no model).\n")
sink()
cat("WS-4 done.\n")
