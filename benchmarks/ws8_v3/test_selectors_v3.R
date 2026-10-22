# WS-8 v3 — unit tests for new selector components.  Fail loudly on violation.
suppressPackageStartupMessages({
  library(stablr); library(glmnet); library(matrixStats); library(sva)
})
v3_dir <- Sys.getenv("WS8_V3_DIR", "benchmarks/ws8_v3")
engine_dir <- Sys.getenv("WS8_ENGINE", "benchmarks/ws8_v2")
source(file.path(engine_dir, "metrics.R"))
source(file.path(engine_dir, "loco.R"))
source(file.path(engine_dir, "selectors.R"))
source(file.path(v3_dir, "selectors_v3.R"))

pass <- function(name) cat(sprintf("PASS  %s\n", name))
fail <- function(name, msg) stop(sprintf("FAIL  %s: %s", name, msg), call. = FALSE)
check <- function(name, cond, msg = "") if (isTRUE(cond)) pass(name) else fail(name, msg)

# ---- Toy design with known cohort structure -----------------------------------
set.seed(1)
n_per <- 60; K <- 4; p <- 12
cohort <- factor(rep(paste0("C", 1:K), each = n_per))
n <- n_per * K
X <- matrix(rnorm(n * p), n, p, dimnames = list(paste0("s", 1:n), paste0("f", 1:p)))
# feature 1: homogeneous age effect; feature 2: cohort-2-only effect;
# feature 3: batch shift in cohort 3 (no age effect)
age <- rnorm(n, 50, 15)
X[, 1] <- X[, 1] + 0.5 * scale(age)[, 1]
X[cohort == "C2", 2] <- X[cohort == "C2", 2] + 0.8 * scale(age[cohort == "C2"])[, 1]
X[cohort == "C3", 3] <- X[cohort == "C3", 3] + 2

# ---- T1: ComBat implementation matches sva::ComBat on training data -----------
mod <- model.matrix(~ age)[, -1, drop = FALSE]
cf <- .combat_fit(X, cohort, mod)
ref <- sva::ComBat(dat = t(X), batch = cohort, mod = model.matrix(~age),
                   par.prior = TRUE, mean.only = FALSE)
ref <- t(ref)
delta <- max(abs(cf$adjusted - ref) / (1 + abs(ref)))
check("T1 combat_matches_sva", delta < 1e-6, sprintf("max rel diff %.2e", delta))

# ---- T2: frozen apply is deterministic and label-free --------------------------
# Covariate-free fit: with mod = NULL the training stand_mean IS the grand
# mean, so eb_newbatch on a training batch must reproduce the training
# adjustment to machine precision (same moments, same per-batch hyperpriors).
cf0 <- .combat_fit(X, cohort, mod = NULL)
a1 <- .combat_apply(cf0, X[cohort == "C1", ], mode = "eb_newbatch")
a2 <- .combat_apply(cf0, X[cohort == "C1", ], mode = "eb_newbatch")
check("T2a apply_deterministic", identical(a1, a2))
a3 <- .combat_apply(cf0, X[cohort == "C3", ], mode = "standardize_only")
a1c3 <- .combat_apply(cf0, X[cohort == "C3", ], mode = "eb_newbatch")
check("T2b modes_differ", max(abs(a1c3 - a3)) > 0,
      "eb_newbatch and standardize_only should differ when batch offset exists")
d_tr <- max(abs(a1 - cf0$adjusted[cohort == "C1", ]) / (1 + abs(cf0$adjusted[cohort == "C1", ])))
check("T2c newbatch_on_train_exact", d_tr < 1e-8, sprintf("rel diff %.3e", d_tr))

# ---- T3: M6 reduces to M0 ranking when filter disabled -------------------------
s6 <- sel_m6_meta_homog(X, age, cohort, q_thresh = 0, bh_thresh = 2)  # pass all
s0 <- sel_m0_meta(X, age, cohort)
check("T3a m6_eq_m0_ranking", identical(s6$ranked, s0$ranked))
check("T3b m6_all_native", length(s6$native) == p)
s6f <- sel_m6_meta_homog(X, age, cohort)
check("T3c m6_finds_homog_signal", "f1" %in% s6f$native,
      paste("native:", paste(s6f$native, collapse = ",")))
check("T3d m6_schema", is.character(s6f$ranked) && is.numeric(s6f$scores) &&
        !is.null(names(s6f$scores)) && identical(names(s6f$scores), colnames(X)))
# heterogeneity stats present and sane
check("T3e m6_extra", all(c("pooled_padj", "q_p", "tau2", "I2") %in% names(s6f$extra)) &&
        all(s6f$extra$q_p >= 0 & s6f$extra$q_p <= 1))

# ---- T4: M8 contract + omega behavior on toy design ----------------------------
s8 <- sel_m8_contrast_stabl(X, age, cohort, n_boot = 20L, seed = 1, omega = 2)
check("T4a m8_schema", is.character(s8$ranked) && length(s8$ranked) == p &&
        identical(names(s8$scores), colnames(X)))
check("T4b m8_finds_common", which.max(s8$scores) == 1L,
      sprintf("top feature = %s (score %.2f)", names(which.max(s8$scores)), max(s8$scores)))
s8w <- sel_m8_contrast_stabl(X, age, cohort, n_boot = 20L, seed = 1, omega = 8)
check("T4c m8_omega_stable", which.max(s8w$scores) == 1L,
      "common signal should survive stronger deviation penalty")
# omega -> Inf approximates pooled lasso on common terms only: cohort-2-only
# signal (f2) should be selected LESS under high omega
check("T4d m8_omega_direction",
      s8w$scores["f2"] <= s8$scores["f2"] + 1e-9,
      sprintf("f2 score omega=2: %.2f, omega=8: %.2f", s8$scores["f2"], s8w$scores["f2"]))

# ---- T5: M9 contract + transform hook ------------------------------------------
s9 <- sel_m9_combat_stabl(X, age, cohort, n_boot = 20L, seed = 1)
check("T5a m9_schema", is.character(s9$ranked) && length(s9$ranked) == p &&
        is.function(s9$test_transform) && !is.null(s9$train_transformed))
xt <- s9$test_transform(X[1:10, ])
check("T5b m9_transform_shape", is.matrix(xt) && nrow(xt) == 10 && ncol(xt) == p &&
        identical(colnames(xt), colnames(X)))
check("T5c m9_finds_signal", which.max(s9$scores) == 1L,
      sprintf("top = %s", names(which.max(s9$scores))))

# ---- T6: deterministic tie-break ------------------------------------------------
s6b <- sel_m6_meta_homog(X, age, cohort)
check("T6 deterministic", identical(s6f$ranked, s6b$ranked))

cat("\nALL TESTS PASSED\n")
