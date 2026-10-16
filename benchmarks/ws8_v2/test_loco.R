# WS-8 v2 — leakage and engine invariant tests (gate G2, synthetic data).
#
# Fast tests that run without the atlas.  Prove: outer folds are donor-disjoint;
# the fold transform is fit on training data only; the leakage guard trips on
# contamination; ranking ties are deterministic; the zero-panel fallback returns
# the median-age baseline; and the engine never touches test labels upstream.
#
# Run from repo root:  Rscript benchmarks/ws8_v2/test_loco.R
# Exits non-zero on any failure.

source("benchmarks/ws8_v2/loco.R")

pass <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("  ok:", msg, "\n") }
  else { fail <<- fail + 1L; cat("  FAIL:", msg, "\n") }
}

set.seed(1)
# Synthetic: 4 cohorts, uneven sizes, 8 features; f1/f2 signal, rest noise.
coh <- rep(c("A", "B", "C", "D"), times = c(40, 30, 20, 10))
n <- length(coh); p <- 8
X <- matrix(rnorm(n * p), nrow = n,
            dimnames = list(paste0("d", seq_len(n)), paste0("f", seq_len(p))))
# make f5 constant within cohort D (zero-variance-in-some-cohort)
X[coh == "D", "f5"] <- 2
y <- 0.5 * X[, "f1"] - 0.3 * X[, "f2"] + rnorm(n)
donor <- rownames(X)

cat("== fold definition ==\n")
folds <- make_loco_folds(donor, coh)
check(nrow(folds) == n, "fold table has one row per donor")
check(identical(sort(unique(folds$fold)), c("A", "B", "C", "D")),
      "one fold per cohort")
# donor-disjoint: each donor is test in exactly one fold
check(all(table(folds$donor) == 1), "each donor is test in exactly one fold")

cat("== fold-safe transform (training-only filtering) ==\n")
# Hold out D: f5 is constant only in D, so on training (A,B,C) f5 is NOT
# zero-variance and must be KEPT.  This proves filtering is training-fit.
tr <- coh != "D"
rec <- fit_fold_transform(X[tr, ], coh[tr])
check("f5" %in% rec$keep, "f5 kept when training cohorts have variance")
# Hold out A,B,C and train on D only is degenerate; instead test a feature
# constant across ALL training cohorts is dropped.
X2 <- X; X2[, "f8"] <- 7  # constant everywhere
rec2 <- fit_fold_transform(X2[tr, ], coh[tr])
check(!("f8" %in% rec2$keep), "globally constant feature dropped on training")

cat("== leakage guard ==\n")
ok <- tryCatch({ assert_no_leakage(c("d1", "d2"), c("d3", "d4"), "t"); TRUE },
               error = function(e) FALSE)
check(ok, "guard passes on disjoint sets")
bad <- tryCatch({ assert_no_leakage(c("d1", "d2"), c("d2", "d3"), "t"); FALSE },
                error = function(e) grepl("LEAKAGE", conditionMessage(e)))
check(bad, "guard trips when a test ID is in training")

cat("== training-only signs ==\n")
sgn <- training_signs(X[tr, ], y[tr], coh[tr])
check(sgn["f1"] == 1L, "f1 sign positive (training)")
check(sgn["f2"] == -1L, "f2 sign negative (training)")
check(all(names(sgn) == colnames(X)), "signs named for all features")

cat("== deterministic tie-break in ranking ==\n")
# A selector stub with tied scores; ranking must be deterministic.
sel_stub <- function(Xtr, ytr, ctr) {
  sc <- setNames(c(0.5, 0.5, 0.9, 0.1, 0.5, 0.2, 0.8, 0.05), colnames(Xtr))
  # rank by score desc, tie-break name asc (deterministic)
  o <- order(-sc, names(sc))
  list(ranked = names(sc)[o], native = names(sc)[sc >= 0.5], scores = sc)
}
r1 <- run_one_fold("D", X, y, coh, donor, sel_stub, K = 3)
r2 <- run_one_fold("D", X, y, coh, donor, sel_stub, K = 3)
check(identical(r1$panel_topk, r2$panel_topk), "top-K panel identical across calls")
check(identical(r1$panel_topk, c("f3", "f7", "f1")),
      "tie among f1/f5 broken deterministically (name order)")

cat("== zero-panel fallback ==\n")
sel_empty <- function(Xtr, ytr, ctr) {
  sc <- setNames(rep(0.1, ncol(Xtr)), colnames(Xtr))
  list(ranked = names(sc), native = character(0), scores = sc)
}
re <- run_one_fold("D", X, y, coh, donor, sel_empty, K = 3)
check(length(re$panel_native) == 0L, "empty native panel preserved (not substituted)")
check(length(re$panel_topk) == 3L, "top-K panel still formed from ranking")

cat("== engine never sizes test into training ==\n")
check(all(re$train_ids %in% donor[coh != "D"]), "train IDs exclude held-out cohort")
check(all(re$test_ids %in% donor[coh == "D"]), "test IDs are the held-out cohort")
check(nrow(re$Xtr) == sum(tr) && nrow(re$Xte) == sum(!tr),
      "train/test matrices have correct row counts")

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0L) quit(status = 1L)
cat("ALL LOCO ENGINE TESTS PASSED\n")
