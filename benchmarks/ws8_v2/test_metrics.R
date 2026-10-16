# WS-8 v2 — fast synthetic tests for metric/export helpers (gate G1 unit tests).
#
# Toy-data invariants that run without the atlas.  Run from repo root:
#   Rscript benchmarks/ws8_v2/test_metrics.R
# Exits non-zero on any failure.

source("benchmarks/ws8_v2/metrics.R")

pass <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("  ok:", msg, "\n") }
  else { fail <<- fail + 1L; cat("  FAIL:", msg, "\n") }
}

cat("== macro vs micro ==\n")
# Unbalanced toy: cohort A (n=3) errors 1,1,1; cohort B (n=1) error 10.
errs <- c(1, 1, 1, 10); grp <- c("A", "A", "A", "B")
mm <- mae_macro_micro(errs, grp)
# macro = mean(1, 10) = 5.5 ; micro = 13/4 = 3.25
check(abs(mm$macro - 5.5) < 1e-9, "macro = 5.5 on unbalanced toy")
check(abs(mm$micro - 3.25) < 1e-9, "micro = 3.25 on unbalanced toy")
check(abs(mm$macro - mm$micro) > 1, "macro != micro when cohort sizes differ")

cat("== safe_mean_scores (fixes A1) ==\n")
es <- matrix(c(0.9, 0.5, 0.3,
               0.8, 0.4, 0.6), nrow = 3, byrow = TRUE,
             dimnames = list(c("f1", "f2", "f3"), c("c1", "c2")))
ms <- safe_mean_scores(es)
check(length(ms) == 3 && all(names(ms) == c("f1", "f2", "f3")),
      "mean score has one value per feature, named")
check(abs(ms["f1"] - 0.70) < 1e-9, "f1 mean = (0.9+0.5)/2 = 0.70")
check(!any(is.na(ms)), "no NA in mean scores")

cat("== safe_match_scores (fixes A2) ==\n")
sc <- c(a = 0.7, b = 0.2, c = 0.9)
out <- safe_match_scores(sc, c("c", "a"))
check(identical(unname(out), c(0.9, 0.7)), "name-matched before numeric coercion")
check(identical(names(out), c("c", "a")), "names preserved")
check(!any(is.na(out)), "no NA after name matching")

cat("== signed_rho_replication ==\n")
feats <- c("x", "y", "z")
tsign <- c(x = 1, y = -1, z = 1)
hrho  <- c(x = 0.5, y = -0.4, z = -0.3)
rr <- signed_rho_replication(feats, tsign, hrho)
# signed: x=+0.5, y=+0.4, z=-0.3 -> R = 0.2
check(abs(rr$R - 0.2) < 1e-9, "R = mean(0.5, 0.4, -0.3) = 0.2")
check(abs(rr$sign_agreement - 2/3) < 1e-9, "sign agreement = 2/3")

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0L) quit(status = 1L)
cat("ALL METRIC/EXPORT TESTS PASSED\n")
