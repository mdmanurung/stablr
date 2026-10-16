# Adversarial Scientific Audit Program — 2026-10-17

Comprehensive audit of stablr per the fifteen-phase adversarial program
commissioned 2026-10-17. This directory is the program's persistent record:
coverage ledger, findings, specification matrix, and validation log.

## Audited revision

- Commit: `5cf603b` (main, PR #21), package version 0.3.0
- Working tree: 40 untracked files (excluded S1 raws) + in-flight v2
  reproduction run writing `benchmarks/ws8_v2/results/` (validation task,
  not audited source changes)
- Environment: R 4.4.2, gcc 13.3.0, Ubuntu 24.04.5, 16 cores / 16 GB
- Key deps: glmnet 5.0, ranger 0.18.0, knockoff 0.3.6, data.table 1.18.6.1,
  ggplot2 4.0.3, Rcpp 1.1.2, RcppArmadillo 15.6.0.1, testthat 3.3.2

## Program structure

| Phase | Scope | Status |
|---|---|---|
| 0 | Reconnaissance + fresh baseline | **complete** (tests 3325 pass/0 fail/1 skip; R CMD check 0E/0W/1NOTE) |
| 1 | Specification reconstruction from primary literature | **complete** (delegated independent pass, integrated + spot-verified; F-001..F-010) |
| 2 | Core algorithm audit (bootstrap, learners, knockoffs, FDP+) | in progress (fdp_control.R, artificial_features.R, knockoff_mvr.R clean; MVR solver numerically verified) |
| 3-5 | Leakage, multi-omic, outcome-family audits | pending |
| 6 | Independent validation program (synthetic DGPs) | pending |
| 7 | Aging case-study adversarial audit | pending (v2 bit-for-bit run in flight as input) |
| 8-9 | Software, compiled code, test-suite audit | compiled-code pass **complete** (F-011..F-018, P0 runtime-confirmed); test-suite audit pending |
| 10-11 | Documentation + distribution standards | pending |
| 12-15 | Publication readiness, remediation, referee report | pending |

## Files

- `coverage-ledger.csv` — 49 R files, 72 exported functions, test mapping,
  per-file audit disposition (updated as the audit proceeds)
- `findings.md` — living findings report (P0-P3, evidence-linked)
- `specification-matrix.md` — Phase 1 output (from the independent pass,
  integrated and verified by the lead auditor)
- `baseline-2026-10-17.md` — fresh baseline evidence

## Standing rules observed

- No silent workarounds; no benchmark-driven overfitting; no fabricated
  evidence. Adverse findings are preserved.
- The AI audit does not satisfy the package's independent human-review
  requirement (SCI-04 confirmation packet remains a human task).
- Frozen benchmark protocols are never retroactively edited; amendments are
  new versions.

## Baseline (fresh, this commit)

Recorded in `baseline-2026-10-17.md` upon completion of the background runs.
