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
| 2 | Core algorithm audit (bootstrap, learners, knockoffs, FDP+) | **complete** (fdp_control.R, artificial_features.R, knockoff_mvr.R clean; MVR solver numerically verified — L3 partially resolved; completion findings F-019..F-023; first-pass leads L1–L5 dispositioned) |
| 3-5 | Leakage, multi-omic, outcome-family audits | **complete** (Phases 3+4 delegated independent passes + Phase 5 lead; no P0/P1; F-CV-01..08, P4-F24..F29, F-024..F-026) |
| 6 | Independent validation program (synthetic DGPs) | **complete** (F-P6-1 headline: min_fdr_ anti-conservative 1.5-2.2× on pure nulls; regime-dependent calibration; env-permutation null clean) |
| 7 | Aging case-study adversarial audit | **complete** (V8/V12/V13 + vignette §4/§10.3/§10.5 re-derived exactly from committed artifacts; saturation mechanism attributed: rf threshold-path median floor + max-over-lambda scoring + fixed 0.5 bar decoupled from FDP+ — F-P7-1; CP convention inconsistency F-P7-2; WS-5 calibration-scope gap F-P7-3; evidence in phase7_validation/) |
| 8-9 | Software, compiled code, test-suite audit | compiled-code pass **complete** (F-011..F-018; F-016 corrected to P3 fragility after runtime probe); test-suite audit **complete** (all ledger gaps refuted as false leads) |
| 10-11 | Documentation + distribution standards | **complete** (F-P10-1 release-contract skew carried) |
| FR | Final pre-publication per-function review (all 74 exports) | **complete** (2026-10-18: 28 OK / 34 OK-with-notes / 12 issue; 1 P1 fixed — rf auto-grid NULL-seed failure B-01; 5 fail-closed robustness fixes + regression tests; 9 doc corrections; advanced-vignette Jaccard-heatmap axis fix; 27 findings deferred to Phase 13 with rationale — see `final-review.md`) |
| 12 | Cross-phase synthesis + residual risk + Phase 13 plan | **complete** (2026-10-18: 90 finding IDs consolidated — 73 open/owed, 15 fixed in `6aa5a6e`, 1 refuted, 1 positive; release blockers F-011/F-012/F-013/F-P6-1/F-001/F-002 verified open at HEAD; merged 6-workstream Phase 13 plan; program verdict: not release-ready pending PR-1..PR-3 + fresh release run + SCI-04 — see `phase12-synthesis.md`) |
| 13 | Remediation execution | **in progress** (PR-1 DONE 2026-10-18: F-011/F-012/F-013 fail-closed fixes + compiled-boundary regression tests, 9 blocks / 25 expectations; suite 0 failed / 0 errors; parity contract intact. PR-2 DONE 2026-10-18: F-CV-01/02/03 + P4-F28 RNG-hygiene fixes — hermetic `stabl_fit`, pre-fork per-fold seeds, per-consumer stream derivation; 2 additional reachable unseeded `cv.glmnet` foldid sites fixed (F-CV-07 reachability corrected); `test-rng-hygiene.R` 9 blocks; suite 602 tests 0 failed / 0 errors, 52 warnings at pre-PR-2 parity; Python-parity contract intact, R-internal multiomic golden values deliberately updated with dated rationale. PR-A DONE 2026-10-18: 1.0.0-spec contract closure — F-CV-09/10/11 RNG residuals fixed (`.with_local_seed` wrappers, outputs identical); `explore`/`n_explore` removed (empty support stays empty, spec §5.4); `python_legacy` + `cooperation_selection="validation"` fenced from the public workflows (F-CV-04/P4-F24 fixed by removal; internal branches retained for characterization tests, char-F1/F2-val goldens preserved); FDP+ presentation relabelled (`plot_fdr_graph` "FDP+ diagnostic", `fdr_target` → `reference_fdp_plus`, `min_fdr_` stochastic-bound docs, strict-`>`/tie-policy contract documented); F-019 fixed (hashed per-environment streams); F-CV-05 verified unreachable with regression-pinned invariant; suite 605 tests 0 failed / 0 errors, 52 warnings at exact parity; R CMD check 0E/0W/1N (pre-existing size NOTE); parity contracts byte-identical — next PR-B vignette bundle per the 1.0.0 spec program) |
| 14-15 | Reconciliation, referee report | pending (release gates unchanged) |

## Files

- `coverage-ledger.csv` — 49 R files, 72 exported functions, test mapping,
  per-file audit disposition (updated as the audit proceeds)
- `findings.md` — living findings report (P0-P3, evidence-linked)
- `final-review.md` — final pre-publication per-function review (2026-10-18):
  verdict rollup for all 74 exports, new findings N-01..N-06, B-01..B-08,
  C-1..C-14, D-1..D-13 with dispositions, ensemble-coherence notes,
  defaults-calibration verdict, publication-readiness verdict
- `phase12-synthesis.md` — Phase 12 cross-phase synthesis (2026-10-18):
  per-phase program status, consolidated 90-ID findings inventory with
  HEAD-verified remediation status, release-blocking set, residual-risk
  statement, merged Phase 13 remediation plan, program-level verdict
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
