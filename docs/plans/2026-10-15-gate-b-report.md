# Gate B report — binary MVP evidence (PR6-B)

Date: 2026-10-15. Branch: docs/pkgdown-0.3.0. Evidence HEAD: 2eb6929
(PR0d). This report collects the predeclared evidence for the binary
gate: the IV interpretation/validated-workflow benchmark, the full test
suite, and the WS-8/WS-9 repairs. Per governance, adverse results are
reported as measured, not suppressed.

## 1. IV benchmark (predeclared 2026-10-15, run after the predesign)

Script: `benchmarks/interpretation_validation.R`. Predesign:
`benchmarks/interpretation_validation_predesign.md` (frozen before
any run; one dated addendum records two post-run repairs with the design
unchanged). Outputs: `iv_results.csv` (110 runs), `iv_effects.csv` (350),
`iv_folds.csv` (190), `iv_importance.csv` (550), `iv_failures.csv` (0),
`iv_summary.md`.

Run integrity: 110 planned, 110 completed, 0 failed; the script exits
non-zero if any planned run is missing (exit code 0 observed).

| Gate | Criterion | Value | Verdict |
|---|---|---|---|
| G-IV-1 | completeness (all planned runs accounted for) | 110/110, 0 failed | MET |
| G-IV-2 | IV-1 effect-direction recovery >= 0.90 | 1.00 | MET |
| G-IV-3 | IV-2 substitutes in one group >= 0.80 | 0.60 | NOT MET |
| G-IV-4 | IV-5 leaky mean AUROC >= nested mean AUROC | +0.0148 (0.830 vs 0.815) | MET |
| G-IV-5 | IV-0 median validated AUROC in [0.35, 0.65] | 0.473 | MET |
| G-IV-6 | IV-3 one-feature folds non-degenerate | 4 folds, 0 degenerate | MET |

### G-IV-3 adverse finding (committed, not suppressed)

The criterion failed on the selection layer, not the grouping layer.
Across the 10 IV-2 seeds, whenever f1a/f1b/f1c were co-selected they were
ALWAYS assigned to a single correlation group (10/10 seeds, no split).
In 4 seeds (5, 6, 7, 10) FDP+ selection did not retain all three
substitutes (seed 5 selected f1a only), so the predeclared
"all three selected AND grouped" operationalization reached only 0.60.
Interpretation: correlation grouping is correct on whatever is selected;
stability selection with FDP+ control prunes redundant substitutes, which
is expected selection behavior, not a grouping defect. Users should not
expect all members of a correlated block to be co-selected.

### Other notable measurements

- IV-0 null validated mode: mean accuracy 0.515, mean AUROC 0.468,
  median AUROC 0.473 — null behaviour is centred near chance as required.
- IV-5: leaky full-data selection is directionally optimistic over
  correctly nested selection on identical outer folds (+0.015 mean AUROC,
  pooled across seeds; not a per-replicate requirement).
- IV-4 study-disjoint validated mode completes under prevalence and
  feature-scale drift (mean AUROC 0.828, 4 leave-one-study-out folds).
- Validated-mode signal designs (IV-1/IV-2) reach mean AUROC 0.815/0.820
  with median 3-4 selected features against 2-3 true signals.

## 2. Full test suite

`testthat::test_local(".")`: 3263 pass / 0 fail / 1 skip at evidence HEAD
(includes the PR0d empty-selection regression test and the IV contract
tests in `tests/testthat/test-interpretation-validation.R`, 48
expectations).

## 3. WS-8 / WS-9 repair status (context for Gate B language)

- WS-9 (PR0b, commit 5939565): corrected 3-feature truth set; 38/40 runs
  completed (2 infeasible-bootstrap failures on seed 8, recorded).
  knockoff_equi mean FDP 0.140 (O1) / 0.156 (O2), TPR 1.0. Graduation
  NOT MET (incomplete grid; O2 FDP > 0.15) — adverse result committed.
  WS-9 is an ordinal (Gate S) diagnostic; it does not block binary Gate B
  but forbids ordinal graduation language.
- WS-8 (PR0c, commit e60bfb7): evidence labels corrected to
  post-selection descriptive; no predictive-validity claim.

## 4. Repairs triggered by this gate

- PR0d (commit 2eb6929): empty-selection folds aborted
  `.analyze_fold_importance()`; fixed with a regression test. Found by
  IV-0 seeds 6/9/10 before any gate claims were made.

## 5. Verdict

Gate B evidence collection is COMPLETE and the binary exploratory /
validated workflow is supported by predeclared synthetic evidence:
5 of 6 IV criteria MET, one adverse criterion (G-IV-3) committed with a
selection-layer explanation. Validated mode remains binary-only;
multiclass and ordinal graduation remain reserved for the Gate S track
(WS-9 graduation NOT MET stands). No application-data claims are made;
all IV evidence is simulation-only.
