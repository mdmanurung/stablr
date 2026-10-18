# Phase 12 — Cross-phase synthesis (2026-10-18)

Consolidates the entire fifteen-phase adversarial audit program into one
program-level record: per-phase status, the complete findings inventory with
current remediation status, the verified release-blocking set, a residual-risk
statement, the Phase 13 remediation plan, and an honest program-level
publication-readiness verdict.

**Audited revision:** `6aa5a6e` (main). **Synthesis date:** 2026-10-18.
**Inputs:** `findings.md`, `final-review.md`, `remediation-matrix.md`,
`coverage-ledger.csv`, `specification-matrix.md`, `baseline-2026-10-17.md`,
`phase6_validation/`, `phase7_validation/`.

**Verification method (this synthesis):** every open/closed status claim below
was re-verified against HEAD `6aa5a6e` on 2026-10-18, not carried forward from
earlier notes: `git log 5cf603b..HEAD -- src/ R/cooperative-multiview.path.R`
is **empty** (no commit since the audit baseline touched the vendored engine),
and the code sites of F-011, F-012, F-013, F-005, F-CV-02, and F-P10-1 were
inspected directly (evidence in §4).

---

## 1. Program status by phase

| Phase | Scope | Status | Headline outcome |
|---|---|---|---|
| 0 | Reconnaissance + fresh baseline | **complete** | Tests 3325 pass / 0 fail; R CMD check 0E/0W/1NOTE at `5cf603b` |
| 1 | Specification reconstruction (primary literature) | **complete** | F-001..F-010; two P1 interpretive/extension gaps (F-001, F-002) |
| 2 | Core algorithm audit (bootstrap, learners, knockoffs, FDP+) | **complete** | fdp_control / artificial_features / knockoff_mvr clean; MVR solver numerically verified (L3 partially resolved); F-019..F-023 |
| 3 | Leakage audit | **complete** | No leakage found; RNG-hygiene cluster F-CV-01..08 (P2: 01–04) |
| 4 | Multi-omic audit | **complete** | Alignment clean; P4-F24..F28 robustness/doc items; P4-F29 refuted |
| 5 | Outcome-family audit | **complete** | F-024..F-026 (P3); ordinal stays experimental (WS-9) |
| 6 | Independent validation (synthetic DGPs) | **complete** | **F-P6-1 (P1): `min_fdr_` anti-conservative 1.5–2.2× on pure nulls** (realized FDP 0.918/0.934 vs reported 0.423/0.662; 18–33 features on no-signal data); F-P6-2 regime dependence; F-P6-3 positive (env-permutation null clean, max min-score 0.168); winner's-curse mechanism quantified (per-threshold mean 1.022 → min-over-100-thresholds 0.674); F-001 bound holds (P(FDP ≥ (1+Δ)q+) = 0.60/0.47/0.17/0.03 at Δ = 0.5/1/2/5 vs bounds 0.67/0.50/0.33/0.17 — within bound, nearly tight at small Δ) |
| 7 | Aging case-study adversarial audit | **complete** | V8/V12/V13 + vignette §4/§10.3/§10.5 re-derived exactly; **F-P7-1 (P2): null-saturation mechanism attributed** (rf threshold-path median floor + max-over-lambda scoring + fixed 0.5 bar that never consults FDP+); F-P7-2 (CP convention inconsistency); F-P7-3 (WS-5 calibration-scope gap) |
| 8 | Software + compiled-code audit | **complete** | **F-011 (P0): sparse-path OOB read/SIGABRT**; F-012/F-013 (P1); F-014..F-018; F-016 corrected P2→P3 (crash claim refuted; `$` partial-matching fragility stands) |
| 9 | Test-suite audit | **complete** | 3325 baseline assertions verified; all ledger "no tests" flags refuted as token-matcher false leads (L5 refuted); real gap: no tests at the compiled boundary (links to F-011..F-013) |
| 10–11 | Documentation + distribution standards | **complete** | Dependency hygiene, NEWS honesty, vignettes clean; F-P10-1 (P2): locked release contract is 0.1.1-era vs DESCRIPTION 0.3.0; F-P10-2 (P3) |
| FR | Final per-function review (all 74 exports) | **complete** (2026-10-18) | 28 OK / 34 OK-with-notes / 12 issue; 42 new findings: **B-01 (P1) fixed in-commit** + 6 code fixes + 8 doc IDs fixed; 27 deferred to Phase 13; suite 54 files 0 failures; check 0E/0W/2N; PASS 3067 in-check |
| **12** | **Cross-phase synthesis (this document)** | **complete** (2026-10-18) | 90 finding IDs consolidated: 73 open/owed, 15 fixed in `6aa5a6e`, 1 refuted, 1 positive; release blockers verified open at HEAD (§4) |
| 13 | Remediation execution | pending | Plan in §7 below |
| 14 | Reconciliation | pending | After remediation |
| 15 | Referee report | pending | Final deliverable |

---

## 2. Consolidated findings inventory

Ninety distinct finding IDs are tracked across the program (48 from Phases
1–11, 42 from the final review), plus the five Phase 2 first-pass leads
(L1–L5, dispositioned below). Status as of HEAD `6aa5a6e`:

| Status | Count | Items |
|---|---|---|
| **Open — P0** | 1 | F-011 |
| **Open — P1** | 5 | F-001, F-002, F-012, F-013, F-P6-1 |
| **Open — P2** | 16 | F-003, F-004, F-005, F-006, F-014, F-015, F-017, F-CV-01..04, P4-F24, P4-F28, F-P6-2, F-P7-1, F-P10-1 |
| **Open — P3** | 27 | F-007, F-008, F-010, F-016, F-018, F-019..F-026, F-CV-05..08, P4-F25..F27, F-P7-2, F-P7-3, F-P10-2, B-04, C-4, D-2, D-4 |
| **Open — P4** | 23 | N-03..N-06, B-06..B-08, C-7..C-14, D-5..D-8, D-10, D-12, D-13, FR-hygiene |
| Open — validation task (not a defect) | 1 | F-009 (cooperative-engine objective vs direct augmented-matrix solve) |
| **Fixed in `6aa5a6e`** | 15 | Code: B-01 (P1), B-02, N-02, D-1, C-3, C-5, D-3. Docs: N-01, B-03, B-05, C-1, C-2, C-6, D-9, D-11. Plus the pre-existing `test-public-interface-audit` correction (72 → 74 exports) |
| Refuted | 1 | P4-F29 (L5, a lead, also refuted) |
| Positive validation result | 1 | F-P6-3 (env-permutation null clean) |

**Total open/owed to Phase 13: 73 items** (72 defects + the F-009 validation
task). Phase 2 lead dispositions: L1 folded into F-005; L2 (PD repair is a
documented approximation) folded into the documentation cluster; L3
**partially resolved** (MVR solver property tests + first-order-condition
re-derivation all passed; knockpy cross-comparison blocked — not installed);
L4 (max-over-lambda + threshold-minimization optimism) resolved
quantitatively by the Phase 6 winner's-curse analysis and F-P7-1; L5 refuted.

### 2.1 What was fixed in the final-review commit (`6aa5a6e`)

Seven code fixes with regression tests (`tests/testthat/test-final-review-regressions.R`,
17 expectations): B-01 (P1 — `auto_rf_threshold_grid` default `seed=NULL`
always failed; fixed via `%||%` draw in `.rf_fit_permutation_importance`),
B-02 (fail-closed on NULL importance), N-02 (`.validate_proportion` at the
`get_support.stabl_fit` choke point), D-1 (`length(y)==nrow(x)` checks in
`boxplot_features`/`scatterplot_features`), C-3 (lone `y_valid` fails
closed), C-5 (`.cooperative_prediction_matrix` row-mismatch fails closed),
D-3 (advanced-vignette Jaccard heatmap: `remove_diag=FALSE` + explicit seed
dimnames). Eight documentation IDs (N-01, B-03, B-05, C-1, C-2, C-6, D-9,
D-11), including hand-edits to the three non-roxygen man pages. Full test
suite 54 files / 0 failures; `R CMD check` 0E/0W/2N; all 9 vignettes re-render.

---

## 3. Scope reconciliation: the final-review verdict vs the program ledger

`final-review.md` states: "After this revision pass there are **no known P1
defects**." That verdict is **scoped to the final review's own new findings**
(the 42 per-function findings): the one P1 it raised (B-01) was fixed in the
same commit and every other item was fixed or explicitly deferred.

It does **not** mean the program-level ledger is clear of P1s. The program-level
open P0/P1 set (F-011, F-012, F-013, F-P6-1, plus the interpretive P1s F-001
and F-002) predates the final review, was deliberately deferred to Phase 13 by
the program's sequencing, and remains open at HEAD — verified in §4. Any
reading of the final-review verdict as "the package has no open P1 defects"
would be wrong. The correct program-level statement is §9.

---

## 4. Open release-blocking findings (verified at HEAD `6aa5a6e`, 2026-10-18)

| ID | Sev | Defect | Verification evidence at HEAD | Remediation class |
|---|---|---|---|---|
| **F-011** | **P0** | Sparse-path out-of-bounds read (SIGABRT observed) when `attr(x,"xm"/"xs")` is missing; silent heap OOB under `-DNDEBUG` | `git log 5cf603b..HEAD -- src/ R/cooperative-multiview.path.R` is **empty**; `R/cooperative-multiview.path.R` ~1071–1072 still passes `as.double(attr(x,"xm"))` / `as.double(attr(x,"xs"))` to `spwls_exp` with no length/positivity validation | Code fix in `elnet.fit` before `.Call` (length-validate `xm`/`xs` vs `ni`, `xs > 0`, else `jerr=10001`) + regression test |
| **F-012** | P1 | NA in `x` silently drops columns (`jerr=0`, no warning) | No `anyNA` guard on `x`/`y`/`weights` in `elnet.fit` (only an unrelated `mu.eta.val` NA check at line 615) | `anyNA` check + regression test |
| **F-013** | P1 | maxit non-convergence returns bare `list(jerr=-1)`, no warning; `jerr.multiview` conflates code ranges | Line 1090: `return(list(jerr = wls_fit$jerr))` unchanged | Warn at the `elnet.fit` boundary; distinct messages per code range |
| **F-P6-1** | P1 | `min_fdr_` anti-conservative 1.5–2.2× on pure nulls (realized FDP 0.918/0.934 vs reported 0.423/0.662; 18–33 features selected on no-signal data) | Documentation-first remediation not yet performed; threshold rule unchanged (parity contract) | Documentation: quantify null behavior in `?stabl_fit`/`?get_fdrs`/vignette; report the F-001 bound alongside `min_fdr_`; consider a null-calibration diagnostic. **Do NOT silently change the threshold rule** |
| **F-001** | P1 (interpretive) | `min_fdr_` alpha semantics: min over a random q+ curve; only the probabilistic F-001 guarantee holds | Phase 6 addendum: bound holds empirically (within bound at every Δ, nearly tight at small Δ); interpretive framing not yet added to docs | Folded into the F-P6-1 documentation fix |
| **F-002** | P1 (unvalidated extension; matrix classes it P2-fence) | Environment-aware min/mean aggregation is novel and unvalidated | No validation study; no "experimental" labeling in help/pkgdown | Validation study (WS-2/WS-5 style) or explicit "experimental, unvalidated" labeling |

The remaining 67 open items (16 P2 minus F-002, 27 P3, 23 P4, F-009) are
robustness, hygiene, documentation, or convention items — they degrade
quality but none is a known correctness or integrity blocker. Two P2s
deserve explicit mention because they are statistical, not cosmetic:
F-P7-1 (the native env-stability panel rule `min_score >= 0.5` never
consults the FDP+ machinery and saturates under the null — remediation is
an FDP+-gated option or an explicit "descriptive filter" contract in the
code, plus a leniency warning on the rf auto-grid) and F-CV-01..03 (RNG
hygiene: forked-worker seed sharing, global `.Random.seed` clobbering,
seed reuse across folds/omics).

---

## 5. Residual-risk statement

The residual risk of releasing stablr at HEAD `6aa5a6e`, stated by category:

1. **Vendored-engine memory/NA/convergence risks (F-011, F-012, F-013).**
   The vendored glmnetpp-based cooperative engine can abort the R process
   (F-011), silently drop columns (F-012), or swallow non-convergence
   (F-013) on inputs that are reachable but rare in the package's own
   benchmarks and vignettes. These are latent, input-dependent hazards —
   the audit reproduced F-011's SIGABRT and F-012/F-013's behavior
   independently, so they are not theoretical. Until PR-1 lands, any user
   passing sparse matrices with missing `xm`/`xs` attributes, NA-bearing
   matrices, or hitting maxit limits is exposed.
2. **Statistical-calibration risks (F-P6-1, F-P6-2, F-001, F-002, F-P7-1).**
   On pure-null data the default threshold rule selects 18–33 features with
   reported `min_fdr_` 0.42–0.66 against realized FDP ≈ 0.92–0.93. The
   F-001 probabilistic bound itself holds (Phase 6 addendum), but the
   point estimate is winner's-curse-biased (per-threshold mean 1.022 →
   min-over-100-thresholds 0.674) and the native env-stability panel rule
   ignores the decoy contrast entirely. These are properties of the
   published method's default rule, not implementation bugs; the parity
   contract forbids silently changing them, so the residual risk after
   documentation is **user misinterpretation**, which only docs and
   diagnostics can mitigate.
3. **Reproducibility risks (F-CV-01..03, P4-F28).** With
   `random_state=NULL` and parallel CV workers, results can depend on
   worker count; `stabl_fit` clobbers global RNG state; one seed is reused
   across folds/omics/cooperative foldid. None affects point estimates
   under fixed seeds; all affect bit-reproducibility guarantees the
   release contract implies.
4. **Audit-methodology limits.** (a) This entire program is an AI audit —
   it does **not** satisfy SCI-04; independent human confirmation remains
   outstanding and cannot be self-granted. (b) The Phase 7 S1 rests on a
   5-of-20-replicate pilot. (c) The v1 pipeline's donor-level derivation
   script was never archived (vignette §2.3); that step cannot be audited,
   only superseded by the v2 rebuild. (d) The coverage ledger's
   token-matcher produced false "no tests" leads (all refuted in Phase 9),
   so coverage claims rest on the corrected manual ledger. (e) L3's
   knockpy cross-comparison remains blocked (not installed); the MVR
   solver is verified by property tests only.
5. **Unremediated tail.** 67 non-blocking open items (P2–P4) are deferred
   with rationale in `remediation-matrix.md`; the four behaviour-design
   P3s (B-04, C-4, D-2, D-4) require user-facing API or convention
   decisions and are documented as known issues, not silent.

---

## 6. What the audit established positively (carried forward)

- **No data leakage** in any CV/train-validation path (Phase 3).
- **Multi-omic alignment is correct** (Phase 4); no P0/P1.
- **Core FDP+ implementation matches the specification** (Phase 2):
  strict-exceedance counting, tie-breaking, empty-selection fail-closed
  behavior, and input contracts all verified; MVR solver numerically sound.
- **The aging case study's numbers are reproducible from committed
  artifacts** (Phase 7): V8/V12/V13 and vignette §4/§10.3/§10.5 re-derived
  exactly; the v2 pipeline reproduces bit-for-bit (submodule pin rationale).
- **The test suite is real**: 3325 baseline assertions, 3067 in-check PASS,
  and every ledger "no tests" flag was a matcher false lead (Phase 9).
- **Distribution hygiene is clean** (Phases 10–11): dependency split
  correct, GPL-2 licensing consistent with the vendored engine, NEWS
  preserves adverse results verbatim.
- **The API surface is fully audited** (FR): 74/74 exports; defaults
  calibrated against the package's own benchmark evidence; the one P1
  found (B-01) fixed with a regression test.

---

## 7. Phase 13 remediation plan (merged sequencing)

Merges the remediation-matrix sequencing proposal with the 27 FR-deferred
items. Every remediation must preserve the Python-parity contract or
document the deliberate divergence (no silent workarounds). Frozen benchmark
protocols are never retroactively edited.

| PR | Content | Items | Est. |
|---|---|---|---|
| **PR-1 (blocking)** | P0+P1 code fixes + regression tests: F-011 (`xm`/`xs` length/positivity validation → `jerr=10001`), F-012 (`anyNA` guard), F-013 (boundary warning + distinct jerr messages). Compiled-boundary tests (Phase 9 gap) | F-011, F-012, F-013 | S |
| **PR-2** | P2 RNG-hygiene cluster: per-worker seed derivation when `random_state=NULL`, save/restore `.Random.seed` in `stabl_fit`, independent streams per fold/omic/cooperative foldid | F-CV-01, F-CV-02, F-CV-03/P4-F28 | S |
| **PR-3** | P2 documentation/fence cluster: null-behavior quantification + F-001 bound in `?stabl_fit`/`?get_fdrs`/vignette; regime-dependence note; env-stability panel rule contract (FDP+-gated option or "descriptive filter" statement) + rf auto-grid leniency warning; RP warning extended to all learners; F-002 "experimental" labeling; F-CV-04 warning; CP-convention statement (F-P7-2); WS-5 scope sentence in vignette §3.1 (F-P7-3); L2 approximation note; F-010 lambda-grid provenance | F-P6-1 (+F-001), F-P6-2, F-P7-1 (doc/warning half), F-005/L1, F-002, F-CV-04/P4-F24, F-P7-2, F-P7-3, L2, F-010 | M |
| **PR-4** | P3 batch: compiled-code hygiene (F-014 copy-guards, F-015 dev.ratio rename/fix, F-016 `warm_fit$a`, F-017 xs zero-guard, F-018 vendoring prune + version pin), seed/error-class/message items (F-019..F-021), doc items (F-022, F-023, F-P10-2), ordinal/RF validation items (F-024..F-026), nested-CV seeds (F-CV-05..08), fused-lambda/hard-threshold/mse items (P4-F25..F27), FR behaviour-design P3s (B-04 family param, C-4 stratified folds, D-2 task_type cross-check, D-4 diagonal convention — each needs a decision, flag in PR description) | 27 items | M |
| **PR-5** | P4 batch: validation guards, doc gaps, hygiene (incl. FR-hygiene `vignettes/figure/` removal → clears the 2nd check NOTE) | 23 items | S–M |
| **Owner action** | Re-ratify the release contract at 0.3.0 or amend the contract version (F-P10-1); decide F-009 validation-task execution (direct augmented-matrix solve check) | F-P10-1, F-009 | M |

**Open question for Phase 13 kickoff:** the matrix batches F-014/F-015/F-017
with the P3 hygiene list although `findings.md` labels them P2 (F-017 is
SUSPECT). Either confirm the re-batch rationale in the matrix or promote
them to PR-3/PR-4 scope explicitly.

**Sequencing gates:** PR-1 → PR-2 → PR-3 must land before the fresh release
run; PR-4/PR-5 may land in parallel; the release run and SCI-04 follow only
after PR-1..PR-3.

---

## 8. Release gates (human tasks, unchanged)

1. `fresh-release-run-pending` — fresh n=150 release run after PR-1..PR-3.
2. `sci04-independent-confirmation-pending` — independent human statistician
   confirmation. **The AI audit does not satisfy SCI-04 and must not be
   represented as doing so.**
3. Phase 15 referee report — final program deliverable, after remediation.

---

## 9. Program-level publication-readiness verdict

**Not release-ready at HEAD `6aa5a6e`.** The blocking set is concrete and
small: one P0 (F-011, vendored-engine OOB), two P1 code defects (F-012,
F-013), one P1 calibration-documentation obligation (F-P6-1 with F-001),
and one P1 unvalidated-extension obligation (F-002) — all verified open at
HEAD — plus the release-contract re-ratification (F-P10-1) and the two
standing human gates.

What the audit program has established is substantial and favorable: the
statistical core matches its specification, there is no leakage, the
multi-omic machinery aligns, the case-study numbers reproduce bit-for-bit
from committed artifacts, the test suite is genuine and green, the
documentation and distribution hygiene are clean, and the full 74-export
API surface has been per-function audited with every finding either fixed
or explicitly deferred with rationale. The remaining work to release is
execution of the §7 plan (PR-1..PR-3 are small and fully specified), the
fresh release run, and independent human confirmation (SCI-04) — the last
of which no part of this program can substitute.
