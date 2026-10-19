# SCI-04 Independent Confirmation Packet — null selected-fraction gate v2 (regenerated 2026-10-19)

Packet status: **drafted for independent confirmation, 2026-10-19**, at stablr
commit `2899a2a` (branch `main`), after the Phase 13 remediation program
(PR-1/PR-2, PR-A through PR-G, and the vignette PR).

This packet **regenerates** the 2026-10-16 confirmation packet
(`SCI-04-independent-confirmation-packet.md`, drafted at commit `6cd041f`)
under that packet's own rule: "any mismatch means the artifact changed after
this packet was drafted and the packet must be regenerated." Two Section 4
artifacts changed after `6cd041f` (Section 4a below); all gate arithmetic,
contract pins, reference tests, and sealed fixtures are unchanged and
verified at this packet's commit. The 2026-10-16 packet remains on disk
unedited; this regeneration supersedes it as the review entry point.

Object under confirmation: the **provisional owner acceptance** recorded in
`inst/analysis/contracts/SCI-04-null-selected-fraction-v2-accepted.md`
(current sha256 `fd55432bdfc290e6b375ce7e88667ef48104162a6e14ce36d1ff616e77092d7e`),
which provisionally accepts the draft decision record
`SCI-04-null-selected-fraction-v2-draft`
(sha256 `ed977c503ada334b9f4467d8663073b36abc206e58c3ce0fd72ef6be6aabaf9d`).

The release contract pins `sci04-independent-confirmation-pending` until this
confirmation completes. Nothing in this packet authorizes promotion; the
contract also pins `fresh-release-run-pending` and `promotable` stays
`false`.

---

## 0. Why this packet exists

The draft record's acceptance boundary requires an independent reviewer to
(1) confirm the estimand, independence assumptions, formula, confidence
level, and multiplicity policy; (2) approve or replace the exact record in a
new immutable decision record; (3) record reviewer identity, decision date,
and the accepted record hash; and (4) approve the deterministic reference
tests independently of the failed run. On 2026-10-16 the project owner — the
package author — recorded a blanket acceptance without itemized review
notes. That decision is recorded verbatim and is not rewritten here; it does
not, however, satisfy the independence requirement as written. This packet
gives an eligible independent reviewer everything needed to confirm or
replace the acceptance item by item, from committed public artifacts alone.

## 1. Reviewer eligibility (independence)

The confirming reviewer must attest that they:

- are **not** the package author and have **no commits** in the repository
  history;
- did **not** participate in the working session that produced the
  provisional acceptance, and had no role in drafting the draft record, the
  provisional acceptance, or the sealed RE-01 v1 baseline;
- have **no co-authorship or current collaborative relationship** with the
  package author within the last 24 months;
- are **not** producing competing software whose value depends on this
  gate's outcome.

Co-authorship of the *data* publication (Terekhova et al. 2026, the aging
atlas) is a relationship to **disclose**, not a disqualifier: the gate
concerns the statistical method under null simulation, not the data. All
relationships are disclosed in the response template and recorded in the
confirmation record.

## 2. Itemized confirmation checklist

Each item gets a verdict: **confirm / amend / replace**, with notes. File
anchors are given as grep-able strings; line numbers refer to this packet's
commit (`2899a2a`) and may drift.

- **C1 — Estimand.** The gate quantity is the replicate-level mean selected
  fraction `E[X]`, one observation `X_i = n_selected_i / p_i` per
  independently seeded Monte Carlo replicate, `0 <= X_i <= 1`, with a
  one-sided upper confidence bound. Verify the implementation computes
  exactly this quantity per null cell
  (`run_methodology_validation.R`, anchor `x <- cell$n_selected / cell$p`,
  line 1012 at this commit).
- **C2 — Independence.** Replicates are independently seeded: root seed
  270627 with disjoint per-replicate derivation
  (`data_seed <- seed + family_i * 1e6 + scenario_i * 1e5 + replicate * 1e3`,
  `fit_seed <- data_seed + 1L`; anchor `data_seed <- seed +`, line 1342).
  Hoeffding's inequality requires independent (not merely exchangeable)
  bounded variables. Confirm the derivation shares no randomness across
  replicates.
- **C3 — Boundedness.** `X_i in [0, 1]` by construction
  (`n_selected_i <= p_i`); the bound is valid for [0,1]-bounded variables
  without distributional assumptions. Verify the implementation clamps at 1
  (`min(1, mean(x) + ...)`) and never at the lower end (one-sided upper).
- **C4 — Formula, confidence level, reference arithmetic.** Accepted bound:
  `U = min(1, mean(X) + sqrt(log(1 / alpha) / (2 * n)))`, `alpha = 0.05`,
  criterion `U <= 0.10`. Verify the implementation line matches the formula
  exactly (anchor `sqrt(log(1 / 0.05)`, line 1013), and run the
  deterministic reference tests (Section 5) pinning `margin(100) = 0.1224`,
  `margin(150) = 0.0999`, `margin(200) = 0.0865`, `margin(250) = 0.0774`
  (tolerance 1e-3), vacuity at `n = 100`, and the feasibility boundary
  (`margin(149) > 0.10 >= margin(150)`).
- **C5 — Multiplicity policy.** Seven methodology gates and four late-fusion
  gates, each with its own predeclared fixed criterion and **no across-gate
  alpha adjustment**. Confirm this per-gate policy or specify a replacement
  in the confirmation record.
- **C6 — One-sidedness.** Only the upper bound is used; a low `mean(X)` is
  never penalized. Confirm.
- **C7 — Fail-closed semantics.** Missing, duplicated, non-finite, or
  out-of-range replicate fractions make the cell incomplete and the gate
  fails (`run_methodology_validation.R`, anchor
  `.gate_cell_status <- function`, line 945; the pass term
  `is.finite(selected_fraction) && selected_fraction <= 0.10` cannot be
  true for an incomplete cell). Verify, including that the normalization
  layer forces `pass` off for any non-accepted decision status
  (`release_evidence.R`, anchor `gate decision status is`, line 2105).
- **C8 — No-tuning and sealed fixtures.** The v2 proposal was selected from
  the sampling model before any corrected bound was computed against the
  failed locked run; the sealed RE-01 v1 baseline
  (`tests/testthat/fixtures/release_evidence_baseline/v1/`, fixture README
  sha256 `2a7d2025d952fe90053194c58fd06a1ef4a1c787783cd43e6f13b2fe4823cf4b`)
  is diagnostic only and excluded from any Release Evidence Packet. Verify
  the fixture README sealing statement and that no gate evaluation reads
  locked-run data. (Fixture integrity re-verified at this commit:
  `sha256sum -c manifest.sha256` reports OK for all entries.)
- **C9 — Reference tests approved independently of the failed run.** Execute
  the tests yourself (Section 5) and record the outcome. They pin the
  arithmetic and the presence of the feasibility clause without reading any
  locked-run data. (Maintainer-run outcome at this commit, for orientation
  only and not a substitute for the reviewer's own run: 30 expectations,
  0 failed / 0 errors / 0 warnings / 0 skipped, 2026-10-19.)
- **C10 — Provisional acceptance chain.** Verify: (a) the draft record hash
  quoted inside the accepted record matches the actual draft file hash;
  (b) the reviewer response
  (`inst/analysis/review/SCI-04-acceptance-record.md`) records the owner's
  verbatim statement, identity, and date, and discloses the blanket nature
  and the independence gap; (c) the accepted record's current hash matches
  Section 4 of this packet and the amended hash line in the acceptance
  record (the acceptance record carries both the at-acceptance hash and the
  current hash after the 2026-10-16 provisional-status amendment);
  (d) the release contract's `null_selected_fraction` entry carries
  `gate_version = "SCI-04-null-selected-fraction-v2-accepted"`,
  `decision_status = "accepted"`, and both blockers.
- **C11 — Operating point.** At the pinned `n = 150` the margin alone is
  0.0999, so the criterion `U <= 0.10` requires `mean(X) <= 0.0001` — a
  near-degenerate, deliberately conservative operating point. The reviewer
  may: confirm `n = 150` as proposed; raise `n` (margin 0.0865 at `n = 200`,
  requiring `mean(X) <= 0.0135`; 0.0774 at `n = 250`, requiring
  `mean(X) <= 0.0226`); or renegotiate the 0.10 criterion. Any change is
  made in the new immutable record with justification and must **never** be
  tuned against the failed locked run.

## 3. Decision semantics

- **Confirm.** Write a new immutable decision record
  `SCI-04-null-selected-fraction-v2-confirmed` (status `accepted`,
  superseding the provisional accepted record, which stays on disk
  unedited) containing: the confirmed bound and settings (identical to the
  accepted record unless C11 amends them), the itemized C1–C11 verdicts,
  the sha256 of the accepted record as confirmed, reviewer identity,
  decision date, and eligibility attestation.
- **Confirm-with-amendments / replace.** Same, with the amended settings
  and justification; the deterministic reference tests must still pass or
  be updated in the same commit with the amended arithmetic (never tuned
  against the locked run).
- The draft record and the provisional accepted record are **never edited**.

## 4. Locked artifact inventory (sha256, at packet commit `2899a2a`)

| Artifact | Role | sha256 |
|---|---|---|
| `inst/analysis/contracts/SCI-04-null-selected-fraction-v2.md` | draft decision record (sealed, superseded) | `ed977c503ada334b9f4467d8663073b36abc206e58c3ce0fd72ef6be6aabaf9d` |
| `inst/analysis/contracts/SCI-04-null-selected-fraction-v2-accepted.md` | provisional owner acceptance — object under confirmation | `fd55432bdfc290e6b375ce7e88667ef48104162a6e14ce36d1ff616e77092d7e` |
| `inst/analysis/review/SCI-04-acceptance-record.md` | owner's reviewer response (verbatim statement, independence disclosure, amended hash line) | `9444739c1cd19d64345553cff43bbb8d5c9ad905980be031b74ab90bbd352c4a` |
| `inst/release/release-contract.json` | release contract (both blockers pinned, `replicates: 150`) | `1b25f043920f59f5849e8001f1f94a4cea23f9d772d5ffd879386e0b7b6c9588` |
| `inst/analysis/release_evidence.R` | contract enforcement + fail-closed normalization | `f7bbf7d36400f870ee63362e844fce36d76e94a76ae28443f9a4f6d89417abf4` |
| `inst/analysis/run_methodology_validation.R` | gate implementation (Hoeffding evaluation) | `9abf6f096ce9e24000748340af40be50d8b07a56ccaccb244493a3a77dda63eb` |
| `tests/testthat/test-sci04-gate-feasibility.R` | deterministic reference tests | `fc7636cd1458fccffb5429a38f7b470140601d1189b61ede9e0039e6f160612b` |
| `tests/testthat/test-methodology-validation-runner.R` | runner pins incl. accepted-record assertions | `ee46978dfd909ed8d10dcbaffab05ff46db38d790a4381a9f8110f2332d5e63f` |
| `tests/testthat/test-release-evidence-contract.R` | contract pins incl. blockers | `f6b99d2c6d594ef90ced56a4b2c13a5dd3cc81daa02d083eca192ea95ce6ec79` |
| `tests/testthat/test-release-evidence-packets.R` | normalization fail-closed pins | `3b3507e66f8f28bb8386c6ae367af739aacd12c0ca51bf7227116cf210f1127d` |
| `tests/testthat/fixtures/release_evidence_baseline/v1/` (11 files) | sealed RE-01 v1 baseline, diagnostic only | per-file hashes in `manifest.sha256` (fixture README `2a7d2025…3cf4b`) |

Verify these hashes before reviewing; any mismatch means the artifact
changed after this packet was drafted and the packet must be regenerated.

### 4a. Delta against the 2026-10-16 packet (commit `6cd041f`)

Two artifacts changed after the previous packet was drafted; both changes
are behavior-preserving with respect to this gate, and every checklist
anchor was re-verified at this commit:

1. `inst/analysis/run_methodology_validation.R`
   (`5a097065…720218` -> `9abf6f09…a63eb`): Phase 13 PR-G migrated the
   script's similarity-metric calls from removed pre-1.0.0 aliases to the
   canonical 1.0.0 names (`jaccard_similarity_matrix()`,
   `adjusted_similarity_pairwise()`, `pearson_similarity_pairwise()`,
   `adjusted_similarity_summary()`, `pearson_similarity_summary()`,
   `get_feature_names()`; internal benchmark metrics now use `stablr:::`).
   The aliases were thin wrappers, so computed values are unchanged; the
   Python-parity math remains pinned by the frozen fixtures, and the gate
   arithmetic lines (C1/C2/C4/C7 anchors) are byte-identical. Full suite
   after the migration: 62 files, 0 failed / 0 errors.
2. `inst/analysis/review/SCI-04-acceptance-record.md`
   (`be696e53…138ba` -> `9444739c…52c4a`): amended in commit `5cf603b`
   (PR #21) to record the accepted record's hash after the 2026-10-16
   provisional-status amendment, alongside the at-acceptance hash — exactly
   the amended-hash-line disclosure item C10(c) checks. No other content
   changed.

No other Section 4 artifact changed. The release contract
(`1b25f043…6c9588`), the enforcement layer, the reference tests, and the
sealed fixtures are untouched since the previous packet.

## 5. Verification commands

```sh
# Artifact hashes (Section 4)
sha256sum inst/analysis/contracts/SCI-04-null-selected-fraction-v2*.md \
          inst/analysis/review/SCI-04-acceptance-record.md \
          inst/release/release-contract.json inst/analysis/release_evidence.R \
          inst/analysis/run_methodology_validation.R tests/testthat/test-*.R

# Sealed-fixture integrity (C8)
(cd tests/testthat/fixtures/release_evidence_baseline/v1 && sha256sum -c manifest.sha256)

# Deterministic reference tests + governance pins (no locked-run data)
R_LIBS=<library-path> Rscript -e 'testthat::test_local(".", filter="sci04-gate-feasibility|release-evidence-contract|release-evidence-packets|methodology-validation-runner", stop_on_failure=FALSE, reporter="summary")'

# Margin arithmetic spot-check
Rscript -e 'for (n in c(100, 149, 150, 200, 250)) cat(n, sprintf("%.4f", sqrt(log(1/0.05)/(2*n))), "\n")'
```

## 6. What may happen only after confirmation

1. The maintainer records the reviewer's confirmation (new immutable record
   per Section 3) and flips the release contract in the same commit or
   series: `null_selected_fraction` gate version -> the confirmed record id,
   `decision_status` stays `accepted`, and the
   `sci04-independent-confirmation-pending` blocker is removed.
2. `release_evidence.R`'s accepted-id check and the runner's
   `.selected_fraction_gate_version` are updated to the confirmed id in the
   same series; the fail-closed override for non-accepted statuses stays.
3. The fresh full-candidate release run at the confirmed `n` is executed and
   a new Release Evidence Packet assembled; `fresh-release-run-pending` is
   removed only then. The sealed RE-01 v1 fixtures remain excluded from any
   packet.
4. Any post-confirmation change to `alpha`, the criterion, or `n` requires a
   new immutable decision record — never an edit to an existing record.

## 7. Reviewer response template

```
Reviewer name / affiliation / role:
Eligibility attestation (Section 1) + disclosures:
C1 estimand (confirm/amend/replace + notes):
C2 independence of replicates (…):
C3 boundedness and clamping (…):
C4 formula, alpha = 0.05, reference arithmetic (…):
C5 multiplicity policy (confirm or replacement):
C6 one-sidedness (…):
C7 fail-closed semantics (…):
C8 no-tuning and sealed fixtures (…):
C9 reference tests independently executed (Y/N + outcome):
C10 provisional acceptance chain verified (Y/N):
C11 operating point (confirm n = 150 / raise n / renegotiate criterion + justification):
Decision (confirm / confirm-with-amendments / replace):
Confirmed (or replacement) record id:
Confirmed record sha256:
Decision date:
Signature / verification method:
```

## 8. Practical notes

- Everything needed is in the public repository at the packet commit. No
  access to the sealed locked-run data, the aging-atlas data, or any
  simulation re-run is required: this is a protocol and implementation
  confirmation, not a data re-analysis.
- Expected effort: the itemized pass is a reading and test-running exercise
  (roughly one to two hours).
- Suitable reviewer profile: a statistician familiar with concentration
  inequalities and Monte Carlo release gating, independent of the package
  author under Section 1.

## 9. Phase 13 context (orientation, not part of the C1–C11 checklist)

Between the previous packet (`6cd041f`) and this one (`2899a2a`), the
repository completed the Phase 13 remediation program. None of it changes
the gate under confirmation; all of it is auditable in the committed
adversarial-audit trail at `docs/audit/2026-10-17-adversarial/`:

- `findings.md` — the full findings ledger (finding IDs F-001…, F-CV-…,
  P4-…, A1–A6, N-01…, B-01…, C-1…, D-1…), including the per-PR Phase 13
  ledger blocks with gate evidence (suite counts, R CMD check status) for
  PR-1, PR-2, and PR-A through PR-G.
- `README.md` — phase index with the Phase 13 row recording each PR's
  scope and gates.
- `final-review.md` — the pre-publication per-function review (FR) of all
  then-exported functions.
- `phase12-synthesis.md` — cross-phase synthesis, residual risk, and the
  merged Phase 13 plan this program executed.
- `remediation-matrix.md`, `specification-matrix.md`,
  `coverage-ledger.csv` — finding-to-fix and spec-to-implementation
  traceability.
- `baseline-2026-10-17.md`, `phase6_validation/`, `phase7_validation/` —
  the independent synthetic-DGP validation phases and their artifacts.

Program summary for orientation: PR-1/PR-2 closed the release-blocking
fail-closed and RNG-hygiene findings; PR-A closed the 1.0.0 contract
residuals; PR-B shipped the Python-parity and plot-reading vignettes;
PR-C/D/E added the reliability toolchain (`stabl_diagnose()`,
`stabl_null_check()`, `stabl_replication_profile()`); PR-F added the
gaussian validated path and portable final models (`stabl_finalize()`,
`predict.stabl_predictor()`); PR-G cut the public API to the lean 1.0.0
surface (54 exports / 43 S3 methods, migration guide
`docs/MIGRATION_0.3_TO_1.0.md`, removal ledger
`inst/release/removal-ledger.csv`) and retired the experimental oblique
random forest to `attic/`; the vignette PR added the `stablr-reliability`
article. DESCRIPTION remains `0.3.0.9000`; the 1.0.0 version bump is part
of release certification, which this confirmation gates.
