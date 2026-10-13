# SCI-04 Independent Review Packet — null selected-fraction gate v2

Packet status: **drafted for external review, 2026-10-13**, at stablr commit
`9c6d8fc8083a9e33a212a618d6741686880b44e6` (branch `improvements/0.3.0`).
Decision record under review: `SCI-04-null-selected-fraction-v2-draft`
(status `proposed_unaccepted`).

This packet is the review entry point for the decision record
`inst/analysis/contracts/SCI-04-null-selected-fraction-v2.md`. It inventories
the locked artifacts, states what the reviewer is asked to confirm or replace,
and defines what may happen only after acceptance. Nothing in this packet
authorizes an implementation or a release.

---

## 1. What the reviewer is asked to decide

Per the acceptance boundary of the v2 record, the reviewer must:

1. confirm the **estimand, independence assumptions, formula, confidence
   level, and multiplicity policy** (Sections 3–4);
2. **approve or replace** the exact record
   `SCI-04-null-selected-fraction-v2-draft` in a **new immutable decision
   record** (the draft itself is never edited into acceptance);
3. record **reviewer identity, decision date, and the accepted record hash**;
4. approve the **deterministic reference tests independently of the failed
   run** (Section 6).

Until all four are complete, `null_selected_fraction` has no numeric release
bound, carries status `proposed_unaccepted`, and fails every promotable
release packet. This is enforced in code, not only by convention:
`inst/analysis/release_evidence.R` errors on any methodology validation whose
SCI-04 gate is not exactly
`gate_version = "SCI-04-null-selected-fraction-v2-draft"`,
`decision_status = "proposed_unaccepted"`.

## 2. Background: why v1 was rejected and sealed

The historical v1 implementation passed `sum(X_i)` — a sum of continuous
fractions — to a binomial Wilson interval. Continuous fractions are not
Bernoulli successes; the resulting "bound" is not a valid confidence bound
for any estimand of this gate. The v1 result is **sealed as diagnostic
evidence only**:

- Sealed packet: RE-01 bounded-behavior baseline v1, normalized from commit
  `c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1` under the pinned R4_51 runtime.
- In-repo location: `tests/testthat/fixtures/release_evidence_baseline/v1/`
  (sha256 manifest: `manifest.sha256`, file
  `2a7d2025d952fe90053194c58fd06a1ef4a1c787783cd43e6f13b2fe4823cf4b`).
- The fixture README states the sealing contract: these files "must never be
  assembled into a Release Evidence Packet".

The v1 numbers (e.g. `null_selected_fraction` bounds of 0.575–0.941 in the
sealed gates table) are artifacts of the invalid v1 statistic and carry no
decision weight. They are retained solely so the reviewer can audit what the
v1 implementation did.

## 3. The proposed decision (v2 record, verbatim substance)

**Estimand.** The gate concerns the mean fraction of selected features under
each null simulation cell. One Monte Carlo replicate contributes one
observation `X_i = n_selected_i / p_i`, with `0 <= X_i <= 1`; the release
quantity is a one-sided upper confidence bound on the replicate-level mean
`E[X]`.

**Proposed bound.** The distribution-free one-sided Hoeffding upper
confidence bound:

```
U = min(1, mean(X) + sqrt(log(1 / alpha) / (2 * n)))
```

**Fixed settings.** `alpha = 0.05`; release criterion `U <= 0.10`;
experimental unit = the independently seeded Monte Carlo replicate.
Missing, duplicated, non-finite, or out-of-range replicate fractions make
the cell incomplete and must fail closed.

**Predeclared replicate count.** `n >= 150`. The feasibility clause (audit
E1, 0.2.2) shows the criterion is achievable only if the margin alone is
below 0.10, since even a perfect run (`mean(X) = 0`) yields `U = margin(n)`:

| replicates `n` | margin | `U <= 0.10` achievable |
|---|---|---|
| 100 | 0.1224 | no (vacuous) |
| 150 | 0.0999 | yes (requires `mean(X) <= 0.0001`) |
| 200 | 0.0865 | yes (requires `mean(X) <= 0.0135`) |
| 250 | 0.0774 | yes (requires `mean(X) <= 0.0226`) |

The previously considered `n = 100` makes the criterion impossible to
satisfy. At the minimum-feasible `n = 150` the observed mean must be near
zero; the accepted record may raise `n` or renegotiate the 0.10 criterion
explicitly, but neither may be tuned against the failed locked run.

## 4. Assumptions the reviewer must explicitly confirm or replace

1. **Independence.** Replicates are independently seeded Monte Carlo runs
   (contract root seed 270627 with per-replicate seed derivation); no
   randomness is shared across replicates. Hoeffding's inequality requires
   independent (not merely exchangeable) bounded variables.
2. **Boundedness.** `X_i in [0, 1]` by construction (`n_selected_i <= p_i`);
   the bound is valid for [0,1]-bounded variables without distributional
   assumptions.
3. **Same-size replicates.** The bound as written assumes a common `p_i`
   within a cell (so `X_i` are comparable draws of one variable). Cells that
   vary `p` across replicates would need a per-replicate normalization
   review before the bound applies.
4. **Confidence level.** One-sided `alpha = 0.05` (95% upper bound).
5. **Multiplicity policy.** SCI-04 is one of seven scientific gates in the
   methodology validation
   (`cell_completeness`, `python_metrics_parity`, `null_select_any`,
   `null_selected_fraction`, `null_select_90pct`, `signal_mean_fdp`,
   `signal_tpr`), each with its own predeclared fixed criterion and no
   across-gate alpha adjustment. The reviewer is asked to confirm this
   per-gate policy or specify a replacement in the accepted record.
6. **One-sidedness.** Only the upper bound is used; low `mean(X)` is never
   penalized.

## 5. No-tuning statement

- The v2 proposal was **selected from the sampling model before any
  corrected bound was computed** against the failed run.
- The **failed locked run** (the run evaluated under the invalid v1 Wilson
  statistic) is sealed in the RE-01 baseline; no number from it entered the
  v2 derivation.
- The **feasibility clause** (Section 3 table) is derived from the formula
  and the predeclared settings alone; the deterministic reference tests
  pin that arithmetic without reading any locked-run data.
- The reviewer must approve the reference tests **independently of the
  failed run**, and any post-acceptance change to `alpha`, the criterion,
  or `n` requires a new immutable decision record — never an edit to the
  draft.

## 6. Reference tests (reviewer approval item 4)

`tests/testthat/test-sci04-gate-feasibility.R` pins, deterministically and
without locked-run data:

- the margin formula against the contract table
  (`margin(100) = 0.1224`, `margin(150) = 0.0999`, `margin(200) = 0.0865`,
  `margin(250) = 0.0774`, tolerance 1e-3);
- vacuity at `n = 100` (`U` of a perfect run exceeds 0.10);
- the feasibility boundary (`margin(149) > 0.10`, `margin(150) <= 0.10`);
- the presence of the feasibility clause in the contract file.

Run:

```sh
R_LIBS=<library-path> Rscript -e 'testthat::test_local(".", filter="sci04-gate-feasibility")'
```

## 7. Locked artifact inventory (sha256, at packet commit `9c6d8fc`)

| Artifact | Role | sha256 |
|---|---|---|
| `inst/analysis/contracts/SCI-04-null-selected-fraction-v2.md` | decision record under review | `ed977c503ada334b9f4467d8663073b36abc206e58c3ce0fd72ef6be6aabaf9d` |
| `inst/release/release-contract.json` | locked release contract (`promotable: false`; blocker `SCI-04-null-selected-fraction-v2-pending-independent-acceptance`) | `db2843fc9f6c3c7925bf8294436c4561191e9a12845699120286b913d572b386` |
| `tests/testthat/test-sci04-gate-feasibility.R` | deterministic reference tests | `fc7636cd1458fccffb5429a38f7b470140601d1189b61ede9e0039e6f160612b` |
| `inst/analysis/release_evidence.R` | gate evaluation + fail-closed enforcement | `59cb6da5e6b3f3e31e51bcffbffff80529d47b42766148110e5077d756df1c45` |
| `tests/testthat/fixtures/release_evidence_baseline/v1/` (11 files) | sealed RE-01 v1 baseline, diagnostic only | per-file hashes in `manifest.sha256` (manifest hash `2a7d2025…3cf4b`) |

The reviewer should verify these hashes before reviewing; any mismatch means
the artifact changed after this packet was drafted and the packet must be
regenerated.

## 8. What may happen only after acceptance

1. The reviewer (or maintainer, recording the reviewer's decision) writes a
   **new immutable decision record** containing: the accepted bound and
   settings, reviewer identity, decision date, and the sha256 of the
   accepted record.
2. Only then may `inst/release/release-contract.json` be updated:
   `null_selected_fraction` gate version → the accepted record id,
   `decision_status` → `accepted`, the blocker removed, and the methodology
   validation `settings.replicates` raised to the accepted `n >= 150`
   (the locked contract currently pins `replicates: 100`, which the v2
   record declares infeasible — this update is part of the acceptance
   flip, not something that may be done before it).
3. `release_evidence.R`'s fail-closed check (SCI-04 must be
   `proposed_unaccepted` under the draft id) is replaced by a check for the
   accepted record id and status, in the same commit or series as the
   contract flip.
4. A fresh release run at the accepted `n` is executed and a new Release
   Evidence Packet assembled; the sealed RE-01 v1 fixtures remain excluded
   from any packet.

## 9. Reviewer response template

```
Reviewer name / affiliation:
Decision (approve / approve-with-amendments / replace):
Estimand confirmed (Y/N, notes):
Independence + boundedness assumptions confirmed (Y/N, notes):
Formula and alpha = 0.05 confirmed (Y/N, notes):
Multiplicity policy confirmed or replaced (notes):
Reference tests approved independently of the failed run (Y/N):
Accepted record id:
Accepted record sha256:
Decision date:
Signature / verification method:
```
