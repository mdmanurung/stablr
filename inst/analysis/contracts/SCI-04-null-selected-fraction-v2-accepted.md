# SCI-04 null selected-fraction gate v2 — accepted decision record

Status: accepted

Decision ID: `SCI-04-null-selected-fraction-v2-accepted`

Supersedes: `SCI-04-null-selected-fraction-v2-draft`
(sha256 `ed977c503ada334b9f4467d8663073b36abc206e58c3ce0fd72ef6be6aabaf9d`),
which remains sealed as drafted and is never edited into acceptance.

## Accepted decision

The bound, confidence level, release criterion, and experimental unit are
accepted exactly as proposed in the superseded draft record:

`U = min(1, mean(X) + sqrt(log(1 / alpha) / (2 * n)))`

with `alpha = 0.05`, release criterion `U <= 0.10`, experimental unit = the
independently seeded Monte Carlo replicate, and the predeclared release
replicate count `n >= 150` (the release contract pins `n = 150`, the
minimum feasible count under the feasibility clause). Missing, duplicated,
non-finite, or out-of-range replicate fractions make the cell incomplete
and fail closed. No setting was tuned against the failed locked run.

## Independent review record

- Reviewer: Mikhael Dito Manurung (project owner, acting as the external
  reviewer required by the acceptance boundary).
- Decision: approve (blanket acceptance of the draft record as proposed).
- Decision date: 2026-10-16.
- Verification method: reviewer's acceptance statement recorded verbatim
  from the project owner during the working session:
  "as external reviewer i accept gate".
- Reviewer response template:
  `inst/analysis/review/SCI-04-acceptance-record.md`.

## Provisional status of this acceptance

The acceptance boundary of the superseded draft requires an **independent**
reviewer to item-by-item confirm the estimand, assumptions, formula,
multiplicity policy, and reference tests. The acceptance recorded here was
provided by the project owner — the package author — as a blanket statement,
with no itemized review notes. It therefore does **not** satisfy the
independence requirement of the review packet as written, and this record is
provisional: an itemized confirmation (or replacement) by a reviewer
independent of the package remains required before any promotable Release
Evidence Packet relies on this gate. The release contract carries
`sci04-independent-confirmation-pending` as a blocker until then.

## What this record authorizes

1. Implementation of the Hoeffding evaluation in
   `inst/analysis/run_methodology_validation.R`: the gate no longer fails
   closed by construction; it evaluates `U` per null cell and passes when
   `U <= 0.10`.
2. The release-contract flip: `null_selected_fraction` gate version ->
   `SCI-04-null-selected-fraction-v2-accepted`, `decision_status` ->
   `accepted`, the pending-acceptance blocker removed, and the methodology
   validation `settings.replicates` raised to 150.
3. `inst/analysis/release_evidence.R` enforces the accepted record id and
   status in place of the draft fail-closed check.

A fresh release run at the accepted replicate count and assembly of a new
Release Evidence Packet remain required before anything is promotable; the
sealed RE-01 v1 fixtures stay excluded from any packet. Any post-acceptance
change to `alpha`, the criterion, or `n` requires a new immutable decision
record.
