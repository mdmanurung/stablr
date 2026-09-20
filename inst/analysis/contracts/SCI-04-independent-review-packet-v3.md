# SCI-04 independent statistical review packet v3

Status: implementation authorized by package owner; independent decision pending

Packet ID: `SCI-04-independent-review-packet-v3`

Implementation authority: `mdm123`

This additive packet supersedes the requested decision content in v2. It does
not replace either historical packet and is not an independent acceptance.
The owner has authorized implementation so that the formula, applicability
rules, malformed-data handling, and deterministic reference values can be
reviewed as executable code.

## Proposed complete decision

For each applicable global-null cell and independently seeded Monte Carlo
replicate `i`, define `X_i = n_selected_i / p_i`. The estimand is `E[X_i]` for
that cell. Each `X_i` must be finite and in `[0, 1]`.

With exactly `n = 600` replicates, the one-sided fixed-sample Hoeffding upper
bound is

`U = min(1, mean(X) + sqrt(log(1 / 0.05) / (2 * 600)))`.

The radius is `0.049964422955689106`. A cell passes when `U <= 0.10`; equality
passes. The complete release decision is an intersection-union test: every
applicable null cell must pass. Per-cell `alpha = 0.05` is therefore unadjusted
for the release decision. This controls false promotion at no more than 0.05
when at least one required cell violates its null requirement. It does not
provide simultaneous 95% confidence coverage across all 44 null cells.

If simultaneous 95% reporting is claimed, it is a separate Bonferroni analysis
over the 44 applicable null cells and is not the release gate. At a target
radius of 0.05, that reporting analysis requires at least 1,356 replicates per
cell; the 600-replicate release run must not be described as simultaneous 95%
coverage.

## Applicability and failure handling

`SCI-04-applicability-v1.csv` declares all 96 cells. The eight combinations of
fixed-X `knockoff` with a high-dimensional scenario are `not_applicable` with
reason `fixed_x_requires_n_greater_than_p`. The remaining 88 cells are
applicable; 44 of those are global-null cells. An applicable cell cannot pass
using a fallback generator.

Missing cells, missing or duplicate replicate identifiers, a count other than
600, non-OK fits, generator mismatch, missing/non-finite/non-integer `p` or
`n_selected`, `p <= 0`, or `n_selected` outside `[0, p]` fail closed. Declared
not-applicable cells are retained in output with `not_applicable` status and
are excluded from gate summaries; they are never counted as passing cells.

Exact deterministic examples and invalid-input expectations are recorded in
`SCI-04-reference-cases-v1.json`. The normative machine-readable proposal is
`SCI-04-statistical-specification-v3.json`.

## Independent-review action

An independent reviewer must accept or reject this complete decision using
`sci-decision-v2.schema.json`. Acceptance must bind a fresh sealed review
bundle, this packet, the applicability matrix, the reference cases, and their
hashes through an authenticated detached signature and external reviewer trust
root. `mdm123` is recorded only as implementation authority and cannot satisfy
the independent-review role.

Until that signature validates, the active embedded release contract remains
v1, non-promotable, and blocked at SCI-04. No passing test result from the
owner-authorized implementation authorizes release.
