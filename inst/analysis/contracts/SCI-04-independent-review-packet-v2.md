# SCI-04 independent statistical review packet v2

Status: awaiting independent decision

Packet ID: `SCI-04-independent-review-packet-v2`

This packet supersedes the review request in v1 without replacing or editing
that historical packet. It requests a decision; it is not an accepted
scientific contract, and the executor cannot act as the independent reviewer.

## Declared validation matrix

The release-profile matrix contains 96 declared cells:

- four outcome families;
- three scenarios per family;
- four artificial-feature generators; and
- two dimensional regimes.

Eight fixed-X/high-dimensional cells are currently proposed as
`not_applicable`, because fixed-X construction requires conditions that the
declared high-dimensional design does not satisfy. The proposal therefore
leaves 88 applicable cells. Of those, 44 are global-null cells to which the
SCI-04 selected-fraction gate would apply. These counts are proposals for the
reviewer to confirm or replace; they are not implemented applicability rules.

The future contract must carry the complete cell-by-cell applicability matrix,
including an explicit reason code for every `not_applicable` cell. An
applicable cell may not fall back to a different artificial-feature generator.

## Estimand and invalid historical calculation

For each applicable global-null validation cell, Monte Carlo replicate `i`
contributes `X_i = n_selected_i / p_i`, where `0 <= X_i <= 1`. The proposed
estimand is the replicate-level mean selected-feature fraction for that cell.
The experimental unit is one independently seeded replicate, not a feature and
not a pooled selected-feature count.

Historical v1 code passed `sum(X_i)` to a binomial Wilson interval. Continuous
replicate fractions are not Bernoulli outcomes. The sealed historical result
must remain negative, non-promotable diagnostic evidence.

## Infeasible draft requiring replacement or explicit confirmation

The unaccepted v2 proposal used a one-sided Hoeffding bound with nominal
`alpha = 0.05`, criterion `U <= 0.10`, and 100 replicates. Its radius alone is
`0.12238734153404082`, so the gate cannot pass even when `mean(X) = 0`. If 44
applicable null cells form one Bonferroni family, the radius at 100 replicates
is approximately `0.1829`; the exact value and required replicate count must be
derived from the reviewer-approved multiplicity policy, not hardcoded here.

## Decisions required from the independent reviewer

The signed decision record must bind this packet's complete review-bundle
manifest and decide all of the following:

1. the estimand and replicate-independence assumptions;
2. the accepted one-sided confidence bound and confidence level;
3. the multiplicity family and adjustment over reviewer-approved applicable
   cells;
4. a jointly feasible threshold and replicate count;
5. the complete 96-cell applicability matrix, including the proposed eight
   fixed-X/high-dimensional exclusions and their reason codes;
6. fail-closed handling of missing, duplicate, non-finite, out-of-range,
   malformed, and incomplete replicate data;
7. deterministic reference inputs and exact expected outputs, including
   boundary equality; and
8. whether any additional scientific gate or reporting field is required.

The decision record must conform to `sci-decision.schema.json`, be signed with
the reviewer's detached signature, and authenticate against an externally
supplied `reviewer-trust-root.schema.json` record. Reviewer identity, trust
root, public key, decision, and signature remain outside the package tarball.

## Stop rule and post-review migration

Until a complete authenticated accepted decision exists, release-contract v1
remains strict, historical, and non-promotable. No formula, replicate count,
multiplicity policy, applicability semantics, methodology schema bump, or
promotable release-contract v2 may be implemented.

The authorized post-review continuation is additive: schema-dispatched v1/v2
validators, the full applicability matrix as contract data, explicit
`not_applicable` reason codes, no fallback in applicable cells, methodology
gate/artifact schema bumps, one canonical source of settings, and
reviewer-supplied deterministic tests. Multiplicity must be counted over the
reviewer-approved applicable cells rather than a hardcoded historical count.
