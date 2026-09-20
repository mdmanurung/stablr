# SCI-04 independent statistical review packet v1

Status: awaiting independent decision

Packet ID: `SCI-04-independent-review-packet-v1`

This packet requests a decision; it does not authorize implementation. The
executor preparing it is not an independent reviewer and cannot accept the
gate.

## Estimand and experimental unit

For each global-null validation cell, Monte Carlo replicate `i` contributes

`X_i = n_selected_i / p_i`, with `0 <= X_i <= 1`.

The estimand is the replicate-level mean selected-feature fraction for that
cell. The experimental unit is one independently seeded Monte Carlo replicate,
not an individual feature and not a selected-feature count pooled across
replicates. The release profile currently declares 100 replicates per cell.

There are 48 global-null cells in the locked matrix: four families, three null
scenarios, and four artificial-feature generators. The reviewer must decide
whether these 48 decisions form one multiplicity family, several predeclared
families, or separately interpreted cells.

## Invalid historical calculation

The v1 code passed `sum(X_i)` to a binomial Wilson interval. Continuous
replicate fractions are not Bernoulli outcomes, so that calculation does not
bound this estimand. The sealed RE-01 fixtures retain the historical result as
diagnostic, explicitly non-promotable evidence.

## Current unaccepted proposal

`SCI-04-null-selected-fraction-v2.md` proposes the one-sided Hoeffding bound

`U = min(1, mean(X) + sqrt(log(1 / alpha) / (2 * n)))`

with nominal `alpha = 0.05`, criterion `U <= 0.10`, and 100 replicates. No
multiplicity policy has been accepted.

## Feasibility defect requiring reviewer resolution

At 100 replicates and nominal `alpha = 0.05`, the Hoeffding radius is
`0.12238734153404082`. Therefore even `mean(X) = 0` gives `U > 0.10`; the
proposed gate cannot pass for any valid data. At least 150 replicates are
needed before the radius alone is at most 0.10.

If all 48 null cells use a Bonferroni allocation (`alpha / 48`), the
100-replicate radius is `0.1852961586820121`, and at least 344 replicates are
needed before the radius alone is at most 0.10. These calculations do not
select a replacement policy. They demonstrate that the draft settings are
internally incompatible and must not be implemented unchanged.

## Required reviewer decisions

The independent reviewer must record:

1. the estimand and independence assumptions;
2. the accepted one-sided bound and confidence level;
3. the multiplicity family and adjustment, if any;
4. a jointly feasible release threshold and replicate count;
5. fail-closed handling of missing, duplicate, non-finite, out-of-range, and
   incomplete replicate values;
6. exact deterministic reference inputs and expected outputs, including
   boundary equality; and
7. reviewer identity, decision date, accepted decision-record path, and
   SHA-256 after the decision record is finalized.

Missing, duplicated, non-finite, out-of-range, or incomplete replicates are
proposed to fail closed and never be silently removed. This handling is also
subject to reviewer acceptance.

## Stop rule

Until an independent reviewer supplies a complete accepted decision, the
embedded release contract remains non-promotable, `null_selected_fraction`
remains `proposed_unaccepted`, and no corrected SCI-04 implementation or
release-profile packet may claim acceptance.
