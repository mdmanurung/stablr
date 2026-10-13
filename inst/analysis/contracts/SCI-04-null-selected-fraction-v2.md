# SCI-04 null selected-fraction gate v2

Status: proposed, not accepted

Decision ID: `SCI-04-null-selected-fraction-v2-draft`

## Scope

The gate concerns the mean fraction of selected features under each null
simulation cell. One Monte Carlo replicate contributes one observation

`X_i = n_selected_i / p_i`, with `0 <= X_i <= 1`.

The historical v1 implementation passed `sum(X_i)` to a binomial Wilson
interval. That treats continuous fractions as Bernoulli successes and is not a
valid bound for this estimand. The sealed RE-01 packet and the committed
`tests/testthat/fixtures/release_evidence_baseline/v1/` files preserve that
historical result; they are diagnostic evidence only and must not be rewritten
or promoted.

## Proposed decision

Subject to independent statistical review, use the distribution-free
one-sided Hoeffding upper confidence bound for the replicate-level mean:

`U = min(1, mean(X) + sqrt(log(1 / alpha) / (2 * n)))`.

The proposed fixed settings are `alpha = 0.05`, a release criterion of
`U <= 0.10`, and the predeclared release replicate count. The experimental
unit is the independently seeded Monte Carlo replicate. Missing, duplicated,
non-finite, or out-of-range replicate fractions make the cell incomplete and
must fail closed.

This proposal was selected from the sampling model before any corrected bound
was computed. It must not be tuned against the failed locked run.

## Feasibility of the release criterion (audit E1, 0.2.2)

The Hoeffding margin depends only on `alpha` and the replicate count `n`:
`margin(n) = sqrt(log(1 / alpha) / (2 * n))`. The release criterion
`U <= 0.10` is achievable only if `margin(n) < 0.10`, because even a perfect
run (`mean(X) = 0`) yields `U = margin(n)`. At `alpha = 0.05`:

| replicates `n` | margin | criterion `U <= 0.10` achievable at all |
|---|---|---|
| 100 | 0.1224 | **no** (vacuous: `U >= 0.1224 > 0.10` always) |
| 150 | 0.0999 | yes (requires `mean(X) <= 0.0001`) |
| 200 | 0.0865 | yes (requires `mean(X) <= 0.0135`) |
| 250 | 0.0774 | yes (requires `mean(X) <= 0.0226`) |

The previously considered replicate count (n = 100) therefore makes the
proposed criterion **impossible to satisfy**: the gate would fail every run
regardless of the method's true null-selection behavior. This is a property
of the formula and the predeclared settings alone (no locked-run data was
used). Accordingly:

- the predeclared release replicate count for this gate must satisfy
  `n >= 150`; a cell designed with `n < 150` is infeasible by construction
  and must fail closed at design time, before any simulation is run;
- the deterministic reference tests pin this arithmetic
  (`tests/testthat/test-sci04-gate-feasibility.R`): `margin(100) > 0.10`,
  `margin(150) <= 0.10`, `margin(200) <= 0.10`;
- the minimum-feasible replicate count remains demanding (at `n = 150` the
  observed mean must be near zero). If the accepted decision record wants a
  less degenerate operating point, it must raise `n` or renegotiate the
  0.10 criterion explicitly; neither may be tuned against the failed
  locked run.

## Acceptance boundary

This record does not authorize an implementation. The proposed bound must not
be implemented or used for a release decision until an independent reviewer:

1. confirms the estimand, independence assumptions, formula, confidence level,
   and multiplicity policy;
2. approves or replaces this exact version in a new immutable decision record;
3. records reviewer identity, decision date, and the accepted record hash; and
4. approves deterministic reference tests independently of the failed run.

Until then, `null_selected_fraction` has no numeric release bound, carries
status `proposed_unaccepted`, and fails every promotable release packet.
