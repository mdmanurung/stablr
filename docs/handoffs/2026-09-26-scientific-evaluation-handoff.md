# stablr scientific evaluation — session handoff

Recorded: 2026-09-26 (Europe/Amsterdam)

Package: `/exports/archive/hg-funcgenom-research/mdmanurung/_misc/stablr`
branch `release/v0.1.1-completion`, HEAD `590cf6b0ab5c03f2af493a80e6426b4c6dd16b9f`

Companion repository: `/exports/para-lipg-hpc/mdmanurung/stablr-experiments` at
`4320b08cb5529ec911dde6732227de5c2550abbc`

Evidence tree: `/exports/para-lipg-hpc/mdmanurung/stablr-validation-artifacts/0.1.1/completion-20260920-fac5a622/`

## Read this first

**Nothing is committed.** The working tree carries 48 modified files and 11 new
ones: **+4662 / -451**. This is verified, coherent work, not a work-in-progress
diff, but it is uncommitted and therefore at risk. The first action on resume
should be to commit it.

The authoritative plan with the full execution log is
`docs/plans/2026-09-26-stablr-scientific-evaluation-plan.md`. Do not reconstruct
status from candidate outputs or from this file alone; use the plan's
"Execution log" section, which records per-phase before/after evidence.

## Environment

No R is on `PATH`. The single runtime is R 4.5.1:

```sh
export PATH=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin:$PATH
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5
RT=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin
```

**`sparsegl` lives only in `/exports/para-lipg-hpc/mdmanurung/R/4.5`, which
`--vanilla` discards.** Every Rscript invocation must prepend it, or the suite
silently skips the sparse-group-lasso assertions and still reports green:

```sh
Rscript --vanilla -e '.libPaths(c("/exports/para-lipg-hpc/mdmanurung/R/4.5", .libPaths())); <code>'
```

This is itself a recorded finding, not just a convenience: the release-evidence
harness mandates `--vanilla` and sets `R_PROFILE_USER=/dev/null`, so a real
release validation run under the sanitized environment would drop that library
and lose the coverage silently. See the plan's Phase 0 section.

## Completed

### Phase 0 — baseline and measurement instrument

Authoritative source-suite baseline, re-measured after all remediation:

| | before this plan | now |
|---|---|---|
| test files | 37 | 42 |
| `test_that` blocks | 418 | 574 |
| assertions | 2586 | 3383+ |
| failed / errors / skipped | 0 / 0 / 0 | 0 / 0 / 0 |

`R CMD check --as-cran --no-manual` improved from a stale recorded 2 WARNING /
5 NOTEs to **1 WARNING / 2 NOTEs**. The WARNING is `'qpdf' is needed for checks
on size reduction of PDFs` (missing host tool); the NOTEs are CRAN incoming
reporting `New submission` and Conda injecting `-march=nocona`. Neither is a
package defect.

Measured, not assumed: the source suite runs 3383 assertions with 0 skips, while
`R CMD check` on the tarball runs fewer because `test-release-hygiene.R` and
`test-vignette-readiness.R` are `.Rbuildignore`d and every CI job checks the
tarball. The OOL claim guard added in Phase 1 lives in one of those files and is
therefore source-only.

### Phase 1 — OOL parity reconciliation

Delivered `docs/release/0.1.1-ool-parity-reconciliation.md`. Verdicts: H1
supported, H2 supported, H3 refuted, H4 partially refuted, H5 clean, H6
supported.

- Spearman **0.808 withdrawn** from all five paper files. It came from a Bounded
  Run with no candidate SHA, no manifest and no seal.
- Jaccard 0.778 re-anchored to the sealed v2 **failed** packet, with its
  `not_ready` status, failed gate values, and the reason it failed (fixed-X
  knockoff infeasible at `p > n`, silent fallback to random permutation, so a
  different FDP+ null than the reference) now stated in the manuscript.
- Regression guard in `test-release-hygiene.R` fails if 0.808 reappears or if the
  manuscript stops presenting OOL as a negative result. Verified to bite.

### Phase 2A–2J — code correctness remediation

All test-first, every fix verified to fail before the change. No pinned fixture
was re-pinned.

| group | files | headline defects fixed |
|---|---|---|
| 2A | `knockoff_mvr.R`, `artificial_features.R` | all-zero `S` returned silently with an empty fallback history; `NaN` from an unguarded Cholesky downdate; divergence read as convergence; unscoped warning handler substituting a non-minimiser; silent `corpcor` fork; overflow input crashing the dispatcher |
| 2B | `learner_adapters.R` | adaptive lasso silently degenerating to plain lasso on small-scale responses |
| 2C | `stabl_fit.R`, `stabl_accessors.R` | `hard_threshold` silently overriding the computed FDP+ threshold; unguarded `family`/`l1_ratio`/`alpha`; zero-column `x` |
| 2D | `stabl_fit.R`, `bootstrap_helpers.R`, `artificial_features.R` | one RNG rule adopted and enforced: a seeded public entry point never advances the caller's `.Random.seed` |
| 2E | `conditions.R` | glmnet patterns 1–2 structurally unreachable; `eigen`/solve failures escaping untyped so the knockoff fallback never engaged |
| 2F | `metrics.R`, `visualization.R` | `adjusted_similarity` documentation contradicting its own formula; unvalidated universe size; FDP labelled FDR |
| 2G | `input_validation.R` | cooperative fusion defaulting to `rho = 0`, i.e. early fusion with no agreement penalty |
| 2H | `multiomic_workflows.R`, `late_fusion.R` | `stacked_multi_omic` returning `-Inf` with an all-`NA` column; AUC `0.5` for an absent class; two contradictory binary-outcome mappers (one produced all-`NA` for character outcomes) |
| 2I | `input_validation.R`, `bootstrap_helpers.R`, `stabl_fit.R` | `data.frame` outcomes resolved IDs from the column axis; silent group undersampling; non-finite training omics reaching glmnet |
| 2J | `nested_cv.R`, `learner_adapters.R` | unpaired inner-CV candidate comparison injecting bootstrap noise into a single-number selection rule; four duplicated lambda-grid resolvers |

### Phase 4 (partial) — release-evidence lifecycle tests

Delivered `tests/testthat/test-release-evidence-lifecycle.R`, 55 assertions,
all passing. Covers the previously untested state machine: tree-record tamper
detection, archive-member safety, canonical-JSON stability, no-replace writes,
explicit-packet-set validation (duplicate/missing/extra/foreign/failed/
unpromotable/unsealed), and the non-promotable-contract block.

**A real maintainability finding came out of this.** `.release_expected_contract()`
is a hand-maintained R literal that *flattens* the nested
`release-contract.json` (hoisting `runtime.dependency_roots` and
`runtime.environment`, projecting each required validation's settings, artifacts
and gates into parallel vectors). Nothing verified the two against each other, so
a contract edit could silently diverge from the code enforcing it. The new test
asserts the projection agrees with the JSON.

### Phase 6 — claim ledger

Delivered `docs/release/0.1.1-claim-ledger.md`: 21 claims, each with its
backing artifact and status (SUPPORTED / QUALIFIED / WITHDRAWN / UNVERIFIED),
plus the manuscript changes made and the six items still required before
submission.

## Key results

### The OOL `>= 0.90` Spearman gate is unreachable by construction

Measured, not inferred. Five independent `stabl_fit` runs on the exact pinned
OOL design (150 x 1317, 10 pinned lambdas, gaussian lasso, 500 bootstraps,
`knockoff_equi`), differing only in `random_state`:

```
CEILING (mean pairwise R-vs-R Spearman) = 0.8385
range = 0.8291 .. 0.8490
frozen gate requires >= 0.90  =>  UNREACHABLE
```

Two runs sharing every modelling choice and differing only in which bootstrap
resamples were drawn measure the noise floor of the estimator. **The gate sits
above that floor, so the protocol's Spearman threshold is a defect in the gate,
not in the package.** v2's 0.7847 is close to the ceiling, not evidence of a
package defect.

The Jaccard gate has the same structural problem in milder form: between two
independent R runs, Jaccard of the selected sets ranges 0.43–0.83 with only 4–6
features selected.

**The FDP-threshold gate is the one worth investigating**, because it is a
systematic difference rather than a noise floor: R gives 0.70–0.77 against the
frozen Python 0.58, a delta of 0.12–0.19 against a 0.05 threshold, exceeding the
seed-to-seed spread of 0.07.

**Consequence:** the full-profile OOL run recommended in the reconciliation
report is **deferred**. It cannot pass the current gates, so running it would
produce another sealed negative packet that says nothing new. The protocol's
Spearman and Jaccard thresholds must be revised first.

### FDP+ calibration study — IN PROGRESS, was running

`/tmp/opencode/p3/fdp-calibration.R`, log at `/tmp/opencode/p3/calibration.log`,
results at `/tmp/opencode/p3/calibration-results.rds`.

Design: 2 families x 4 artificial-feature generators x 4 scenarios = 32 cells x
30 replicates = 960 fits, at `n_bootstraps = 200`, `n_lambda = 30` — production-like
settings, deliberately unlike the committed RE-01 baseline which used
`n_bootstraps = 2, n_lambda = 2` and is therefore uninformative about calibration.
Cells mirror `inst/analysis/run_methodology_validation.R`: n = 60, p = 30, AR(1)
correlation 0 / 0.7, 5 planted signals for the signal scenarios.

It reports, per cell, the mean selected fraction with a one-sided Hoeffding upper
bound (the SCI-04 statistic) under the null, and mean TPR / FDP under signal.

**Status at interruption: 26 of 32 cells complete, still running.** The slow
cells are the fixed-X `knockoff` generator, which takes ~8–11 minutes per cell
because `knockoff::create.fixed` is called once per injected decoy. That cost is
itself worth recording: 30 decoys per replicate at ~1.5 s each.

**To resume the calibration:** if the process is gone, re-run with a reduced scope
that skips the fixed-X `knockoff` generator, which costs ~50 of the ~60 minutes
and is the slowest path by an order of magnitude:

```sh
N_BOOT=200 N_LAM=30 REPS=30 Rscript --vanilla \
  -e '.libPaths(c("/exports/para-lipg-hpc/mdmanurung/R/4.5", .libPaths())); source("/tmp/opencode/p3/fdp-calibration.R")'
```

Edit the `GENERATORS` vector in that script to drop `"knockoff"` first. Note that
this changes the cell set, so the applicability accounting must be restated
rather than merged with the completed run.

Also in flight at interruption: `R CMD check --as-cran --no-manual` on the
rebuilt tarball, log at `/tmp/opencode/check4.log`, started but not finished.

## Not yet done

### Phase 3 — finish the calibration study and write it up

1. Complete or re-run the study.
2. Decide whether the mean selected fraction under the null is compatible with
   the SCI-04 threshold of 0.10. If it is not, the gate is measuring a defect in
   the method rather than validating it, and the protocol needs revision — the
   same conclusion the OOL Spearman gate reached.
3. Write the result to `docs/release/`, with the explicit statement that the
   committed RE-01 baseline is not evidence in either direction.

### Phase 3 — SCI-04 v3 independent re-derivation

Partially done by inspection; not yet written up as a decision. The statistics
check out:

- Hoeffding radius `sqrt(log(20)/1200) = 0.0499644` at n = 600 — correct for a
  0.05 shift at alpha 0.05.
- Intersection–union over the 44 applicable null cells with **no multiplicity
  adjustment** — correct; the IUT type-I error is bounded by `max_i alpha_i` =
  0.05, not by `44 * 0.05`.
- `n = log(20) / (2 * 0.05^2) = 599.1`, so 600 replicates is the right size.

Two substantive reservations a reviewer must resolve, both independent of the
arithmetic:

1. **The estimand is a proxy.** `mean_replicate_selected_feature_fraction` is not
   the FDP. The release claim would be "FDP+ controls the false discovery
   proportion" while the gate measures "what fraction of features get selected".
   The linkage needs to be argued explicitly, or a second estimand based on the
   decoy-estimated FDP added.
2. **Attainability.** If the calibration study shows the null selected fraction is
   far above 0.05, then no number of replicates makes the gate pass, and
   accepting the proposal would authorise a release gate guaranteed to fail for
   reasons unrelated to statistical rigor.

### Phase 7 — decision records and contract v2

Not started. Required artifacts:

- An SCI-04 decision record conforming to
  `inst/analysis/contracts/sci-decision-v2.schema.json`. Binding values already
  computed:
  - review bundle: `/exports/para-lipg-hpc/mdmanurung/stablr-validation-artifacts/0.1.1/completion-20260920-fac5a622/review/sci04-v3-590cf6b-20260920/SCI-04-review-bundle-v3/`
    (sealed, `sha256sum -c manifest.sha256` verifies; manifest SHA-256
    `765c279bb59a9a987e66291cb42c95d47617ff5e862f76ca6639d374ed1f1d33`)
  - target: `SCI-04-independent-review-packet-v3`, SHA-256
    `3b984f6ae5c741f1593d0934c877c9a21da5f9e5fcfb1b9154edc044020af8a3`
  - applicability matrix SHA-256
    `594a1c870116dc0614adc9f30091c3114f1b7fd16e83dfe60e7e252925851009`
  - reference cases SHA-256
    `60091c815feb184752590940931fa9f3315e85c24c1483d597c3a0a2479022f5`
  - `reviewer_id`: the user chose an agent audit in place of an independent human
    statistician. **The schema has no independence field**, so a conforming
    record alone would misrepresent the review as independent. The caveat must
    live in the contract, `NEWS.md`, `README.md`, the manuscript, and a
    hash-bound sidecar provenance file that the contract references. Do not let
    the signature mechanism imply a human signatory.
  - The decision is likely **rejected** on the estimand-linkage and
    attainability points above, with a precise statement of what would make it
    acceptable. Do not accept a gate that cannot pass.
- An OOL decision record conforming to
  `stablr-experiments/analysis/contracts/ool-decision.schema.json`. Blocked on
  the protocol revision, since the current Spearman and Jaccard thresholds are
  above the measured noise floor.
- Contract v2, additive, leaving v1 immutable, removing the SCI-04 blocker only
  once a decision record verifies, and carrying an explicit `claim_scope` string
  so no consumer can read `promotable: true` as a universal FDP+ guarantee.
- Detached-signature verification. Not implemented anywhere in the package; the
  three decision schemas and the reviewer trust root are currently never loaded
  by any code or test.

### Phase 4 (remainder)

- The 7 `inst/release/schemas/*.json` are still never used to validate anything.
  Either wire a minimal in-house Draft 2020-12 subset validator into the sealing
  path with a negative test per schema, or state plainly in `inst/release/` that
  they are documentation only.
- `run_required_validation`, `assemble_release_evidence`,
  `.release_finalize_validation_attempt`, `.release_adapter_child_main`,
  `.release_build_runtime_artifact`, `.release_install_candidate` and
  `inst/analysis/release_validation_child.R` still have **no** coverage. They
  shell out to a real Rscript against an installed tarball, so they need an
  integration environment, not a unit test. The subprocess-free state machine is
  now covered.
- Record the effective `.libPaths()` and the skipped-test set in every Validation
  Run Packet, and treat an unexpected skip as a gate failure.
- Remove the empty, untracked `inst/analysis/adapters/` and
  `tests/testthat/fixtures/release_evidence/`.

### Phase 5 — CRAN readiness

- Re-run the package gate once the in-flight `check4` finishes and resolve
  anything it reports.
- Relocate the two `.Rbuildignore`d test files' assertions into a file that ships,
  or add a required CI job that runs `devtools::test()` from source. The OOL claim
  guard is currently source-only for this reason.
- Measure and record an actual `covr` number. None exists anywhere in the
  repository, the upload has `fail_ci_if_error: false`, and the coverage job is
  excluded from the evidence chain.
- Exact-candidate matrix receipts: Linux release/devel/oldrel-1, macOS, Windows,
  R 4.1 minimal, plus Win-builder and r-hub with sanitizers. The sanitizers are
  expected to exercise the `NaN` propagation and unguarded-`chol` paths fixed in
  Phase 2A.
- Only then de-draft `cran-comments.md`. `test-release-hygiene.R:167` requires it
  **not** to contain "0 errors | 0 warnings" until the real receipts exist.

## Open item for the author

**MVR fidelity above ~300 features.** `knockoff::create.gaussian` builds a
*diagonal* `S` from `diag_s`, which requires
`Sigma ⪯ diag(diag(S_mvr)) ⪯ 2*Sigma`. That sandwich does not hold in general for
the diagonal of the MVR solution, so the MVR path can reach `create.gaussian`
with an infeasible `diag_s` and silently degrade. Reproduced at p = 301, where
the required interval for `s_j` is empty. Phase 2A added a cost-based warning
ceiling and recorded it in provenance, but did **not** paper over this, because
whether the vendored port's diagonal extraction is faithful is a scientific
judgement, not a numerical guard. It also decides whether MVR belongs in the
Phase 3 calibration grid at all.

## Corrections made to my own earlier claims

Recorded so later phases do not act on the wrong version:

1. `.is_glmnet_numerical_error` pattern 3 was **not** broken where first
   reported. glmnet 5.0 emits the single string
   `one multinomial or binomial class has 1 or 0 observations; not allowed` for
   both arities, and the pattern matches it. The real defects were that patterns
   1 and 2 are structurally unreachable (glmnet signals negative C++ codes as
   warnings; the wrapper only installs an `error=` handler).
2. `stabl_fit`'s `lambda_grid` has **no default**. It is a required argument;
   `auto_lambda_grid()` runs only when a caller passes the string `"auto"`.
3. `fdr_similarity` / `tpr_similarity` / `fscore_similarity` already deduplicate,
   because `intersect()` and `setdiff()` are set-based. No parity divergence
   exists there. A narrower divergence does exist in `adjusted_similarity`, now
   documented and pinned.
4. The Phase 2A brief for the silent-solve defect was wrong: a triangular matrix
   with non-zero diagonal is never singular, so `forwardsolve`/`backsolve` cannot
   silently fail that way. The real route is a Cholesky factor that no longer
   factors `2*Sigma - S`, guarded by a reconstruction check.
5. `.release_expected_contract()` is not a duplicate of the JSON; it is a
   *flattened projection* of it. Whole-object identity is the wrong assertion.

## Commit guidance

The diff is large and coherent. Suggested split, matching the phases:

1. `docs:` plan, OOL reconciliation, claim ledger
2. `fix(conditions):` reachability and coverage of the numerical-error classifiers
3. `fix(metrics):` documentation, universe validation, FDP labels
4. `fix(knockoff):` MVR and artificial-feature numerics
5. `fix(stabl-fit):` contract validation, threshold suppression, RNG rule, bootstrap guards
6. `fix(multiomic):` cooperative defaults, stacking guards, binary outcome mapper
7. `test:` release-evidence lifecycle coverage
8. `docs(paper):` withdraw the unsupported OOL figure, re-anchor to the sealed packet

Check `git status` before staging; `stablr_0.1.1.tar.gz`, `stablr.Rcheck/` and
`src/*.o`/`*.so` are build byproducts and must not be staged.
