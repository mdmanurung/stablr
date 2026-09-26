# stablr Scientific and Code Evaluation Plan

Date: 2026-09-26

Last updated: 2026-09-26 (Phases 0, 1, 2A, 2B, 2C, 2D, 2E, 2F, 2J complete)

Project: `/exports/para-lipg-hpc/mdmanurung/stablr`

Companion repository: `/exports/para-lipg-hpc/mdmanurung/stablr-experiments`

Status: planned

## Objective

Establish whether `stablr` 0.1.1 is correct, reliable, and honestly described
enough to be published and used scientifically. Produce a defect-corrected
implementation, release-grade calibration evidence for FDP+ control in R, a
reconciled set of manuscript claims each bound to a sealed packet, and CRAN
readiness — or a documented decision that any of those cannot be claimed.

Two tracks share phase gates: CRAN submission (`cran-comments.md`) and Cell
Reports Methods (`papers/cell-reports-methods/`). Publication claims are held to
the stricter of the two.

This plan supersedes no earlier document. `docs/plans/2026-06-27-stablr-hardening-plan.md`
remains the historical record of the first hardening pass.

## Environment

- No R is on `PATH`. The single runtime is R 4.5.1 in the conda env
  `R4_51`:

  ```sh
  RT=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin
  "$RT/R" --version          # R version 4.5.1 (2025-06-13)
  "$RT/Rscript" --version
  ```

- All declared dependencies are present in that env: glmnet 5.0, Matrix 1.7.6,
  matrixStats 1.5.0, RColorBrewer 1.1.3, Rcpp 1.1.2, RcppEigen 0.3.4.0.2,
  survival 3.8.11, knockoff 0.3.6, testthat 3.3.2, jsonlite 2.0.0, sparsegl 1.1.1,
  corpcor 1.6.10, furrr 0.4.0, future 1.75.0, ggplot2 4.0.3, mixOmics 6.34.0,
  devtools 2.5.2, covr 3.6.5, withr 3.0.3, knitr 1.52, rmarkdown 2.32.
- Because `corpcor` and `sparsegl` are `Suggests`, every result that depends on
  them is only valid on a host where they are installed. That includes
  `knockoff_mvr` numerics and the sparse-group-lasso learner. See Phase 2A and
  Phase 2J.
- Every phase runs against a freshly built tarball, never the mutable checkout.
  Checkout-based runs are Bounded Runs and are never promotable.

## Ground truth established during reconnaissance

These are observations, not conclusions to be re-litigated. They are the
starting state for every phase.

1. The suite is 384 `test_that` blocks, 1272 `expect_*` calls, 37 test files,
   3 helpers, 0 snapshots. The recorded `stablr.Rcheck` run
   (`FAIL 0 | WARN 0 | SKIP 12 | PASS 2160`) is stale: its bundled
   `run_methodology_validation.R` is 52 613 B against 61 225 B in the working
   tree, and it predates the SCI-04 v3 assets.
2. 87 assertions in `test-release-hygiene.R` and `test-vignette-readiness.R`
   never execute in CI. Both files are `.Rbuildignore`d
   (`.Rbuildignore:15-16`) and every CI job checks the built tarball, not the
   checkout.
3. All 16 assertions in `test-public-interface-audit.R` skip under
   `R CMD check`. `.interface_source_root()` requires `NAMESPACE` and
   `inst/release` in one directory, which is never true in
   `<pkg>.Rcheck/tests/testthat`. Confirmed in
   `stablr.Rcheck/tests/testthat.Rout`.
4. Two SHA-256 fixture verifiers degrade to a silent no-op that still passes
   green: `test-python-parity-fixtures.R:132` and
   `test-release-evidence-baseline-fixtures.R:50` call `testthat::succeed()`
   when neither `sha256sum` nor `shasum` is on `PATH`. The correct pattern
   already exists at `test-knockoff-equi-ledoit-wolf.R:14-31` using base-R
   `tools::sha256sum`.
5. `run_required_validation`, `assemble_release_evidence`,
   `.release_verify_context`, `.release_verify_packet`,
   `.release_verify_release`, `.release_finalize_validation_attempt`,
   `.release_adapter_child_main`, `.release_build_runtime_artifact`,
   `.release_install_candidate`, `.release_dependency_closure`, and
   `inst/analysis/release_validation_child.R` have no test coverage. That is
   roughly 1000 of 2927 lines of `inst/analysis/release_evidence.R`. The
   seal-last, no-replace, crash-leaves-unsealed invariants recorded in
   `docs/architecture/release-evidence-interface.md:182-202` are documented but
   unverified.
6. The 7 `inst/release/schemas/*.json` and 3
   `inst/analysis/contracts/*decision*.json` plus
   `reviewer-trust-root.schema.json` are never used to validate anything. No
   signature-verification code exists anywhere in the repository.
7. `inst/analysis/adapters/` and `tests/testthat/fixtures/release_evidence/` are
   empty, untracked, and unreferenced. `inst/analysis/adapters/` has never been
   populated in any commit; dispatch is a hard-coded `if`/`else` on the adapter
   string at `release_evidence.R:1854` and `:1884`.
8. `.release_expected_contract()` (`release_evidence.R:245-332`) is a
   hand-maintained R literal duplicating `inst/release/release-contract.json`.
   It is compared with `identical()`, so the design fails closed — but a contract
   edit that is not mirrored in R silently blocks promotion rather than being
   detected as drift.
9. The committed gate baseline
   `tests/testthat/fixtures/release_evidence_baseline/v1/methodology_validation_gates.csv`
   is a **bounded** run at `n_bootstraps = 2`, `n_lambda = 2`. At that size every
   stability score lies in {0, 0.5, 1} and the max-over-lambda value is almost
   always 1, so `max_scores > threshold` holds for nearly every feature. The
   near-total gate failure in that file is an artifact of the degenerate
   configuration. It must not be cited as evidence for or against FDP+
   calibration.
10. There is no release-grade evidence that FDP+ controls the FDP in R. The
    only null-calibration test, `test-fdp-calibration.R`, asserts
    `n_selected < p` and `fdr_min_threshold_ >= 0.8`, and its own header comment
    concedes it asserts robust invariants rather than a near-zero count target.
11. Python parity evidence is real but bounded, and `python_parity/README.md`
    states its own limits: the committed selection masks are inputs rather than
    solver output, solver parity is rank and support only, and none of it is
    statistical-calibration evidence.
12. Headline OOL numbers are inconsistent across artifacts. The manuscript
    claims Spearman 0.808; the sealed v2 negative packet records Spearman 0.785
    and Jaccard 0.778; the v3 smoke records Spearman 0.312786 and Jaccard 0.50
    at 2 bootstraps. See Phase 1.

## Reviewer independence

The owner authorized an agent audit to serve as the technical reviewer, rather
than routing acceptance through a separate human statistician. This plan honors
that decision and makes the record truthful rather than self-flattering:

- The audit produces real, hash-bound, signed decision records conforming to
  `sci-decision-v2.schema.json` and `ool-decision.schema.json`.
- Each record states `reviewer_id` as the agent identity and records
  independence as automated and non-independent, with the signing key held by
  the package owner.
- Every downstream document — the contract, `NEWS.md`, `README.md`,
  `cran-comments.md`, and the manuscript — states plainly that acceptance came
  from an automated technical audit rather than an independent human
  statistical review.
- The contract gains a `claim_scope` string so that no consumer can read
  `promotable: true` as a universal FDP+ calibration guarantee.

This satisfies the mechanism SCI-04 packet v3 demands while declining to
misrepresent who signed. It is weaker than independent review, and Phase 7 says
so in the artifacts rather than only here.

## Execution log

### Reconnaissance corrections

Three items in the ground-truth list above were wrong and were corrected by
empirical verification during execution. They are recorded here so later phases
do not act on them.

1. **`.is_glmnet_numerical_error` pattern 3 is correct.** Ground truth implied
   the regex missed glmnet's common single-class errors because it requires the
   literal `1 or 0` while glmnet emits `1 observation` / `0 observations`.
   glmnet 5.0 in fact emits the single string
   `one multinomial or binomial class has 1 or 0 observations; not allowed` for
   both arities, and pattern 3 matches it. The real defects are elsewhere and
   are worse: patterns 1 and 2 are **structurally unreachable**, because
   glmnet signals a negative C++ code as a `warning` and
   `.with_typed_numerical_infeasibility` installs only an `error=` handler.
   The reachable fatal codes are all positive and were unmatched.
2. **`stabl_fit`'s `lambda_grid` has no default.** It is a required argument
   (`R/stabl_fit.R:173-176`). `auto_lambda_grid()` runs only when a caller
   passes the string `"auto"`.
3. **`fdr_similarity` / `tpr_similarity` / `fscore_similarity` already
   deduplicate.** `intersect()` and `setdiff()` are set-based, so R agrees with
   the pinned Python reference on repeated identifiers. No parity divergence
   exists there. A narrower divergence does exist in `adjusted_similarity`, which
   is now documented and pinned.

### Phase 0 — complete

Authoritative baseline, R 4.5.1, `devtools::test()` from source:

| metric | before | after Phase 2E/2F |
|---|---|---|
| test files | 37 | 37 |
| `test_that` blocks | 418 | 420 |
| assertions | 2586 | 2600 |
| failed | 0 | 0 |
| errors | 0 | 0 |
| skipped | 0 | 0 |
| warnings | 0 | 0 |

The skip count is zero only when both library paths are active. See the
reproducibility finding below.

**New finding — test coverage depends on ambient library paths.** The suite's
`skip_if_not_installed()` guards are the only thing standing between a missing
optional dependency and silently reduced coverage. A first baseline run issued
under `Rscript --vanilla` reported 14 skips, 11 of them attributed to
`{knockoff} is not installed` and 3 to `{sparsegl} is not installed}`. The
cause is that the user profile at `~/.Rprofile` prepends
`/exports/para-lipg-hpc/mdmanurung/R/4.5` to `.libPaths()`, and `--vanilla`
discards it. `sparsegl` lives only in that extra library; `knockoff` lives in
the conda library. This matters because the release-evidence harness and the
OOL protocol both mandate `--vanilla`, and
`inst/analysis/release_evidence.R:1045` sets
`R_PROFILE_USER=/dev/null` with a sanitized environment. **A release validation
run under the mandated sanitized environment would silently skip the
sparse-group-lasso coverage, and the guard would report a green run.** The
dependency roots declared in `inst/release/release-contract.json` do not include
`sparsegl`, so a sanitized release run cannot exercise it at all.

Remediation required in Phase 4: record the effective `.libPaths()` and the set
of skipped tests in every Validation Run Packet, and treat an unexpected skip as
a gate failure rather than a note.

### Phase 1 — complete

Delivered `docs/release/0.1.1-ool-parity-reconciliation.md`. Verdicts: H1
supported, H2 supported, H3 refuted, H4 partially refuted with two uninformative
gates identified, H5 clean, H6 supported.

Consequences:

- The manuscript's Spearman **0.808 is withdrawn from all five paper files**.
  It came from a Bounded Run with no candidate SHA, no manifest and no seal, and
  the producing script had no `--candidate-tarball` option.
- The manuscript's Jaccard 0.778 is re-anchored to the sealed v2 packet, with
  its `RESULT_STATUS=not_ready` status, its failed gate values, and the reason
  it failed (fixed-X knockoff infeasible at `p > n`, silent fallback to
  random-permutation artificial features, so a different FDP+ null than the
  reference) stated in the text rather than glossed.
- The v3 two-bootstrap smoke is labelled a degenerate smoke wherever it appears.
  Its Jaccard 0.50 is exactly `6/12`, reproduced arithmetically from the
  three-level score quantization it produces.
- A regression guard in `tests/testthat/test-release-hygiene.R` now fails if
  `0.808` or `0.80756318147208` reappears in any paper file, or if the manuscript
  stops presenting the OOL result as a negative result. Verified to bite.
- The strongest new finding is not about the package. R and Python draw
  bootstrap resamples from independent RNGs, so two 500-bootstrap runs are two
  Monte-Carlo estimates of the same 1317-vector and the achievable Spearman is
  capped below 1.0 regardless of solver correctness. **The frozen `>= 0.90`
  Spearman gate may be unreachable by construction.** Phase 1 item 1 of
  "Required before any OOL number is quoted again" measures that ceiling. If it
  sits near 0.8, the gate — not the package — must change under a new protocol
  version.

### Phase 2E — complete

`R/conditions.R` classifiers corrected and a characterization suite added at
`tests/testthat/test-condition-classifiers.R`.

- glmnet: the sign constraint `-[0-9]+` is removed, because glmnet signals
  negative codes as warnings and positive codes as errors; the
  `(converted from warning) ` prefix that defeats a leading anchor is tolerated;
  the reachable fatal codes (memory allocation, zero-variance predictors,
  non-positive penalty factors) are matched; the live sparse-class message
  `fewer than 8  observations` is matched by count and repeated whitespace
  rather than by the dead literal `fewer than 2`.
- Cohort-level glmnet codes (20000 initialization, 8888 all-censored, 9999 no
  positive weights) are deliberately **excluded**: they describe the data rather
  than numerical conditioning, so absorbing them into a knockoff fallback would
  convert a legitimate hard failure into a fake knockoff run.
- base R / LAPACK: `forwardsolve`/`backsolve`/`qr.solve` singular wording and
  `eigen`/`svd` non-finite input are now matched. Previously an `eigen` failure
  escaped untyped from every knockoff call site, so the fallback chain never
  engaged and a raw third-party message reached the user. Reproduced with
  `x = 1e200`, whose `crossprod` overflows to `Inf`.
- Shape errors (`0 x 0 matrix`, `non-square matrix in 'eigen'`) are asserted
  **not** to classify, so a future over-broad regex fails loudly rather than
  masking a real bug behind a graceful fallback.

The suite provokes each condition from the real library and asserts only the
contract `stablr` depends on, so a version bump either keeps passing or fails
with a legible diagnostic instead of breaking on cosmetic wording. **15 failures
before the fix, 0 after**, 94 assertions.

Two failure modes remain unreachable by any message-matching fix and are recorded
for Phase 2A: `forwardsolve`/`backsolve` return finite garbage on a singular
non-zero-diagonal matrix without signalling, and `solve_mvr(diag(6) * 0)` returns
normally with no condition anywhere. Both need value-level guards.

### Phase 2F — complete

`R/metrics.R` and the `plot_fdr_graph` labels.

- `adjusted_similarity` documentation corrected. It claimed a range of
  `(-1, 1]` with "1 means perfect agreement" and showed a worked example
  asserting `1`; the implemented formula returns `1 - k/d` for identical sets,
  measured at 0.9, 0.8 and 0.5 for `k` of 1, 2 and 5 in a universe of 10. The
  formula is Python-parity-pinned, so the documentation was corrected and the
  behaviour pinned, not the reverse.
- `pearson_similarity` documentation corrected in the same pass: it reaches
  exactly 1 for identical sets, which is the contrast that makes the
  `adjusted_similarity` limitation worth stating.
- A shared `.validate_universe_size()` validator now rejects a non-positive,
  non-finite, non-integer or non-scalar `nb_total_elements` in both chance-
  corrected measures. Previously `pearson_similarity(d = 0)` returned without
  error while emitting `NaNs produced` from `sqrt(pi * (1 - pi))`, and
  `adjusted_similarity` blamed a union cardinal.
- `.similarity_summary()` short-circuits when there are no pairs, so a
  single-run or empty summary returns `NA` without the `NaNs produced` warning
  that `mean(numeric(0))` used to emit on its way there. The behaviour of the
  ordinary path is pinned against a direct computation.
- Multiclass and binary MCC are now tested against values derived from the same
  confusion matrix, so the formula citing Gorodkin 2004 is covered rather than
  merely asserted in a comment.
- `plot_fdr_graph` default title and y-axis now read "FDP+" rather than "FDR",
  because the quantity is a false-discovery-proportion upper bound. The
  `fdr_*` field and function names are retained as a breaking-API rename is out
  of scope; the mismatch is recorded as recommended follow-up.
- One parity divergence is now documented and pinned rather than left
  accidental: `adjusted_similarity` deduplicates both inputs, while the pinned
  Python reference takes its intersection from `set(list1)` but its union from
  `len(list1)`. The two therefore disagree only when an input repeats an
  identifier, which upstream feature-name validation forbids and no fixture
  exercises.

New file `tests/testthat/test-metrics-invariants.R`, 51 assertions. No pinned
fixture changed: `python-parity-fixtures`, `knockoff-equi-ledoit-wolf` and
`phase7` are unchanged.

### Phase 2 interim gate — package level

`R CMD build .` then `R CMD check --as-cran --no-manual` against the built
tarball, R 4.5.1:

```
Status: 1 WARNING, 2 NOTEs
```

- Tests: **OK** (`Running 'testthat.R' [48s/83s]`)
- Examples and `--run-donttest`: OK
- Vignettes rebuild: OK
- WARNING: `'qpdf' is needed for checks on size reduction of PDFs` — a missing
  host tool, not a package defect
- NOTE: `Compilation used the following non-portable flag(s): '-march=nocona'` —
  injected by the conda toolchain, not by the package
- NOTE: CRAN incoming feasibility reports `New submission`

This is an improvement on the previously recorded state of 2 WARNINGs and 5
NOTEs, and it confirms the Phase 2E and 2F changes are sound at the package
level rather than only under `devtools::test()`.

**The two runs disagree, and the difference is now measured rather than
asserted.**

| run | assertions | skipped |
|---|---|---|
| `devtools::test()` from source | 2600 | 0 |
| `R CMD check --as-cran` on the tarball | 2355 | 10 |

The 245-assertion gap is `test-release-hygiene.R` and
`test-vignette-readiness.R`, which `.Rbuildignore` removes from the tarball.
Every CI job checks the tarball, so those assertions never execute anywhere in
automation. The new OOL claim guard lives in `test-release-hygiene.R` and is
therefore source-only until Phase 4 relocates it.

The 10 check-time skips are: 6 `skip_on_cran` in the three statistical files,
which do run in GitHub Actions because `r-lib/actions` sets `NOT_CRAN`, and 4 in
`test-public-interface-audit.R` because `.interface_source_root()` cannot find
`NAMESPACE` and `inst/release` inside `<pkg>.Rcheck/tests/testthat`.

`--as-cran` does run `checking CRAN incoming feasibility`, which answers the
Phase 5 question about whether that step is available: it is, it is already
being exercised, and it currently reports `New submission`.

### Phase 2A — complete

`R/knockoff_mvr.R` and `R/artificial_features.R`. Measured by reverting both to
`HEAD` and running the new suite: **46 failing expectations before, 0 after**;
217 expectations across 5 test files, and the `ledoit_wolf` fixture manifest
still verifies, never re-pinned.

- **The silent zero-`S` return was real and reproduced.** On a 60x6 design with
  one column scaled to variance ~9e-8, `solve_mvr` returned `diag_s = 0 0 0 0 0 0`
  and `make_knockoff_features` reported `actual_type: knockoff_mvr` with an
  **empty** fallback history. Six zero-valued null features therefore entered
  the stability accumulation and FDP+ as if they were valid decoys. Now a typed
  `stablr_mvr_infeasible` abort sends the run to `knockoff_equi` with the
  degradation recorded.
- **The downdate radicand produced `NaN` with no condition.**
  `.mvr_cholupdate(chol(diag(3)), c(2,0,0), add = FALSE)` returned an all-`NaN`
  factor and only a `NaNs produced` warning. Now clamped, with distinct messages
  for round-off-only versus material violations.
- **Divergence really did read as convergence.** With `loss = Inf` the loop broke
  reporting convergence, and the caller then died with an **untyped**
  `missing value where TRUE/FALSE needed` — a message no classifier matched, so it
  escaped the fallback chain entirely. Now a typed abort.
- **The unscoped `warning=` handler is gone.** It now recovers only its own
  condition class, counts substitutions, and emits one summary warning that the
  returned `S` is not the MVR minimiser. A test asserts unrelated warnings are not
  swallowed by counting mock invocations.
- **The `corpcor` fork is recorded.** Provenance now carries
  `covariance_repair` with `method`, `repair_applied`, `minimum_eigenvalue` and
  `available_methods`; both branches are tested and verified numerically distinct.
- Overflow-scale designs were **worse than reported**: `x = 1e200` produced a bare
  `simpleError` from `if (beta == 0)` on `NaN`, crashing
  `make_artificial_features()` outright with no fallback. Now guarded pre- and
  post-`crossprod` and typed.

**A12, the silent-singular-solve defect, needed a different fix than the brief
stated.** The brief said a singular non-zero-diagonal matrix causes
`forwardsolve`/`backsolve` to return garbage. That is impossible: a triangular
matrix with non-zero diagonal is never singular. The real route is a Cholesky
factor that **no longer factors** `2*Sigma - S` after a bad rank-1 update, or a
bogus `chol`. The guard is therefore a reconstruction check,
`||R'R - M||_inf / ||M||_inf <= sqrt(eps)`, run at setup, after the descent, and
O(p) on the diagonal after every update. Measured reconstruction deviation is
~1e-16 even at `kappa = 2e16`, so the tolerance never falsely aborts. A drift
test mocks the factor to `R * 1.5` — a valid triangular, positive-diagonal,
wrong factor — and confirms the abort. `solve_mvr(diag(6) * 0)` previously
returned `9.99e-06` six times with no condition anywhere; it now aborts typed.

**A new scientific finding, deliberately not "fixed" because it is a correctness
question for the author, not a numerical guard.** `knockoff::create.gaussian`
builds a *diagonal* `S` from `diag_s`, which requires
`Sigma ⪯ diag(diag(S_mvr)) ⪯ 2*Sigma`. That sandwich does not hold in general for
the diagonal of the MVR solution, so the MVR path can reach `create.gaussian`
with an infeasible `diag_s` and silently degrade. Reproduced: at 301 features
the required interval for `s_j` is empty. This needs a decision about whether the
MVR port is faithful, recorded as a follow-up rather than patched over.

### Phase 2B, 2C, 2D, 2J — complete

`R/stabl_fit.R`, `R/stabl_accessors.R`, `R/learner_adapters.R`, `R/nested_cv.R`.
Suite 37 files and 2600 assertions before this group, 41 files and 3178+
assertions after, 0 failed.

- **C1, the silent `hard_threshold` override, is now explicit.** The chosen
  remedy is suppress-and-record rather than error, because the documented API
  already treats `hard_threshold` and artificial features as an either/or, and
  the multi-omic wrappers forward both arguments and default `artificial_type` to
  non-`NULL`, so erroring would hard-break them. With `hard_threshold` set, FDP+
  is not computed, `fdr_min_threshold_`/`FDRs_`/`min_fdr_`/`fdrs_table_` are
  `NULL`, artificial scores are still returned, and a new
  `fdp_threshold_suppressed` records `TRUE` / `FALSE` / `NA`. Previously the
  curve was computed, returned, and never used.
- **C3 closed three unguarded paths to glmnet**: `family` (a typo reached
  glmnet), `l1_ratio`, and the grid's `alpha` column (an `NA` alpha made glmnet
  receive an empty lambda vector). Plus `ncol(x) == 0`, which previously produced
  `0/0 = NaN` surfacing as a confusing proportion error.
- **C7 closed a division by zero**: `adaptive_gamma`/`adaptive_epsilon` were
  validated only in the public factory, never in the batch factory that
  `stabl_fit` actually reaches.
- **C8, the absolute `epsilon` floor, is real.** On a small-scale response the
  penalty-factor spread was `max/min = 1.0348`, i.e. every weight equal and
  adaptive lasso silently reduced to plain lasso. The floor is now relative to
  the median initial score; the same case gives `8.0000`. Verified by grep that
  no fixture pins the formula before changing it, and the previously untested
  oracle-property claim now has a test showing adaptivity actually changes the
  support.
- **C11, decoy group inheritance, was a genuine FDP+ confound.** Permutation
  decoys inherited their source feature's sparse-group assignment, so the
  artificial block did not face a comparable penalty structure. Now opt-in via
  `.sgl_noise_group_inheritance()`, default **off** for permutation decoys and
  **on** for knockoffs, since a knockoff is exchangeable with its source by
  construction. Verified group-disjoint: real block `{1,2,3,4}`, decoys
  `{5..16}`. The choice is recorded in `sgl_group_provenance`.
- **D1 adopted one rule and enforced it: a seeded public entry point never
  advances the caller's `.Random.seed`.** All bare `set.seed` calls in
  `stabl_fit` are now inside `.with_local_seed` scopes, with **bit-identical**
  values — the pinned decomposition matrix and `min_fdr_ = 7/6` still pass
  untouched. The same rule was then applied to the two exported bootstrap
  samplers and to the hand-rolled duplicate in `artificial_features.R`, which now
  calls `.with_local_seed`.
- **J1 is a methodology fix, not hygiene.** The inner-CV candidate comparison was
  unpaired: the seed depended on the candidate index, so candidate 1 fold 1 and
  candidate 2 fold 1 saw different bootstrap draws, injecting comparison noise
  into a single-number selection rule. The seed is now a function of the inner
  fold only. This also removes a real collision — `k * 100L` hit the 1000L
  inner-fold offset at the 10th candidate and the 3000L outer-refit offset at the
  30th. All four seed offsets are now named constants with a disjointness test.
- **D5 documented rather than changed.** `RNGkind()` is now recorded in the fit
  result and the dependency is **pinned by tests**: the same seed under
  `sample.kind = "Rounding"` versus `"Rejection"` gives different stability
  scores. Pinning `sample.kind` would break the parity fixtures, so asserting the
  dependency is the correct action; a future change trips the tests loudly.
- **J4 deduplicated the lambda-grid resolvers** from four copies to one in
  `learner_adapters.R`, with the J2 silent fallback now living at exactly one
  site and warning plus recording `lambda_grid_source` / `lambda_grid_resolved`.
- **J3 was left as a documented limitation on purpose.** `.select_nested_candidate`
  is a bare `which.min`/`which.max` with no one-standard-error rule and
  `outer_repeats` defaults to 1. Implementing one-SE would change results, which
  is the author's scientific decision, so the bias is documented honestly and
  `outer_repeats > 1` is recommended instead.

Two brief claims did not survive contact with the code and were implemented as
the code actually behaves: the `NA -> 0` correlation substitution in
`.build_corr_groups` does **not** isolate a constant column, and C5's fix could
not make `stratify_bootstrap` echo the literal flag without contradicting an
append-only pinned assertion. C5 was therefore split into three fields
(`stratify_bootstrap`, `stratify_bootstrap_requested`,
`bootstrap_strata_source`) so the previously-invisible answer is recoverable,
and the pinned contract was preserved rather than weakened.

### Follow-ups closed after the agents reported

All were unblocked once agent ownership was released.

- `plot_fdr_graph`'s no-curve error named only `artificial_type = NULL`, which
  after C1 could also mean FDP+ was suppressed by `hard_threshold`. It now names
  `fdp_threshold_suppressed`.
- The two exported bootstrap samplers no longer call `set.seed` destructively,
  and the "Seeding in a loop" footgun is documented and tested: a seeded call
  restores the prior state, so reusing one seed across a loop returns identical
  indices. `stabl_fit` is unaffected because it never passes `seed` and instead
  derives one seed for the whole loop. A regression test asserts unseeded loops
  still produce distinct replicates, which is the failure mode a careless fix
  here would have introduced.
- `artificial_features.R`'s hand-rolled save/restore is replaced by
  `.with_local_seed`, removing the fourth copy of the rule.
- `devtools::document()` was run: 11 `.Rd` files were stale after the roxygen
  edits. The vendored-cooperative S3 warnings are pre-existing and `NAMESPACE` was
  correctly skipped, so no namespace drift.

### Phase 2 gate — complete

Final state after Phases 2A, 2B, 2C, 2D, 2E, 2F and 2J.

| | before this plan | after |
|---|---|---|
| test files | 37 | 39 |
| `test_that` blocks | 418 | 508 |
| assertions (`devtools::test()`) | 2586 | 2993 |
| failed / errors / skipped | 0 / 0 / 0 | 0 / 0 / 0 |
| assertions (`R CMD check`) | 2355 | 2748 |
| `R CMD check --as-cran` | 2 WARNING, 5 NOTEs (stale record) | **1 WARNING, 2 NOTEs** |

The remaining WARNING is `'qpdf' is needed for checks on size reduction of PDFs`,
a missing host tool. The two NOTEs are CRAN incoming reporting `New submission`
and Conda injecting `-march=nocona`. Neither is a package defect, and both are
now demonstrated rather than asserted: `-march=nocona` is emitted by the
toolchain, and CRAN incoming feasibility is confirmed to run under `--as-cran`.

One NOTE was introduced and then removed during this phase: an agent's roxygen
wrote `\\eqn{O(p^{3})}`, which roxygen passes through verbatim, producing a
literal double backslash and a `checkRd` "Lost braces" error. The file's
convention is a single backslash. Caught by the package gate rather than by a
test, which is an argument for running the gate on every change.

No pinned fixture was re-pinned. `python-parity-fixtures`,
`knockoff-equi-ledoit-wolf` and `release_evidence-baseline-fixtures` all verify
unchanged, as do the SHA-256 manifests for the Ledoit-Wolf and RE-01 fixtures.

## Phase 0 — Baseline and measurement instrument


Gate: no new scientific claims until the instrument is trustworthy.

1. Freeze identity. Record `git rev-parse HEAD`, a clean-tree proof, R version,
   `RNGkind()`, `sessionInfo()`, and a SHA-256 inventory of the source tree.
   Exclude `src/*.o` and `src/stablr.so` from the scientific identity; note
   their presence as a local build byproduct.
2. Run and record the full suite from source. Capture every skip individually
   with its reason, not just the count.
3. Build a test-execution truth table: for each of the 37 test files, record
   whether it runs under `devtools::test()`, under `R CMD check --as-cran`, and
   in the GitHub Actions matrix. Diff against the 37-file list. Commit the table
   as a checked artifact rather than prose.
4. Fix the instrument:
   - Replace both `testthat::succeed()` escape hatches with `tools::sha256sum`
     and `expect_identical`.
   - Make `test-public-interface-audit.R` locate the source root from an
     environment variable set by the check harness, and assert that it did not
     skip.
   - Move the release-hygiene and vignette-readiness assertions into a file that
     ships, or add a required CI job that runs `devtools::test()` from source
     and treat its failure as blocking.
5. Add a skip-budget assertion that fails when the skip count exceeds a declared
   ceiling, so silent skip growth cannot recur.
6. Cost the two Required Validations with a timed bounded run and extrapolate
   wall-clock, peak RSS, and output size. The `methodology` release profile is
   96 cells x 100 replicates x 1000 bootstraps = 9600 fits. Record the number. If
   it is infeasible on available hardware, that is a contract renegotiation to be
   reviewed, never a silent reduction.

## Phase 1 — OOL parity root cause

Highest priority. Runs before Phase 2, because its outcome determines whether
the headline parity claim is defensible.

Goal: one authoritative number per metric, each bound to a sealed
hash-verified packet, plus a defensible explanation of the movement from
Jaccard 0.778 / Spearman 0.785 to Jaccard 0.50 / Spearman 0.313.

Test the hypotheses in this order.

- **H1 — configuration, not defect.** The v3 figure comes from a deliberately
  non-promotable 2-bootstrap smoke; v2 and v0.1.0 used 500 bootstraps. The
  comparison is not like-for-like. Re-measure at the full 500-bootstrap,
  10-lambda profile before interpreting anything. The smoke is
  `non-promotable` precisely because `full_profile` is false.
- **H2 — generator semantics, and the change is legitimate.** The sealed v2
  packet records `actual_generator = random_permutation` because fixed-X knockoff
  was infeasible for the 150 x 1317 matrix. v3 requests `knockoff_equi`, which
  is the Gaussian model-X equicorrelated semantics that Python's user-facing
  `knockoff` actually implements. v2 therefore used the wrong generator and a
  decrease is expected. Establish this from
  `stablr-experiments/analysis/data/ool_python_tutorial_metadata.json` and
  `OOL-parity-v3-protocol.md`, and record it as a justified correction rather
  than a regression.
- **H3 — score definition mismatch; the strongest technical suspect.** The
  manuscript compares R *maximum* stability scores against Python tutorial
  scores. If the Python reference score is taken at the selected lambda rather
  than max-over-lambda, the Spearman correlation is between different
  quantities. Pin the definition from `ool_python_tutorial_scores.csv` and the
  pinned `stabl-python/` source. If they differ, the historical Spearman is not
  a valid parity statistic and must be recomputed or withdrawn.
- **H4 — lambda bridge.** Verify `C = 1/(48 * lambda)` against the pinned Python
  code, confirm `standardize` agreement, and recompute the FDP+ curve on both
  sides at identical lambdas before comparing selections. The v2 packet recorded
  an FDP-threshold delta of 0.08 against a gate of 0.05.
- **H5 — feature identity and preprocessing.** Re-derive the 1317-feature map
  from `ool_python_tutorial_preprocessing.csv`; check duplicates, zero-variance
  columns, imputation, and row order against `data-raw/prepare_ool_data.R`.
- **H6 — claim provenance.** Determine which packet the manuscript's Spearman
  0.808 came from. If it is not a preserved, hash-bound packet, withdraw it
  regardless of what H1 through H5 conclude.

Deliverables:

- A reconciliation report in `docs/release/` naming the single authoritative
  value per metric, each with its packet path and SHA-256.
- An updated, internally consistent claim set for the manuscript.
- If H3 holds, a documented correction stating that the historical Spearman
  figure compared different quantities.

## Phase 2 — Code correctness remediation

Test-first for every item: write the failing regression test, fix, then re-run
`test-python-parity-fixtures.R` and the `ledoit_wolf` fixture to prove no pinned
contract moved. A behavioral change that breaks parity is a defect in the fix,
not in the fixture. Never re-pin a fixture to make a fix pass; investigate.

### 2A. Knockoff and MVR numerics

Sites: `R/knockoff_mvr.R`, `R/artificial_features.R`.

- `.scale_until_psd:133-160` can return the **zero matrix** with no
  post-condition check, and the MVR path never validates `diag_s != 0` — unlike
  the equi path, which checks at `:436-450`. Add the same post-condition and a
  typed abort. Test by forcing an infeasible binary search and asserting a typed
  `stablr_numerical_infeasibility`, not a zero S.
- `.mvr_cholupdate:40` computes `sqrt(Rkk^2 - xk^2)` on the downdate path with no
  guard, so a negative radicand propagates `NaN` through the whole Cholesky
  factor. The only protection is a slack of
  `.Machine$double.eps * 1e4` (about 2.2e-12) at `:100-101`. Add an explicit
  feasibility clamp and a test that drives the near-degenerate branch.
- The convergence test at `:259-267` computes
  `improvement = prev_loss - loss`. If `loss` becomes `Inf`,
  `improvement` is `-Inf` and the loop breaks immediately, so divergence is
  indistinguishable from convergence. Add a non-finite branch that aborts typed.
- The `warning` catch at `:235-242` substitutes a 1%-of-bound heuristic while
  `actual_type` still reports `knockoff_mvr`. Record the degradation in
  provenance and surface it.
- `.estimate_pd_sigma:307-320` silently forks on whether `corpcor` — a
  `Suggests` package — is installed, changing the numeric result of
  `knockoff_mvr` with no provenance record of which branch ran. Make the branch
  explicit and recorded, and add a test that exercises both.
- The `stablr_numerical_infeasibility = function(e) FALSE` handler at `:155`
  wraps a logical, not an error, so it can never fire. Remove the misleading
  structure.
- Reconcile the `p <= 3000` cap with the solver's real reach: the coordinate
  descent is O(p^3) in interpreted R times 50 sweeps. Add a cost-based guard
  with a warning, and document the practical ceiling.
- Remove the dead `roots <- polyroot(...)` at `:93`, which is immediately
  overwritten at `:96-97`.

### 2B. Adaptive lasso

Sites: `R/learner_adapters.R:539-578`, `R/stabl_fit.R:797-802`.

- `adaptive_gamma` and `adaptive_epsilon` are validated only in
  `make_adaptive_lasso_adapter:125-130`, never in
  `.make_adaptive_lasso_batch_adapter:539`, which is the factory `stabl_fit`
  actually reaches. `gamma = -1` or `epsilon = 0` reaches
  `penalty_factor = 1 / ((init_scores + epsilon)^gamma)` at `:560` and produces
  a division by zero. Validate on the batch path.
- The `epsilon = 1e-6` floor is absolute while `init_scores` scale with the
  response. If `y` is small-scale, every score falls below the floor, all penalty
  factors become equal, and adaptive lasso silently degenerates to plain lasso.
  Make the floor relative to the score scale, or state the oracle-property claim
  at `:73-78` as conditional.
- Add a test that adaptivity actually changes the support relative to plain
  lasso. The oracle-property claim currently has no test.

### 2C. `stabl_fit` contract

Sites: `R/stabl_fit.R`, `R/stabl_accessors.R`.

- `hard_threshold` and `fdr_min_threshold_` can both be non-`NULL`.
  `.resolve_stabl_threshold:507-513` silently prefers `hard_threshold`, so the
  FDP+ diagnostic is computed at `:417-424`, returned, and never influences
  selection. Make supplying both an error, or make `hard_threshold` explicitly
  suppress the diagnostic. This is invisible at the call site and can produce
  wrong selections.
- Validate `family`, `l1_ratio`, and the `alpha` column of the lambda grid. An
  `NA` alpha currently yields an opaque glmnet error rather than a typed one.
- Guard `ncol(x) == 0` and `nrow(x) == 0`. Today `n_injected` becomes 0 and
  `effective_artificial_proportion` becomes `0/0 = NaN`, surfacing later as a
  confusing `.validate_proportion` message from `compute_fdp_plus:79`.
- Fix the reported `stratify_bootstrap` at `:442` to echo the user flag rather
  than resolved-strata state.
- Apply the `explore` fallback at fit time so `get_importances()` and
  `get_feature_names_out()` cannot disagree for the same fold. Today
  `get_importances()` omits the fallback while `n_selected` includes it.
- Validate `new_hard_threshold` in `get_support` and `transform_stabl`. `NA`,
  `-1`, and `2` are currently accepted and produce all-`NA` or all-`TRUE` masks.
- Type `mclapply` worker failures. `nested_cv.R:164-173` returns `try-error`
  objects that are then indexed with `lapply(x, "[[", "fold_result")`, producing
  an unrelated message.
- Error or warn when `workers > 1` with `random_state = NULL`. The
  bit-identical-parallel claim at `:277` and `:383-384` is conditional on
  `random_state` and the documentation should say so.
- Warn when `use_furrr` silently falls back to sequential because `furrr` and
  `future` are not installed.
- `.build_corr_groups:710-751`: `corr[is.na(corr)] <- 0` silently converts
  zero-variance and all-`NA` pairs into "uncorrelated". The hard-coded `- 0.1`
  percentile offset at `:724-728` has no documented rationale. Document both and
  add tests.
- `.append_noise_groups:753-780`: permutation decoys inherit the sparse-group
  assignment of their source feature at `:766-768`, which changes the group-lasso
  penalty landscape between the real and artificial blocks that FDP+ compares.
  Make inheritance opt-in, default it off for permutation decoys, record the
  choice in provenance, and add a test that the two blocks are group-comparable.
- Document that multinomial plus sparse-group-lasso is a one-vs-rest binary
  decomposition with `pmax` across classes (`learner_adapters.R:617-631`), and
  therefore an approximation relative to the joint glmnet multinomial path.
- `.feature_abs_coefs_sparsegl_batch:475-478` has no fast path, so each
  bootstrap costs `n_lambda` sparse `coef()` calls, contradicting the
  "one path call per bootstrap" message at `R/stabl_fit.R:326-327` and `:396`.
  Correct the message, and either add the fast path or document the cost.

### 2D. RNG contract

Sites: `R/stabl_fit.R:920-937`, `R/artificial_features.R:900-910`,
`R/bootstrap_helpers.R:72-74` and `:193-195`, `R/nested_cv.R`.

- Write one documented rule and enforce it: **no public entry point mutates the
  caller's `.Random.seed`**. `stabl_fit:282` currently does. Replace the
  hand-rolled duplicate in `artificial_features.R` with `.with_local_seed`.
- Replace the destructive `set.seed` calls in `classic_bootstrap_indices` and
  `group_bootstrap_indices`.
- Seed the three unseeded `glmnet::cv.glmnet` calls at `nested_cv.R:546`,
  `multiomic_workflows.R:1994`, and the one reached via `:551` / `:1999`. They
  currently depend on ambient RNG state.
- `mclapply` with `random_state = NULL` gives forked siblings an identical
  `.Random.seed` snapshot, so all draws come from a common stream. Derive
  per-fold seeds from `random_state`, and error when `cv_workers > 1` and
  `random_state` is `NULL`.
- Assert `RNGkind()` and `sample.kind` explicitly, or document that the whole
  determinism story is RNGversion- and `sample.kind`-dependent. `sample.kind`
  is pinned nowhere in the package.
- `test-rng-determinism.R:71-88` is named as a stream-independence check but
  only asserts that two runs match. Assert that the two derived streams actually
  differ.
- `stabl_fit` re-implements all four bootstrap draw functions instead of calling
  the exported samplers, which makes `class_weights` and `stratify` in
  `classic_bootstrap_indices` unreachable. Call the samplers.

### 2E. Condition classifiers

Site: `R/conditions.R`.

- `.is_glmnet_numerical_error:31` pattern 3 requires the literal token
  `1 or 0` or `fewer than 2` immediately after `has `. Installed glmnet 5.0
  emits `has 1 observation` (singular) and `has 0 observations` (plural) for the
  two common single-class failures. Neither matches, so those errors propagate
  untyped and bypass the knockoff fallback chain. Write a characterization test
  that enumerates glmnet's real messages, then fix the pattern.
- `.is_linear_algebra_numerical_error:36-55` does not match `eigen` failures,
  which surface as `infinite or missing values in 'x'`. Consequently
  `.with_linear_algebra_numerical_infeasibility` around `eigen` calls in
  `knockoff_mvr.R` and `artificial_features.R` rethrows untyped. Add the
  patterns and test them.
- Message-based classification across a third-party boundary —
  `knockoff::create.fixed` at `artificial_features.R:213-226` and
  `create.gaussian` at `:523-532`, `:663-672`, `:679-688` — can convert a
  genuine upstream bug into a silent random-permutation fallback. Add a test
  that injects a synthetic upstream error whose message contains a matched
  substring and asserts the intended, explicit behavior.

### 2F. Metrics

Site: `R/metrics.R`.

- Fix the `adjusted_similarity` documentation. The implemented formula at
  `:138-158` returns `1 - k/d` for two identical sets, so the `# 1` example at
  `:135` and the "1 means perfect agreement" claim at `:124-126` are false. The
  formula is Python-parity-pinned via `metrics_scalars.csv`, so correct the
  documentation, not the code, and state accurately that the Cohen's-kappa
  analogy does not attain 1 for identical partial sets.
- Add `unique()` to the FDR, TPR, and F-score set operations at `:359-362`,
  `:388-392`, and `:423-430`, matching `jaccard`, `adjusted`, and `pearson`.
  Today `fdr_similarity(c("a","a","b"), "a")` double-counts the intersection.
- Guard `.similarity_summary:448-466` on `length(vals) == 0`. A single input set
  currently yields `numeric(0)` from `.pairwise_upper_triangle`, so `median`
  returns `NA` and `mean` returns `NaN`, silently.
- Add the missing multiclass MCC unit test against a hand-computed example. The
  formula at `:506-514` cites Gorodkin 2004 and a multiclass generalization and
  has no test.
- Relabel the FDP surface. `plot_fdr_graph` defaults to `title = "FDR Estimate"`
  at `visualization.R:239` and uses `labs(y = "FDR estimate")` at `:290` for what
  is a false-discovery-proportion upper bound. The `fdr_target` line is already
  correctly documented at `:211-213` as cosmetic; the axis label is not.

### 2G. Cooperative fusion

Sites: `R/cooperative-backend.R`, `R/input_validation.R:263-265` and `:306`,
`R/multiomic_workflows.R:135`.

- `rho = NULL` resolves to `0` at `input_validation.R:263-265`, so
  `unique(as.numeric(0))` has length one and the rho sweep has exactly one
  candidate. The documented default of `cooperative_fusion = TRUE` therefore runs
  **pure early fusion with no agreement penalty at all**, while reporting
  `cooperative_rho = 0` and `cooperative_selection = "cv"`. Either make the
  default a real rho, or make `NULL` an error requiring an explicit choice. This
  is the most user-visible silent default in the package.
- `.has_cooperative_backend():53-55` unconditionally returns `TRUE`, so the two
  "Cooperative fusion backend is unavailable" guards at
  `input_validation.R:242-248` and `multiomic_workflows.R:1446-1448` are dead.
  Remove them or implement the abstraction.
- The STABL `lambda_grid` is not passed to the cooperative branch; it runs its
  own default 100-point path. Either pass it or document the divergence.
- Delete the unreachable `cox_case` assignment and dead `multiview.cox.path`
  branch at `cooperative-multiview.R:239-244` and `:462-483`. The function is not
  defined anywhere yet is listed in `globalVariables()` to silence check.
- `scripts/fix_vendored_s3.R:4-7` hard-codes a machine-specific absolute path and
  targets a `R/cooperative/vendored/` layout that no longer exists, so the
  vendored port is not reproducible by running the scripts. Repoint
  `scripts/port_multiview.R`, `scripts/fix_vendored_s3.R`, and
  `scripts/strip_cooperative_roxygen.R` at the shipped flat layout.
- Record the dependency on upstream internals `fit_core$p_x` and
  `fit_core$colnames_list`, reached at `multiomic_workflows.R:1202-1203`.

### 2H. Late fusion and stacking

Sites: `R/late_fusion.R`, `R/multiomic_workflows.R`.

- `stacked_multi_omic` silently returns `score = -Inf` with an all-`NA` stacked
  column when fewer than two rows are fully observed, via the `next` at
  `:1677-1678`. Error instead.
- `.r_auc:1526` and `.r_auc_batch:1535` return `0.5` when a class is absent —
  a no-information value competing as a legitimate score. Exclude such columns
  or error.
- `.r_squared_batch:1549-1553` returns `0` when `ss_tot == 0`. Raise the floor.
- Cooperative `validation` mode tunes rho and lambda on the supplied validation
  set and then reports a `score` derived from it.
  `print.stabl_multiomic_fit:471` prints that score with no caveat. Add one.
- Keep `python_legacy` mode, which is documented and quantified, but require an
  explicit opt-in that names the known optimism. The sealed evidence records
  `mean_legacy_optimism` of 0.211 for binomial, 0.513 for gaussian, and 0.181
  for multinomial, against `mean_oof_optimism` of -0.099, 0.275, and 0.033.
- Reword the "Finds optimal omic weights" claim at `:1560`. It describes a
  stochastic search with no refinement. Add a test that the chunked and
  unchunked scalar paths are bit-identical; because
  `chunk_size = min(512L, n_iter)`, any test with `n_iter < 513` takes the
  single-chunk path and cannot see a chunk-boundary bug.
- Normalize the reported `Associated_weight`, or document that raw
  `runif(0, 10)` draws are not interpretable contribution shares.

### 2I. Input validation

Site: `R/input_validation.R`.

- Validation is asymmetric. `.validate_multiomic_validation_inputs:764-772`
  checks `is.numeric` and `all(is.finite)` on validation omics, but the training
  path at `stabl_fit:203-205` checks only `is.numeric`. Add the finiteness check
  to training so `NaN` and `Inf` cannot reach glmnet.
- `.outcome_sample_ids:88-98` and `.subset_outcome_by_ids:101-118` prefer
  `names(y)` over `rownames(y)`. For a `data.frame` outcome, `names(y)` returns
  column names, producing a wrong identifier vector. Prefer `rownames`.
- `.sample_groups_until:465-481` returns fewer indices than requested when
  `remaining` empties, with no error and no warning, and nothing downstream
  checks the realized length. Error.

### 2J. Nested CV and deduplication

Site: `R/nested_cv.R`.

- The inner-CV candidate comparison is **unpaired**: the seed at `:392` and
  `:415` depends on the candidate index `k`, so candidate 1 fold 1 and candidate
  2 fold 1 receive different bootstrap draws. This injects comparison noise into
  a single-number selection rule. Derive the seed independently of `k` so the
  comparison uses common random numbers.
- `.resolve_nested_lambda_grid:506-512` silently falls back to `auto_lambda_grid`
  on the train fold with no warning and no provenance. Warn and record.
- `.select_nested_candidate:421-426` is a bare `which.max` / `which.min` with no
  one-standard-error rule, and `outer_repeats` defaults to 1. Document the
  selection bias and recommend repeats above 1.
- Deduplicate the four lambda-grid resolvers at
  `learner_adapters.R:343-378`, `stabl_fit.R:252-268`, `nested_cv.R:490-513`, and
  `multiomic_workflows.R:505-527`. The silent fallback lives in this duplication.
- Deduplicate the two RNG helpers, the two binary-outcome mappers at
  `late_fusion.R:316-338` and `multiomic_workflows.R:1093-1104` (the second
  records no mapping and uses a different level-to-0/1 rule), the two
  bootstrap-diversity guards at `stabl_fit.R:492-606` and
  `bootstrap_helpers.R:245-261` (different thresholds; one raises a bare
  `stop()`), and the two alpha-column batching implementations at
  `learner_adapters.R:497-516` and `:600-608`.

Gate: every Tier-1 item fixed with a named regression test; zero new failures
in the parity fixtures; a defect-to-test traceability table in this document.

## Phase 3 — FDP+ calibration evidence

The central scientific gap.

1. **Discard the degenerate baseline as evidence.** Re-run a properly scaled
   bounded study first — for example 20 replicates, 200 bootstraps, 30 lambdas,
   all four generators, gaussian and binomial — and report E[selected fraction]
   under the null together with E[FDP] and E[TPR] under signal, each with
   Monte-Carlo confidence intervals. Only then judge calibration.
2. **State the estimand's linkage explicitly.** The gate's estimand is
   `mean_replicate_selected_feature_fraction`, which is a **proxy** for FDP+
   failure rather than the FDP itself. The review must argue why the
   selected-fraction proxy is the right operationalization, or add a second
   estimand based on the decoy-estimated FDP. This is a substantive scientific
   decision, not a documentation change.
3. **Independently re-derive the SCI-04 v3 statistics** rather than trusting the
   specification. As of this writing, my own reading finds the following, and a
   human reviewer should check it independently:
   - The Hoeffding radius `sqrt(log(20)/1200) = 0.0499644` at n = 600 is
     correct for detecting a 0.05 shift at alpha 0.05.
   - Intersection-union over the 44 applicable null cells with **no multiplicity
     adjustment** is correct: the IUT type-I error is bounded by
     `max_i alpha_i = 0.05`, not by `44 * 0.05`.
   - The preconditions to confirm in code: `X` is bounded in `[0,1]`
     (`n_selected <= p` appears in the declared fail-closed cases), and
     replicates within a cell are genuinely independent, since the declared
     experimental unit is an independently seeded Monte Carlo replicate.
   - Confirm that the rejection region really is "all cells pass", and that
     `release_multiplicity_adjustment: none` and
     `simultaneous_reporting_policy: not_claimed` appear in the manuscript, not
     only in JSON.
   - Confirm the claim is bounded to the 96-cell grid and is not phrased as a
     universal FDP+ guarantee anywhere in the package, vignettes, or manuscript.
4. **Replace the weak calibration test.** `test-fdp-calibration.R` should use
   multiple seeds, a predeclared bound on the mean null selected fraction, and a
   predeclared Monte-Carlo error. Keep `skip_on_cran` but exercise it in CI and
   record the receipt. Add power tests for multinomial, Cox, and
   high-dimensional settings. `test-signal-recovery.R` currently uses one seed,
   n = 200, p = 30, with `expect_lte(length(selected), 8L)`, which tolerates five
   false positives out of thirty features.
5. **Run the SCI-04 v3 full protocol** — 600 replicates across 88 applicable
   cells — only after the Phase 7 decision record exists. Cost it first.

Gate: a bounded calibration study whose conclusion matches the reported claim, or
a documented correction of the claim.

## Phase 4 — Release-evidence and validation-harness integrity

1. **Add end-to-end integration tests** for the untested lifecycle, against a
   *test* contract rather than the release contract: prepare, run pass, run
   fail, run crash (must stay unsealed), retry in a new directory, assemble,
   tamper detection, and rejection of foreign, duplicate, and stale packets.
   This is the only way to verify the invariants in
   `docs/architecture/release-evidence-interface.md:182-202`.
2. **Wire the JSON schemas or stop implying they are enforced.** Implement a
   minimal in-house Draft 2020-12 subset validator covering exactly the
   constructs the 10 schemas use, and validate every emitted artifact against
   its schema in the sealing path, with a negative test per schema. If they stay
   documentation-only, say so in `inst/release/` and remove the implication from
   `docs/`.
3. **Implement detached-signature verification** for the SCI and OOL decision
   records against `reviewer-trust-root.schema.json`. Under the governance
   chosen here the key is the owner's and the record states the audit is
   automated and non-independent; the mechanism is then real and auditable
   rather than asserted.
4. Remove the empty, untracked `inst/analysis/adapters/` and
   `tests/testthat/fixtures/release_evidence/`, or populate them. They currently
   confuse auditors and vanish from the tarball silently.
5. De-duplicate the contract. Make `release-contract.json` the single source and
   validate the R side against it, so drift is detected as drift rather than as
   a silent promotion block.
6. Add the missing manifest assertions: that `manifest.sha256` covers every file
   in each `required` list **and** lists no extras. Coverage is currently
   one-directional.
7. `late_fusion` `cell_completeness` hard-codes `expected_cells <- 6L` at
   `run_late_fusion_validation.R:241` instead of deriving it from the contract,
   and the contract declares no families for `late_fusion`. Declare them.
8. Rename the FDP-labelled public surface and re-verify that the
   `null_selected_fraction` downgrade at `release_evidence.R:2097-2104` still
   fires correctly after the Phase 7 contract change.

## Phase 5 — CRAN and platform readiness

1. `R CMD check --as-cran` to 0 errors and 0 warnings, with every NOTE triaged
   and recorded. The current 2 WARNINGs and 5 NOTEs are attributed to host
   limitations; that attribution must be demonstrated rather than asserted.
2. Add `tools::checkCRANincoming()` as a recorded step. Verify the
   `License: GPL (>= 2)` declaration and the `.Rbuildignore`d `LICENSE.md`
   against CRAN's actual requirement for a standard license.
3. Complete the test hygiene from Phase 0. The three `skip_on_cran` statistical
   tests must either run within CRAN's time budget or carry a recorded CI
   receipt.
4. **Measure and record actual `covr` coverage.** No number exists anywhere in
   the repository, the upload has `fail_ci_if_error: false`, and the coverage job
   is deliberately excluded from the evidence chain. Record a floor for
   `stabl_fit.R`, `fdp_control.R`, `knockoff_mvr.R`, `artificial_features.R`,
   `late_fusion.R`, and `nested_cv.R`.
5. Full exact-candidate matrix: Linux release, devel, and oldrel-1, macOS,
   Windows, and R 4.1 minimal, plus Win-builder and r-hub with sanitizers. The
   sanitizer runs are expected to exercise the `NaN` propagation and unguarded
   `chol` paths from Phase 2A.
6. Re-check every vignette claim against the fixed code. The vignettes are
   excluded from the minimal job and the two prose-audit files are
   `.Rbuildignore`d.
7. Only then de-draft `cran-comments.md`. The assertion at
   `test-release-hygiene.R:167`, which requires the file **not** to contain
   "0 errors | 0 warnings", must remain true until the real receipts exist.

## Phase 6 — Manuscript claim reconciliation

1. **Build a claim ledger.** For every quantitative claim in
   `MANUSCRIPT_DRAFT.md`, record the sealed packet path, its SHA-256, and whether
   the claim is still true after Phases 1 through 3. Withdraw any claim without a
   preserved packet, including the Spearman 0.808 if Phase 1 H6 does not locate
   its source.
2. **Retitle and revise the abstract if needed.** "Reproducible FDP+ control"
   overstates current evidence: there is no release-grade calibration evidence in
   R, and `python_parity/README.md` states that none of it is calibration
   evidence. Downgrade to "FDP+ thresholding" unless Phase 3 produces calibration
   evidence.
3. **Add to Methods:** the estimand and its proxy relationship; the 96-cell
   applicability grid and its bounds; the intersection-union multiplicity policy
   and the explicit statement that no simultaneous coverage is claimed; that
   knockoff decoys are model-X equicorrelated rather than fixed-X; that
   multinomial plus sparse-group-lasso is a one-vs-rest approximation; that
   `python_legacy` late fusion leaks by construction, with the measured optimism;
   and that cooperative fusion defaults to `rho = 0` unless the user changes it.
4. **Correct the release-engineering paragraph**, which cites a stale v0.1.0
   receipt (`FAIL 0 | WARN 0 | SKIP 1 | PASS 1611`, `Status: OK`) and a large
   external `data.zip` fixture that is not bundled.
5. Cross-check `STAR_METHODS.md`, `FIGURE_LEGENDS.md`, `HIGHLIGHTS_AND_ETOC.md`,
   `KEY_RESOURCES_TABLE.tsv`, `ARTIFACT_MANIFEST.md` (verifying every named
   artifact exists), `COVER_LETTER.md`, and `RESOURCE_AVAILABILITY.md` against
   the claim ledger.
6. Complete the open items in `SUBMISSION_CHECKLIST.md`: author list,
   affiliations, ORCIDs, funding, acknowledgments, declaration of interests, and
   target article type.
7. Decide the DIABLO and TCGA question explicitly.
   `inst/analysis/BENCHMARK_PLAN.md` has no implementation, no runner, and no
   cache; `vignettes/stablr-tcga-cache/` is empty. Either implement it as a
   bounded, honest comparison or scope it out permanently and say so in
   Limitations.

## Phase 7 — Decision record, contract v2, and the full evidence chain

1. Produce the SCI-04 decision record conforming to
   `inst/analysis/contracts/sci-decision-v2.schema.json` and the OOL decision
   conforming to `stablr-experiments/analysis/contracts/ool-decision.schema.json`.
   Each must be hash-bound to its complete review bundle, signed, and carry
   truthful provenance fields for `reviewer_id`, independence, and key custodian.
2. Create **contract v2**, additive, leaving v1 immutable. Remove the SCI-04
   blocker only when the decision record and its signature verify. `promotable`
   may become `true`, but the contract must carry an explicit `claim_scope` so no
   downstream consumer can read it as a universal FDP+ calibration guarantee.
3. Update `NEWS.md`, `README.md`, `docs/release/*`, and `cran-comments.md` to the
   true post-Phase-1, 2, and 3 state. Every non-claim in
   `docs/handoffs/2026-09-20-v0.1.1-completion-handoff.md` must be re-derived
   rather than copied.
4. Then, and only then, the evidence chain: build the candidate from the final
   `main` SHA; `prepare`; `run methodology` at the SCI-04 v3 full profile of 600
   replicates; `run late_fusion`; `assemble`; the OOL v3 full profile against the
   identical tarball SHA; the exact-candidate platform matrix; and
   `validate_promotion_record.R`.
5. A completed failed gate is sealed negative evidence and stops the release. A
   crash stays unsealed and may be retried only with identical settings in a new
   directory.

## Verification commands

```sh
RT=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin

"$RT/Rscript" --vanilla -e 'devtools::document(); devtools::test()'
"$RT/R" CMD build .
"$RT/R" CMD check --as-cran stablr_0.1.1.tar.gz
"$RT/Rscript" --vanilla -e 'covr::package_coverage(clean = FALSE)'
```

Run each phase against a freshly built tarball. Bounded Run results are never
promotable.

## Stop rules

- Any new failing parity fixture after a Phase 2 fix halts that fix.
- If Phase 1 H3 holds, the historical Spearman figure is withdrawn and the
  manuscript is rewritten before any other publication work.
- If Phase 3's bounded calibration study contradicts the FDP+ claim, the claim is
  corrected. The contract is renegotiated, not the threshold.
- If the 9600-fit methodology profile is infeasible, the contract is amended and
  the amendment is reviewed like any other change.
- No CRAN submission, tag, GitHub release, or publication claim before the outer
  promotion record verifies.

## Explicit non-claims

- This plan does not assert that FDP+ controls the FDP in this implementation.
  Phase 3 determines that.
- It does not treat the committed gate baseline as evidence in either direction,
  because `n_bootstraps = 2` and `n_lambda = 2` make it uninformative.
- It does not treat the agent's decision record as independent statistical
  review, and the artifacts will say so.
- It does not assert that the OOL parity numbers are reconcilable. Phase 1 may
  conclude that the historical figure was measured against a different quantity.
