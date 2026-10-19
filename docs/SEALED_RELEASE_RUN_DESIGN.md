# Sealed full-candidate validation run — execution design

Status: **ready for execution once the owner gates it** (2026-10-19, repo at
`2899a2a`). This document walks through the fresh full-candidate release run
end to end: what is sealed, what runs, which gates evaluate, and what happens
on failure. It is a design walkthrough, not an execution record — nothing
here has been run, and running it is a human-gated step.

Scope note: "full OOL" refers to the predeclared full Onset-of-Labor parity
analysis executed from the `stablr-experiments` companion repository
(`analysis/generate_publication_parity_table.R --profile full`), distinct
from the bundled `load_ool_data()` quick-build tutorial subset. The release
evidence tracks (methodology + late-fusion simulations) execute through this
package's own sealed machinery (`inst/analysis/release_evidence_cli.R`).

---

## 1. What the run must produce

A **Release Evidence Packet** assembled from content-addressed,
hash-verified validation packets, normalized fail-closed by
`inst/analysis/release_evidence.R`. Only a packet in which every required
validation passed under the exact contract settings lets the release
contract drop `fresh-release-run-pending`. The contract
(`inst/release/release-contract.json`) currently pins `promotable: false`
with two blockers: `sci04-independent-confirmation-pending` and
`fresh-release-run-pending`.

## 2. Prerequisites — all must complete before the run

1. **SCI-04 independent confirmation.** The itemized confirmation per
   `inst/analysis/review/SCI-04-independent-confirmation-packet-2026-10-19.md`
   (regenerated 2026-10-19 at `2899a2a`). The confirmation flip series
   updates: the release contract (`null_selected_fraction` gate version ->
   the confirmed record id; blocker removed), `release_evidence.R`'s
   accepted-id check, and the runner's `.selected_fraction_gate_version`.
2. **Pre-flight fix — release-profile replicate count (found during this
   design review).** `.release_profile_settings()` in
   `inst/analysis/run_methodology_validation.R` (line ~63) still locks
   `replicates = 100L`, and `test-methodology-validation-runner.R`
   (line ~142) pins `100L`, while the contract and
   `.release_expected_contract()` pin **150** (raised at the 2026-10-16
   provisional acceptance). The release-profile branch overrides the
   adapter-supplied value with the locked one, and the candidate-only
   guard checks seed/target_fdp/workers but **not** replicates — so as the
   repo stands, the sealed run would execute at n = 100, where the SCI-04
   Hoeffding margin alone is 0.1224 > 0.10: the `null_selected_fraction`
   gate fails vacuously and the packet is non-promotable. This is
   fail-closed (no false pass is possible), but the run cannot succeed
   until the locked profile and its test pin are raised to the confirmed
   `n` in the confirmation-flip commit series. The late-fusion runner has
   no such gap (candidate-only enforces exactly 50/20/500/220711).
3. **Contract candidate identity.** The embedded contract names candidate
   `stablr` version `0.1.1` (`contract_id: stablr-0.1.1-release-v1`). The
   1.0.0 candidate requires a contract revision — new `contract_id`,
   `candidate.version = "1.0.0"`, settings re-pinned at the confirmed
   values — recorded as an immutable governance change in the same series
   as the DESCRIPTION bump to `1.0.0`. The DESCRIPTION bump happens only
   here, immediately before cutting the candidate tarball.
4. **Runtime pin.** An exact, executable `Rscript` path. Preparation
   refuses to run under any other runtime than the supplied one
   (`prepare_release_context()` compares against `R.home("bin")`), and the
   contract's `runtime.environment` pins apply: `TZ=UTC`, `LC_ALL=C`,
   `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`,
   `R_PROFILE_USER=R_ENVIRON_USER=/dev/null`.
5. **Dependency closure.** The contract's `runtime.dependency_roots`
   (jsonlite, glmnet, Matrix, matrixStats, RColorBrewer, Rcpp, RcppEigen,
   survival, knockoff) must be installed under the pinned runtime;
   preparation copies the closure into the sealed store library.
6. **`stablr-experiments` checkout** at its pinned revision for the OOL
   parity runner, on a machine with the pinned runtime.

## 3. Execution sequence

### Step 0 — cut the candidate

From a clean, tagged tree at the release commit (after the Section 2
flips): `R CMD build .` -> `stablr_1.0.0.tar.gz`; record its sha256. The
tarball — not a mutable checkout — is the object under validation
("mutable-checkout execution does not qualify",
`docs/RELEASE_STATUS.md`).

### Step 1 — prepare the sealed context

```sh
<pinned-Rscript> --vanilla inst/analysis/release_evidence_cli.R prepare \
  --candidate /absolute/path/stablr_1.0.0.tar.gz \
  --runtime <pinned-Rscript> \
  --store /absolute/new/path/release-store \
  --out-ref /absolute/new/path/release-context.json
```

Preparation runs in a fresh R process with no preloaded stablr namespace;
hashes the tarball before/during/after adoption into the
content-addressed store; extracts and hash-verifies the embedded release
contract; copies the dependency closure; sanitizes the environment; and
installs the candidate into the store's private library. Any identity
mismatch aborts with a typed error and discards the staged attempt.

### Step 2 — methodology validation (release profile)

```sh
<pinned-Rscript> --vanilla inst/analysis/release_evidence_cli.R run \
  --context-ref /absolute/path/release-context.json \
  --validation-id methodology \
  --out-ref /absolute/new/path/packet-methodology.json
```

Executes the candidate's own `analysis/run_methodology_validation.R`
under the contract settings: 150 replicates x {gaussian, binomial,
multinomial, cox} x {random_permutation, knockoff, knockoff_equi,
knockoff_mvr} x 6 scenarios, 1,000 bootstraps, 30 lambdas, seed 270627,
32 workers. This is the large track — plan for a 32-worker HPC node and
hours-scale wall time; the paired design gives every artificial-feature
strategy the same data and fit seed within a replicate. Seven scientific
gates evaluate (all currently `accepted` in the contract):
`cell_completeness`, `python_metrics_parity`, `null_select_any`,
`null_selected_fraction` (SCI-04, Hoeffding upper bound <= 0.10 at the
confirmed n), `null_select_90pct`, `signal_mean_fdp`, `signal_tpr`.

### Step 3 — late-fusion validation

```sh
<pinned-Rscript> --vanilla inst/analysis/release_evidence_cli.R run \
  --context-ref /absolute/path/release-context.json \
  --validation-id late_fusion \
  --out-ref /absolute/new/path/packet-late-fusion.json
```

50 replicates, 20 bootstraps, 500 stacking iterations, seed 220711
(candidate-only enforced). Four gates: `cell_completeness`,
`reduced_optimism`, `noninferiority` (0.02 margin), `fallback_rate`
(0.05).

### Step 4 — full OOL parity (from `stablr-experiments`)

```sh
# From the pinned stablr-experiments checkout; the output directory must
# not pre-exist.
<pinned-Rscript> --vanilla analysis/generate_publication_parity_table.R \
  --candidate-tarball /absolute/path/stablr_1.0.0.tar.gz \
  --runtime <pinned-Rscript> \
  --outdir /absolute/new/path/ool-candidate-evidence \
  --profile full
```

Predeclared gates (the 0.1.1 candidate's sealed negative result is the
reference: support Jaccard 0.778, score Spearman 0.785, FDP-threshold
delta 0.08): support Jaccard >= 0.80, score Spearman >= 0.90,
FDP-threshold delta <= 0.05, plus the lambda and training-only
preprocessing gates and the actual-generator gate. Note for the 1.0.0
candidate: the 0.1.1 run's requested fixed knockoff construction was
infeasible at 150 x 1,317 (`n > p` required) and fell back, typed and
recorded, to `random_permutation`, failing the actual-generator gate —
the 1.0.0 predeclaration should state the expected generator for the
full OOL dimensions explicitly so the gate evaluates the intended
construction.

### Step 5 — assemble the Release Evidence Packet

```sh
<pinned-Rscript> --vanilla inst/analysis/release_evidence_cli.R assemble \
  --context-ref /absolute/path/release-context.json \
  --packet-ref /absolute/path/packet-methodology.json \
  --packet-ref /absolute/path/packet-late-fusion.json \
  --out-ref /absolute/new/path/release-evidence-packet.json
```

Assembly normalizes fail-closed: any incomplete cell, any gate whose
decision status is not exactly the accepted record, any schema or hash
mismatch -> the packet records failure and `promotable` stays `false`.

### Step 6 — record

On pass: commit the packet references and remove
`fresh-release-run-pending` from the contract in the same series. On any
failure: retain the output directory untouched as the sealed negative
artifact (the 0.1.1 precedent), record the failing gates in
`docs/RELEASE_STATUS.md`, and keep both blockers.

## 4. What is verified already vs. what the run decides

Already verified statically (this design review, repo at `2899a2a`):
contract/enforcement/runner settings agreement (except the Section 2.2
replicate gap), fail-closed paths in preparation and normalization, the
sealed RE-01 v1 fixture exclusion, and the SCI-04 reference arithmetic
(30/30 reference-test expectations pass at HEAD).

Only the sealed run decides: the seven methodology gates and four
late-fusion gates at release scale under the candidate tarball, and the
full OOL parity gates against the predeclared thresholds. No number from
a mutable-checkout run, a bounded profile, or the sealed RE-01 v1
baseline may enter the packet.
