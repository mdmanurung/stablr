# stablr Artifact-First Release Evidence Module Implementation Plan

Date: 2026-07-31

Last updated: 2026-07-31

Project: `/exports/para-lipg-hpc/mdmanurung/stablr`

Status: in progress

## Objective

Implement the accepted three-call Release Evidence Module Interface for local,
artifact-first scientific validation:

```r
context <- prepare_release_context(candidate_tarball, runtime, store)
packet  <- run_required_validation(context, validation_id)
release <- assemble_release_evidence(context, packet_refs)
```

The work centralizes candidate/runtime identity, immutable packet lifecycle,
fail-closed sealing, and exact aggregation while preserving the methodology and
late-fusion Adapters' ownership of scientific calculations and Scientific
Gates.

This plan is architecture-scoped. Completion will make release evidence
trustworthy; it will not make the current branch publication-ready. The known
scientific defects and failed locked methodology gates remain independent
release blockers and must not be waived by this refactor.

Accepted design: [Artifact-First Release Evidence Interface](../architecture/release-evidence-interface.md)

## Review Findings

| Finding | Resolution |
|---|---|
| The current runners duplicate Git discovery, package loading, hashing, manifest writing, and release enforcement. | Move those responsibilities behind one Release Evidence Module and delete the duplicate release paths after both Adapters migrate. |
| Artifact-first execution has a bootstrap risk if checkout code launches the validator. | Invoke the thin CLI from a fresh extraction of the Candidate Source Artifact; install that same tarball into an isolated library and verify the loaded namespace path. |
| A caller can currently select release settings and obtains different failure semantics from direct and CLI entry points. | Put all promotable settings and Required Validations in the embedded Release Contract; expose only `validation_id` at run time and use one typed result path. |
| Scientific failure must remain diagnosable while still stopping release. | Seal a failed Validation Run Packet on normal completion, then signal `stablr_validation_failed` with its Packet Reference. |
| Crashes and kills cannot be mistaken for failed scientific evidence. | Keep abnormal attempts unsealed and ineligible; never synthesize a Packet Reference. |
| Independent Slurm jobs cannot share in-memory state safely. | Make Release Contexts and Packet References immutable, versioned JSON records and require explicit references during assembly. |
| Runtime identity was underspecified. | Define a v1 Runtime Artifact covering the `Rscript` and companion R executables, R home/version, platform, locale, BLAS/LAPACK, isolated dependency closure, installed candidate tree, and resolved native libraries; re-verify it before and after a run. |
| Directory scanning can silently select stale evidence. | Make packet discovery impossible through the Interface; assembly accepts only explicit Packet References and rejects duplicates, extras, and foreign identities. |
| The real release profiles are too expensive for routine tests. | Use a fixed full-branch bounded methodology matrix and a six-cell tiny late-fusion run, repeat them for determinism, and reserve full locked runs for a later exact-candidate release stage. |
| The late-fusion runner rejects the tracked-dirty main checkout before producing scientific artifacts. | Assert that production inputs match HEAD, clone that exact HEAD locally into a fresh clean temporary repository, and run both characterizations there without weakening the guard. |
| Refactoring the runners could change scientific results or RNG behavior. | Freeze normalized golden artifacts, repeat serial A/B runs byte-for-byte, compare serial and parallel scientific tables without sorting, and enumerate every normalized field. |
| Temporary characterization output is too easy to overwrite or lose. | Retain raw outputs, logs, runtime identity, comparisons, and manifests in a unique ignored local baseline directory; publish its SHA-256 manifest marker-last and label it non-promotable. |
| JSON support is not currently declared. | Add `jsonlite` as a release-tooling `Suggests` dependency, require it through the Release Contract, record its exact Runtime Artifact identity, and fail closed when unavailable. |
| Existing dirty and untracked work belongs to the user. | Hash the tracked diff and every protected untracked/ignored file before RE-01, compare them afterward, and limit mutations to an explicit allowlist. |
| Architecture completion could be mistaken for scientific acceptance. | End with a publication stop rule and keep the current scientific blockers visible in final evidence and documentation. |

All plan-review findings are addressed. No unresolved design question blocks
implementation.

## Progress Summary

- [ ] RE-01: Freeze current behavior and safety baseline
- [ ] RE-02: Add the versioned Release Contract and serializable Interface types
- [ ] RE-03: Implement Candidate Source Artifact and Runtime Artifact preparation
- [ ] RE-04: Implement the immutable Validation Attempt and packet state machine
- [ ] RE-05: Migrate methodology validation into its scientific Adapter
- [ ] RE-06: Migrate late-fusion validation into its scientific Adapter
- [ ] RE-07: Implement fresh-process Required Validation execution
- [ ] RE-08: Implement explicit-reference Release Evidence Packet assembly
- [ ] RE-09: Replace legacy release entry points with the thin CLI and documentation
- [ ] RE-10: Validate the complete Module without claiming scientific release success

## Implementation Steps

### [x] RE-01: Freeze current behavior and safety baseline

- **Status:** complete — sealed non-promotable baseline at the reviewed source
  identity with an isolated `Rdsdp` 1.0.6 plus `knockoff` 0.3.6 closure
- **Outcome:** A reproducible pre-refactor baseline records current tests,
  bounded Adapter schemas, deterministic outputs, known failures, and the dirty
  worktree preservation set.
- **Actions:**
  1. Require HEAD
     `c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1`. Before any run, prove that
     `DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, and `inst/analysis/` have no
     tracked difference from that HEAD.
  2. Create the new, no-replace raw baseline directory
     `.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/`. Record porcelain-v2
     status, the tracked binary diff, source hashes, and a path/type/mode/size/
     SHA-256 manifest for the preservation-only set: `.gitignore`, `CONTEXT.md`,
     `docs/architecture/`, `docs/release/0.1.1-next-tasks.md`, `papers/`, and
     `.claude/reviews/stabl_lw/summary.md`.
  3. Limit RE-01 writes to the raw baseline directory,
     `tests/testthat/helper-release-evidence-baseline.R`,
     `tests/testthat/test-release-evidence-baseline-fixtures.R`,
     `tests/testthat/fixtures/release_evidence_baseline/v1/`, and RE-01 status/
     evidence entries in this plan. The execution task list and handoff are
     frozen packet files, not implementation output.
  4. Resolve the R4_51 `Rscript` real path and SHA-256. Record the companion R
     identity, `R.version.string`, `sessionInfo()`, locale, timezone,
     BLAS/LAPACK, `.libPaths()`, and exact versions and library paths of loaded
     scientific/test dependencies. Require `pkgload`, `survival`, and `knockoff`.
     For authorized continuation 1 only, preserve the stopped runtime record
     and append a new epoch beneath
     `.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/continuations/01-cran-knockoff/`.
     Retrieve only CRAN `PACKAGES.gz`, `Rdsdp_1.0.6.tar.gz`, and
     `knockoff_0.3.6.tar.gz`; verify index MD5, pinned SHA-256, and internal
     package identity; then install only those two packages with the companion
     R4_51 `R CMD INSTALL` into a new attempt-local staging library. Require
     SHA-256 values
     `dd462374f38ae45791c024805582e1275936cb96c2127fd2dceeb3df0f18bfdd`
     and
     `afec4c2bdb8aad4ab4a7f33527a546deda80229588ce7a4d60d37114d3b046ba`,
     respectively. Promote the staging library only after provenance, native
     library, import-resolution, SDP smoke, and unchanged-R4_51-library checks
     pass; then set its absolute path first in `R_LIBS_USER` for every resumed
     R process. No conda/home library mutation, automatic dependency install,
     other network resource, historical-evidence rewrite, cleanup, or retry is
     authorized.
  5. Run commands under `TZ=UTC`, `LC_ALL=C`, `OMP_NUM_THREADS=1`,
     `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
     `R_PROFILE_USER=/dev/null`, and `R_ENVIRON_USER=/dev/null` with
     `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla`.
  6. Run the full source test suite in the main checkout before fixture changes,
     capturing the complete log and exit status.
  7. Create a fresh local temporary clone at the exact required HEAD, verify it
     is clean, and run methodology serial A, serial B, and parallel C into three
     distinct new output directories with these fixed settings:
     `profile="bounded"`, `replicates=2L`, `n_bootstraps=2L`, `n_lambda=2L`,
     `families=c("gaussian","binomial","multinomial","cox")`,
     `artificial_types=c("random_permutation","knockoff","knockoff_equi","knockoff_mvr")`,
     all six default scenarios, `seed=270627L`, `target_fdp=0.1`, and workers
     `1L`, `1L`, and `2L`, respectively. Preserve fit errors and gate failures
     as characterization.
  8. In the same clean clone, run late-fusion A and B through the function
     Interface with `replicates=1L`, `n_bootstraps=2L`, `n_iter=2L`, and
     `seed=220711L`. Require all four artifacts. Accept only the exact expected
     post-artifact condition `Late-fusion release gates failed; release must
     stop.`; any load, fit, provenance, or completeness failure blocks RE-01.
  9. Normalize without sorting or changing row/column order: remove only
     `elapsed_sec` from methodology replicates and `mean_elapsed_sec` from its
     summary; copy methodology warnings, parity, and gates and every late-fusion
     CSV byte-for-byte. In methodology manifests normalize only `created_at`,
     `package_mode`, Git start/end/stability fields, artifact path values, and
     artifact hash values. In late-fusion manifests normalize only
     `package_mode`, Git start/end/stability fields, and artifact hash values.
     Preserve seeds, NA encoding, numeric text, conditions, warnings, gate
     booleans, and ordering.
  10. Require normalized serial A/B artifacts to be byte-identical. Require
      serial/parallel methodology scientific CSVs to be identical after the two
      elapsed columns are removed; record the intentional manifest `workers`
      difference. Do not weaken a mismatch.
  11. Save the normalized serial-A golden artifacts, normalization helper,
      fixture integrity test, README, and SHA-256 manifest under
      `tests/testthat/fixtures/release_evidence_baseline/v1/`. Each golden file
      references the ignored raw baseline identity and is explicitly
      non-promotable.
  12. Run the full source test suite again after fixture creation. Then compare
      the protected-file manifest and tracked diff with the pre-run records.
  13. Write a raw SHA-256 closure manifest, create
      `SEALED.NONPROMOTABLE` last, and make the raw baseline read-only. Record
      exact commands, logs, exits, comparisons, and fixture hashes in this
      plan's evidence log.
  14. Cite `.claude/reviews/stabl_lw/summary.md` by date, reviewed SHA, and
      checksum. Record stable blocker IDs for every unchanged scientific defect,
      including validation-label-controlled prediction shape. Record provenance
      binding and split entry-point enforcement separately as architecture work
      for RE-02 through RE-09.
- **Dependencies:** Accepted architecture record; R4_51; current checkout.
- **Affected areas:** the explicit RE-01 write allowlist above; all production
  code, metadata, runners, other documentation, and user-owned paths are
  preservation-only.
- **Validation:**
  `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla -e 'devtools::test()'`
  before and after fixture creation, normalized A/B byte comparisons,
  serial/parallel scientific-table comparisons, fixture SHA-256 verification,
  and protected pre/post manifest comparison.
- **Acceptance:** Both full suites pass; every required raw and golden artifact
  exists; deterministic comparisons pass; the late-fusion condition is exactly
  the expected post-artifact gate failure; the raw baseline is sealed and
  non-promotable; protected content is unchanged; and no production code,
  metadata, runner, or documentation outside the explicit allowlist changes.
- **Evidence:** The historical stopped epochs and their exact source archives
  remain unchanged. Continuation `04-r4-toolchain` used the R4_51 compiler and
  binutils to install the verified `Rdsdp` 1.0.6 and `knockoff` 0.3.6 sources
  only into its attempt-local library. Package/import identity, native linkage,
  and the two-value SDP smoke passed. Complete 62,996-path R4_51 manifests were
  byte-identical before installation, immediately after installation, and at
  closure. Runtime identity was captured under the pinned environment. The
  pre-fixture suite passed with 1,941 passes, zero failures or warnings, and the
  same eight environment skips; the post-fixture suite passed with 1,966
  passes, zero failures or warnings, and those same skips. Methodology serial A,
  serial B, and two-worker C each produced all six bounded artifacts. After
  removing only the two elapsed fields and normalizing only the declared
  manifest values, A/B were byte-identical, serial/parallel scientific CSVs
  were byte-identical, and the sole manifest difference was `workers: 1`
  versus `workers: 2`. Late-fusion A/B each produced all four artifacts and the
  exact expected condition `Late-fusion release gates failed; release must
  stop.`; normalized A/B artifacts were byte-identical. Ten normalized serial-A
  and late-fusion-A fixtures, their SHA-256 manifest, helper, README, and focused
  integrity tests are retained under
  `tests/testthat/fixtures/release_evidence_baseline/v1/` and explicitly marked
  non-promotable. The ignored raw packet is sealed marker-last at
  `.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/SEALED.NONPROMOTABLE` after
  production-input, protected-file, repository-state, dependency-library, and
  historical-evidence preservation checks. Linus review SHA-256 remains
  `849e3ff94adde90648f85e0a36283d91c158f8564f0dbd3c698a34890c05462c`.
  `SCI-01` through `SCI-10`, `ARCH-01`, and `ARCH-02` remain open after this
  characterization baseline.

### [ ] RE-02: Add the versioned Release Contract and serializable Interface types

- **Status:** not started
- **Outcome:** The Candidate Source Artifact contains one strict Release
  Contract, and Release Contexts and references have validated versioned JSON
  schemas and typed conditions.
- **Actions:**
  1. Add `jsonlite` to `Suggests` and declare it as a required release-tooling
     runtime dependency in the contract.
  2. Add `inst/release/release-contract.json` with schema version, candidate
     package/version constraints, Required Validation IDs, fixed settings,
     Adapter keys, artifact schemas, Scientific Gate versions, and Runtime
     Artifact requirements.
  3. Add JSON schemas under `inst/release/schemas/` for the Release Contract,
     Release Context, Runtime Artifact, Validation Run Packet, Packet Reference,
     and Release Evidence Packet.
  4. Start `inst/analysis/release_evidence.R` with strict parsers, canonical JSON
     writing, SHA-256 helpers, immutable reference classes, and the accepted
     typed conditions.
  5. Reject unknown schema versions, missing/extra Required Validations,
     duplicate IDs, caller settings, unknown Adapter keys, relative runtime
     paths, missing hashes, and unsupported algorithms.
- **Dependencies:** RE-01.
- **Affected areas:** `DESCRIPTION`, `inst/release/`,
  `inst/analysis/release_evidence.R`,
  `tests/testthat/test-release-evidence-contract.R`.
- **Validation:** Unit tests round-trip every record, mutate each required
  field, attempt to weaken settings, and verify the exact typed error class.
  `R CMD build` inspection confirms the contract and schemas are inside the
  tarball.
- **Acceptance:** The caller cannot add, remove, replace, or downscale a
  Required Validation through any Interface argument or serialized record.
- **Evidence:** Pending

### [ ] RE-03: Implement Candidate Source Artifact and Runtime Artifact preparation

- **Status:** not started
- **Outcome:** `prepare_release_context()` returns a re-verifiable context bound
  to a store-owned tarball, isolated installed package, and complete Runtime
  Artifact.
- **Actions:**
  1. Require an absolute regular tarball path, compute SHA-256, copy it into a
     content-addressed store location, and re-hash the copy.
  2. Inspect archive members before extraction; reject absolute paths, `..`,
     links escaping the archive root, multiple package roots, wrong package
     names, or a missing Release Contract.
  3. Resolve every contract-declared release-tooling and scientific dependency
     from the pinned runtime, copy its exact installed tree into a new
     context-local library, and hash the copy. Perform no network resolution or
     install from an unrecorded library.
  4. Install only the store-owned tarball into that context-local library using
     the companion R executable derived from the supplied absolute `Rscript`.
  5. Execute with `--vanilla` and sanitized `R_PROFILE*`, `R_ENVIRON*`, and
     library variables. Verify `system.file(package = "stablr")` resolves inside
     the isolated library.
  6. Construct the Runtime Artifact from the `Rscript` and companion R
     executable hashes, R home/version, platform, locale, BLAS/LAPACK, installed
     candidate tree, copied R package closure, and resolved native libraries.
     Fail on unresolved or unhashable members.
  7. Probe the store's same-filesystem atomic rename and exclusive marker
     creation behavior before accepting it.
  8. Persist a canonical context, make prepared inputs read-only where possible,
     and return only a verified `stablr_release_context` reference.
- **Dependencies:** RE-02.
- **Affected areas:** `inst/analysis/release_evidence.R`, a thin bootstrap path
  in `inst/analysis/run_release_evidence.R`, and
  `tests/testthat/test-release-evidence-prepare.R`.
- **Validation:** Tests cover tarball mutation after preparation, traversal
  entries, missing contract, wrong package, wrong/preloaded namespace,
  unavailable hash tools, missing declared dependencies, dependency drift,
  attempted network resolution, native-library drift, unsupported platform
  discovery, and an atomicity-probe failure.
- **Acceptance:** No promotable context can refer to checkout code, a caller-owned
  mutable tarball, a global user library, an unidentified runtime member, or a
  preloaded namespace.
- **Evidence:** Pending

### [ ] RE-04: Implement the immutable Validation Attempt and packet state machine

- **Status:** not started
- **Outcome:** Normal completion seals exactly one passing or failed Validation
  Run Packet; abnormal completion remains unsealed; retries never mutate prior
  evidence.
- **Actions:**
  1. Allocate unique Validation Attempts under the store without deriving
     identity from timestamps or output-directory discovery.
  2. Restrict Adapter outputs to declared relative regular files inside the
     attempt; reject symlinks, escapes, duplicate names, undeclared files, and
     post-enumeration mutation.
  3. Define the common gate table and derive status only after artifact schema,
     completeness, and checksum verification.
  4. Write packet metadata and the complete checksum manifest, atomically move
     the finalized directory, then create the seal marker last with no-replace
     semantics.
  5. Make sealed packets read-only where possible and verify them again before
     returning a Packet Reference.
  6. On complete failed gates, seal the failed packet and signal
     `stablr_validation_failed` with its Packet Reference. On interruption or
     unexpected error, retain an unsealed attempt and signal
     `stablr_validation_incomplete` without a Packet Reference.
  7. Ensure a retry always receives a new attempt and packet identity.
- **Dependencies:** RE-03.
- **Affected areas:** `inst/analysis/release_evidence.R`,
  `tests/testthat/test-release-evidence-packets.R`.
- **Validation:** State-transition, kill/crash, failed-gate, tampering,
  retry-immutability, duplicate-output, symlink, and concurrent-writer tests.
- **Acceptance:** Only a seal-last, checksum-valid, schema-valid packet can
  produce a Packet Reference; failed packets remain inspectable but ineligible
  for assembly.
- **Evidence:** Pending

### [ ] RE-05: Migrate methodology validation into its scientific Adapter

- **Status:** not started
- **Outcome:** Methodology scientific behavior lives in one Adapter; generic
  provenance, package loading, hashing, sealing, and release exit logic are
  deleted from its Implementation.
- **Actions:**
  1. Move simulation, paired task execution, summaries, warning artifacts,
     parity observations, and Scientific Gate calculations into
     `inst/analysis/adapters/methodology_validation.R`.
  2. Make the Adapter accept only the contract-derived settings and attempt
     path supplied privately by the Module.
  3. Return declared artifact paths, Adapter metadata, and one complete common
     gate table. Represent completeness and Python parity as explicit gates so
     the Module sees one status source.
  4. Keep scientific formulas, thresholds, seed construction, row order, and
     RNG behavior unchanged during extraction. Do not bless the known invalid
     selected-fraction Wilson gate; retain it as a separately tracked scientific
     blocker until a reviewed methodology change replaces its gate version.
  5. Retain a clearly labeled Bounded Run entry for development overrides; it
     must set `promotable = false` and cannot call packet sealing.
  6. Remove Git discovery, `pkgload::load_all()`, preloaded-namespace reuse,
     generic SHA helpers, manifest code, and the split direct-versus-CLI release
     decision from `run_methodology_validation.R`.
- **Dependencies:** RE-01, RE-02, RE-04.
- **Affected areas:** `inst/analysis/run_methodology_validation.R`, new
  `inst/analysis/adapters/methodology_validation.R`,
  `tests/testthat/test-methodology-validation-runner.R`, and Adapter schema
  fixtures.
- **Validation:** The RE-01 bounded fixture matches exactly after excluding only
  lifecycle-owned fields. Tests prove failed gates return normally to the
  Module, parity participates in the gate table, release settings cannot be
  overridden, and Bounded Runs cannot create Packet References.
- **Acceptance:** Deleting the Release Evidence Module would recreate generic
  lifecycle complexity in this Adapter; deleting the Adapter would remove only
  methodology-specific science. This passes the deletion test for a deep Module.
- **Evidence:** Pending

### [ ] RE-06: Migrate late-fusion validation into its scientific Adapter

- **Status:** not started
- **Outcome:** Late-fusion scientific behavior becomes the second real Adapter
  at the Required Validation Seam and shares the same lifecycle semantics.
- **Actions:**
  1. Move simulations, paired legacy/OOF fits, summaries, warnings, and the
     reduced-optimism, noninferiority, and fallback Scientific Gates into
     `inst/analysis/adapters/late_fusion_validation.R`.
  2. Return declared artifacts, metadata, and complete common gate rows rather
     than writing a generic manifest or stopping on failed gates.
  3. Remove checkout discovery, Git provenance, `pkgload::load_all()`, preloaded
     namespace reuse, hashing, and sealing logic from
     `run_late_fusion_validation.R`.
  4. Put the 50-replicate release settings in the Release Contract. Preserve a
     separate non-promotable Bounded Run path for tiny development runs.
  5. Preserve scientific values, seeds, row ordering, warnings, and current
     formulas during extraction; track known late-fusion scientific defects
     separately rather than hiding them in the architecture change.
- **Dependencies:** RE-01, RE-02, RE-04.
- **Affected areas:** `inst/analysis/run_late_fusion_validation.R`, new
  `inst/analysis/adapters/late_fusion_validation.R`,
  `tests/testthat/test-methodology-validation-runner.R`, and a new focused
  late-fusion Adapter test file.
- **Validation:** The RE-01 bounded fixture matches after excluding only
  lifecycle fields. Tests prove all three gate families reach the common table,
  failed gates do not bypass packet sealing, release settings cannot be
  overridden, and the Adapter loads the isolated candidate namespace.
- **Acceptance:** Both real Adapters satisfy one private Interface without
  duplicating lifecycle logic.
- **Evidence:** Pending

### [ ] RE-07: Implement fresh-process Required Validation execution

- **Status:** not started
- **Outcome:** `run_required_validation()` selects only the contract-declared
  Adapter, executes it in a fresh artifact-bound process, and feeds its normal
  completion into the packet state machine.
- **Actions:**
  1. Re-verify the Release Context, candidate, runtime, installed tree, and
     contract immediately before allocating the attempt.
  2. Resolve the Adapter through a private fixed registry keyed only by the
     Release Contract; expose no caller registration or dependency injection.
  3. Spawn the context runtime with `Rscript --vanilla`, sanitized profiles and
     libraries, an explicit serialized job descriptor, and the isolated
     candidate library.
  4. In the child, verify the candidate namespace and Adapter paths, execute the
     Adapter, and write a completion record only after its artifacts and gate
     table are closed.
  5. In the parent, distinguish normal passing/failed completion from abnormal
     exit, re-verify runtime identity, then call the RE-04 state transition.
  6. Return a serializable Packet Reference for passing runs; preserve the
     failed reference on `stablr_validation_failed`.
- **Dependencies:** RE-03 through RE-06.
- **Affected areas:** `inst/analysis/release_evidence.R`, Adapter child runner,
  `tests/testthat/test-release-evidence-run.R`.
- **Validation:** Tests cover wrong Adapter IDs, caller injection attempts,
  polluted parent namespaces, user/site profile contamination, child nonzero
  exit, kill without completion, runtime mutation, pass, scientific failure,
  and serialization across independent R processes.
- **Acceptance:** Every promotable validation process is fresh, candidate-bound,
  runtime-bound, contract-selected, and packetized through one Interface.
- **Evidence:** Pending

### [ ] RE-08: Implement explicit-reference Release Evidence Packet assembly

- **Status:** not started
- **Outcome:** `assemble_release_evidence()` creates one immutable aggregate only
  from the exact passing packet set required by the Release Contract.
- **Actions:**
  1. Accept only explicit Packet References; do not list or scan packet
     directories.
  2. Verify reference and manifest hashes, seals, schemas, status, candidate,
     runtime, contract, gate-table completeness, and Required Validation ID.
  3. Require exactly one packet per Required Validation and reject duplicates,
     extras, missing IDs, failed packets, unsealed attempts, Bounded Runs, and
     foreign candidate/runtime/contract identities.
  4. Build the aggregate manifest from accepted references without
     recalculating Adapter formulas.
  5. Publish the Release Evidence Packet with the same checksum, atomic move,
     no-replace, and marker-last rules; return `stablr_release_ref`.
- **Dependencies:** RE-04 and RE-07.
- **Affected areas:** `inst/analysis/release_evidence.R`,
  `tests/testthat/test-release-evidence-assemble.R`.
- **Validation:** Exhaustive packet-set tests plus post-seal tampering,
  cross-candidate, cross-runtime, stale-contract, duplicate-ID, and directory
  bait tests.
- **Acceptance:** No implicit, failed, incomplete, bounded, stale, or foreign
  evidence can enter a Release Evidence Packet.
- **Evidence:** Pending

### [ ] RE-09: Replace legacy release entry points with the thin CLI and documentation

- **Status:** not started
- **Outcome:** One CLI mirrors the Interface; legacy runner commands are either
  explicit Bounded Runs or route release execution through the new Module.
- **Actions:**
  1. Complete `inst/analysis/run_release_evidence.R` with `prepare`, `run`, and
     `assemble` subcommands. Write machine-readable references to stdout and
     diagnostics to stderr.
  2. Document the artifact-derived bootstrap: freshly extract the Candidate
     Source Artifact and invoke the CLI contained in that extraction, which
     installs and verifies the same store-owned tarball.
  3. Remove or hard-error the old `--profile release` execution path unless it
     delegates to the contract-bound Module. Keep bounded overrides visibly
     non-promotable.
  4. Update `inst/analysis/README.md` and
     `docs/release/0.1.1-validation.md` with packet semantics, typed failures,
     explicit assembly, and the distinction between Bounded Runs and release
     evidence.
  5. Reconcile `docs/release/0.1.1-next-tasks.md` only after re-reading its
     user-owned untracked content; patch relevant commands without overwriting
     unrelated restart history.
  6. Add CI tests for schemas and the tiny lifecycle fixture, not the long
     scientific release profiles.
- **Dependencies:** RE-05 through RE-08.
- **Affected areas:** `inst/analysis/run_release_evidence.R`, both legacy
  runners, `inst/analysis/README.md`, release documentation, and CI test
  selection.
- **Validation:** CLI function-parity tests compare R-call and CLI references
  and typed failures. Search confirms no remaining promotable Git-checkout,
  `pkgload`, preloaded-namespace, directory-discovery, or split-gate path.
- **Acceptance:** There is one release lifecycle Implementation and one thin
  caller Interface; old commands cannot emit promotable evidence independently.
- **Evidence:** Pending

### [ ] RE-10: Validate the complete Module without claiming scientific release success

- **Status:** not started
- **Outcome:** The architecture is verified through unit, adversarial,
  artifact-build, and small end-to-end tests while full scientific release
  validation remains blocked until its separate defects are repaired.
- **Actions:**
  1. Run focused contract, preparation, packet, Adapter, execution, assembly,
     and CLI tests.
  2. Run the full source test suite with R4_51.
  3. Build the source tarball in an isolated temporary directory and inspect it
     for the Release Contract, schemas, Module, CLI, and both Adapters; verify it
     still excludes `docs/`, `papers/`, build products, and user artifacts.
  4. Run installed-package tests and `R CMD check --as-cran`, separating package
     failures from documented host limitations.
  5. Execute tiny pass, failed-gate, and killed-process fixture candidates from
     their built tarballs. Prove pass assembly succeeds, failed evidence seals
     but cannot assemble, and killed evidence remains unsealed.
  6. Re-run bounded Adapter equivalence comparisons from RE-01.
  7. Audit the final diff for duplicate lifecycle logic, accidental public
     Interface growth, unrelated dirty-file changes, and unsupported readiness
     claims.
  8. Record that the real multi-day Required Validations were not run and that
     a Release Evidence Packet for `stablr` does not yet exist.
- **Dependencies:** RE-01 through RE-09.
- **Affected areas:** all planned implementation and test files; no publication,
  tag, remote, CRAN, or release state.
- **Validation:**
  - `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript -e 'devtools::test()'`
  - isolated `R CMD build .`
  - isolated `R CMD check --as-cran <built-tarball>`
  - fixture lifecycle pass/fail/kill commands recorded in the evidence log
- **Acceptance:** The three-call Interface passes adversarial lifecycle tests,
  both scientific Adapters preserve bounded outputs, the built artifact is the
  only executed package source, and no Release Evidence Packet is claimed for
  the current scientifically blocked candidate.
- **Evidence:** Pending

## Final Validation

Implementation is complete only when all RE-01 through RE-10 acceptance checks
pass and the evidence log records exact commands, R identity, candidate hashes,
test/check outcomes, and any host limitations.

The publication stop rule remains absolute:

1. Architecture tests and package checks are not scientific validation.
2. Bounded Runs are not promotable evidence.
3. The current locked methodology result failed and cannot be reinterpreted.
4. No `stablr` Release Evidence Packet may be assembled until separate reviewed
   fixes close the known scientific defects and every Required Validation passes
   against one exact Candidate Source Artifact and Runtime Artifact.
5. CRAN, CI, Win-builder, R-hub, tag, and release work remains outside this
   local architecture plan.

## Decision and Evidence Log

- 2026-07-31 — Artifact-first, two-layer, fail-closed evidence lifecycle
  accepted.
- 2026-07-31 — Three-call Interface accepted; explicit packet references and
  contract-selected Adapters required.
- 2026-07-31 — Plan reviewed against both current validator Implementations,
  existing tests, release documentation, the Linus review, dirty-worktree
  constraints, and R4_51. Draft patched and initialized; implementation not
  started.
- 2026-07-31 — RE-01 execution review added exact bounded matrices, clean-clone
  execution, A/B determinism, preservation hashes, golden fixtures, and a
  sealed non-promotable raw baseline. RE-02 through RE-10 remain unauthorized
  for this execution.
- 2026-07-31 — RE-01 execution stopped before tests and characterizations at
  the fail-closed runtime gate: required package `knockoff` was unavailable in
  the mandated R4_51 library. No installation or download was authorized; the
  unique raw attempt remains incomplete and unsealed.
- 2026-07-31 — The user authorized a same-attempt, append-only continuation to
  retrieve exactly CRAN `PACKAGES.gz`, `Rdsdp` 1.0.6, and `knockoff` 0.3.6;
  verify their pinned MD5/SHA-256 and internal identities; install the two
  packages into an ignored attempt-local library without mutating R4_51; and
  resume RE-01 only after dependency, native-library, smoke, library-drift, and
  historical-preservation gates pass. All remote publication and RE-02 or later
  work remains unauthorized.

## Risks, Assumptions, and Blockers

- **Scientific blockers:** `SCI-01` classification bootstrap feasibility;
  `SCI-02` named stacking alignment; `SCI-03` fold allocation; `SCI-04`
  selected-fraction gate estimation; `SCI-05` actual artificial-feature
  provenance; `SCI-06` high-dimensional construction; `SCI-07` catch-all
  fallbacks; `SCI-08` OOF provenance field access; `SCI-09` late-fusion result
  semantics; and `SCI-10` validation-label-controlled prediction shape remain
  unresolved. They block publication but do not block this architecture
  implementation. Source: `.claude/reviews/stabl_lw/summary.md`, dated
  2026-07-31 and reviewed at
  `c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1`; RE-01 records its SHA-256.
- **Architecture blockers:** `ARCH-01` executed-code provenance binding and
  `ARCH-02` split programmatic/CLI release enforcement are addressed only by
  RE-02 through RE-09, not by RE-01 characterization.
- **Failed evidence:** The observed 52/240 passing methodology gates and
  1,212/9,600 fit errors remain failed historical evidence. Architecture work
  must preserve, not overwrite or relabel, those artifacts.
- **Scientific invariance:** Moving code between files can perturb RNG or
  warning order. RE-01 fixtures and bounded before/after comparisons are
  mandatory before deleting legacy logic.
- **Runtime closure:** Native-library discovery can be platform-specific. The v1
  Implementation supports the current Linux/POSIX runtime and fails closed
  elsewhere until a second real runtime Adapter justifies a new Seam.
- **Filesystem semantics:** The evidence store is assumed to provide same-mount
  atomic rename and exclusive creation. Preparation must probe and reject a
  store that cannot demonstrate them.
- **Long-running validation:** Full methodology and late-fusion release runs are
  deliberately excluded from routine implementation validation. They occur
  only after separate scientific remediation and an exact candidate freeze.
- **JSON dependency:** Release tooling requires `jsonlite`; normal package use
  should not. The Release Contract and Runtime Artifact make that requirement
  explicit.
- **Dirty worktree:** `.gitignore`, `CONTEXT.md`,
  `docs/release/0.1.1-next-tasks.md`, and `papers/` are pre-existing or
  user-owned work. Every implementation stage must compare status before and
  after and avoid broad staging or cleanup.
- **No remote authority:** This plan is local only. It authorizes no push,
  branch promotion, scheduler submission, CRAN action, tag, or publication.
