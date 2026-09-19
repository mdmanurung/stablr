# No-Parent-History Executor Handoff: Freeze RE-01 Baselines

## Objective

Complete RE-01 only: freeze reproducible, deterministic, non-promotable
methodology and late-fusion characterization fixtures plus exact local evidence,
without changing production code or any protected user-owned file.

## Repository

`/exports/para-lipg-hpc/mdmanurung/stablr`

Revised Plan Path:
`/exports/para-lipg-hpc/mdmanurung/stablr/docs/plans/2026-07-31-release-evidence-module.md`

Task List Path:
`/exports/para-lipg-hpc/mdmanurung/stablr/docs/plans/2026-07-31-release-evidence-re01-tasks.md`

## Frozen Plan

1. Require HEAD
   `c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1`; stop if tracked production
   inputs under `DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, or `inst/analysis/`
   differ from it.
2. Fail closed if
   `.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/` already exists. Create
   it as the unique raw attempt and retain status, binary diff, protected-file
   inventory, commands, logs, exits, raw outputs, normalized outputs,
   comparisons, source hashes, and runtime identity there.
3. Before execution, manifest path, type, mode, size, SHA-256, or symlink target
   for the preservation-only set: `.gitignore`, `CONTEXT.md`,
   `docs/architecture/`, `docs/release/0.1.1-next-tasks.md`, `papers/`, and
   `.claude/reviews/stabl_lw/summary.md`. Compare after all runs.
4. Resolve and record the R4_51 `Rscript` real path and SHA-256, companion R,
   `R.version.string`, `sessionInfo()`, locale, timezone, BLAS/LAPACK,
   `.libPaths()`, and versions/library paths for loaded scientific and test
   dependencies. Require `pkgload`, `survival`, and `knockoff`; do not install
   or download anything.
5. Use this process environment for every R run:
   `TZ=UTC`, `LC_ALL=C`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`,
   `MKL_NUM_THREADS=1`, `R_PROFILE_USER=/dev/null`, and
   `R_ENVIRON_USER=/dev/null`. Use
   `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla`.
6. Capture a pre-fixture full-suite log and exit status from
   `devtools::test()`. Any failure blocks RE-01.
7. Make a fresh local temporary clone of this repository at the exact required
   HEAD, verify its status is clean, and use only that clone for both validator
   characterizations. Do not weaken or monkey-patch provenance checks.
8. Run methodology serial A, serial B, and parallel C into distinct new raw
   output directories. Fix all settings to:
   `profile="bounded"`, `replicates=2`, `n_bootstraps=2`, `n_lambda=2`,
   `families=gaussian,binomial,multinomial,cox`,
   `artificial_types=random_permutation,knockoff,knockoff_equi,knockoff_mvr`,
   `scenarios=null_independent,null_correlated,signal_independent,signal_correlated,null_high_dim,signal_high_dim`,
   `seed=270627`, `target_fdp=0.1`, and workers 1, 1, and 2.
9. Run late-fusion A and B through `run_late_fusion_validation()` with
   `replicates=1`, `n_bootstraps=2`, `n_iter=2`, and `seed=220711`. Require
   `late_fusion_replicates.csv`, `late_fusion_summary.csv`,
   `late_fusion_warnings.csv`, and `late_fusion_manifest.txt`. The only accepted
   nonzero result is the exact post-artifact condition
   `Late-fusion release gates failed; release must stop.` Any earlier or
   different error blocks RE-01.
10. Implement the normalization helper with an explicit allowlist. For
    methodology, remove only `elapsed_sec` from replicates and
    `mean_elapsed_sec` from summary; copy warnings, parity, and gates unchanged.
    Normalize only `created_at`, `package_mode`, Git start/end/stability fields,
    artifact path values, and artifact hash values in its manifest. Copy every
    late-fusion CSV unchanged; normalize only `package_mode`, Git
    start/end/stability fields, and artifact hash values in its manifest. Never
    sort rows or columns. Preserve seeds, NA encoding, numeric text, warnings,
    condition messages, and gate booleans.
11. Require normalized serial A/B byte identity for all methodology and
    late-fusion artifacts. Require serial/parallel methodology scientific CSV
    identity after the two explicit elapsed columns are removed. Record the
    intentional parallel manifest `workers` difference. Stop on any other
    mismatch.
12. Save normalized serial-A golden files under
    `tests/testthat/fixtures/release_evidence_baseline/v1/methodology/` and
    `tests/testthat/fixtures/release_evidence_baseline/v1/late_fusion/`, plus a
    README and fixture SHA-256 manifest. Add
    `tests/testthat/helper-release-evidence-baseline.R` and
    `tests/testthat/test-release-evidence-baseline-fixtures.R`. The integrity
    test must verify the manifest and expected schemas without regenerating the
    scientific runs.
13. Capture a post-fixture full-suite log and exit status. Recheck production
    paths, status, binary diff, and the protected-file manifest. Any unexpected
    change blocks RE-01.
14. Hash the complete raw baseline closure, excluding the closure manifest and
    seal marker themselves; create `SEALED.NONPROMOTABLE` last; then make the
    raw directory read-only. The marker and fixture README must state that this
    is bounded characterization, not release or publication evidence.
15. Update only RE-01 checkbox/status/evidence, the overall progress summary,
    dated evidence log, and blocker checksum/reference in the revised plan.
    Leave RE-02 through RE-10 `not started`.

## Authorization

Allowed:

- Read repository files and ignored review evidence.
- Create local temporary directories and one clean local clone under `/tmp`.
- Run local R4_51 tests and the fixed bounded/tiny characterizations above.
- Write only the raw baseline directory, the two named test files, the v1
  fixture directory, this task list's checkbox statuses, and RE-01 evidence/
  status entries in the revised plan.
- Create hashes, logs, manifests, and a marker-last non-promotable seal.
- Make only the completed raw baseline directory read-only.

Approval required:

- Any write outside the explicit allowlist.
- Any destructive cleanup, deletion, overwrite, or replacement.
- Any dependency installation, download, network access, scheduler submission,
  credential access, remote write, or external publication action.
- Continuing into RE-02 or later.

Forbidden or out of scope:

- Production/runtime source, metadata, runner, public Interface, or scientific
  formula changes.
- Editing `.gitignore`, `CONTEXT.md`, architecture documentation,
  `docs/release/0.1.1-next-tasks.md`, `papers/`, or the Linus review.
- Repairing, renaming, deleting, or relabeling historical validation artifacts.
- Running a full locked release profile or claiming release/publication
  evidence.
- Git commit, push, branch promotion, tag, CRAN, CI dispatch, or release work.

## Validation Commands

- Pre/post full suite, with the fixed environment above:
  `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla -e 'devtools::test()'`
- Methodology serial A and B, each in a distinct output directory:
  `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla inst/analysis/run_methodology_validation.R --profile bounded --replicates 2 --n-bootstraps 2 --n-lambda 2 --workers 1 --families gaussian,binomial,multinomial,cox --artificial-types random_permutation,knockoff,knockoff_equi,knockoff_mvr --scenarios null_independent,null_correlated,signal_independent,signal_correlated,null_high_dim,signal_high_dim --seed 270627 --target-fdp 0.1 --out <new-output>`
- Methodology parallel C: the same command with `--workers 2` and another new
  output directory.
- Late fusion A/B: source `inst/analysis/run_late_fusion_validation.R` into a
  fresh environment and call `run_late_fusion_validation(<new-output>,
  replicates=1L, n_bootstraps=2L, n_iter=2L, seed=220711L)`.
- Fixture integrity:
  `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript --vanilla -e 'testthat::test_file("tests/testthat/test-release-evidence-baseline-fixtures.R")'`

## Stop Conditions

- Stop if the required HEAD, clean production-input assumption, or raw-target
  nonexistence check fails.
- Stop if a protected file changes, disappears, changes type, or changes hash.
- Stop if a required dependency, runtime identity, permission, or input is
  missing.
- Stop if either full suite fails.
- Stop if a validator omits a required artifact, or late fusion exits for any
  reason other than the exact expected post-artifact gate failure.
- Stop on serial A/B drift or unexpected serial/parallel scientific drift.
- Stop before any write or action outside the authorization ledger.

## Required Artifacts

- Completed task list and RE-01 plan evidence/status update.
- Ignored raw baseline directory with initial/final preservation evidence,
  runtime record, commands, logs, exit statuses, raw and normalized outputs,
  comparisons, closure manifest, and marker-last non-promotable seal.
- Ten normalized golden validator artifacts: six methodology files and four
  late-fusion files, plus README and SHA-256 manifest.
- Normalization helper and fixture-integrity test.
- Exact Linus-review checksum and stable SCI-01 through SCI-10 plus ARCH-01 and
  ARCH-02 blocker references.
- Final report of changed paths, test counts, comparisons, protected-state
  equality, and residual publication blockers.

## Continuation 1 Authorization (2026-07-31)

The user authorized one narrowly bounded, local dependency continuation of the
same unsealed attempt after the fail-closed `knockoff` stop. This section
supersedes only the following earlier clauses:

1. The prohibition on dependency download, installation, and network access is
   replaced by permission to retrieve exactly the three HTTPS resources listed
   below and to install exactly `Rdsdp` 1.0.6 and `knockoff` 0.3.6 into an
   attempt-local R library.
2. The raw-target nonexistence requirement is replaced by a requirement that
   the existing raw attempt is exactly
   `.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/`, is incomplete and
   unsealed, and has no existing `continuations/01-cran-knockoff/` directory.
3. The missing-`knockoff` stop is replaced by permission to close only its
   verified two-source dependency closure and then resume RE-01.

All other authorization boundaries, fixed scientific settings, preservation
requirements, stop conditions, and the prohibition on RE-02 or later remain in
force. In particular, do not modify the R4_51 conda environment, its installed
library, a home/user library, production source, or the historical stopped
record.

### Continuation Boundary and Append-Only Record

1. Before network access or any other continuation write, recheck exact HEAD,
   unchanged production inputs, absence of any seal, and the stopped-attempt
   comparisons and hashes. Fail if any check is not exact.
2. Capture a continuation-boundary manifest for every one of the 49 files then
   present in the stopped attempt, including path, type, mode, size, SHA-256, or
   symlink target. Record explicitly that historical
   `status/raw-files-at-blocked-stop.txt` contains 48 paths because it was
   written before `exits/blocked-state-recheck.exit.txt`; do not repair,
   rewrite, or relabel that historical inventory.
3. Append the authorization record and all subsequent files only beneath the
   new no-replace directory
   `.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/continuations/01-cran-knockoff/`.
   Use new continuation-specific script, command, log, exit, manifest, runtime,
   and status names. No continuation command may overwrite any earlier file.
4. Retain the original exit-1 runtime evidence byte-for-byte. Later successful
   runtime identity is a new continuation epoch, not a correction of the
   stopped epoch.

### Exact Source Closure and Local Installation

Only these resources may be retrieved:

- `https://cran.r-project.org/src/contrib/PACKAGES.gz`
- `https://cran.r-project.org/src/contrib/Rdsdp_1.0.6.tar.gz`
- `https://cran.r-project.org/src/contrib/knockoff_0.3.6.tar.gz`

Retain each response header, effective HTTPS URL, retrieval timestamp, and
downloaded source file. Verify the package-index MD5 values, the SHA-256 values
below, and the internal `DESCRIPTION` package name and version before install:

- `Rdsdp_1.0.6.tar.gz`: MD5 `8f9101ffe5982d629e92396e2ff1fce0`;
  SHA-256 `dd462374f38ae45791c024805582e1275936cb96c2127fd2dceeb3df0f18bfdd`.
- `knockoff_0.3.6.tar.gz`: MD5 `202203a05fa2d3f65266795d63ab47df`;
  SHA-256 `afec4c2bdb8aad4ab4a7f33527a546deda80229588ce7a4d60d37114d3b046ba`.

Before installation, recursively manifest the complete R4_51 installed library.
Create a new empty continuation-local `runtime/library.staging/`, then use the
R4_51 companion `R CMD INSTALL` with an explicit library path to install the
verified `Rdsdp` source first and verified `knockoff` source second. Disable
automatic dependency installation and use the default staged-install locking;
do not use `install.packages()`, `--no-lock`, another repository, or another
download. Set `R_LIBS_USER` to the absolute staging path for installation and
verification. Promote staging to the final `runtime/library/` only after both
installs and every verification below pass.

Require that the final continuation environment has the local library first in
`.libPaths()`, contains locally only `Rdsdp` and `knockoff`, resolves every
other imported dependency from R4_51, contains and hashes
`Rdsdp/libs/Rdsdp.so`, reports no `not found` entry from `ldd`, and returns two
finite values from `knockoff::create.solve_sdp(diag(2))`. Capture the complete
R4_51 library manifest again and require byte identity with the pre-install
manifest. Any additional dependency, source mismatch, package written outside
the continuation library, lock/install/smoke failure, or R4_51 library drift
stops the attempt unsealed without cleanup or retry.

### Resumed RE-01 and Finalization

For every resumed R process, retain the original fixed timezone, locale, thread,
and profile environment and additionally set `R_LIBS_USER` to the absolute
final continuation `runtime/library/`. Record a new complete runtime identity,
then execute the still-unfinished RE-01 task list exactly. Use `final-*` names
for the fresh terminal epoch; do not reuse the historical `blocked-final-*`
names. Repeat the final production-input, protected-file, tracked-diff, and
R4_51 library preservation checks even where the historical task list contains
a checked stop-era item.

Seal only after the complete closure—including the stopped epoch, continuation
authority, retained sources, local installed library, resumed commands and
outputs, comparisons, and final preservation evidence—passes. The closure
manifest excludes only itself and `SEALED.NONPROMOTABLE`; create the marker
last, make the complete attempt read-only, and add nothing afterward. Any new
failure produces a uniquely named continuation stop record and leaves the
attempt incomplete, writable, unsealed, and non-promotable.
