# RE-01 Execution Task List

Date: 2026-07-31

Repository: `/exports/para-lipg-hpc/mdmanurung/stablr`

Scope: RE-01 only. RE-02 through RE-10 are not authorized.

- [x] Read the executor handoff and complete implementation plan; confirm the
  authorization and stop conditions.
- [x] Fail if HEAD is not
  `c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1`, if the raw baseline target
  already exists, or if production inputs differ from HEAD.
- [x] Create the unique ignored raw baseline attempt and record initial status,
  binary diff, source hashes, and protected-file manifest.
- [ ] Record the complete pinned R4_51 runtime and dependency identity; fail if
  `pkgload`, `survival`, or `knockoff` is unavailable.
- [ ] Run and capture the pre-fixture full `devtools::test()` suite.
- [ ] Create and verify a clean local temporary clone at the exact required
  HEAD.
- [ ] Run the fixed full-branch bounded methodology matrix as serial A, serial
  B, and parallel C in distinct output directories.
- [ ] Run the six-cell tiny late-fusion characterization as A and B; require all
  four artifacts and only the expected post-artifact gate-failure condition.
- [ ] Add the explicit normalization helper and normalize all A/B/C artifacts
  without sorting or broad field removal.
- [ ] Prove serial A/B normalized byte identity and serial/parallel scientific
  table identity; save comparison reports.
- [ ] Save serial-A golden fixtures, README, and SHA-256 manifest under the v1
  fixture directory.
- [ ] Add and run the fixture-integrity test.
- [ ] Run and capture the post-fixture full `devtools::test()` suite.
- [x] Recheck production inputs and prove the protected-file manifest and
  tracked diff are unchanged.
- [ ] Write the raw closure manifest, create `SEALED.NONPROMOTABLE` last, and
  make the ignored raw baseline read-only.
- [x] Update only RE-01 status/evidence and blocker references in the main plan;
  keep later steps not started.
- [x] Record final `git status --short`, changed paths, validation results, and
  residual blockers for the parent.

## Continuation 1: Isolated CRAN Closure

The checked items above describe the preserved stopped epoch. They do not waive
fresh terminal checks for this continuation.

- [x] Recheck exact HEAD, unchanged production inputs, stopped-attempt hashes
  and comparisons, absence of a seal, and absence of
  `continuations/01-cran-knockoff/`; stop on any mismatch.
- [x] Create a no-replace continuation directory and, before network access,
  manifest all 49 stopped-attempt files. Document without altering it that the
  historical 48-path blocked-stop inventory preceded its final exit record.
- [x] Append the dated user authority and narrowly superseded clauses to the
  continuation record; preserve every stopped-epoch file byte-for-byte.
- [x] Retrieve and retain only CRAN `PACKAGES.gz`, `Rdsdp_1.0.6.tar.gz`, and
  `knockoff_0.3.6.tar.gz`, with headers, effective HTTPS URLs, timestamps, and
  command exits.
- [x] Verify the package-index MD5 values, pinned SHA-256 values, and internal
  package names/versions before installation; stop on any mismatch.
- [x] Capture a recursive pre-install manifest of the complete R4_51 installed
  library.
- [ ] Install verified `Rdsdp` then `knockoff` with R4_51 `R CMD INSTALL` into
  a new empty continuation-local staging library, with no automatic dependency
  installation and no writes to the conda or home library.
- [ ] Verify package provenance, resolved import locations, `Rdsdp.so` hash and
  `ldd`, and a finite two-value `create.solve_sdp(diag(2))` smoke; promote the
  staging library to its final attempt-local path only after all checks pass.
- [ ] Re-manifest the complete R4_51 installed library and require byte identity
  with its pre-install manifest.
- [ ] Record a new complete runtime identity with the absolute final local
  library deliberately first in `R_LIBS_USER`; retain the original runtime
  failure record.
- [ ] Resume every still-unfinished original RE-01 task under that exact local
  library and the frozen runtime environment.
- [ ] Repeat final production-input, protected-file, tracked-diff, R4_51
  library, and changed-path preservation checks, regardless of stop-era
  checkmarks.
- [ ] Write new `final-*` terminal records and the complete closure manifest;
  create `SEALED.NONPROMOTABLE` last, make the attempt read-only, and perform
  no later write.
