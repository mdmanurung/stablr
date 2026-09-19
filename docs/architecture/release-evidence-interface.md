# Artifact-First Release Evidence Interface

Status: Accepted

Decision date: 2026-07-31

Scope: local `stablr` scientific release validation

## Outcome

`stablr` will use one deep **Release Evidence Module** for candidate identity,
runtime identity, packet lifecycle, sealing, and aggregation. Methodology and
late-fusion validation remain two real **Adapters** at the Required Validation
**Seam**. Each Adapter owns its scientific calculations and Scientific Gates;
the Module owns whether the resulting evidence is complete, immutable, and
eligible for release promotion.

This decision locks the evidence architecture. It does not repair or waive any
currently known scientific defect, and it does not make the current branch
publication-ready.

## External Interface

The canonical R Interface is deliberately limited to three calls:

```r
context <- prepare_release_context(candidate_tarball, runtime, store)
packet  <- run_required_validation(context, validation_id)
release <- assemble_release_evidence(context, packet_refs)
```

A thin command-line wrapper mirrors `prepare`, `run`, and `assemble` without
adding release logic.

### `prepare_release_context(candidate_tarball, runtime, store)`

The caller supplies:

- `candidate_tarball`: an existing `R CMD build` source tarball;
- `runtime`: the absolute path to the pinned `Rscript` executable; and
- `store`: a local POSIX evidence-store path.

The Module safely inspects and copies the tarball into the store, computes its
SHA-256, reads and validates the embedded Release Contract, installs only that
tarball into an isolated library containing a store-owned copy of the
contract-declared dependency closure, constructs the Runtime Artifact, and
returns an immutable `stablr_release_context` reference. Preparation performs
no dependency download. The stored context is the source of truth; modified
caller-side list fields are never trusted.

Preparation fails closed if the tarball is malformed, the Release Contract is
missing or invalid, any required hash cannot be resolved, the isolated install
fails, the runtime is unsupported, or the evidence store cannot demonstrate
the required local atomic-write behavior.

### `run_required_validation(context, validation_id)`

`validation_id` must name a Required Validation declared by the embedded
Release Contract. The contract selects its Adapter, settings, schemas, and gate
versions. The caller cannot inject an Adapter or override release settings.

The Module re-verifies the context, Candidate Source Artifact, Runtime
Artifact, and isolated library, then executes the selected Adapter in a fresh
`Rscript --vanilla` process. A passing run returns an immutable
`stablr_packet_ref`.

If the Adapter completes normally but any required artifact or Scientific Gate
is missing or failed, the Module seals a failed Validation Run Packet and then
signals `stablr_validation_failed`. That condition carries the failed Packet
Reference so diagnostics remain reachable while the CLI exits nonzero.

If the process crashes, is killed, or raises an unexpected operational error,
the Validation Attempt remains unsealed and ineligible. No Packet Reference is
fabricated.

### `assemble_release_evidence(context, packet_refs)`

`packet_refs` is an explicit list of Packet References. The Module never scans
the store to discover evidence. It verifies every packet's seal, manifest hash,
schema, status, Candidate Source Artifact, Runtime Artifact, Release Contract,
Required Validation identity, and gate-table completeness.

Assembly requires exactly one passing Validation Run Packet for every Required
Validation, with no duplicate, missing, extra, failed, unsealed, bounded, or
foreign packet. Success returns an immutable `stablr_release_ref` for the
sealed Release Evidence Packet. Assembly never recalculates a scientific
formula.

## Lifecycle

```text
Candidate Source Artifact + pinned runtime
                    |
                    v
              Release Context
                    |
                    v
            Validation Attempt
              /             \
 normal completion           crash / kill / unexpected error
        |                              |
        v                              v
 Module derives status          unsealed and ineligible
      /       \
 passing       failed
   |             |
   v             v
sealed passing  sealed failed ----> typed failure with Packet Reference
packet          packet
   |
   +---- explicit passing Packet References ----> Release Evidence Packet
```

Retries always create a new Validation Attempt and never mutate or replace a
prior passing or failed packet.

## Module and Adapter Ownership

### Release Evidence Module

The Module owns the following Implementation:

- Candidate Source Artifact hashing, safe extraction, and store-owned copying;
- Release Contract parsing, schema checks, and enforcement;
- isolated installation and fresh-process execution;
- Runtime Artifact construction and before/after identity verification;
- Validation Attempt allocation;
- artifact-closure validation, SHA-256 manifests, and seal-last publication;
- status derivation from complete gate tables;
- immutable Packet References and Release Context serialization;
- exact packet-set verification and Release Evidence Packet assembly; and
- typed lifecycle errors.

This concentration provides **Locality** for provenance and failure semantics
and **Leverage** across both validators and future Required Validations.

### Required Validation Adapters

The methodology and late-fusion Adapters own:

- scientific simulation or benchmark execution;
- result and warning artifacts;
- Scientific Gate formulas, versions, scopes, observations, and criteria; and
- Adapter-specific scientific metadata.

On normal completion, an Adapter returns named relative artifact paths, a
complete gate table, and metadata. The common gate-table Interface contains at
least `gate_id`, `gate_version`, `scope`, `observed`, `criterion`, `pass`, and
`reason`. An Adapter does not seal packets, hash the evidence closure, choose
the overall packet status, or exit the process merely because a gate failed.

The two existing scientific Adapters make this a real Seam. Filesystem, hashing,
subprocess, and runtime-discovery seams remain private to the Module
Implementation until more than one production Adapter genuinely exists.

### Command-Line Wrapper

The command-line wrapper only parses arguments, invokes the same three
functions, serializes references, and maps typed conditions to nonzero exit
codes. It contains no Scientific Gate, identity, sealing, or aggregation logic.

## Identity and Serialization

- Release Contract, Release Context, Runtime Artifact, packet manifest, Packet
  Reference, and release manifest use versioned JSON schemas.
- JSON bytes and every packet artifact are SHA-256 bound. Hash failure is an
  error, never `NA` provenance.
- Contexts and references are serializable so independent local or Slurm jobs
  can run Required Validations without shared R memory.
- Store layout is private Implementation. Callers pass explicit references and
  must not derive meaning from directory names.
- The v1 Runtime Artifact records and hashes the `Rscript` and companion R
  executables, R version and home, platform, locale, BLAS/LAPACK identity, the
  isolated store-owned R package closure, the installed `stablr` tree, and
  resolved native libraries.
- Every validation re-verifies runtime and candidate identity before execution
  and before sealing. Runtime drift leaves the attempt unsealed.
- The first Implementation supports the current local Linux/POSIX runtime and
  fails closed on unsupported runtime discovery. Remote stores, containers, and
  scheduler submission are not part of this Interface.

## Execution Invariants

1. Promotable validation executes package and Adapter code only from the built
   Candidate Source Artifact, never from a checkout, `pkgload::load_all()`, or
   a preloaded namespace.
2. The fresh child process uses `--vanilla`, a sanitized R environment, and the
   context's isolated library. It verifies that the `stablr` namespace path is
   inside that library.
3. The Release Contract is inside the Candidate Source Artifact and is itself
   covered by the candidate hash.
4. The Module derives packet status only after the Adapter returns a complete
   gate table and the full artifact closure passes schema and hash checks. A
   packet passes only when every required gate value is exactly `TRUE`; `NA`,
   missing, duplicated, and unknown gates fail.
5. Sealing is marker-last and no-replace. Only a verified sealed location is
   eligible for a Packet Reference.
6. Symlinks, paths outside the Validation Attempt, missing hashes, duplicate
   artifact names, and undeclared artifacts fail closed.
7. A Bounded Run is always marked non-promotable and cannot be converted into a
   Validation Run Packet or Release Evidence Packet.

## Typed Error Modes

- `stablr_release_contract_error`: missing, invalid, unknown, or weakened
  Release Contract data;
- `stablr_release_identity_error`: candidate, manifest, packet, or namespace
  identity mismatch;
- `stablr_release_runtime_error`: unsupported, incomplete, or changed Runtime
  Artifact;
- `stablr_validation_incomplete`: the Adapter did not complete normally, so no
  sealed Validation Run Packet exists;
- `stablr_validation_failed`: a failed Validation Run Packet was sealed; the
  condition includes its Packet Reference; and
- `stablr_release_assembly_error`: the explicit packet set is invalid or
  incomplete.

## Rejected Shapes

- A general plug-in protocol engine was rejected for now because callers need
  only two real Adapters. Exposing registry and storage seams would make the
  Interface shallower without current Leverage.
- A mutable stateful release handle was rejected because Required Validations
  must run independently and serialize cleanly across Slurm jobs.
- Directory-discovery assembly was rejected because stale or foreign packets
  could be selected implicitly.
- Git-checkout provenance was rejected because it does not identify the code
  loaded into R and cannot bind untracked package inputs.

## Non-Goals

- repairing the known bootstrap, sample-alignment, fold-allocation, fallback,
  provenance-field, high-dimensional, or gate-estimator defects;
- changing Scientific Gate thresholds or formulas as part of the architecture
  migration;
- running the multi-day locked methodology validation during routine tests;
- replacing CRAN, CI, Win-builder, R-hub, or platform checks; or
- claiming publication readiness before all scientific and external gates pass.
