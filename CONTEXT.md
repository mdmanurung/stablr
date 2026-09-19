# stablr Release Evidence

This context defines the evidence that may support a `stablr` release and keeps
development diagnostics separate from publication-grade validation.

## Language

**Candidate Source Artifact**:
The immutable `stablr` source tarball whose SHA-256 identity is executed by a release validation run.
_Avoid_: release checkout, clean source, candidate tree

**Runtime Artifact**:
The immutable hashed manifest identifying the R executable, dependency set, native libraries, and platform used for release validation.
_Avoid_: session info, environment dump, R environment

**Release Evidence Packet**:
An immutable aggregate that binds one Candidate Source Artifact and one Runtime Artifact to the passing Validation Run Packets required for release promotion.
_Avoid_: results folder, validation output, release evidence directory

**Validation Run Packet**:
An immutable record of one validation attempt against one Candidate Source Artifact and one Runtime Artifact, including settings, results, Scientific Gates, manifest, and checksums.
_Avoid_: job output, results directory, validation run

**Validation Attempt**:
A mutable temporary workspace used while one required validation is executing and not eligible for release promotion.
_Avoid_: partial packet, running packet

**Bounded Run**:
A development validation run that may execute checkout source and is never promotable release evidence.
_Avoid_: smoke evidence, preliminary release run

**Scientific Gate**:
A predeclared pass-or-fail criterion evaluated from validation results before release promotion.
_Avoid_: quality check, threshold check

**Required Validation**:
A named publication validation for which one passing Validation Run Packet is mandatory before release promotion.
_Avoid_: validation script, release job

**Release Contract**:
The versioned declaration inside a Candidate Source Artifact that defines its Required Validations, evidence schemas, gate versions, and Runtime Artifact requirements.
_Avoid_: validator list, release configuration, workflow settings

**Release Context**:
An immutable, serializable descriptor that binds one Candidate Source Artifact, one Runtime Artifact, one Release Contract, and one local evidence store before any Required Validation starts.
_Avoid_: release workspace, validation environment, mutable run state

**Packet Reference**:
An immutable, serializable locator for one sealed Validation Run Packet, including the packet identity and manifest SHA-256 needed for independent verification.
_Avoid_: output directory, results path, discovered packet

## Relationships

- Every **Validation Run Packet** is bound to exactly one **Candidate Source Artifact**
- Every **Validation Run Packet** is bound to exactly one **Runtime Artifact**
- Every **Candidate Source Artifact** contains exactly one **Release Contract**
- Every **Required Validation** owns its predeclared **Scientific Gates** and the calculations used to evaluate them
- A **Release Contract** declares every **Required Validation** needed by its **Release Evidence Packet**
- A **Release Contract** selects the Adapter and fixed release settings for each **Required Validation**; a caller cannot replace, omit, or weaken them
- A **Release Context** binds exactly one **Candidate Source Artifact**, one **Runtime Artifact**, one **Release Contract**, and one local evidence store
- A **Packet Reference** identifies exactly one sealed **Validation Run Packet** and is invalid if its recorded manifest SHA-256 does not match
- Every **Release Evidence Packet** is bound to exactly one **Candidate Source Artifact** and one **Runtime Artifact**, and references exactly one passing **Validation Run Packet** for each **Required Validation**
- A normally completed **Validation Attempt** is sealed exactly once as a passing or failed **Validation Run Packet**
- A crashed, killed, or otherwise incomplete **Validation Attempt** remains unsealed and cannot be referenced by a **Release Evidence Packet**
- A **Release Evidence Packet** may reference only passing **Validation Run Packets**
- A **Validation Run Packet** passes only when its evidence is complete and every **Scientific Gate** passes
- Failed **Validation Run Packets** are retained and are never replaced or mutated by retries
- Every **Validation Run Packet** contains one or more **Scientific Gates**
- A **Release Evidence Packet** is assembled only from explicit **Packet References** and never by scanning an evidence store
- A **Bounded Run** cannot produce a **Release Evidence Packet**

## Example dialogue

> **Developer:** "The checkout-based **Bounded Run** passed; can this release be promoted?"
> **Maintainer:** "No. Build the **Candidate Source Artifact**, seal each **Validation Run Packet**, and assemble its **Release Evidence Packet** only after every required **Scientific Gate** passes."

## Flagged ambiguities

- "clean source" previously meant either a Git checkout or the code actually executed; release validation now means the hashed **Candidate Source Artifact** only.
