# WS-8 v2 — Reproducible Cross-Cohort Immune-Aging Predictors

Prospective revision of the exploratory WS-8 aging case study. This directory holds a
leakage-proof, matched-budget benchmark of environment-aware STABL against simpler
selectors, with association replication as the primary endpoint and predictive
transportability as a secondary endpoint.

**Read `protocol.md` first.** It is the locked scientific contract. `claims_ledger.md`
traces every claim to a versioned artifact.

## Run order

1. `data_prepare.R` — rebuild donor-level counts/metadata from the raw portal file;
   enforce the `Bad quality` exclusion; log missingness. (Requires the raw file; see
   below.)
2. `data_checks.R` — validate IDs, denominators, cohort sizes, ranges, hashes.
3. `loco.R` — define the 8 fixed outer folds and fold-safe transforms.
4. `selectors.R` — M0–M5 training-only rankings/selections per fold.
5. `run.R` — thin driver: primary replication benchmark + secondary prediction.
6. `null_probe.R` — within-cohort age-permutation null (S0).
7. `simulation.R` — known-signal DGPs (S1). *(created later)*
8. `biology.R` — reference populations, ratios, adjusted associations. *(later)*
9. `figures.R` — regenerate all plots from saved, validated results.

## Inputs

- Raw cell-level file `all_cells_metadata.csv.gz` (MD5
  `50651968855bf6ca3efcdf199e2a9fd5`). Set its directory via env var `WS8_DATA_DIR`
  (default `.`). Do NOT commit the raw per-cell file.
- R package `stablr` (this checkout), plus `glmnet`, `ranger`, `knockoff`.

## Outputs

All outputs under `results/<run_id>/` (versioned). See `protocol.md` §10 for the
required tables. Historical v1 artifacts in `inst/benchmarks/` are preserved
byte-for-byte and never overwritten.

## Limitations

Application case study, not calibration evidence. Exact source-paper clock comparator
is unavailable (Table S1 paywalled); a conceptual comparator is used instead. See
`protocol.md` §8 and `claims_ledger.md`.
