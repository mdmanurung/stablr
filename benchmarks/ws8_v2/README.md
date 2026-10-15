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
   Run once per method group via `WS8_METHODS` (e.g. `M0,M1,M2`, `M3`, `M4`, `M5`),
   writing `results/run_<id>/`.
6. `assemble_primary.R` — merge the per-method replication tables into
   `results/run/` and compute the predeclared contrasts.
7. `null_probe.R` — within-cohort age-permutation null (S0).
8. `biology.R` — reference populations, GZMK/GZMB ratio, heterogeneity (secondary).
9. `figures.R` — regenerate all plots from saved, validated results.

S1 (known-signal DGP simulations) and S2 (threshold/learner sensitivity probes) are
specified in `protocol.md` §9 but were not run in this milestone: the S0 null already
answered the decisive H5 question (selection-count saturation), and the primary
endpoint is ranking-based. They remain available for a future amendment.

## Inputs

- Raw cell-level file `all_cells_metadata.csv.gz` (MD5
  `50651968855bf6ca3efcdf199e2a9fd5`). Set its directory via env var `WS8_DATA_DIR`
  (default `/workspace`). Do NOT commit the raw per-cell file.
- R package `stablr` (this checkout), plus `glmnet`, `ranger`, `knockoff`,
  `data.table`, `ggplot2` (figures only).

## Environment variables

- `WS8_DATA_DIR` — directory containing the raw portal file (default `/workspace`).
- `WS8_V2_IN` — input results directory for downstream scripts
  (default `benchmarks/ws8_v2/results/m1`).
- `WS8_V2_OUT` — output directory for `run.R` (default `results/run` under the
  benchmark); set per method group.
- `WS8_METHODS` — comma-separated selector subset for `run.R` (default all M0–M5).
- `WS8_N_PERM`, `WS8_PERM_START` — permutation count and start index for
  `null_probe.R` (defaults 10 and 1), enabling chunked parallel null runs.

## Outputs

All outputs under `results/<run_id>/` (versioned). See `protocol.md` §10 for the
required tables. Historical v1 artifacts in `benchmarks/` are preserved
byte-for-byte and never overwritten.

## Limitations

Application case study, not calibration evidence. Exact source-paper clock comparator
is unavailable (Table S1 paywalled); a conceptual comparator is used instead. See
`protocol.md` §8 and `claims_ledger.md`.
