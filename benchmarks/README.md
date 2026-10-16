# benchmarks/ — experiment and analysis scripts (not package source)

This directory holds the development benchmarks, case studies, validation
experiments, and their result artifacts for `stablr`. It is **not part of the
package**: nothing here is installed, tested by `R CMD check`, or required to
use the package. It is excluded from the source tarball via `.Rbuildignore`.

## What lives where

| Location | Contents |
|---|---|
| `R/`, `src/` | Package source code |
| `tests/testthat/` | Package unit/regression tests (run by `R CMD check`) |
| `vignettes/` | Package vignettes (self-contained; figures committed) |
| `benchmarks/` | Experiment/analysis scripts, predesigns, protocols, result CSVs |
| `benchmarks/ws8_v2/` | The WS-8 v2 prospective benchmark package (locked protocol, claims ledger, run manifests) |
| `benchmarks/archive/` | Superseded artifacts kept for provenance |
| `scripts/` | Package-maintenance utilities (codegen, vendoring) |

## Conventions

- Workstream scripts carry their study prefix (`b1_*`, `ws2_*`, …, `ws9_*`).
- Every study keeps its own predesign/protocol and a summary; result CSVs are
  committed so the vignettes and summaries are auditable without re-running.
- Scripts assume the repository root as working directory and load the package
  with `pkgload::load_all(".")`.
- Benchmark-local test scripts (e.g. `ws8_v2/test_loco.R`) are experiment
  harnesses, not package tests; they run standalone with `Rscript`.
