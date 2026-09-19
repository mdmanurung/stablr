# Artifact manifest

## Include with submission or repository archive

| Artifact | Path | Status |
|---|---|---|
| Source release tarball | GitHub release `v0.1.0` | Available |
| Package documentation | `https://mdmanurung.github.io/stablr/` | Live |
| Figure 1 | `papers/application-note/artifacts/figure1_workflow_schematic.pdf` | Available |
| Figure 2 | `papers/application-note/artifacts/figure2_ool_validation.pdf` | Available |
| Table 1 | `papers/application-note/artifacts/table1_ool_parity_metrics.csv` | Available |
| Parity metrics | `papers/application-note/artifacts/ool_publication_parity_metrics.csv` | Available |
| Selected features | `papers/application-note/artifacts/ool_publication_selected_features.csv` | Available |
| Tutorial overlap | `papers/application-note/artifacts/ool_publication_tutorial_overlap.csv` | Available |
| Session info | `papers/application-note/artifacts/ool_publication_sessionInfo.txt` | Available |
| Generation log | `papers/application-note/artifacts/generation.log` | Available |

## Current release evidence

| Check | Result |
|---|---|
| `devtools::test()` | `FAIL 0 | WARN 0 | SKIP 1 | PASS 1611` |
| `R CMD check --no-manual --ignore-vignettes --no-build-vignettes` | `Status: OK` |
| Vignette render | Six HTML files rendered |
| pkgdown | Built and deployed |
| OOL recall | `7/7` |
| OOL Jaccard | `0.778` |
| OOL Spearman | `0.808` |

## Archive decision

Zenodo is optional unless the journal requests a DOI-bearing archive at
submission. If used, archive the GitHub release source plus the artifact files
listed above and update `RESOURCE_AVAILABILITY.md`.
