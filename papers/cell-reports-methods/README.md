# Cell Reports Methods submission package

Target journal: Cell Reports Methods

Working title: stablr: stability-based biomarker discovery with reproducible FDP+ control in R

This folder is a Cell Press-style publication package for `stablr` v0.1.0. It
reuses the validated OOL parity artifacts from
`papers/application-note/artifacts/` and reframes the work as a methods/software
resource rather than a short application note.

## Package contents

| File | Purpose |
|---|---|
| `MANUSCRIPT_DRAFT.md` | Main manuscript draft with Cell Press-style sections |
| `COVER_LETTER.md` | Cover letter draft for editorial triage |
| `HIGHLIGHTS_AND_ETOC.md` | Highlights and eTOC blurb |
| `RESOURCE_AVAILABILITY.md` | Lead contact, materials availability, data/code availability |
| `STAR_METHODS.md` | Detailed methods and reproducibility details |
| `KEY_RESOURCES_TABLE.tsv` | Key resources table in tab-separated format |
| `FIGURE_LEGENDS.md` | Figure and table legends |
| `ARTIFACT_MANIFEST.md` | Submission/reproducibility artifact checklist |
| `SUBMISSION_CHECKLIST.md` | Remaining actions before upload |

## Current validated release evidence

- GitHub release: `https://github.com/mdmanurung/stablr/releases/tag/v0.1.0`
- Package documentation: `https://mdmanurung.github.io/stablr/`
- Local package tests: `FAIL 0 | WARN 0 | SKIP 1 | PASS 1611`
- `R CMD check --no-manual --ignore-vignettes --no-build-vignettes`: `Status: OK`
- Six pkgdown articles live on the package site
- OOL parity metrics: `n=150`, `p=1317`, 500 bootstraps, 10 lambda values,
  7/7 Python tutorial selected-feature recall, Jaccard `0.778`, Spearman
  `0.808`

## Scope lock

This package should not claim a new STABL method or a predictive-performance
advance over DIABLO. The supported claim is that `stablr` provides a documented,
tested, R-native implementation of STABL-style stability selection with FDP+
control and practical single-omic/multi-omic workflows.

TCGA nested-CV vs DIABLO remains deferred unless the full outer-fold benchmark
is completed and reported with paired uncertainty intervals.
