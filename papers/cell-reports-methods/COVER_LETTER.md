# Cover letter

Dear Cell Reports Methods editorial team,

We are pleased to submit "`stablr`: stability-based biomarker discovery with
reproducible FDP+ control in R" for consideration as a methods/software article
in Cell Reports Methods.

High-dimensional omic studies need biomarker signatures that are sparse,
interpretable, and reproducible across resampled cohorts. The recently published
STABL framework addresses this need through bootstrap stability selection and
FDP+ control with artificial features, but the reference implementation is
Python-based. `stablr` provides an R-native implementation and release package
for the R/`glmnet` ecosystem, including documented APIs, S3 objects, FDP+
calibration, random-permutation and knockoff artificial features, multi-omic
workflows, diagnostics, export utilities, vignettes, and package-level
verification.

The submitted package is framed as a software and reproducibility resource, not
as a new statistical method. We validate the v0.1.0 release on the Onset of
Labor proteomics tutorial dataset using publication-scale settings (150
samples, 1,317 features, 500 bootstraps, knockoff null features). The R release
recovers all seven Python tutorial selected features, with selected-set Jaccard
similarity of 0.778 and feature-score Spearman correlation of 0.808. The release
passes local R package tests, no-manual `R CMD check`, vignette rendering, and
pkgdown documentation deployment.

We believe this manuscript fits Cell Reports Methods because it provides a
practical, validated, open-source computational method implementation for
reproducible omics biomarker discovery and documents the limits of its current
benchmark claims. The package, release, documentation, and reproduction
artifacts are available at:

- Source and release: `https://github.com/mdmanurung/stablr/releases/tag/v0.1.0`
- Documentation: `https://mdmanurung.github.io/stablr/`
- Companion artifacts: `papers/application-note/artifacts/` in the
  `stablr-experiments` repository

This manuscript is original, is not under consideration elsewhere, and all
authors have approved submission. We declare no competing interests.

Sincerely,

Mikael Manurung
