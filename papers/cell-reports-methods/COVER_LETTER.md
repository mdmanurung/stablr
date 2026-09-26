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
as a new statistical method. We validate the implementation against the pinned
Python STABL reference on fixed synthetic fixtures, where selection masks, FDP+
curves, stacking behaviour, metric scalars and signal-feature ranking all agree.
We additionally report a full-profile Onset of Labor proteomics run as a negative
result: at publication scale (150 samples, 1,317 features, 500 bootstraps) the R
release recovers all seven reference selected features but reaches only Jaccard
0.778 and Spearman 0.785 against frozen thresholds of 0.80 and 0.90, because the
artificial-feature generator it used was calibrated against a different null
distribution than the reference. A corrected run is pending, so we make no OOL
parity claim. The release passes local R package tests, `R CMD check --as-cran`,
vignette rendering, and pkgdown documentation deployment.

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
