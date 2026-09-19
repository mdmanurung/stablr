# STAR Methods

## Method details

### Software implementation

`stablr` is implemented as an R package. The core fitting function,
`stabl_fit()`, accepts a feature matrix, outcome, regularization grid, family,
artificial-feature mode, bootstrap count, and random seed. The package validates
sample alignment by identifiers, fits sparse models across bootstrap samples,
computes feature recurrence scores, estimates FDP+ using artificial null
features, and returns selected features through S3 accessors.

The v0.1.0 release supports gaussian, binomial, multinomial, and Cox families
where the backend learner supports the family. Learner adapters include lasso,
elastic net, adaptive lasso, and optional sparse group lasso. Artificial
features can be generated with random permutation or knockoff-based strategies.

### Multi-omic workflows

For named lists of omic matrices, `stabl_multiomic_train_validate()` fits
per-view models and optional early-fusion, late-fusion, and cooperative-fusion
branches. `stabl_multiomic_cv()` supports outer cross-validation, and
`stabl_multiomic_nested_cv()` provides infrastructure for nested evaluation
against external comparators. Cooperative fusion is exposed as an R workflow
extension and is not claimed as part of the original STABL method.

### OOL parity validation

The validation analysis used the Onset of Labor proteomics tutorial data with
150 samples and 1,317 features. The release was run with gaussian lasso models,
knockoff artificial features, 500 bootstraps, 10 lambda values, and fixed seeds.
The Python tutorial selected-feature set contained seven proteins. `stablr`
selected nine proteins and recovered all seven reference selected features.

The primary validation metrics were:

- selected-feature recall relative to the Python tutorial set;
- Jaccard similarity between selected-feature sets;
- Spearman correlation between R maximum stability scores and Python tutorial
  scores;
- runtime and package/session metadata.

### Divergence between R and Python implementations

The validation target is feature-selection and ranking concordance, not
bit-identical coefficients. R `glmnet` and Python elastic-net implementations
may differ in coordinate descent schedules, convergence thresholds, and
coefficient paths. These implementation-level differences are expected and do
not undermine selection-level parity when matched preprocessing, seeds, lambda
grids, and artificial-feature settings produce concordant selected features and
rankings.

## Quantification and statistical analysis

Jaccard similarity was computed as the size of the intersection divided by the
size of the union of the R and Python selected-feature sets. Spearman
correlation was computed between maximum feature stability/importances for
features present in the tutorial score table. No hypothesis test was used for
the software release claim. The validation is descriptive and release-oriented.

## Reproducibility

Run from the companion experiments repository:

```bash
STABLR_PKG=/exports/para-lipg-hpc/mdmanurung/stablr \
  /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  analysis/generate_publication_parity_table.R
```

Package release verification:

```bash
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript -e "devtools::test()"
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/R CMD build --no-build-vignettes .
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/R CMD check --no-manual --ignore-vignettes --no-build-vignettes stablr_0.1.0.tar.gz
```
