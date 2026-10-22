# stablr 0.3.0.9000 (development)

## Documentation and guidance (2026-10-21, WS-8 v3 candidate benchmark)

- The WS-8 v3 candidate benchmark (`benchmarks/ws8_v3/`) extended the
  WS-8 v2 six-selector matrix with four candidates — M6
  heterogeneity-filtered random-effects meta-analysis, M7 MINT-STABL
  (with tuned control M7c), M8 contrast-penalized STABL, and M9
  ComBat-merge-then-select (with transform controls M9b/M9s) — on the
  eight-cohort aging atlas, 8-fold leave-one-cohort-out, package pinned
  at one commit. No candidate beat the simple random-effects
  meta-analysis ranking (M0, macro R@10 0.395); M7, M8, and M9 were
  significantly worse (sign-flip p = 0.008); the M4 − M0 contrast
  replicated attenuated (−0.068, p = 0.047).
- Repositioned roles on this evidence: `stabl_environment_stability()`
  is the transportability *gate* and batch-dependence diagnostic
  (the membership question); per-cohort association with
  random-effects meta-analysis is the *ranking* procedure. The
  "Environment Stability" and "aging case study" vignettes now
  cross-reference the benchmark (Sections 6 and 10.9).
  ComBat-merge-then-select is documented as an anti-pattern when batch
  structure is associated with the outcome.

## API changes (2026-10-21, Phase 14 PR-A: internal machinery hidden)

- Nine exported functions are internalized — they remain in the package
  namespace (reachable as `stablr:::name`, with `\keyword{internal}` man
  pages retained) but no longer appear on the public reference page, in
  the mixOmics style of keeping machinery behind a small set of fit
  functions and verb generics:
  - input validation: `validate_sample_alignment()`,
    `validate_multiomic_inputs()`, `validate_newdata_schema()` — called
    automatically by `stabl_fit()`, the `stabl_multiomic_*()` workflows,
    and `predict.stabl_predictor()` respectively;
  - artificial features: `make_artificial_features()` — dispatched by the
    `artificial_type` argument of the fit functions;
  - learner-specific grid helpers: `auto_rf_threshold_grid()`,
    `auto_ordinal_lambda_grid()`, `auto_splsda_keepx_grid()` —
    `auto_lambda_grid()` stays exported as the single public grid helper,
    and `lambda_grid = "auto"` covers the other learners internally;
  - bootstrap index generators: `sample_bootstrap_indices()`,
    `group_bootstrap_indices()` — resampling is controlled through
    `stabl_fit()`'s `n_bootstraps`, `stratify_bootstrap`,
    `bootstrap_strata`, and `sample_groups` arguments.
- The public export count drops from 54 to 45.
  `inst/release/public-interface-audit.csv` and `_pkgdown.yml` are
  updated to match; `docs/API_REFERENCE.md`,
  `docs/MIGRATION_0.3_TO_1.0.md`, and `docs/PYTHON_TO_R_MAPPING.md`
  record the new internal status; the grouped-bootstrap demonstration in
  `vignette("stablr-advanced")` is inlined to plain base R.
- Existing code that calls the nine functions directly keeps working via
  `stablr:::`; no behavior changes.

## API changes (2026-10-21, Phase 14 PR-B: redundant families merged)

- Three redundant function families are merged into single entry points in
  the mixOmics style of a small verb-based API.  The ten historical names
  remain exported as deprecated shims — one `.Deprecated()` warning per
  call, legacy argument names mapped silently — and will be removed in
  1.0.0.  The shims are excluded from the pkgdown reference index; their
  man pages remain reachable via `?name`.
- Similarity (8 -> 2): `stabl_similarity()` unifies scalar and pairwise
  similarity with `method = "jaccard"` (default), `"adjusted"`, or
  `"pearson"`: two vectors give a scalar, a list of sets gives the flat
  row-major upper-triangle vector (Python convention), and
  `as_matrix = TRUE` returns the full symmetric matrix (`remove_diag = TRUE`
  drops the diagonal row-wise).  `stabl_similarity_summary()` reduces
  pairwise values to a median/mean statistic with a spread measure.  The
  shims `jaccard_similarity()`, `jaccard_similarity_matrix()`,
  `adjusted_similarity()`, `adjusted_similarity_pairwise()`,
  `adjusted_similarity_summary()`, `pearson_similarity()`,
  `pearson_similarity_pairwise()`, and `pearson_similarity_summary()`
  forward to the two new functions.  One deliberate change: the matrix
  diagonal for `method = "adjusted"` is now the honest computed
  self-similarity (1 - k/d for a set of size k <= d/2) instead of a forced
  1; the jaccard and pearson diagonals remain 1.
- Export (2 -> 1): `save_stabl_results()` gains `what = c("all", "csv")`;
  `what = "csv"` writes only the score tables and requires neither `x`/`y`
  nor ggplot2 (the previous `export_stabl_to_csv()` semantics, including
  lenient reuse of an existing directory).  `export_stabl_to_csv()` is a
  deprecated shim.
- Preprocessing (2 -> 1): `stabl_preprocess_fit()` now accepts an unprepped
  `recipes::recipe` as `transforms` directly (prepped on training data
  only); `stabl_preprocess_fit_recipes()` is a deprecated shim.
- The public export count rises from 45 to 47 (the two new canonical
  similarity functions); the pkgdown reference index lists 37 exported
  functions with the ten shims quarantined.
  `inst/release/public-interface-audit.csv`, `_pkgdown.yml`,
  `docs/API_REFERENCE.md`, `docs/MIGRATION_0.3_TO_1.0.md`,
  `docs/PYTHON_TO_R_MAPPING.md`, the README, and the affected vignettes are
  updated to match; the new `test-api-shims.R` pins the
  deprecation-plus-forwarding contract of all ten shims, and the
  legacy-argument mapping tests live in `test-api-deprecation.R`.

## API changes (2026-10-21, Phase 14 PR-C: S3 and multi-omic audit)

- The streamlining audit was extended to the 43 registered S3 methods and
  the eight exported `stabl_multiomic_*` functions.  All 43 S3
  registrations are required for dispatch and every class is still live,
  so no S3 method could be removed or hidden; the reference index already
  lists only the 18 methods that have their own pages, with the accessor
  methods documented on their generic/class pages (the mixOmics verb
  pattern).  One genuine gap is fixed: `print.stabl_interpretation()` was
  undocumented and is now aliased on the `stabl_interpret()` help page.
- `stabl_multiomic_study_split()` is internalized — study-disjoint
  splitting is subsumed by the `study_labels` arguments of
  `stabl_multiomic_nested_cv()` and `stabl_multiomic_scheme_comparison()`
  (the PR-A machinery precedent).  The implementation, tests, and a
  `\keyword{internal}` man page with a `stablr:::` example are retained;
  existing code keeps working via `stablr:::stabl_multiomic_study_split`.
- `stabl_multiomic_train_validate()` and `stabl_multiomic_cv()` stay
  separate: a merge into one function with a validation-mode argument was
  considered and rejected (they are the two most-used multi-omic entry
  points, and the shared accessor verbs already provide the streamlined
  pattern).
- The public export count drops from 47 to 46.
  `inst/release/public-interface-audit.csv`,
  `inst/release/removal-ledger.csv`, and `_pkgdown.yml` are updated to
  match.

## API changes (2026-10-21, Phase 14 PR-D: the 1.0.0 API cut)

- The ten deprecated shims whose removal was announced for 1.0.0 are now
  removed: the eight 0.3.0 similarity names (`jaccard_similarity()`,
  `jaccard_similarity_matrix()`, `adjusted_similarity()`,
  `adjusted_similarity_pairwise()`, `adjusted_similarity_summary()`,
  `pearson_similarity()`, `pearson_similarity_pairwise()`,
  `pearson_similarity_summary()`), `export_stabl_to_csv()`, and
  `stabl_preprocess_fit_recipes()`.  This is a breaking change: 0.3.0
  code calling these names now fails with "could not find function".
  Use the canonical entry points instead -- `stabl_similarity()` /
  `stabl_similarity_summary()` for selection similarity,
  `save_stabl_results(what = "csv")` for CSV export, and
  `stabl_preprocess_fit()` (which accepts a `recipes::recipe` directly)
  for preprocessing; the full name mapping is in
  `docs/MIGRATION_0.3_TO_1.0.md`.
- The two zero-documentation multi-omic diagnostics,
  `stabl_multiomic_study_diagnostics()` and
  `stabl_multiomic_selection_stability()`, stay exported and are now
  demonstrated in a new "Study-diagnostics utilities" section of the
  multi-omic vignette.
- The public export count drops from 46 to 36 -- the final 1.0.0 API
  surface (36 exported functions plus 43 registered S3 methods).
  `inst/release/public-interface-audit.csv` (79 rows),
  `inst/release/removal-ledger.csv` (the ten shim rows flip to
  `removed`), `NAMESPACE`, `_pkgdown.yml`, the migration, API-reference,
  and Python-mapping docs, and the getting-started vignette are updated
  to match; `test-api-shims.R` is removed and the shim legacy-argument
  tests in `test-api-deprecation.R` are retired with the shims.

## Tests (2026-10-21, functional verification audit)

- New test file `test-functional-audit.R` closes the gaps found by the
  covr line-coverage baseline (87.2%) and the integration-chain matrix:
  - end-to-end chains previously covered only piecewise: study tools
    (`stabl_multiomic_study_split()` -> `stabl_multiomic_train_validate()`
    -> `stabl_multiomic_study_diagnostics()`) and artificial features
    (`make_artificial_features()` -> bootstrap scoring ->
    `compute_fdp_plus()`, with a hand-checked FDP+ formula and a
    planted-signal separation check);
  - the unexported cooperative `"validation"` selection branch retained
    for characterization tests (F-CV-04): binomial deviance/auc/class and
    gaussian mae exercised end to end through the internal path, poisson
    deviance verified directly against a dpois-based computation with the
    public family fence asserted explicitly;
  - multinomial late fusion with a validation split
    (`.apply_multiclass_stack_weights()`);
  - registered S3 methods never executed before:
    `family.cv.multiview()` and `plot.multiview()`;
  - `sample_bootstrap_indices()` `class_weights` happy path and error
    branches, and ordinal-learner error branches.
- `inst/release/public-interface-audit.csv` records the new evidence file
  for the eight interfaces it covers.

## Documentation (2026-10-20, reproducibility pointers)

- README Resources and `benchmarks/README.md` now point to the
  `stablr-reproducibility` compendium (pinned submodule + MD5-pinned raw
  atlas + v1/v2 pipeline drivers + committed results + bit-for-bit
  verification evidence) as the way to re-run the aging case study end to
  end; the WS-8 raw-file requirement is documented where the benchmark
  scripts live.

## New features (2026-10-20, single-cell vignettes)

- New vignette `stablr-cebmf`: STABL selection on cEBMF-denoised
  single-cell expression (10x Genomics PBMC 3k fixture; the cEBMF fit
  runs offline with the Python `cebmf_torch` package and is shown as
  unevaluated code).  Includes a raw-vs-denoised selection comparison
  and a `stabl_diagnose()` demonstration of why FDP+ calibration is
  unreliable at tiny feature counts.
- New vignette `stablr-metacells`: STABL selection on metacell-pooled
  expression, with a transparent dependency-free greedy kNN metacell
  partition, majority-label purity QC, and a stratified bootstrap for
  the class-imbalanced outcome.
- New bundled fixtures under `inst/extdata/pbmc3k/` (raw HVG counts,
  Louvain cluster labels, cEBMF factor scores and loadings) with a
  SHA-256 manifest and provenance metadata; regeneration scripts ship as
  `data-raw/prepare_pbmc3k_fixtures.py` and
  `data-raw/prepare_pbmc3k_fixtures.R`.
- The getting-started vignette now cross-links the *Single-cell
  integrations* articles (`stablr-cebmf`, `stablr-metacells`) from its
  expert-interface pointer paragraph.

## Bug fixes (implementation review, 2026-10-20, batch F-8)

- `plot_fdr_graph()`: the FDR threshold reference line now marks the
  stored `fdr_min_threshold_` contract (falling back to the grid recompute
  only when the stored value is absent) — previously it always recomputed
  `thresh_grid[which.min(FDRs_)]`, which can disagree with the threshold
  the fit actually selected.
- `boxplot_features()`: the jitter is now hermetic and deterministic —
  `geom_point(position = position_jitter(width, height, seed = 1L))`
  instead of `geom_jitter()`, which drew from the ambient stream at build
  time and advanced the caller's global `.Random.seed`.
- `scatterplot_features()`: the LOESS smooth is now skipped per panel when
  a panel has fewer than 4 complete points (with a warning naming the
  affected panels); previously `loess` failed or warned opaquely at
  ggplot build time on tiny panels.
- `scatterplot_features()` now fails closed on a non-numeric `y` —
  previously a factor `y` was silently coerced with `as.numeric()`
  (1..nlevels), scrambling the scatter. The `boxplot_features()` /
  `scatterplot_features()` docs also no longer contradict themselves on
  missing features: requested features absent from `x` are silently
  skipped, and an error is raised only when NONE are present.
- Docs: `summary.stabl_interpretation()` and `plot.stabl_interpretation()`
  now carry full roxygen documentation (previously `@export`-only),
  including the empty-selection error; the plot's `aes()` now uses `.data$`
  pronouns instead of `ev$` (the pronoun is registered with
  `utils::globalVariables()` so `rlang` stays a soft Suggests dependency).
- Docs: `plot.stabl_diagnostics()`'s `which` entry for the resample
  distribution curves now correctly describes them as survival curves
  (1 - ECDF), not ECDFs.

## Bug fixes (implementation review, 2026-10-20, batch F-7)

- `stabl_multiomic_stacked()`: a never-observed omic (all-NA column) is
  masked out of every row, so the stacking score does not depend on its
  weight and the reported `Associated_weight` was an arbitrary `runif`
  draw. Such weights are now reported as `NA` with a warning.
  `.weighted_masked_mean()` treats NA weights as zero (a literal NA would
  poison every row sum, since `0 * NA` is `NA`).
- `stabl_multiomic_stacked()` now warns when fewer than 2 prediction rows
  are scorable — previously the loop never updated and the function
  silently returned a degenerate stack (`score = -Inf`, uniform weights,
  all-NA stacked predictions).
- Late-fusion results (`stabl_multiomic_train_validate(late_fusion =
  TRUE)`$late_fusion) now carry `task_type` for ALL modes, not just
  multiclass (return-contract fix).
- Docs: `stabl_multiomic_train_validate()` `@return` now discloses that
  with `late_fusion = TRUE` the top-level `fits`/`selected_*` come from
  the late-fusion full refit (whose `full_per_omic` is stripped from the
  returned `late_fusion` list) and can differ from the
  `late_fusion = FALSE` result on identical input and seed; the
  `lambda_grid` "computed on full data leaks mildly" caveat is extended to
  the OOF stacking path.

## Bug fixes (implementation review, 2026-10-20, batch F-6)

- `.study_outcome_summary()` now fails closed on a dimensioned (matrix or
  data.frame) `y` — previously it passed the `is.numeric()` check and was
  subset as a flat column-major vector, silently scrambling per-study
  summaries.
- `save_stabl_results()` now checks for ggplot2 BEFORE creating the output
  directory, so a missing backend no longer leaves an empty directory
  behind.
- `.lambda_grid_row_labels()` formats mixed-type grids per column (signif
  for numeric, format otherwise) instead of `as.numeric()` on the whole
  row, which coerced non-numeric columns to NA with a warning. All-numeric
  output is unchanged.
- `stabl_multiomic_scheme_comparison()`: ambiguous calls mixing named
  reserved formals (`x_list`/`y`/`study_labels`) with enough unnamed
  positional arguments to shift positional matching now error up front —
  e.g. `(x, y, sl, y = y2)` silently re-bound the positional `y` to
  `study_labels` (empirical correction: the review's duplicate-via-`...`
  scenario is unreachable because R binds exact-named duplicates to the
  formals first). `x_list`/`y` are also validated before `study_labels`,
  so a malformed `x_list` errors informatively.
- `.study_feature_shift()` now counts NA-data SMDs (e.g. an all-missing
  feature in one partition) separately from the masked 0/0 zero-SD case:
  new `n_missing` column alongside `n_zero_sd` (additive schema change).
- `stabl_environment_stability()` validates x/y sample alignment at entry
  (`validate_sample_alignment()`), so a mismatched `y` errors with
  "Sample mismatch" instead of failing deep in the per-environment subset.
- `load_ool_data()` now fails closed on duplicate sample IDs (previously
  `match()` silently took the first) and on non-numeric feature columns.
- Docs: `stabl_multiomic_study_split()`'s balance criterion is
  `factor(as.character(y))`, degenerate for continuous `y` — documented,
  with a recommendation to pass a categorical (or binned) `y`.

## Bug fixes (implementation review, 2026-10-20, batch F-5)

- Nested-CV gaussian downstream regressor: outer-training splits with fewer
  than 3 samples now route to the intercept-only fallback instead of
  crashing in `cv.glmnet` (the `max(3L, ...)` nfolds floor exceeded n;
  empirically glmnet errors uncaught — "y is constant" from 1-sample
  in-fold training — rather than returning the NaN the review predicted).
  Latent: `stabl_fit` itself fails on such tiny splits first.
- `.evaluate_stabl_candidates_inner()` loses its dead `train_ids` formal
  (never referenced; the inner folds carry their own `train_ids`).
- Numeric `strata` whose quantile breaks collapse under heavy ties now warn
  that stratification is effectively disabled (previously silent).
- New test coverage: unstratified (`stratified = FALSE`) nested CV is now
  exercised end-to-end under class skew — both the feasible completion path
  and the typed `stablr_bootstrap_infeasible` fail-closed error when an
  unstratified draw leaves a class bootstrap-infeasible.

## Bug fixes (implementation review, 2026-10-20, batch F-4)

- `stabl_preprocess_fit(transforms = "zv")` now fails closed when the
  zero-variance filter removes ALL features (previously every downstream
  transform fitted on a 0-column matrix and `stabl_fit()` failed deep).
- `make_artificial_features()` now requires `nrow(x) >= 2`: with a single
  sample, `make_rp_features()`'s within-column `sample(col)` would hit the
  length-1 pitfall and silently replace a value >= 1 with a draw from `1:x`.
- `auto_splsda_keepx_grid()` no longer advertises `keepX > p` (with `p = 1`
  the grid spuriously included `keepX = 2`, invalid for mixOmics).
- Fractional `keepX` values in a sPLS-DA `lambda_grid` now error instead of
  being silently truncated by `as.integer()`.
- The MVR knockoff solver now warns when its PSD scaling search finds no
  feasible factor and returns a degenerate all-zero `S` (behavior inherited
  from knockpy; previously silent).
- `make_ranger_adapter()` gains a `family` argument (default `"gaussian"`)
  passed through to the batch adapter — previously the family could not be
  set through this factory.
- Two comment corrections (no behavior change): the MVR degenerate-nudge
  comment claimed "toward upper bound" for a step that is negative when
  `s_j > 1/cd`; the PD-covariance fallback comment claimed Ledoit-Wolf
  shrinkage where corpcor performs eigenvalue clipping.
- Documentation repairs: the README's eight reference links were relative
  (`reference/*.html`) and 404'd on the GitHub repo page — they now point
  at the pkgdown site absolutely; and the hand-written
  `make_ranger_adapter.Rd` / `make_ordinal_adapter.Rd` (stale since the
  learner-factory internalization, invisible to roxygen2) were regenerated
  from source — the staleness both caused a codoc WARNING and broke the
  pkgdown workflow's missing-topics check, which had left the deployed
  site stale.

## Bug fixes (implementation review, 2026-10-20, batch F-3)

- Accessor and input-validation errors are now typed `stablr_error`
  conditions (`stablr_accessor_error` / `stablr_input_error`
  subclasses), consistent with the package's conditions machinery;
  messages are unchanged.
- `.check_fitted_stabl()` now also requires `feature_names` to be a
  character vector matching `nrow(stabl_scores_)` — a malformed object
  previously made `get_feature_names()` return NULL, masquerading as an
  empty support.  A stored out-of-range `hard_threshold` is now
  re-validated on use.
- `get_cooperative_diagnostics()` on a malformed multiomic CV object
  now errors informatively (missing columns named) instead of
  "undefined columns selected".
- `print.stabl_fit()` no longer crashes when `min_fdr_` is set but
  `fdr_min_threshold_` is NULL (`round(NULL)` errored); the two lines
  are guarded separately.
- `stabl()` now errors with guidance when `y` is a *named* length-1
  string naming a column of `x` (previously treated as data, failing
  deep); the column-selector contract is an unnamed string.
- The ranger learner's seed-retry arithmetic no longer overflows to NA
  near `.Machine$integer.max` (which crashed ranger on a retry); the
  retry seed wraps within the valid range and attempt 1 is bit-identical
  to before.  The internal NULL-seed draw is also bounded to the valid
  range.
- `inst/COPYING.cooperative` now discloses all four silent
  upstream-bug fixes in the vendored cooperative backend (IRLS
  step-halving warm-residual form, `lower.limits` propagation, the
  `length(rho)` typo fix, and the `build_predmat` offset crash fix), so
  future upstream diffs are not misread as regressions.
- Removed the dead vendored `coef_ordered` generic and its two
  unregistered, uncalled methods (disclosed in COPYING.cooperative).

## Bug fixes (implementation review, 2026-10-20, batch F-2)

- `stabl_interpret()` permutation importance no longer uses
  `sample(x)` on a column vector: with a single observation a
  positive-integer-valued column would draw from `1:x`.  The permutation
  is now index-based (RNG-identical for n > 1).
- `stabl_interpret()` now validates `loss` against `family` at entry:
  an invalid `loss` previously errored only when the selection was
  non-empty (silently accepted otherwise), and `family = "gaussian"`
  silently treated `"brier"`/`"log_loss"` as squared error — both now
  error with the allowed values.
- The `stabl_interpret()` man page now documents the actual single
  `evidence_type` column (`post_selection_descriptive` /
  `cross_validated`) instead of promising per-dimension labels that do
  not exist; the exported values are unchanged.
- `stabl_replication_profile()` now requires at least 2 samples per
  study at entry (a 1-sample study previously crashed deep inside
  `stabl_fit()`), and the NULL-seed draw can no longer hit
  `.Machine$integer.max` and fail validation.
- The replication held-out refit uses a collision-proof response name,
  so a selected feature literally named `"y"` is no longer silently
  renamed to `y.1` in the model frame.
- `.validate_effects()` duplicate (feature, study) detection no longer
  false-positives on names containing spaces.

## Bug fixes (implementation review, 2026-10-20, batch F-1)

- `stabl_diagnose()` evidence matching now honours its documented
  contract for fits with unknown metadata: a NULL dimension matches
  only "any" rows instead of collapsing the evidence table to zero
  rows (latent defensive-path logic error).
- Two `assumption_notes` ("no threshold available" and "boundary
  minimum") are now single strings; they were split mid-sentence into
  two vector elements by a stray `c()`.
- The FDP+ minima diagnostic no longer crashes on hand-built fixtures
  whose `min_fdr` is not exactly an element of `FDRs`; it reports zero
  minima instead (real fits were never affected).
- `adjusted_similarity_summary()` and `pearson_similarity_summary()`
  now error when fewer than two sets are supplied (previously returned
  NA/NaN silently).
- `jaccard_similarity_matrix()` errors informatively on empty input
  (previously a cryptic `seq_len` coercion error).
- Binary probability metrics (brier/log_loss/AUROC internals) now
  reject NA truth values, consistent with the classification-metrics
  path.

## Bug fixes (implementation review, 2026-10-20, batch E-4)

- The sPLS-DA learner no longer discards the selection mask for ALL
  features when a bootstrap subsample contains one constant column:
  constant columns are dropped for that bootstrap's fit and the mask is
  mapped back, so only the dropped columns are all-FALSE.  Previously
  rare binary features (common in omics) silently biased every stability
  frequency downward.
- The ordinal learner now validates OBSERVED classes (with
  `droplevels()`) instead of declared factor levels, so a bootstrap
  subsample missing a class fails with a clear error instead of
  reaching ordinalNet with unobserved levels.

## Bug fixes (implementation review, 2026-10-20, batch E-3)

- `stabl_multiomic_cv()` now stratifies classification folds on the
  outcome (the ungrouped path was unstratified), and the single-class
  assessment-fold guard applies to ungrouped CV too: infeasible class
  counts abort with typed `stablr_cv_infeasible` at fold construction
  instead of a deep bootstrap error mid-run.
- `stabl_multiomic_train_validate()` now errors when train and
  validation sets share samples (previously silently accepted — a
  leakage foot-gun in the direct API).
- `validate_sample_alignment()` and the outcome-alignment machinery now
  handle data.frame outcomes correctly: row names are the sample IDs.
  Previously `names()` returned column names, making the matrix-like
  branches unreachable and producing a spurious mismatch error.
- The late-fusion downstream glm/lm uses a collision-proof response
  name, so a selected feature named `.y` is no longer silently
  overwritten and dropped from the model.
- Early-fusion cross-view duplicate feature names now fail with an error
  naming the colliding views (previously a generic duplicate-names error
  from deep inside `stabl_fit`).

## Bug fixes (implementation review, 2026-10-20, batch E-2)

- `stabl_fit()` now fails fast when `sample_fraction` implies
  `n_subsamples = 0` (previously surfaced only after the diversity
  guard's redraws with a message blaming class balance);
  `sample_bootstrap_indices()` and `group_bootstrap_indices()` validate
  `n_subsamples >= 1`.
- Stratified bootstrap samplers no longer hit R's `sample(x, 1)`
  length-1 pitfall: a single-member stratum's draws could come from
  `1:x` instead of being that member's row (wrong-stratum indices).
  All shuffles now use index-then-subset via `sample.int`.
- `stabl_multiomic_nested_cv()`'s stratification guard now accounts for
  inner folds being built on outer-*training* splits: configurations
  such as 5 samples/class with `outer_v = inner_v = 5` are rejected at
  entry instead of crashing mid-run.
- Nested-CV candidate selection now fails closed with an informative
  error when every candidate is NA on the selection metric (e.g. R^2
  with a constant outcome) instead of an opaque subscript error.
- `plot_stabl_path()` renders single-lambda grids as points (previously
  `geom_line` dropped the singleton groups, leaving an empty path).

## Bug fixes (implementation review, 2026-10-20, batch E-1)

- `stabl_preprocess_apply()` now fails closed when a feature had no
  finite training values under `center`/`scale` (previously a NaN mean
  was stored and apply silently emitted an all-NA column unless
  `impute_median` was requested).
- `stabl_multiomic_study_diagnostics()`: a feature that is constant
  within both a study and its complement but at *different* values (true
  SMD = ±Inf — a maximal, perfectly separated study shift) now surfaces
  as `Inf` in `smd_max_abs` instead of being masked as NA and counted in
  `n_zero_sd`.
- `stabl_replication_profile()`: held-out metrics are now computed on
  rows complete in the outcome and selected predictors; a held-out study
  with no scorable rows gets `reason = "incomplete_predictors"` instead
  of a silently misaligned or NA metric.
- `stabl_interpret()` now errors when `x` or `y` contain missing values
  (previously `glm` na-omitted while the importance loss used the full
  vector, so all permutation importances were silently NA).

## Bug fixes (implementation review, 2026-10-19)

- **Critical:** `stabl_null_check(design = "independent")` (the default)
  generated no-op nulls: `sample(y)` carried each value's name along, and
  `stabl_fit()`'s name-based outcome alignment exactly undid the
  permutation, so every null replicate fitted the observed data. The
  permutation now keeps names at their original positions. **Any
  null-check result computed with the default design should be re-run.**
  Null-generation failures also now count toward the `max_failures`
  fail-closed budget (previously a systematically failing `null_simulator`
  returned `n_successful = 0` instead of aborting).
- `stabl_interpret()` reported `coefficient = NA`
  (`"effect_not_estimable"`) for every selected feature with a
  non-syntactic name (e.g. `"feat 1"`); coefficients are now mapped
  positionally. A feature named `.y` is no longer silently dropped from
  the effect model.
- `stabl_preprocess_fit()` / `stabl_preprocess_apply()` now reject
  duplicate column names (previously the second duplicate's data was
  silently replaced by the first's at apply time).
- Late fusion (`stabl_multiomic_train_validate(late_fusion = TRUE)`): the
  full-refit per-omic fits now use the same per-omic RNG streams as the
  `late_fusion = FALSE` path (toggling late fusion no longer changes the
  per-omic panels), the per-fold selection fits and downstream predictors
  use per-omic derived seeds, and the recorded provenance seeds now
  replay the fits exactly.

## Breaking changes — lean public API (1.0.0 contract, Phase 13 PR-G)

- The eight pre-0.3.0 deprecation shims are **removed** (they now error):
  `classic_bootstrap_indices()` -> `sample_bootstrap_indices()`;
  `get_feature_names_out()` -> `get_feature_names()`;
  `jaccard_matrix()` -> `jaccard_similarity_matrix()`;
  `adjusted_similarity_values()` -> `adjusted_similarity_pairwise()`;
  `adjusted_similarity_measure()` -> `adjusted_similarity_summary()`;
  `pearson_similarity_values()` -> `pearson_similarity_pairwise()`;
  `pearson_similarity_measure()` -> `pearson_similarity_summary()`;
  `stacked_multi_omic()` -> `stabl_multiomic_stacked()`.
- The duplicate beginner generics `support()` and `feature_names()` (and
  their S3 methods) are removed; use `get_support()` /
  `get_feature_names()`.
- The single-row learner factories `make_glmnet_adapter()`,
  `make_adaptive_lasso_adapter()`, `make_sgl_adapter()`,
  `make_ranger_adapter()`, `make_ordinal_adapter()`, and
  `make_splsda_adapter()` are internalized: select learners via the
  `base_learner` argument of `stabl_fit()`. The grid helpers
  `auto_lambda_grid()`, `auto_rf_threshold_grid()`,
  `auto_ordinal_lambda_grid()`, and `auto_splsda_keepx_grid()` remain
  exported.
- The experimental oblique random forest is **retired**:
  `base_learner = "oblique_rf"` fails closed with an explanatory error, and
  `make_aorsf_adapter()` / `auto_oblique_threshold_grid()` are removed
  (WS-3 failed the predeclared FDP+ calibration bar; archived source and
  tests live in `attic/`). Use `base_learner = "rf"`.
- `plot_roc()`, `plot_prc()`, `boxplot_features()`, and
  `scatterplot_features()` are internalized (use general plotting packages;
  the TCGA vignette shows inline ggplot2 replacements), as are the
  benchmark metrics `fdr_similarity()`, `tpr_similarity()`, and
  `fscore_similarity()` (realized proportions against known truth, not
  general-purpose similarities).
- See `docs/MIGRATION_0.3_TO_1.0.md` for the full migration guide and
  `inst/release/removal-ledger.csv` for the machine-readable removal
  record. Public interface: 54 exports and 43 S3 methods (from 79 and 47).

## Breaking changes (1.0.0 contract, Phase 13 PR-A)
- `stabl_fit()` loses the `explore`/`n_explore` arguments and
  `get_support()` loses its top-N fallback: an empty selection now stays
  empty instead of being silently relabelled as an "explore" set (spec
  §5.4). Feature ranking remains available via `get_importances()`.
- `late_fusion_training = "python_legacy"` is removed from the public
  multi-omic workflows; only the leakage-safe `"oof"` stacking remains.
  The historical in-sample algorithm survives unexported
  (`.late_fusion_multiomic_fit()`) for parity characterization tests.
- `cooperation_selection = "validation"` is removed from
  `stabl_multiomic_train_validate()` and `stabl_multiomic_cv()`: tuning
  on the assessment/validation set biased performance estimates
  (F-CV-04/P4-F24). Only `"cv"` (shared inner folds) is supported.
- `plot_fdr_graph()` argument `fdr_target` is renamed to
  `reference_fdp_plus` (a purely visual reference line), and the plot's
  default title/axis now say "FDP+ diagnostic" instead of "FDR Estimate":
  the FDP+ curve is a decoy-based stochastic-bound diagnostic, not a
  controlled false-discovery rate.

## Fixes
- Evidence-profile provenance corrections in `stabl_diagnose()` (Phase 13
  PR-C follow-up, post-merge provenance review): the adaptive_lasso row no
  longer cites Python-parity coverage that does not exist (the pinned
  solver cases cover lasso and elastic_net only); the sPLS-DA rows now
  reflect the MINT benchmark's FDP+ calibration results (random
  permutation arm calibrated, ratio 1.01; knockoff arms conservative),
  replacing the incorrect "no calibration artifact" status; the lasso
  block-correlated row's regime is narrowed to the evaluated Phase 6 D3
  conditions with the cross-benchmark tension noted; a new row surfaces
  the MINT parent benchmark's severely anti-conservative lasso +
  knockoff_equi result under cohort heterogeneity (calibration 8.0); and
  the stale "68 byte-pinned" fixture count is replaced by the durable
  SHA-256-manifest reference.
- RNG hygiene residuals (Phase 13 PR-A; F-CV-09/10/11). The three
  remaining bare-`set.seed` sites outside the core fit paths —
  `make_artificial_features()`, the permutation-importance engine in
  `stabl_interpret()`, and the one in `stabl_analyze()` — now wrap their
  seeded regions in `.with_local_seed()`, leaving the caller's global RNG
  untouched. Outputs are identical to before.
- `stabl_environment_stability()` derives per-environment seeds via the
  package's hashed multi-stream derivation instead of consecutive
  `random_state + i` integers (F-019); results under a given
  `random_state` change accordingly.
- Nested-CV seed-offset audit (F-CV-05): the offset-collision risk is
  verified unreachable under the production seed-derivation structure
  (proof in `docs/audit/2026-10-17-adversarial/findings.md`); the
  invariant is pinned by a regression test and a call-site comment. No
  numerical change.
- Compiled-boundary fail-closed validation (Phase 13 PR-1). `elnet.fit`
  now length-validates the sparse standardization attributes and rejects
  non-finite `xs` before `.Call` (previously a missing/short attribute
  reached the compiled engine and read out of bounds), rejects `NA` in
  `x` (dense and stored sparse entries), `y`, and weights with an
  informative error (previously affected columns were dropped silently),
  and warns at the boundary with code-specific messages before returning
  `list(jerr=...)` on non-fatal codes. Behavior on well-formed inputs is
  unchanged.
- RNG hygiene (Phase 13 PR-2). `stabl_fit()` no longer permanently
  clobbers the global RNG: `.Random.seed` is saved on entry and restored
  on exit, so later unseeded randomness in the caller's session is no
  longer a deterministic function of the last `random_state`. Nested CV
  results no longer depend on `cv_workers` (per-fold seeds are drawn
  before workers fork), and all multi-omic fan-out sites (per-fold,
  per-omic, early/late/cooperative fusion, downstream predictors) derive
  independent RNG streams instead of reusing one `random_state`. Two
  downstream `cv.glmnet` fits that drew their internal fold assignment
  from the ambient stream are now seeded. Multi-omic numerical outputs
  under explicit seeds change accordingly (per-omic feature selections
  are unchanged); the Python-parity contract for the core STABL
  algorithm is unaffected and its fixtures pass byte-identically.

## Documentation
- New vignette "Reliability diagnostics, validation, and final models"
  (`stablr-reliability`, Phase 13): a runnable end-to-end tour of the
  reliability toolchain — `stabl_diagnose()` warning codes and threshold
  path, `stabl_null_check()` exchangeability designs including a live
  fail-closed refusal, `stabl_replication_profile()` cohort replication,
  and the validated-gaussian development path
  (`stabl_analyze(mode = "validated")`) finalized into a portable
  `"stabl_predictor"` with `stabl_finalize()`/`predict()`. Registered in
  the `_pkgdown.yml` navbar and articles index.
- Two new vignettes (Phase 13 PR-B, 1.0.0 spec §8): "Python STABL parity:
  precisely what matches" (`stablr-python-parity`) replays the shipped
  byte-identical subset of the frozen Python-reference fixtures
  (`inst/extdata/python_parity_vignette/`, SHA-256 manifest) and states
  explicitly what parity does and does not establish; "Reading stability
  paths and FDP+ plots" (`stablr-reading-plots`) is a plot-interpretation
  guide. Both render without Python, using only installed-package
  resources. `_pkgdown.yml` gains navbar and articles entries, the README
  links both articles, and `docs/PYTHON_TO_R_MAPPING.md` is corrected to
  the canonical post-rename accessors.
- `?stabl_fit` documents `min_fdr_` as a minimum-FDP+ diagnostic with the
  stochastic-bound caveat; `?compute_fdp_plus` gains a contract section
  pinning the strict `>` exceedance comparator and first-index tie policy
  (Python `argmin` parity). The intro, multi-omic, and cooperative
  vignettes and `docs/API_REFERENCE.md` are updated for the 1.0.0
  contract changes above.
- All nine vignettes audited against the current API and re-rendered.
  Replaced 25 deprecated `get_feature_names_out()` calls with
  `get_feature_names()`, `jaccard_matrix()` with
  `jaccard_similarity_matrix()`, `adjusted_similarity(nb_total_elements=)`
  with `n_universe=`, and deprecated `groups=` arguments to `stabl_fit()`
  with `sample_groups=` (advanced, environments). Corrected the
  `stabl_multiomic_train_validate()` group argument to `groups_train`
  (multiomic, cooperative) and the stale `groups` note in the
  environments vignette.
- Added seven figure chunks so key conclusions are visually supported:
  environment-wise minimum-score heatmap (environments), Cox stability
  path, Jaccard similarity heatmap and interpretation scatter (advanced),
  late-fusion prediction scatter (multiomic), cooperative rho surface
  (cooperative), and FDR+ regularisation curve (intro).
- Advanced vignette: the grouped-bootstrap and reproducibility examples
  now pass `hard_threshold = 0.5` so the demonstration selects a
  non-empty set at p = 10 (the default data-driven threshold legitimately
  selects nothing there, which made the Jaccard comparison vacuous); the
  Cox example text now states that the FDR threshold saturates at 1, so
  the selection set is empty and the planted features are recovered in
  the importances ranking and the sparse end of the stability path.
- Multiomic vignette: the late-fusion scatter is accompanied by an honest
  note that the compact build settings (10 bootstraps, 6-lambda grids)
  demonstrate workflow shape, not competitive prediction.
- Intro vignette: section 7.3 benchmark table replaced with the current
  learner-benchmark calibration ratios (benchmarks/learner_comparison/),
  and next-steps pointers added to the getting-started vignette and the
  sPLS-DA / random-forest material.
- Getting-started vignette: corrected the `knockoff_equi` "Why" table
  cell -- decoys are conservative for tree learners, and measured
  calibration is learner- and regime-dependent (see the pairing table in
  `vignette("stablr-intro")`).
- sPLS-DA coverage extended: new section in the advanced vignette
  (multinomial data via the binomial adapter), topic-table row, a bullet
  in the TCGA learner notes, and an adapter pointer in the cooperative
  DIABLO paragraph.

## New features
- Gaussian validated analysis (Phase 13 PR-F, 1.0.0 spec A5):
  `stabl_multiomic_nested_cv()` and `stabl_analyze(mode = "validated")`
  now support `family = "gaussian"` end to end.  The complete procedure
  (fold-fitted preprocessing, STABL selection, candidate tuning,
  downstream predictor) is isolated inside each training partition;
  held-out predictions are numeric and performance is reported as MAE,
  RMSE and R-squared with explicit definitions (R-squared is computed
  about the held-out mean, is reported as-is when negative, and is `NA`
  for constant held-out truths).  The candidate-selection `metric`
  resolves per family (`"rmse"` default, also `"mae"`/`"r2"`; binomial
  keeps `"ber"`/`"accuracy"`), stratification uses quantile bins of the
  numeric outcome, study-disjoint folds are supported, and an empty
  selected panel falls back to deterministic intercept-only
  (training-mean) predictions.  `event_class` and
  `return_probabilities` are rejected for regression.
- `stabl_finalize()` (Phase 13 PR-F, 1.0.0 spec A6): finalize a
  completed development result (a `"stabl_fit"`, or a `"stabl_analysis"`
  -- validated analyses must have been produced with
  `final_refit = TRUE`) into a portable `"stabl_predictor"`.  The object
  stores the selected panel, the training-fitted preprocessing transform,
  the event/level mapping, the downstream engine (intercept-only
  fallback, one-feature glm stored as bare coefficients, or seeded
  `cv.glmnet`), the learner specification, feature provenance, and a
  development summary explicitly labelled as a full-data refit distinct
  from any cross-validated performance estimate.  No training
  observations are serialized (one-feature glm objects are reduced to
  coefficients; prepped recipes have their retained training `template`
  stripped).  `predict.stabl_predictor()` applies the stored transform
  without refitting, matches columns by name (order-invariant), and
  supports `type = "response"` (gaussian) or `"class"` / `"probability"`
  (binomial, respecting the stored event mapping).
  `validate_newdata_schema()` enforces the recorded schema: unique,
  named, numeric columns; missing, renamed or duplicated required columns
  are always an error; extra columns are an error unless
  `extra_columns = "drop"` is given explicitly.
- `stabl_replication_profile()` (Phase 13 PR-E, 1.0.0 spec A3): cohort
  replication and heterogeneity profile for multi-study data (>= 3
  studies). Distinguishes four estimands explicitly: (i) pooled
  selection stability, (ii) training-cohort replication across
  leave-one-study-out (LOSO) fits, conditional on per-cohort
  eligibility (measured and non-constant; ineligibility reasons
  reported), (iii) direction/effect replication -- sign concordance,
  inverse-variance pooled estimate/CI, and Cochran's Q heterogeneity,
  computed only from user-supplied per-study effect estimates (the
  package never invents effect estimates), and (iv) held-out cohort
  performance via a documented downstream `stats::glm()` refit of each
  LOSO panel (gaussian MSE / binomial balanced error rate; empty panels
  contribute `NA` with reason). Every per-study and LOSO fit is a fresh
  `stabl_fit()` call on the relevant partition, so preprocessing and
  tuning are recomputed inside training data by construction, and no
  held-out cohort ever chooses a threshold. Uncertainty is reported as
  cohort-level observations over dependent folds -- no binomial CIs.
  Per-fit seeds derive from the recorded base seed (fully replayable).
  S3 `print()` for class `stabl_replication_profile`.
- `stabl_null_check()` (Phase 13 PR-D, 1.0.0 spec A2): outcome-null
  calibration workbench. Repeatedly reruns the complete `stabl_fit()`
  pipeline on outcome-null versions of the data and summarizes the
  selected-panel-size distribution, P(any null selection), per-feature
  null selection frequencies, reported-FDP+ versus by-construction
  realized FDP, null-saturation proportions, and per-replicate
  decoy-generator provenance. Null generation is prescribed by an
  explicit exchangeability unit: `"independent"` (global permutation —
  refused when `blocks`/`clusters` are supplied), `"within_block"`,
  `"whole_cluster"` (equal-sized clusters), or `"simulated"` (an explicit
  user-supplied `null_simulator`; the package never invents a null DGP).
  Unsupported designs are refused, never silently replaced by a global
  permutation. Per-replicate permutation and fit seeds are derived from
  the recorded base seed (fully replayable); `decoy_reps > 1` separates
  sensitivity to random decoy realizations. Failures are recorded with
  replicate ids, seeds, stage, and message, and the run aborts
  (fail-closed) when they exceed `max_failures`. The result
  (`stabl_null_validation`, with `print()`/`plot()` methods) carries
  `truth_status` distinguishing synthetic global nulls from permutation
  nulls of observational data — the latter quantify null selection
  behavior, not FDP validity for partly non-null data.
- `stabl_diagnose()` (Phase 13 PR-C, 1.0.0 spec A1): selection-reliability
  diagnostics for any `stabl_fit()`/`stabl()` result. Reports real and
  artificial score distributions, decoy-generator provenance
  (requested/actual/fallback), real/decoy selection counts along the full
  threshold path and per lambda, the selection rule in effect (user
  `hard_threshold` vs minimized-FDP+), the strict `>` comparator and the
  FDP+ denominator convention, FDP+ curve features (all minimizing
  thresholds, ties, flat minima, boundary solutions, low-decoy-count
  instability), a null-saturation probe (the RF/oblique-RF
  importance-leakage signature, evaluated at the selected threshold —
  including off-grid user thresholds), a score-dependence assumption
  probe (explicitly not an exchangeability or knockoff-validity test),
  and a learner × decoy × outcome × dimensionality evidence profile with
  benchmark provenance (`validated_in_stablr` / `empirically_concerning`
  / `not_evaluated` / `unsupported`; deliberately no binary safe/unsafe
  score). Returns machine-readable `warning_codes`, `assumption_notes`,
  and `calibration_scope`; `print()` and `plot()` S3 methods included.
  The diagnostic is read-only: it never changes the fit or its support.
  `stabl_fit()` objects gain `base_learner` and `family` fields to
  support the evidence profile (additive; also fixes a latent
  `print.stabl()` fallback that always displayed "gaussian"/"lasso").
- New experimental base learner `base_learner = "splsda"`: sparse PLS
  discriminant analysis via `mixOmics` for binary classification
  (`family = "binomial"` only). The STABL lambda grid sweeps sPLS-DA's
  per-component `keepX`; a feature counts as selected when any retained
  component loading exceeds `splsda_ncomp`-independent magnitude
  `1e-10`. New exports: `make_splsda_adapter()` (standalone adapter
  factory, mirroring `make_glmnet_adapter()`) and
  `auto_splsda_keepx_grid()` (geometric keepX grid used by
  `lambda_grid = "auto"`), plus the `stabl_fit(splsda_ncomp = 2L)`
  argument. Empirical calibration on a multi-cohort benchmark
  (benchmarks/mint_comparison/): sPLS-DA + random permutation was the
  best-calibrated configuration (FDP/FDP+ ratio 1.01); knockoff pairings
  should be validated per dataset.

# stablr 0.3.0

## Maturity benchmarks
- Interpretation/validated-workflow benchmark committed (IV, 110 runs,
  predeclared gates): effect-direction recovery 1.00, null validated AUROC
  centred at chance (median 0.473), leaky full-data selection directionally
  optimistic over nested (+0.015 mean AUROC), one-feature folds
  non-degenerate. Redundancy co-selection gate NOT met (0.60): correlated
  substitutes are always grouped when co-selected (10/10), but FDP+
  selection prunes members in 4/10 seeds. Adverse result committed
  (benchmarks/iv_summary.md; Gate B report in docs/plans/).
- Ordinal regression learner (`base_learner = "ordinal"`, `ordinalNet`
  proportional-odds): simulated calibration (10 seeds x 2 designs) gives
  knockoff_equi mean FDP+ 0.00 (nominal 0.1) but mean TPR 0.37 -- strongly
  conservative; graduation bar (FDP+ <= 0.15, TPR >= 0.8) NOT met, so the
  learner stays experimental (benchmarks/ws9_summary.md). (simulation evidence; summaries in benchmarks/)

- Random-forest adapter graduated within a measured scope (WS-2, 60 runs):
  rf + knockoff_equi meets the predeclared criterion (mean FDP <= 0.15,
  TPR >= 0.8) on linear and low-dimensional nonlinear designs; NOT met at
  p = 50 with strong correlation or for binary outcomes. Man pages and the
  pkgdown section document the scope.
- Oblique random-forest adapter fenced (WS-3, 40 runs): no predeclared
  remedy reached the criterion; knockoff decoy importance is inflated.
  Man pages recommend base_learner = "rf" for knockoff-calibrated selection.
- Cooperative-fusion benchmark committed (WS-4, 60 runs): no fusion mode
  dominates; cooperative fusion is competitive; the late-fusion gaussian
  weakness is documented as a caveat. Cox cooperative fusion is permanently
  deferred (no CV-loss precedent; fail-closed stops guard the engine).
- Environment-stability workflow promoted with a documented scope (WS-5,
  60 runs): min-aggregation controls FDP+ across location-shift and
  prevalence-drift designs; per-environment permutation rows are documented
  anti-conservative.
- Knockoff high-dimensional boundary recorded as permanent (WS-7 decision
  record, docs/plans/).
- Learner-regime guidance table shipped in the intro vignette, section 7.3
  (WS-6).
- plot_stabl_path now draws the path in canonical decreasing-lambda order
  regardless of grid row order, aligns artificial-feature lines to the same
  order, and shows x-axis tick labels.

## Release governance
- SCI-04 null-selected-fraction gate accepted on 2026-10-16 by the project
  owner acting in the reviewer role (accepted record
  `SCI-04-null-selected-fraction-v2-accepted`, reviewer response in
  `inst/analysis/review/`). This is a provisional owner acceptance: the
  review packet's itemized independent confirmation remains outstanding.
  The release contract pins the Hoeffding criterion (alpha = 0.05,
  U <= 0.10, n = 150 replicates) and the methodology validation now
  evaluates the bound per null cell instead of failing closed. Remaining
  release blockers: the fresh full-candidate run at the accepted replicate
  count and the independent confirmation.
- IV-benchmark CI smoke workflow (weekly schedule + manual dispatch) runs
  the interpretation/validation benchmark end-to-end at shrunken settings
  and asserts run completeness.

## Documentation

- Documentation website restyled after the dreamlet site (Flatly theme,
  KaTeX math rendering, maintainer links on the authors page); the pkgdown
  workflow now deploys the site on every push to main.
- Homepage restructured to lead with resources, motivation, the analysis
  workflow diagram, and installation; release-governance detail moved to
  `docs/RELEASE_STATUS.md`.
- All eight vignettes rewritten in the formal register of the limma User's
  Guide (numbered task-oriented sections, declarative prose, output
  interpretation); code chunks are unchanged.
- `inst/CITATION` derives the package version dynamically instead of
  hardcoding it.
- pkgdown reference index extended with the knockoff feature constructors.

## Cross-study reliability: study-disjoint resampling (B1–B3)

- **New `study_labels` argument on `stabl_multiomic_nested_cv()`** declares
  the sampling unit of the cross-validation.  When supplied (a named vector
  covering exactly the sample IDs), whole studies become the CV unit: outer
  folds assign whole studies (`outer_v` equal to the number of studies gives
  leave-one-study-out), inner candidate-selection folds are study-disjoint
  within each outer training split, and `stratified = TRUE` is reinterpreted
  at the study level through a deterministic largest-first greedy pass that
  minimises class-proportion imbalance.  The default `study_labels = NULL`
  preserves the sample-i.i.d. scheme byte-identically (`fold_scheme =
  "sample_iid"`).  Study disjointness concerns the CV folds only —
  `bootstrap_strata` (sampling design) and SGL `feature_groups` (features)
  remain orthogonal.
- **Fail-closed validation** for the study-disjoint scheme: label names must
  cover exactly the sample IDs; at least `max(outer_v, 2)` studies; each
  outer-training split needs at least `inner_v + 1` studies; every
  validation fold (outer and inner) must hold both outcome classes; and
  1-sample studies are rejected.  `stabl_multiomic_nested_cv()` also gained
  a majority-class fallback when a fold's selection leaves fewer than two
  features (previously a glmnet error).
- **Provenance and reporting:** the nested-CV object records `fold_scheme`
  (`"sample_iid"` / `"study_disjoint"`) and the validated `study_labels`;
  under the study-disjoint scheme each fold carries `study_ids` and held-out
  `class_counts`, and `print()` shows the scheme.
- **New `stabl_multiomic_study_diagnostics()`** (inferred scope, flagged):
  per-study outcome composition and per-feature absolute standardized mean
  differences versus the pooled rest, quantifying the study-level shift that
  motivates study-disjoint folds.
- **New `stabl_multiomic_study_split()`** (inferred scope, flagged): a
  deterministic study-disjoint train/validation split (greedy, seeded,
  class-balanced) feeding `stabl_multiomic_train_validate()`.
  `stabl_multiomic_cv()` already supports group-disjoint assessment folds
  through its `groups` argument and now records `fold_scheme` =
  `"group_disjoint"` for provenance.
- **Benchmark evidence:** `benchmarks/b1_study_disjoint_cv.R` compares
  held-out BER under sample-i.i.d. versus study-disjoint outer folds on
  multi-study simulations, isolating two batch mechanisms: class-independent
  location shifts on null features (a deliberate null probe — the gap is
  approximately zero there, as expected) and outcome-prevalence drift across
  studies (gap ≈ +0.04 BER over 10 seeds, positive in 9 of 10 seeds — the
  sample-i.i.d. estimate is optimistic precisely when study-level outcome
  composition differs).  Per-fold selection sets are persisted
  (`b1_selected_features.csv`) for the Round C stability-of-stability
  analysis.

## Cross-study reliability: reporting, comparison, and guardrails (C1–C5)

All Round C items are inferred scope (the approved plan flags them for
audit reconciliation):

- **Per-study performance (C1):** `stabl_multiomic_nested_cv()` now returns
  `performance_by_study` when `study_labels` is supplied — per-study outer
  held-out accuracy and balanced error rate pooled across folds and repeats
  (`NA` balanced error rate for studies whose pooled held-out predictions
  contain a single outcome class).  `NULL` under the sample-i.i.d. scheme.
- **Scheme comparison helper (C2):** new
  `stabl_multiomic_scheme_comparison()` runs the same nested CV under both
  sampling units (sample-i.i.d. and study-disjoint) with matched arguments
  and returns both objects plus a comparison table and the BER gap
  (positive = the sample-i.i.d. estimate was optimistic).
- **Harmonization guardrail (C3, documentation):** the CV docs now warn that
  batch harmonization (e.g. ComBat) fitted on the full dataset before
  cross-validation leaks study-level information across folds; fit
  harmonization within training folds or interpret metrics as conditional
  on that step.  Study-disjoint folds do not repair leakage from pooled
  preprocessing.
- **Stability of stability (C4):** new
  `stabl_multiomic_selection_stability()` computes fold x fold Jaccard
  similarity of the selected feature sets across the outer folds of a
  study-disjoint nested CV (the deferred selection-stability heatmap
  analysis; render with `stats::heatmap()`).
- **Grouped-CV guardrail (C5):** `stabl_multiomic_cv()` now fails closed
  when a grouped assessment fold holds a single outcome class under a
  classification family (mirroring the study-disjoint validation);
  regression is unaffected.

## Two-tier API: beginner `stabl()` + expert `stabl_fit()`

- **New `stabl()`:** a one-stop-shop entry point for first-time and
  application users.  Pass a `data.frame` (or a matrix plus an outcome
  vector; the outcome may also be named as a column of the `data.frame`)
  and it runs STABL with statistically safe defaults — model-X knockoff
  artificial features (`artificial_type = "knockoff_equi"`), FDP+
  calibration, strict knockoff fallback (`knockoff_fallback = "stop"`),
  and an automatic lambda grid.  Roughly eight arguments cover the common
  cases; everything else is forwarded to `stabl_fit()` through `...`.
- **New S3 methods for the `stabl` class:** `print.stabl()` (plain-language
  answer: what was fitted, which features were selected, the achieved FDP+
  estimate), `summary.stabl()`, and `plot.stabl()` (stability paths plus
  the FDP+ curve).  `stabl` objects also remain full `stabl_fit` objects,
  so every accessor and workflow applies unchanged.
- **New beginner generics `support()` and `feature_names()`** with methods
  for `stabl` and `stabl_fit`; they delegate to `get_support()` and
  `get_feature_names()`.

## Renamed functions and arguments (old names kept as deprecation shims)

Each old name still works: it warns once per session and delegates to the
new name.  Positional calls are unaffected.

| Old name | New name | Why |
|---|---|---|
| `stabl_fit(groups =)` | `stabl_fit(sample_groups =)` | `groups` collided with `feature_groups`; sample- vs feature-grouping is now unambiguous |
| `get_feature_names_out()` | `get_feature_names()` | drops the sklearn-ism; reads naturally in R |
| `stacked_multi_omic()` | `stabl_multiomic_stacked()` | joins the `stabl_multiomic_*` family and fixes the multi_omic/multiomic casing split |
| `jaccard_matrix()` | `jaccard_similarity_matrix()` | states what the matrix contains |
| `adjusted_similarity_values()` | `adjusted_similarity_pairwise()` | the triplet now reads pair / all-pairs matrix / summary consistently |
| `adjusted_similarity_measure()` | `adjusted_similarity_summary()` | as above |
| `pearson_similarity_values()` | `pearson_similarity_pairwise()` | as above |
| `pearson_similarity_measure()` | `pearson_similarity_summary()` | as above |
| `classic_bootstrap_indices()` | `sample_bootstrap_indices()` | states what is resampled (samples, not features or groups) |
| `pearson_similarity_pairwise(d =)` / `_summary(d =)` | `nb_total_elements =` | aligns the element-count argument with the adjusted-similarity family (`d =` still accepted with a deprecation warning) |

## Documentation

- pkgdown reference index reorganized into tiers: Getting started, Core,
  Workflows, Expert, Metrics, Plotting, Data, and a dedicated Deprecated
  section.
- New vignette "Getting started with stabl" — a beginner walkthrough from
  a `data.frame` to a plain-language selection.
- README gains a two-tier quick start (beginner `stabl()` first, expert
  `stabl_fit()` second).
- All help pages regenerated; every shim's help page points to its
  replacement.

# stablr 0.2.2

## Scientific hardening (audit A1-A3, E1-E3)

- **Strict fallback policy (flagged behavior change):** when a knockoff
  constructor fails on a supported-dimensional input, the previous behavior
  was to warn and silently fall back to random-permutation features — which
  do not support FDP+ calibration, so the fit's FDP+ estimates silently
  changed target.  `make_artificial_features()` and `stabl_fit()` now take
  `knockoff_fallback = "stop"` (the default): the fallback is an error with
  an actionable message.  `knockoff_fallback = "allow"` restores the
  warn-and-fall-back behavior (recorded in the provenance
  `fallback_history`).  The benign `knockoff_mvr` -> `knockoff_equi`
  fallback is unaffected (both are model-X constructions with valid FDP+
  calibration) and remains a warning under either policy.  Related fix: the
  linear-algebra infeasibility classifier now recognizes the hyphenated
  "not positive-definite" and "cannot solve SDP" messages emitted by the
  knockoff package's SDP solver, so those failures are typed and handled
  consistently instead of escaping as raw errors.
- **RF/oblique-RF confirmatory guard extended:** the runtime guardrail for
  multivariate RF learners now also fires for the fixed-X `knockoff`
  pairing (fixed-X exchangeability does not survive bootstrap row
  subsampling), not only for `random_permutation`.  The model-X pairings
  (`knockoff_equi`, `knockoff_mvr`) remain warning-free.
- **Covariance terminology:** the fixed-X knockoff documentation now states
  the actual contract (second-moment/Gram structure on the input scale)
  rather than "covariance structure".
- **SCI-04 feasibility clause (E1):** the draft
  `SCI-04-null-selected-fraction-v2` contract now documents that the
  proposed release criterion `U <= 0.10` (Hoeffding, `alpha = 0.05`) is
  vacuous at the previously considered n = 100 replicates (margin 0.1224 >
  0.10) and predeclares `n >= 150` as a design-time fail-closed
  feasibility requirement, with deterministic reference tests pinning the
  arithmetic.  The gate remains `proposed_unaccepted` pending independent
  review; no locked-run data was used.
- **F5 benchmark expansion (E2/E3):** the RF validation suite now includes
  the oblique-RF adapter (`aorsf`) when installed, giving the adapter its
  first F5 FDP evidence, and the start-of-run knockoff probe verifies the
  0.2.2 shipped-pair scale contract on a heterogeneous-scale design
  (column sds 0.01-100) instead of a standardized one.

## Audit hotfixes (source audit of 0.2.1, PR #3 follow-up)

- **Fixed-X shipped-pair scale contract completed (flagged behavior
  change):** 0.2.1 rescaled the fixed-X knockoffs to the *standardized*
  design but still shipped them next to the caller's *raw* matrix — a
  mismatched pair whenever input columns were heterogeneous.  The
  construction now runs on the raw design (create.fixed normalizes
  internally, scale only) and the knockoffs are back-transformed by the
  input column norms, so the shipped pair `(x, xk)` satisfies the fixed-X
  Gram condition on the input's own scale: `t(xk) xk = t(x) x` exactly and
  `t(x) xk = t(x) x - diag(s)` with `s_j >= 0`.  The real block is
  untouched; the pre-construction centering/scaling step is removed.
  Gram/cross-Gram identities and paired-block scale matching are now
  asserted on the actual `x_augmented` matrix under heterogeneous scales
  (column sds 0.01-100).  Provenance field
  `knockoff_rescaled_to_standardized_scale` is replaced by
  `knockoff_rescaled_to_input_scale`.
- **Windows CI fixtures:** `.gitattributes` now enforces LF (and marks
  `.rds` binary) for the `multiview_parity` fixtures; the Windows
  `R-CMD-check` manifest failures were checkout line-ending conversion,
  not hash drift.
- **pkgdown index:** the `stablr-environments` vignette and the
  `print.stabl_env_stability` topic are indexed; the site builds clean.
- **Oblique-RF tests:** a stale probe (measuring `create.fixed`'s internal
  unit-norm scale against a sd ~ 1 expectation, and gating a *model-X*
  assertion on a *fixed-X* build check) was removed; the `knockoff_equi`
  signal-recovery assertion now runs wherever the knockoff package is
  installed.
- **0.2.1 claim corrections:** the F5 measurement of empirical FDP 0.665
  was for lasso + *model-X equicorrelated* knockoffs (not fixed-X), and the
  "invisible to the learner" causal statement is moderated: for the lasso
  backend, glmnet's internal standardization moderates the selection-level
  impact of a scale mismatch.  The contract violation remains the reason
  for the fix; the 0.2.1 notes above are corrected in place.

# stablr 0.2.1

## Fixed-X knockoff scale and alignment fix (flagged behavior change)

- **Scale contract:** `knockoff::create.fixed()` column-normalizes its input
  and returns knockoffs on that internal scale (entries ~ 1/√n).  The fixed-X
  path previously paired these knockoffs with the sd-1 standardized design,
  violating the fixed-X construction contract (the paired blocks lived on
  different scales).  Knockoffs are now rescaled so that `t(Xk) Xk = t(X) X`
  holds on the paired scale.  (0.2.2 correction: the F5 measurement of
  empirical FDP 0.665 cited here was for lasso + *model-X equicorrelated*
  knockoffs, not fixed-X; and for the lasso backend the selection-level
  impact of a scale mismatch is moderated by `glmnet`'s internal
  standardization — the contract violation, not a measured FDP figure, is
  the reason for this fix.)  Fits using fixed-X knockoffs
  (`artificial_type = "knockoff"`) may select different features, and FDP+
  estimates on such fits change.
- **Column alignment:** artificial features are now shipped in aligned order
  (injected column k is the artificial counterpart of feature k's selected
  subset).  The previous random column permutation broke the knockoff
  swap-symmetry in the bootstrap design, because stability scores are
  accumulated positionally.  Applies to the fixed-X and model-X knockoff
  types; random-permutation nulls are unaffected.
- **Probes and tests:** the knockoff sanity probe (tests and F5 benchmark)
  now probes the shipped contract (rescaled scale + Gram identity) instead of
  the internal normalized scale; the formerly CI-gated knockoff quality and
  calibration assertions now run wherever the knockoff package is available.
  The end-to-end calibration test was redesigned to the bootstrap-valid
  model-X regime (`knockoff_equi`) with an explicit lambda grid and powered
  signals; fixed-X knockoffs under row subsampling remain documented as an
  invalid combination and are not asserted.

## New experimental oblique random-forest learner (plan item D2)

- `base_learner = "oblique_rf"` via `make_aorsf_adapter()` /
  `auto_oblique_threshold_grid()` (aorsf backend, negation importance),
  mirroring the random-forest adapter: importance-threshold penalty path,
  per-tree seeds drawn from the caller's seeded RNG stream, and the same
  model-X-knockoff-pairing guardrail.  As with `rf`, selection frequencies
  are noisier than the glmnet family's; the adapter is experimental and
  FDP+ calibration is empirical (F5-style benchmarks).  aorsf is a Suggests
  dependency; the feature is guarded by `.require_pkg()`.

## Multiview-parity fixtures (plan item F4)

- SHA-256-pinned fixtures (`tests/testthat/fixtures/multiview_parity/`) lock
  the vendored cooperative-learning backend's outputs to CRAN multiview 1.0
  on gaussian, binomial, and weighted fits, plus the documented B1/B2/B3
  divergence sites (offset handling in `build_predmat`, warm-start residual
  weighting, `lower.limits` clamping).  Upstream `build_predmat` errors with
  an offset (divergence D4 confirmed live); the vendored fix is pinned by
  fixture.  Makes the PATCHES.cooperative.md divergence log testable.

# stablr 0.2.0

Breaking renames ship with deprecation shims that warn for one minor release
before removal.

## Scientific correctness

- **Behavior change (flagged):** the fixed-X knockoff path now standardizes
  the input matrix before calling `knockoff::create.fixed()`, matching that
  function's documented contract. Fits on heterogeneous feature scales may
  select different features. Standardization details (means, sds, zero-
  variance handling) are recorded in the fit provenance, and the fixed-X
  documentation now states the bootstrap-reuse caveat.
- **FDP+ communication:** new "What FDP+ is and is not" vignette section and
  documented interpretation caveats on `get_fdrs()`.

## Vendored cooperative-backend fixes

- **Offset subsetting in `build_predmat`:** the cross-validation prediction
  path now subsets the `offset` argument per fold instead of referencing an
  undefined symbol. The path was previously unreachable because offsets are
  not exposed through the public API.
- **Halved-step residual weighting:** after step-halving, the warm-start
  residual now weights the offset consistently with the non-halved branch.
- **`lower.limits` symmetry:** the cooperative path now passes
  `lower.limits` to `elnet.fit()` alongside `upper.limits`.
- **Divergence log:** `inst/PATCHES.cooperative.md` records every deviation
  from CRAN multiview 1.0.

All backend fixes are inherited from CRAN multiview 1.0 and are
behavior-neutral under unit weights, the only weighting stablr exercises
today.

## API consistency

- **`seed` deprecated:** `classic_bootstrap_indices()` and
  `group_bootstrap_indices()` now take `random_state`; `seed` still works but
  warns, and will be removed in the following minor release. The helpers no
  longer clobber the global RNG state (`withr` promoted from Suggests to
  Imports).
- **Similarity-metric argument renames:** `adjusted_similarity()` and
  `pearson_similarity()` (and their `_values`/`_measure` variants) rename
  their arguments to `set1`/`set2`/`n_universe`; old names keep working
  through deprecation shims.
- **New accessors:** `get_fdrs()` and `get_lambda_grid()` expose fitted FDP+
  estimates and the lambda grid, with documented interpretation caveats.
- **New S3 methods:** `plot.stabl_fit()` and `summary.stabl_fit()` (with
  `print.summary.stabl_fit()`) for fitted STABL objects.
- **Internal:** late-fusion internal formals are dot-prefixed so user
  arguments passed through `...` cannot collide with them; no user-visible
  impact.

## New capabilities (extensions beyond Python STABL, labeled experimental)

- **Random-forest base learner:** `base_learner = "rf"` via
  `make_ranger_adapter()` and `auto_rf_threshold_grid()` — permutation-
  importance threshold paths so the full STABL machinery runs unchanged on
  non-linear signals. Random forests must be paired with model-X knockoff
  artificial features; the anti-conservative random-permutation pairing now
  warns. Calibration is validated empirically by the F5 benchmark
  (`benchmarks/f5_rf_validation.R`).
- **Environment-stability workflow:** `stabl_environment_stability()` fits
  STABL per environment (batch, cohort, site) and aggregates selection
  frequencies by their minimum across environments — an ICP-inspired
  transportability heuristic without causal assumptions or new dependencies —
  with `get_min_scores()`, a print method, and a new vignette.

## Tests and documentation

- Golden-path rendering tests for all four plotting functions (plot-data
  correctness: rows, mappings, threshold geometries, FDP+ curve values).
- Lambda-grid leakage guidance for cross-validated multi-omic workflows, a
  note on the classification-bootstrap subsample deviation, dead-code removal
  in the MVR knockoff solver, and the release-evidence glossary moved from
  the repository root to `docs/CONTEXT.md`.

# stablr 0.1.1

## Release correctness

- Late-fusion weights now default to leakage-safe out-of-fold training. Each
  stacking fold reruns STABL selection and its downstream learner using only
  fold-training samples; `late_fusion_training = "python_legacy"` preserves
  the historical Python-compatible in-sample algorithm.
- Named outcomes, groups, bootstrap strata, and validation outcomes are
  canonically aligned to omic sample identifiers before fusion.
- OOF provenance records named fold assignments, seeds, selected features,
  artificial-feature metadata, warnings, fallbacks, and full-refit diagnostics.
- Scalar stacking now rejects malformed outcomes and infinite observed
  predictions instead of allowing invalid metric inputs.
- Binary event levels and multiclass probability columns are now validated
  explicitly; probability rows must be complete, finite, and normalized, and
  metric calculations no longer silently discard malformed observations.

## Architecture, profiling, and calibration hardening

- **Fixed-X knockoff correctness:** fixed-X construction now rejects row-
  augmented knockoffs and uses the documented random-permutation fallback.
  Truncating only the knockoff matrix to the original row count violated the
  fixed-X Gram identities when `p < n < 2p`; fallback provenance records the
  reason.
- **FDP+ performance and contracts:** threshold counts now sort once per score
  vector and use interval lookup, preserving strict `>` tie semantics and
  first-minimum threshold selection while avoiding feature-by-threshold logical
  allocations. Public inputs now reject malformed, non-finite, out-of-range,
  or column-incompatible score matrices.
- **Learner-module depth:** public single-lambda glmnet and adaptive-lasso
  adapters now delegate to the same path-fitting implementations used by the
  batched bootstrap interface, removing duplicate fitting and coefficient-
  extraction logic without changing seeded fitted outputs.
- **Methodology runner:** artificial-feature strategies are compared on paired
  simulated data and fit seeds. Replicate artifacts record both seeds, and
  summaries report SD and Monte Carlo SE for selection count, empirical FDP,
  and TPR, plus Monte Carlo SE for FDP-exceedance and fallback rates.
- **ROC/PRC diagnostics:** `plot_roc()` and `plot_prc()` reject ambiguous or
  incomplete labels and probabilities, and aggregate tied scores so metrics
  are row-order invariant.
- **Python parity evidence:** reproducible Python 3.11 fixtures pin the public
  STABL reference commit and include frozen matrices, bootstrap schedules,
  artificial features, lambda grids, masks, FDP+ curves, stacking candidates,
  solver-ranking contracts, provenance, and a SHA-256 manifest.
- **MVR knockoffs:** solver checks now cover PSD feasibility, objective
  behaviour, determinism, and caller RNG isolation. Fallbacks and approximate
  high-dimensional chunking carry explicit provenance.
- **Release infrastructure:** expanded the R CMD check platform matrix,
  installed-package methodology runner resolution, always-uploaded check
  artifacts, one-pass Cobertura/Codecov reporting, release-gated pkgdown
  deployment, the canonical manual R-hub v2 workflow for sanitizer and
  specialist-platform checks, GPL and vendored-code provenance, and
  reproducible bundled OOL data preparation.
- **Metadata and documentation:** set the release version to 0.1.1, refreshed
  the citation and method DOIs, credited the STABL and multiview contributors,
  recorded the bundled OOL data source and CC0 terms, and replaced "Pure-R"
  with the accurate "R-native without a Python runtime" description.
  Documentation now describes FDP+ as a data-dependent diagnostic minimizer,
  states its assumptions and fallback limitations, and distinguishes frozen-
  score parity from end-to-end solver ranking parity.
- **Release evidence:** added deterministic tests for sample-order invariance,
  grouped OOF separation, validation-outcome independence, fold perturbations,
  sequential/parallel agreement, legacy stacking parity, MVR feasibility, and
  locked statistical release gates. Validation commands, versions, hashes,
  warnings, profiling limitations, and platform results are recorded in
  `docs/release/0.1.1-validation.md`.
  Locked long-running validators now require a clean Git source tree, record
  both start and end provenance, and reject evidence if the tracked tree moves
  during computation.

# stablr 0.1.0

## Publication-readiness blocker remediation

- **License:** corrected `License:` field from `MIT + file LICENSE` to `GPL (>= 2)`
  and removed the MIT `LICENSE` file. The vendored cooperative-learning engine
  under `src/glmnetpp/` is derived from glmnet/multiview (GPL-2); the compiled
  combined work must therefore be GPL. `inst/COPYING.cooperative` already
  documented the GPL-2 component.
- **Repository URLs:** updated `URL:` and `BugReports:` in DESCRIPTION and the
  package URL in `inst/CITATION` from the upstream Python repo
  (`gregbellan/Stabl`) to the canonical R package repo (`mdmanurung/stablr`).
- **Contributors:** added the original STABL method authors (Hédou, Marić, Bellan,
  Gfeller — Nat Biotechnol 2024) and the cooperative-learning/multiview authors
  (Ding, Li, Narasimhan, Tibshirani — PNAS 2022) as `ctb` in `Authors@R`,
  consistent with the already-correct `inst/CITATION` entries.
- **Cox in late fusion:** `late_fusion = TRUE` with `family = "cox"` now throws
  a clear error instead of silently scoring with R-squared and OLS predictors
  (both invalid for censored data). Early-fusion and per-omic cox remain valid.
  Mirrors the existing cooperative-fusion cox guard.
- **README:** documented the C++17 compiler / toolchain requirement (Rtools on
  Windows, Xcode CLT on macOS, gcc/clang on Linux) and added a canonical
  `remotes::install_github()` install path.
- **Vignettes:** added reader-facing `eval = FALSE` callout blocks to the three
  heavy vignettes (`stablr-multiomic`, `stablr-tcga`, `stablr-tcga-nestedcv`)
  explaining that code is illustrative and how to run it.
- **MCC metric:** added Matthews Correlation Coefficient (Gorodkin 2004 multiclass
  formula) to the list returned by `.classification_metrics()`. MCC is a standard
  metric for imbalanced classification, complementing the existing accuracy, BER,
  and macro-F1. Returns 0 for degenerate (all-one-class) predictions.
- **Elastic-net parity disclosure:** the Python-parity fixtures for elastic-net
  models assert *signal ranking* (which features are selected), not bit-identical
  coefficients. R and Python elastic-net implementations use different coordinate
  descent schedules; ranking parity is the appropriate criterion and is explicitly
  documented in the test file.
- **Verification:** `mixOmics` (Suggests) is used only in TCGA vignettes and is
  already guarded by `requireNamespace()` + global `eval = FALSE`; no change
  needed.
- **Release metadata:** aligned pkgdown URLs, article listings, and the public
  reference index with the current `mdmanurung/stablr` repository. Removed the
  stale Python-parity article link; the Python-to-R mapping now lives in
  `docs/PYTHON_TO_R_MAPPING.md`.
- **Parallel tests:** the parallel RNG determinism test now probes
  `future::multisession`/PSOCK availability and skips with a clear reason in
  socket-restricted check sandboxes.
- **Release hardening:** added shared scalar validators, duplicate sample/feature
  name guards, effective FDP+ decoy scaling, `transform_stabl()` for selected
  feature subsetting, `matrixStats`-backed row maxima, bundled Python-parity
  fixture checks, fixed-seed methodology/stacking scripts, validator rollout
  across workflow/plot/export helpers, and parity-preserving chunked
  `stacked_multi_omic()` evaluation for binary/regression tasks. The multiclass
  stacking path remains scalar by design to preserve probability-normalization
  parity.
## Quality and publication-readiness sweep

- **Infrastructure:** corrected `License:` field; lowered `Depends: R (>=
  4.1.0)` from 4.4.0; moved package-level doc to `R/stablr-package.R`; guarded
  `save_stabl_results()` against a missing ggplot2 Suggests dep; removed stray
  `Rplots.pdf` build artifacts.
- **Documentation:** documented all four `print.*` S3 methods; added `@examples`
  to all 54 exported symbols; replaced duplicated `@param` blocks with
  `@inheritParams stabl_fit`; added `@seealso` cross-links across the metrics
  family and `compute_fdp_plus`.
- **Correctness:** guarded cooperative lambda-index lookup against float
  mismatch; fixed integer overflow in `.derive_nested_seed` for large fold
  counts; added explicit `stop()` default to the cooperation `type_measure`
  switch; replaced `seq(start+1, start+n)` with `seq_len(n)` idiom at all index
  ranges; added tie-policy comments at every `which.min`/`which.max` site.
- **Performance:** centralised `rowMaxs` helper with edge-case fast paths for
  0-row and 1-column matrices; eliminated quadratic `c(accumulator, chunk)`
  vector growth in three bootstrap helpers; extracted
  `.BOOTSTRAP_COEF_THRESHOLD` named constant; halved pairwise-similarity call
  count in `jaccard_matrix`, `adjusted_similarity_values`, and
  `pearson_similarity_values` via upper-triangle loops.
- **Internal quality:** extracted `.require_ggplot2()` to replace six identical
  inline guards; collapsed `.cooperative_selected_features` two-pass loop to a
  single `lapply`; vectorised the finite-check in `.multiclass_log_loss`;
  extracted `.scores_to_long_df` to deduplicate the real- and
  artificial-feature long-data-frame loops in `plot_stabl_path`.
- **Tests:** added tests for `get_stabl_scores` and `load_ool_data`; modernised
  `skip()` to `skip_if()` idiom. Current local no-vignette release check:
  FAIL 0 | WARN 0 | SKIP 1 | PASS 1611.
  No-manual/no-vignette `R CMD check` is `Status: OK`; all six source
  vignettes rendered to `/tmp/stablr-vignette-review`; pkgdown built to
  `/tmp/stablr-pkgdown`.
- **CI:** added `.github/workflows/R-CMD-check.yaml` (ubuntu release+devel,
  macOS release) with a separate covr/Codecov coverage job; gated heavy
  vignettes (`stablr-multiomic`, `stablr-tcga`, `stablr-tcga-nestedcv`) with
  `eval = FALSE` in their setup chunks.

## Native cooperative learning (multiview engine)

- Vendored CRAN `multiview` v1.0 cooperative-learning engine (GPL-2 submodule in
  `inst/COPYING.cooperative`) for gaussian and binomial cooperative fusion.
- No runtime dependency on the external `multiview` package for cooperative
  workflows.
- Compiled C++ core (`wls_exp`, glmnetpp headers) with `NeedsCompilation: yes`.
- New `vignette("stablr-advanced")` documents Cox, knockoffs, grouped bootstrap,
  reproducibility metrics, export helpers, and outer CV.

## Application-note release candidate

- R-native STABL implementation without a Python runtime, with glmnet backends (lasso, elastic net,
  adaptive lasso, optional sparse group lasso).
- FDP+ thresholding with random-permutation and knockoff artificial features.
- Multi-omic train/validation, outer CV, and nested CV workflows.
- Early, late, and cooperative fusion integration paths.
- Bundled Onset of Labor (OOL) tutorial subset in `inst/extdata/`.
- Python parity fixtures and vignettes for cross-language validation.
- Publication reproduction script: `inst/analysis/generate_publication_parity_table.R`.

## Naming

The package is **`stablr`**. Cooperative multi-omic fusion is a workflow inside
`stablr`, not a separate package (`costablr` is not used).
