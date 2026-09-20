# Methodology Validation Experiments

`run_methodology_validation.R` runs deterministic simulation experiments that
evaluate methodological concerns without changing package runtime behavior.
Within each scenario and replicate, every artificial-feature strategy receives
the same simulated data and fit seed. This paired design isolates strategy
differences from Monte Carlo differences in the generated data or bootstrap
stream.

## Hypotheses

1. FDP+ thresholding should avoid all-feature collapse under null data, including
   correlated predictor settings.
2. In signal settings, empirical FDP and TPR should quantify the calibration and
   power tradeoff for each artificial-feature strategy.
3. Knockoff fallback frequency should be observable from captured warnings and
   from the additive provenance now stored on fitted objects.
4. Bundled Python metric fixtures should remain numerically identical to exported
   R metric functions.

## Default Run

```sh
Rscript inst/analysis/run_methodology_validation.R --out /tmp/stablr-methodology
```

Defaults are bounded for local development: 3 replicates, 12 bootstraps per fit,
6 auto-lambda values, the Gaussian family, and
`random_permutation,knockoff_equi` artificial features across six scenarios:

- `null_independent`: low-dimensional independent null predictors.
- `null_correlated`: low-dimensional correlated null predictors.
- `signal_independent`: low-dimensional independent predictors with 5 strong
  signals.
- `signal_correlated`: low-dimensional correlated predictors with 5 signals.
- `null_high_dim`: high-dimensional correlated null predictors.
- `signal_high_dim`: high-dimensional correlated predictors with 5 signals.

Scale up with CLI flags such as `--replicates`, `--n-bootstraps`, `--n-lambda`,
`--workers` (Unix-alike independent-replicate process workers),
`--families gaussian,binomial,multinomial,cox`,
`--artificial-types random_permutation,knockoff,knockoff_equi,knockoff_mvr`,
and `--scenarios`.

The locked release profile cannot be downscaled through CLI overrides. It runs
100 paired Monte Carlo replicates per cell, 1,000 bootstraps, and 30 lambdas for
all four advertised outcome families and all four artificial-feature strategies:

```sh
Rscript inst/analysis/run_methodology_validation.R \
  --profile release --workers 32 --out /tmp/stablr-v0.1.1-methodology
```

The release profile requires an identifiable, clean Git source tree before it
starts. The runner records both start and end commit/tree state and rejects the
artifacts if tracked provenance changes during the run. The late-fusion release
runner applies the same invariant.

Release gates are evaluated separately for every
family/scenario/artificial-strategy cell. Missing, incomplete, skipped, or
errored cells fail. The accepted `null_select_any` gate uses its one-sided
Wilson bound. The distinct `null_selected_fraction` gate is SCI-04: its
historical Wilson calculation was invalid for a continuous replicate-level
fraction. The owner-authorized `release_v2` implementation uses 600 replicates,
a one-sided per-cell Hoeffding upper bound at alpha 0.05, a threshold of 0.10,
and an intersection-union decision over all applicable null cells. Eight
fixed-X/high-dimensional cells are declared not applicable. The implementation
remains pending an independent authenticated signature and is not wired into
the active release contract; the active v1 gate therefore still fails every
promotable packet. No replicate may select at least 90% of features. For each
strong-signal cell, the one-sided 95% upper bound for mean empirical FDP must
be at most 0.12 and the corresponding lower bound for mean TPR must be at least
0.50.

## Artifacts

All artifacts are written under the user-specified `--out` directory.

- `methodology_validation_replicates.csv`: one row per scenario, replicate, and
  artificial-feature strategy. Columns include the paired data and fit seeds,
  scenario metadata, warning counts, selected-feature counts, true positives,
  false positives, empirical FDP, TPR, FDP+ minimum and threshold,
  real/artificial stability-score summaries, selected feature names, and any
  error text.
- `methodology_validation_summary.csv`: grouped means by scenario and artificial
  type, including SD and Monte Carlo SE for selected-feature count, empirical
  FDP, and TPR; empirical FDP exceedance rate relative to `--target-fdp` with
  its Monte Carlo SE; fallback rates with Monte Carlo SE; and elapsed time.
- `methodology_validation_warnings.csv`: captured warning messages with scenario,
  replicate, artificial type, and warning index. Random-permutation fallback and
  MVR-to-equi fallback rates are derived from this table; fitted objects also
  expose `artificial_provenance` for direct per-fit reporting.
- `python_metrics_parity.csv`: observed R metric values compared against bundled
  Python parity fixtures with absolute error and `ok`/`mismatch` status.
- `methodology_validation_manifest.txt`: run settings, hypotheses, start/end Git
  provenance, source-stability status, artifact paths, and artifact SHA-256
  hashes.
- `methodology_validation_gates.csv`: per-cell gate bounds, completeness status,
  criteria, and pass/fail results.

## Non-Additions

The runner deliberately does not change FDP+ thresholding or promote any new
artificial-feature strategy to a default. The current fitting rule chooses the
minimum of the diagnostic FDP+ curve. A plotted `fdr_target` is a visual
reference only: it is not used by fitting and is not a universal false-discovery
guarantee. Methodological changes should be evaluated with these artifacts
before they are added to package runtime paths.
