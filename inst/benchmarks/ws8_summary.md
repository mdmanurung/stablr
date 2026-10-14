# WS-8 — Aging-atlas application showcase summary

- Data: 2609 donors, 60 features (59 cell types + covariates), 8 cohorts.
- Primary (rf + knockoff_equi, 200 bootstraps, classic bootstrap): 24/60 features pass min-agg (>= 0.5).
- Paper-consistent features recovered by min-agg (7/8): CT:MAIT, CT:CD4+ RTE, CT:CD8+ RTE, CT:CD8+ T Naive, CT:CD8+ Temra, CT:CD4+ Temra, CT:CD8+ Tem GZMK+.
- LOCO age MAE: all-feature ridge 12.14 y (genuine held-out evaluation); selected-feature model 11.72 y (POST-SELECTION estimate: the 24 features were chosen with all eight cohorts, so this is apparent for the selection procedure, not unbiased complete-pipeline LOCO). Paper clock: 8.2 y full / 8.9 y reduced; protocols differ -- ours is donor-level ridge on cohort-held-out data, not their trained clock.
- Lasso sensitivity min-agg selected: see ws8_results.csv (selected_lasso column).

## Per-cohort LOCO MAE

| cohort | n | MAE all (held-out) | MAE selected (post-selection) |
|---|---|---|---|
| ABF300 | 166 | 9.30 | 9.43 |
| AIDA | 625 | 14.60 | 14.11 |
| CIMA | 421 | 15.60 | 7.62 |
| KCL | 188 | 9.67 | 10.27 |
| OneK1K | 981 | 15.47 | 17.65 |
| SoundLife | 96 | 6.75 | 9.56 |
| SPAC | 33 | 12.33 | 14.32 |
| UCSF | 99 | 13.41 | 10.83 |

Application case study — NOT calibration evidence.

Evidence labels (added 2026-10-15): the full-atlas min-agg selection is
EXPLORATORY. The all-feature LOCO column is a genuine held-out evaluation.
The selected-feature LOCO column is a POST-SELECTION estimate: the panel
was chosen with all eight cohorts, so its error is optimistic for the full
workflow. A selection-isolated LOCO that re-selects features inside each
seven-cohort training partition is a separate correction (PR5). Historical
MAE values, figures, and datasets are preserved unchanged; no rerun is
implied here.

