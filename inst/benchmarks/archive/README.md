# Archived benchmark outputs

Files here are preserved historical outputs that were superseded by a
corrected rerun. They are kept for provenance and must not be cited as
valid evidence.

- `ws9_results_precorrection.csv`, `ws9_summary_precorrection.md`
  (archived 2026-10-15): WS-9 ordinal benchmark outputs generated with the
  incorrect truth set (`paste0("f", seq_along(beta))`, all 10 predictors
  instead of the 3 nonzero coefficients) and an incomplete run grid (38/40;
  seed 8 missing in both arms of the correlated design). Empirical FDP/TPR
  in these files are invalid. See the dated erratum in
  `../MATURITY_PREDESIGN.md` and the corrected `../ws9_results.csv`.
