# WS-8 v2 — Claims Ledger

Every narrative claim, its source table/artifact, the exact comparison, status
(`supported` / `tentative` / `unsupported`), and limitations. Updated at each
milestone. Initialized at M0 (2026-10-15); no v2 outcome claims exist yet.

## M0 audit findings (source-code and historical-artifact claims)

| # | Claim | Evidence | Status | Limitation |
|---|---|---|---|---|
| A1 | v1 `mean_score` is entirely missing | `ws8_results.csv`: `mean_score` NA for 60/60 rows; `colMeans(env_scores)` returns 8 cohort means, name-indexed by 60 features → all NA | supported (artifact) | none |
| A2 | v1 `selected_lasso` is entirely missing | `ws8_results.csv`: `selected_lasso` NA for 60/60 rows; `as.numeric()` strips names before name-indexing | supported (artifact) | v1 vignette "lasso selects 6, only CD8+ T Naive overlaps" is NOT verifiable from the committed CSV |
| A3 | v1 conflates macro and micro MAE | script `mean(loco_all$mae)` = 14.03 (micro, donor-weighted); vignette 12.14 = macro (unweighted cohort mean) | supported (recomputed) | two estimands presented as one |
| A4 | SPAC knockoff covariance is a PD projection | `.estimate_pd_sigma()` uses `corpcor::make.positive.definite`; SPAC n=33 < p=60 rank-deficient | supported (code) | PD repair ≠ knockoff validity |
| A5 | Paper design matches | Terekhova 2026 Immunity: 2,609 donors, 8 cohorts, 59 cell types, ~12M cells; clock = 55 full / 20 robust; association rule = per-cohort `% ~ age + sex + ancestry`, BH, consistent ≥5/8 | supported [1,2] | exact 20/55 lists in Table S1, paywalled, not retrieved |
| A6 | Protocol §3.2 "5/8 frequency vote → 59/59" | inconsistent with `env_scores` (many features < 0.5 in individual cohorts, e.g. KCL 0.265) | unsupported (do not propagate) | the 59/59 rule is unspecified; treated as unverified historical claim |

## v2 outcome claims

(none yet — populated at M3/M4/M5/M6)

## References

[1] Terekhova et al. (2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012
    (cell.com fulltext S1074-7613(26)00387-0).
[2] KCL Pure record, https://kclpure.kcl.ac.uk/portal/en/publications/single-cell-analyses-of-global-cohorts-outline-determinants-of-im/
