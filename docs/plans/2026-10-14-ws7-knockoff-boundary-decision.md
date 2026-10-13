# WS-7 decision record: the 3,000-predictor knockoff boundary is permanent

Date: 2026-10-14. Workstream: WS-7 of the 0.3.0 feature-maturity roadmap.
Status: DECIDED — the boundary stays; no chunked construction will be
reintroduced.

## Decision

`stabl_fit()` rejects knockoff artificial features above
`.maximum_supported_knockoff_features <- 3000L` predictors with the
fail-closed error class `stablr_high_dimensional_knockoff_unsupported`.
This boundary is permanent. The former chunked knockoff construction was
removed because it did not establish a valid global knockoff design, and no
replacement construction is planned for 0.3.0.

## Grounding

- Model-X knockoffs require a feature distribution model; validity of the
  selected-set guarantees depends on the knockoff design as a whole, not
  per chunk (Candès, Fan, Janson, Lv 2018, doi:10.1111/rssb.12265).
- Tractable exact constructions at scale are limited: metropolized knockoff
  sampling (Bates et al. 2020, doi:10.1080/01621459.2020.1729163) and MRC
  knockoffs (Spector & Janson 2020) target specific structures; neither
  gives a globally valid cheap construction for arbitrary p > 3000 in a
  bootstrap loop.
- Robustness of inference degrades when the knockoff construction is only
  approximate (Barber, Candès, Samworth 2020, doi:10.1214/19-AOS1852).
- Blain et al. 2024 ("When Knockoffs fail") diagnose non-exchangeability
  under high equicorrelation with p > n — the regime a chunked shortcut
  would silently enter.
- SDP-based constructions are computationally infeasible at p > 3000 inside
  a bootstrap loop that refits knockoffs per resample.

## Consequences

1. Users above the boundary keep two calibrated-in-scope paths:
   `random_permutation` decoys (no dimensional restriction; anti-conservative
   by construction and documented as a guard) and feature pre-screening
   before the stabl run (documented in the intro vignette).
2. The error is fail-closed with a named class, so downstream automation can
   detect it.
3. If a literature-backed globally valid construction for large p appears,
   reopening this decision requires a new predeclared benchmark cycle.

## Disposition ladder application

The boundary does not worsen default behaviour (it prevents an invalid
silent construction), is backed by the literature above, and passed its
packaged tests. Disposition: KEEP as a permanent documented boundary.
