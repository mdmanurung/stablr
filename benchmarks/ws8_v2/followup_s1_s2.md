# WS-8 v2 — Deferred follow-up plan: S1 simulations and S2 sensitivity probes

Status: deferred, non-blocking. The v2 package (protocol, amended results,
claims ledger V1–V10, vignette Section 10) is complete and merged (PR #15).
This note scopes the two probes that were explicitly postponed so they can be
picked up without re-deriving context. Neither probe changes the merged
claims unless its results contradict them; the falsification hooks in
section 4 say what to look for.

## 1. S1 — Synthetic calibration simulations

**Why.** The real atlas has no gold standard. S1 asks whether the R(K)
endpoint and the M0–M5 selector matrix behave correctly when the truth is
known, and whether M4's loss to M0/M1 is specific to this atlas or
structural (e.g., the min-aggregation veto under cross-cohort heterogeneity).

**Design.** Simulate 8 synthetic cohorts mirroring the real design: donor
counts as observed (166, 625, 421, 188, 981, 96, 33, 99), 59 features,
ages 18–97. Cell-proportion features are compositional — generate on a
log-ratio scale and map back so features sum to 1 per donor. Four scenarios:

| Scenario | True structure | Question |
|---|---|---|
| S1a homogeneous | 10 true age effects, same sign and size in all cohorts | Do all selectors recover them? Does R(K) track truth? |
| S1b heterogeneous | 10 true effects with random-effects heterogeneity (tau² > 0), 2 sign-inconsistent | Does the min-aggregation veto punish M4 vs M0/M1 as on the real atlas? |
| S1c correlated null | no true effects, realistic block correlation | Selection-count saturation (H5) and FDP+ under knockoffs |
| S1d sparse strong | 3 strong consistent effects only | Ceiling check: can any selector separate rank signal from noise? |

**Execution.** Reuse `selectors.R`, `loco.R`, `metrics.R` unchanged; only a
data generator is new (`s1_simulate.R`). 20 replicates per scenario, 50
bootstraps (not 200) per STABL fit, same seed discipline (base seed
20261013 + replicate index). Estimated runtime: one M0–M5 pass over 8 folds
took ~3–4 h at 200 bootstraps; at 50 bootstraps one replicate is ~1 h, so
80 replicates ≈ 20–25 h — run as chunked background jobs, 4–5 replicates
per chunk, checkpointing per-replicate CSVs.

**Readout.** Per scenario: R@10 per selector vs truth, sign agreement,
selection count, FDP+. Success = endpoint ranks selectors sensibly under
known truth and the null scenario reproduces the saturation seen on real
data.

## 2. S2 — Real-atlas sensitivity probes

| Probe | Change | Cost | Question |
|---|---|---|---|
| S2a no-SPAC | drop SPAC; 7-cohort LOCO, full M0–M5 | ~2–3 h | SPAC (n = 33 < 59 features) needs a PD-projected knockoff covariance and has veto power under min-aggregation. How much of the saturation (V8) and the M3/M4 gap (V1) does it drive? |
| S2b aggregation quantile | parameterize min-aggregation: 10% and 25% cross-cohort quantiles between min (M4) and majority (M5) | small `selectors.R` patch + 2 runs (~3–4 h each) | Map the aggregation-strength curve; M5 > M4 already suggests strict min is too conservative (V5) |
| S2c knockoff variant (optional) | knockoff_mvr instead of knockoff_equi for M4 only | ~2 h | Robustness of the primary method to the knockoff construction |

K sensitivity (K = 5, 10, 15, 20) is already covered by the committed
`results/run/heldout_replication.csv`; do not redo it.

## 3. Order and effort

1. S2a first (cheapest, directly probes the biggest caveat: SPAC/PD veto).
2. S2b second (answers the open M4-vs-M5 aggregation question).
3. S1 if time allows (largest cost, calibrates the endpoint itself).
4. S2c optional.

Total: ~1–2 days of background compute on one 16-core worker, no new
infrastructure.

## 4. Falsification hooks — what would change the merged claims

- S2a: if removing SPAC removes the null saturation (selection count well
  below 59/59 in most permutations), V8's framing must be revised from
  "selection count carries no signal" to "saturation is a small-cohort
  artifact".
- S2b: if a 10–25% quantile beats both M4 and M5 on R@10, the locked
  min-aggregation choice is suboptimal; V5's "tentative" status upgrades
  and the vignette's aggregation discussion needs amending.
- S1b: if M4 loses to M0/M1 even under homogeneous truth (S1a), the
  deficit is structural to min-aggregation, not atlas-specific; if M4 wins
  under S1a but loses under S1b, heterogeneity is the driver — either way
  the H2 headline interpretation sharpens.
- S1c: if the null scenario does NOT saturate synthetically, the real-data
  saturation is driven by the atlas's correlation structure, not the
  selector alone — narrow V8 accordingly.

Any claim change goes through the claims ledger with a new amendment to
`benchmarks/ws8_showcase_predesign.md`, same as Amendment 1.
