# Pushing the ws8-v3-candidates branch

The module in this archive is complete and self-contained. From a local clone
of mdmanurung/stablr:

```bash
git checkout -b ws8-v3-candidates
# copy the archive contents over the repo root:
tar -xzf ws8_v3_module.tar.gz -C /path/to/stablr   # creates benchmarks/ws8_v3/
git add benchmarks/ws8_v3
git commit -m "benchmarks: WS-8 v3 candidate mechanisms for cross-cohort reproducible biomarker discovery

M6 het-filtered RE-meta, M7 MINT-STABL (+M7c), M8 contrast-STABL,
M9/M9s/M9b ComBat+STABL baseline, benchmarked against the v2 matrix
(M0-M5 rerun like-for-like on pinned a8f3a58) on the Terekhova 2026
immune-aging atlas. Includes LOCO driver, null probes, full-atlas fits,
assembly/statistics, figures, unit tests, and assembled results."
git push -u origin ws8-v3-candidates
```

Notes:
- The module reuses the v2 engine (`benchmarks/ws8_v2/{metrics,loco,selectors}.R`)
  and data (`benchmarks/ws8_v2/results/m1/donor_tables.rds`) already in the repo.
- All scripts take repo-relative default paths; run from the repository root.
  See benchmarks/ws8_v3/README.md for details and the package-version policy.
- Results were produced with stablr pinned at a8f3a58d74f17658fd56467e0355aade1ad667ab.
