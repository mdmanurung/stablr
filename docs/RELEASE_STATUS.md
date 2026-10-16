# Release governance and review status

The bundled OOL vignette is a workflow tutorial. It is not application-note
evidence. Release claims require a separate run from an exact, hashed candidate
source tarball and the predeclared full OOL inputs and gates; mutable-checkout
execution does not qualify.

```bash
# From stablr-experiments/ after building the candidate tarball. The output
# directory must not exist.
RUNTIME=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/lib/R/bin/Rscript
"$RUNTIME" --vanilla analysis/generate_publication_parity_table.R \
  --candidate-tarball /absolute/path/to/stablr_0.3.0.tar.gz \
  --runtime "$RUNTIME" \
  --outdir /absolute/new/path/ool-candidate-evidence \
  --profile full
```

For the 0.1.1 candidate, the frozen full-candidate protocol produced a
negative result: support Jaccard 0.778 (required 0.80), score Spearman 0.785
(required 0.90), and FDP-threshold delta 0.08 (maximum 0.05). The lambda and
training-only preprocessing gates passed. The requested fixed knockoff
construction was infeasible for the 150-by-1,317 matrix (`n > p` is
required), so the typed, recorded fallback was `random_permutation` and the
actual-generator gate also failed. The sealed negative artifact is retained.

Release status for the current candidate: the SCI-04 null-selected-fraction
gate was accepted on 2026-10-16 by the project owner acting in the reviewer
role (accepted record `SCI-04-null-selected-fraction-v2-accepted`; reviewer
response in `inst/analysis/review/SCI-04-acceptance-record.md`). This is a
**provisional owner acceptance**: the review packet's acceptance boundary
requires an itemized review by a reviewer independent of the package, which
is still outstanding (`sci04-independent-confirmation-pending` in the release
contract). The contract pins the accepted Hoeffding criterion
(`alpha = 0.05`, `U <= 0.10`, `n = 150`) and the gate evaluates instead of
failing closed. The entry point for the itemized independent confirmation is
`inst/analysis/review/SCI-04-independent-confirmation-packet.md`. Remaining
release blockers: the fresh full-candidate release run at the accepted
replicate count, the itemized independent confirmation, and assembly of the
new Release Evidence Packet.
