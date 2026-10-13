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

Release status for the current candidate is governed by the SCI-04
null-selected-fraction gate, which is fail-closed at `proposed_unaccepted`
until an external reviewer accepts the predeclared criterion
(`inst/analysis/review/SCI-04-independent-review-packet.md`); the fresh
full-candidate release run then executes at the accepted sample size.
