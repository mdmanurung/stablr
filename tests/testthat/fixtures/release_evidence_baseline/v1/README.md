# RE-01 bounded behavior baseline v1

These ten normalized artifacts characterize `stablr` commit
`c7bd2f8b8ade95e1050d15db6e0f6443f21ce6a1` under the pinned R4_51
runtime. They are bounded diagnostic fixtures, not promotable release evidence.

The ignored raw identity is
`.claude/reviews/stabl_lw/re01-baseline-c7bd2f8-v1/continuations/04-r4-toolchain/`.
Methodology files come from normalized `methodology-serial-a-v3`; late-fusion
files come from normalized `late-fusion-a`. Serial A/B and late-fusion A/B
were byte-identical after normalization. Serial/parallel scientific tables
were byte-identical; only the manifest worker count differed.

Methodology settings: bounded profile, two replicates, two bootstraps, two
lambdas, all four families, all four artificial generators, all six default
scenarios, seed 270627, target FDP 0.1, and one worker. Late-fusion settings:
one replicate per cell, two bootstraps, two iterations, and seed 220711.
The late-fusion run produced all four artifacts and then emitted the expected
condition `Late-fusion release gates failed; release must stop.`

Normalization is deliberately narrow: methodology `elapsed_sec` and
`mean_elapsed_sec`; specified methodology timestamp, package-mode, Git,
artifact-path, and artifact-hash values; and specified late-fusion
package-mode, Git, and artifact-hash values. No row, column, warning, seed,
numeric text, condition, gate, or scientific result was reordered or removed.

Known SCI-01 through SCI-10 defects remain represented by this pre-fix
baseline. These files must never be assembled into a Release Evidence Packet.
