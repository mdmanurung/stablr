#!/usr/bin/env bash
# WS-8 v2 S1 — driver: run the standard benchmark driver (run.R) over every
# simulated replicate at WS8_N_BOOT=50 (pilot bootstrap count), 8 concurrent.
# Each replicate writes results/s1/<scenario>_rep<r>/ and a log in
# /workspace/ws8_out/s1_run_logs/.  Re-runnable: existing outputs are skipped.
set -u
cd "$(dirname "$0")/../.."   # repo root

IN_ROOT="${WS8_S1_OUT:-/workspace/ws8_out/s1_inputs}"
RES_ROOT="benchmarks/ws8_v2/results/s1"
LOG_ROOT="/workspace/ws8_out/s1_run_logs"
mkdir -p "$RES_ROOT" "$LOG_ROOT"

jobs=()
for d in "$IN_ROOT"/*_rep*/; do
  tag=$(basename "$d")
  out="$RES_ROOT/$tag"
  if [ -f "$out/heldout_replication.csv" ]; then
    echo "skip $tag (done)"; continue
  fi
  jobs+=("$tag")
done
echo "to run: ${#jobs[@]} replicates"

printf '%s\n' "${jobs[@]}" | xargs -P 4 -I TAG bash -c '
  tag=TAG
  R_LIBS=/workspace/Rlib-shared \
  WS8_V2_IN="'"$IN_ROOT"'/$tag" \
  WS8_V2_OUT="'"$RES_ROOT"'/$tag" \
  WS8_N_BOOT=50 \
  Rscript benchmarks/ws8_v2/run.R > "'"$LOG_ROOT"'/$tag.log" 2>&1
  echo "done $tag (exit $?)"
'
echo "S1 driver complete."
