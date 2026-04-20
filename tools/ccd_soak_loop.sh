#!/usr/bin/env bash
# Run the CCD soak in a loop, restarting on crash so a known BVH churn
# crash can't halt an overnight run. Each inner run uses a different RNG
# seed so we explore new projectile sequences on every restart. All
# violations accumulate into ccd_soak_violations.log (append-mode).
#
# Usage: tools/ccd_soak_loop.sh [duration_seconds]
#   duration_seconds: total wall-clock budget. Default: 8 hours.
#   Each inner run is capped at 30000 frames (~8 minutes of sim time).
#
# Output:
#   ccd_soak_violations.log - every VIOLATION line logged by the soak
#   soak_runs.log           - per-run exit code + seed + violation count
#
# To triage an overnight run the morning after:
#   grep VIOLATION ccd_soak_violations.log | sort -u | head
# Each line has the spawn_rng / kind / scale / vel / omega needed to
# reproduce via `nudge_tests.exe --ccd-repro <args>`.

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
EXE="$ROOT/build/Release/nudge_tests.exe"
if [ ! -x "$EXE" ]; then
  echo "nudge_tests.exe not found at $EXE; build Release first." >&2
  exit 1
fi

DURATION=${1:-28800} # 8 hours default
INNER_FRAMES=30000   # ~8 min per inner run
RUNS_LOG="$ROOT/soak_runs.log"
VIOL_LOG="$ROOT/ccd_soak_violations.log"

# Append session header, keep history across multiple overnight runs.
echo "=== soak session start $(date -Is) budget=${DURATION}s ===" >> "$RUNS_LOG"

end_ts=$(( $(date +%s) + DURATION ))
run=0
total_viol=0
total_crash=0
while [ "$(date +%s)" -lt "$end_ts" ]; do
  run=$((run + 1))
  seed=$(printf '%08x' $((RANDOM * RANDOM + $(date +%N 2>/dev/null || echo $$))))
  pre_lines=$(wc -l <"$VIOL_LOG" 2>/dev/null || echo 0)
  CCD_SOAK_SEED="0x$seed" "$EXE" --bench-ccd-soak "$INNER_FRAMES" >/dev/null 2>>"$RUNS_LOG"
  exit_code=$?
  post_lines=$(wc -l <"$VIOL_LOG" 2>/dev/null || echo 0)
  v=$((post_lines - pre_lines))
  total_viol=$((total_viol + v))
  if [ "$exit_code" -ne 0 ]; then total_crash=$((total_crash + 1)); fi
  echo "run=$run seed=0x$seed exit=$exit_code new_violations=$v" >> "$RUNS_LOG"
  echo "[$run] seed=0x$seed exit=$exit_code +${v} violations (running total: $total_viol)"
done

echo "=== soak session end $(date -Is) runs=$run total_violations=$total_viol crashes=$total_crash ===" >> "$RUNS_LOG"
echo "Soak done. $run runs, $total_viol violations, $total_crash crashes."
echo "Full log: $VIOL_LOG"
echo "Run log:  $RUNS_LOG"
