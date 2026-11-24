#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
JOBS=${JOBS:-100}
if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || (( JOBS < 1 || JOBS > 100 )); then
  printf 'JOBS must be an integer from 1 through 100; received %s\n' "$JOBS" >&2
  exit 2
fi
CONFIG=${CONFIG:-config/meqtl_comprehensive_conditioning.yaml}
SNAKEFILE=workflow/meqtl_comprehensive_conditioning/Snakefile
PROFILE=workflow/meqtl_comprehensive_conditioning/slurm_profile
OUTPUT_ROOT=$(uv run python -c 'import sys,yaml; print(yaml.safe_load(open(sys.argv[1]))["workflow_output"])' "$CONFIG")
CONTROLLER_DIR="$OUTPUT_ROOT/controller"
SLURM_LOG_DIR="$OUTPUT_ROOT/slurm_logs"
mkdir -p "$CONTROLLER_DIR" "$SLURM_LOG_DIR"
STAMP=$(date +%Y%m%d_%H%M%S)
LOG="$CONTROLLER_DIR/controller.$STAMP.log"
PID_FILE="$CONTROLLER_DIR/controller.$STAMP.pid"
export MEQTL_SLURM_LOG_DIR="$SLURM_LOG_DIR"
command=(uv run snakemake -s "$SNAKEFILE" --configfile "$CONFIG" --profile "$PROFILE" --jobs "$JOBS" "$@")
setsid "${command[@]}" >"$LOG" 2>&1 < /dev/null &
PID=$!
printf '%s\n' "$PID" > "$PID_FILE"
printf 'controller_pid=%s\ncontroller_log=%s\npid_file=%s\n' "$PID" "$LOG" "$PID_FILE"
