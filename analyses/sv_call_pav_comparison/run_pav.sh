#!/bin/bash
#SBATCH --job-name=pav_smoke
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/slurm-%j.out
#SBATCH --error=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/slurm-%j.err
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
ROOT=${CENSV_ROOT}
cd "$ROOT/analyses/sv_call_pav_comparison"
export OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$ROOT/software/pav:$ROOT/software/pav/dep/svpop${PYTHONPATH:+:$PYTHONPATH}"
echo "=== ENV ==="
echo "host: $(hostname)"; echo "pwd: $(pwd)"; echo "date: $(date)"
echo "=== run PAV (Logsdon-style: defaults + reference + aligner + assembly_table) ==="
"$HOME/.pixi/bin/pixi" run --manifest-path "$ROOT/pixi.toml" \
  snakemake -s "$ROOT/software/pav/Snakefile" -j 32 --rerun-incomplete \
  2>&1
echo "=== DONE ==="; date
