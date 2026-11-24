#!/bin/bash
#SBATCH -J a2_len
#SBATCH -c 1
#SBATCH --mem 4G
#SBATCH -t 4:00:00
#SBATCH --array=0-15
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
D=${CENSV_ROOT}/analyses/sv_call_reliability
S=${CENSV_COLLAB_ROOT}/exploration/simulate_svim_asm_accuracy/scripts/measure_array_lengths.py
PY=${CENSV_ROOT}/.pixi/envs/default/bin/python
mkdir -p $D/output/a2_netlength/lengths
CH=$(for i in $(seq 1 22) X; do printf "chr%s " $i; done)
$PY $S --chroms $CH --shard $SLURM_ARRAY_TASK_ID --n-shards 16 \
  --output $D/output/a2_netlength/lengths/shard_${SLURM_ARRAY_TASK_ID}.tsv 2> >(grep -v 'index file is older' >&2)
echo SHARD_DONE $SLURM_ARRAY_TASK_ID
