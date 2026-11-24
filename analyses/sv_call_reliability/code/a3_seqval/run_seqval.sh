#!/bin/bash
#SBATCH -J a3_seqval
#SBATCH -c 1
#SBATCH --mem 12G
#SBATCH -t 6:00:00
#SBATCH --array=1-6
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
D=${CENSV_ROOT}/analyses/sv_call_reliability
PY=${CENSV_ROOT}/.pixi/envs/default/bin/python
O=$D/output/a3_seqval
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $O/targets.tsv)
[ -z "$line" ] && { echo "no target for task $SLURM_ARRAY_TASK_ID"; exit 0; }
IFS=$'\t' read -r tid chrom pos svtype len src freq <<< "$line"
$PY $D/code/a3_seqval/seqval.py "$tid" "$chrom" "$pos" "$svtype" "$len" $O 2> >(grep -v 'index file is older' >&2)
echo A3_DONE $tid
