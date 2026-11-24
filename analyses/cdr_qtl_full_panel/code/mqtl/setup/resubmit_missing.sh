#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
P1=${CENSV_ROOT}
WU=$P1/analyses/cdr_qtl_full_panel
M=$WU/output/mqtl

cd "$M"
sort taglist.txt > /tmp/all_tags.txt
ls per_hap/ 2>/dev/null | sed 's/\.tsv\.gz$//' | sort > /tmp/done_tags.txt
comm -23 /tmp/all_tags.txt /tmp/done_tags.txt > missing_tags.txt
N=$(wc -l < missing_tags.txt)
echo "done: $(wc -l < /tmp/done_tags.txt)   missing: $N"
if [ "$N" -eq 0 ]; then echo "nothing to resubmit"; exit 0; fi
mkdir -p "$WU/logs/mqtl_lift"
JID=$(sbatch --parsable --nodes=1 --mem=12G -c 1 -t 0:30:00 --array=1-"$N"%50 \
  --output="$WU/logs/mqtl_lift/%A_%a.log" \
  --wrap="cd $P1 && TAG=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $M/missing_tags.txt) && ~/.pixi/bin/pixi run python -u $WU/code/mqtl/cis_mqtl_lift.py \$TAG")
echo "resubmitted JID=$JID for $N tags"
