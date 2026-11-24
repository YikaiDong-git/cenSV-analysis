#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
P1=${CENSV_ROOT}
WU=$P1/analyses/cdr_qtl_full_panel
M=$WU/output/mqtl

cd "$M"
if [ -d per_hap_bam ] && [ ! -d per_hap_bam_mirrorbug ]; then
  mv per_hap_bam per_hap_bam_mirrorbug
  echo "preserved per_hap_bam -> per_hap_bam_mirrorbug ($(ls per_hap_bam_mirrorbug | wc -l) files)"
fi
mkdir -p per_hap_bam
for d in matrices scans; do
  if [ -d "$d/bam" ] && [ ! -d "$d/bam_mirrorbug" ]; then
    cp -r "$d/bam" "$d/bam_mirrorbug"
    echo "preserved $d/bam -> $d/bam_mirrorbug ($(ls $d/bam_mirrorbug | wc -l) files)"
  fi
done

mkdir -p "$WU/logs/mqtl_liftbam_fix"

N=$(wc -l < taglist.txt)
JID=$(sbatch --parsable --nodes=1 --mem=12G -c 1 -t 1:00:00 --array=1-"$N"%60 \
  --output="$WU/logs/mqtl_liftbam_fix/%A_%a.log" \
  --wrap="cd $P1 && TAG=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $M/taglist.txt) && ~/.pixi/bin/pixi run python -u $WU/code/mqtl/lift/cis_mqtl_lift_bam.py \$TAG")
echo "JID=$JID  tasks=$N  throttle=60"
