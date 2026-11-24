#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${HOME:?set HOME (see README)}"
set -euo pipefail
P=${CENSV_ROOT}
DISC=$P/analyses/ld_maps
OUT=$DISC/output/_censv_frames
PIXI="${HOME}/.pixi/bin/pixi run --manifest-path $P/pixi.toml python"
mkdir -p "$OUT" "$DISC/logs"

CH=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 \
    chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
printf "%s\n" "${CH[@]}" > "$OUT/chroms.txt"

sbatch --job-name=censvld --output="$DISC/logs/censvld_%A_%a.log" \
  --array=1-${#CH[@]}%8 --mem=8G -c 2 --time=00:30:00 \
  --wrap="c=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $OUT/chroms.txt); \
          $PIXI $DISC/code/extract_censv_ld.py --chrom \$c --out $OUT"
