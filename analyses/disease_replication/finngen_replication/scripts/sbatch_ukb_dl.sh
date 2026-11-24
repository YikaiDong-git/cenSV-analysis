#!/bin/bash
#SBATCH --job-name=ukb_ear_v2
#SBATCH --time=3:00:00
#SBATCH --mem=16G
#SBATCH -c 4
: "${CENSV_LAB:?set CENSV_LAB (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
REPL=${CENSV_ROOT}/analyses/disease_replication/finngen_replication
TOOLDIR=${CENSV_WORKSPACE}/_shared/tools/process_neale_lab_GWAS/v_master_83578d1_20260514
MAN=$TOOLDIR/neale_lab_manifest_20180731.tsv
DBSNP=${CENSV_LAB}/projects/resources/dbSNP_to_parquet
VARFILE=${CENSV_LAB}/projects/PRS/LD_REFERENCE_PANEL/output/tsv/variant_list_EUR.tsv
LIST=$REPL/manifests/ukb_phenos.tsv
OUTBASE=$REPL/ukb_neale
ROW=$SLURM_ARRAY_TASK_ID
CODE=$(awk -F'\t' -v r="$ROW" 'NR==r{print $1}' "$LIST")
[ -z "$CODE" ] && { echo "no code at row $ROW"; exit 2; }
echo "[task $ROW] UKB $CODE"
MROW=$(awk -F'\t' -v c="$CODE" '$1==c && $4=="both_sexes"{print; exit}' "$MAN")
[ -z "$MROW" ] && { echo "no both_sexes manifest row for $CODE"; exit 3; }
FILE=$(printf '%s' "$MROW" | awk -F'\t' '{print $5}')
URL=$(printf '%s' "$MROW" | awk -F'\t' '{print $6}' | sed 's/^wget //; s/ -O .*$//')
echo "  FILE=$FILE"; echo "  URL=$URL"
mkdir -p "$OUTBASE/raw/$CODE" "$OUTBASE/processed/$CODE"
BGZ="$OUTBASE/raw/$CODE/$FILE"
if [ ! -s "$BGZ" ]; then
  echo "  [dl] wget $CODE"
  wget -q --tries=3 "$URL" -O "$BGZ.part"
  mv -f "$BGZ.part" "$BGZ"
fi
echo "  size: $(stat -c %s "$BGZ")"
OUTFILE="$OUTBASE/processed/$CODE/${CODE}.hg38.tsv"
cd "$TOOLDIR"
~/.pixi/bin/pixi run python convert.py -s "$BGZ" -o "$OUTFILE" -d "$DBSNP" -ib hg19 -vf "$VARFILE" -cf hg19ToHg38.over.chain -lt liftOver
[ -s "$OUTFILE" ] || { echo "ERROR: empty $OUTFILE"; exit 4; }
touch "$OUTBASE/processed/$CODE/.process_ok"
echo "  DONE $CODE -> $OUTFILE ($(wc -l < "$OUTFILE") rows)"
