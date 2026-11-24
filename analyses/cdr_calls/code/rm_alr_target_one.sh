#!/usr/bin/env bash
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
SAMPLE=$1; HAP=$2; TAG="${SAMPLE}_${HAP}"
case "$HAP" in hap1) H=1;; hap2) H=2;; *) H="${HAP#hap}";; esac
ROOT=${CENSV_WORKSPACE}
SH=$ROOT/_shared/datasets/centromere_globe_logsdon2025/v_20260611
B=$ROOT/tools/CDR-Finder
D=$ROOT/Project1_Centromere/analyses/cdr_calls
SEQTK=$B/env/.pixi/envs/default/bin/seqtk
SAM=$B/env/.pixi/envs/default/bin/samtools
BT=$B/env/.pixi/envs/default/bin/bedtools
RM=$B/env_rm/.pixi/envs/default/bin/RepeatMasker
THREADS=${SLURM_CPUS_PER_TASK:-12}
WD=${CENSV_SCRATCH}/cdr_rmalr/$TAG
OUTBED=$D/inputs/rm_alr/${TAG}.bed; mkdir -p "$(dirname "$OUTBED")"
rm -rf "$WD"; mkdir -p "$WD/arrays"

unzip -o -j "$SH/cen_seq_HPRC.zip" "cen_array/${SAMPLE}_*${SAMPLE}#${H}#*" -d "$WD/arrays" >/dev/null 2>&1 || true
NFA=$(ls "$WD/arrays"/*.fa 2>/dev/null | wc -l)
[ "$NFA" -gt 0 ] || { echo "[FAIL] $TAG no cen_array fastas in globe"; exit 2; }
cat "$WD/arrays"/*.fa > "$WD/cat.fasta"
$SAM faidx "$WD/cat.fasta"
$SEQTK rename "$WD/cat.fasta" "seq" > "$WD/renamed.fasta"
$SAM faidx "$WD/renamed.fasta"
paste -d$'\t' <(cut -f1 "$WD/renamed.fasta.fai") <(cut -f1 "$WD/cat.fasta.fai") > "$WD/key.tsv"
export PATH=$B/env_rm/.pixi/envs/default/bin:$PATH
RepeatMasker -engine rmblast -species human -dir "$WD/rm" -qq -pa "$THREADS" "$WD/renamed.fasta" > "$WD/rm.log" 2>&1 || true
RMOUT="$WD/rm/renamed.fasta.out"
[ -s "$RMOUT" ] || { echo "[FAIL] $TAG RM produced no .out"; tail -20 "$WD/rm.log" 2>/dev/null; exit 3; }
awk -v OFS='\t' 'FNR==NR{kv[$1]=$2; next}
  FNR>3 && $10 ~ "ALR" {
    orig=kv[$5]; ci=index(orig,":"); name=substr(orig,1,ci-1); coords=substr(orig,ci+1);
    split(coords,c,"-"); off=c[1]+0;
    n=split(name,p,"_"); contig=p[3]; for(i=4;i<=n;i++) contig=contig"_"p[i];
    print contig, off+$6, off+$7
  }' "$WD/key.tsv" "$RMOUT" | sort -k1,1 -k2,2n \
  | $BT merge -d 1000 -i - | awk '($3-$2) > 100000' > "$OUTBED"
NREG=$(wc -l < "$OUTBED"); MB=$(awk '{s+=$3-$2}END{printf "%.1f",s/1e6}' "$OUTBED")
echo "[OK] $TAG RM-ALR target: $NREG regions, ${MB}Mb (from $NFA arrays)"
rm -rf "$WD"
