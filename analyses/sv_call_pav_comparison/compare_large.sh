#!/bin/bash
#SBATCH --job-name=pav_vs_svim_large
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/compare-large-%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
ROOT=${CENSV_ROOT}
cd "$ROOT"
PAV=analyses/sv_call_pav_comparison/HG002.vcf.gz
SVIM_H1=$(ls variant_calling/svim_asm/HG002_hap1_mapq20_*.vcf.gz | head -1)
SVIM_H2=$(ls variant_calling/svim_asm/HG002_hap2_mapq20_*.vcf.gz | head -1)
CEN_BED=resources/chm13v2.0_censat_v2.1.active.bed
PIX="$HOME/.pixi/bin/pixi run --manifest-path $ROOT/pixi.toml"

echo "=== Q1: GENOME-WIDE large SVs (no centromere restriction) ==="
echo "--- PAV genome-wide ---"
$PIX bcftools view -H "$PAV" 2>/dev/null \
  | awk -F"\t" "{ for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){sv=a[1]+0;if(sv<0)sv=-sv; for(c=0;c<=8;c++){ cut[c]; if(sv>=int(thr[c]))n[c]++ }; break } } } BEGIN{thr[0]=50;thr[1]=100;thr[2]=500;thr[3]=1000;thr[4]=5000;thr[5]=10000;thr[6]=50000;thr[7]=100000;thr[8]=500000} END{for(c=0;c<=8;c++) printf \"  PAV >= %dbp: %d\n\", thr[c], n[c]+0}"

echo "--- svim-asm genome-wide (h1 + h2 union by chr,pos,end,svtype) ---"
TMP=$(mktemp -d)
for h in h1 h2; do
  if [ "$h" = "h1" ]; then F=$SVIM_H1; else F=$SVIM_H2; fi
  $PIX bcftools view -H "$F" 2>/dev/null \
    | awk -F"\t" "BEGIN{OFS=\"\t\"} {svlen=0;svtype=\"\"; for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){svlen=a[1]+0;if(svlen<0)svlen=-svlen} if(match(\$i,/SVTYPE=([^;]+)/,b))svtype=b[1] } if(svtype!=\"\"){end=\$2+(svtype==\"DEL\"?svlen:1); print \$1,\$2,end,svtype,svlen}}"
done | sort -u -k1,1 -k2,2n -k3,3n -k4,4 > "$TMP/svim_gw.tsv"
awk "BEGIN{thr[0]=50;thr[1]=100;thr[2]=500;thr[3]=1000;thr[4]=5000;thr[5]=10000;thr[6]=50000;thr[7]=100000;thr[8]=500000} {sv=\$5+0; for(c=0;c<=8;c++) if(sv>=thr[c]) n[c]++} END{for(c=0;c<=8;c++) printf \"  svim >= %dbp: %d\n\", thr[c], n[c]+0}" "$TMP/svim_gw.tsv"

echo
echo "=== Q2: TOP 20 LARGEST SVs from each ==="
echo "--- PAV top-20 ---"
$PIX bcftools view -H "$PAV" 2>/dev/null \
  | awk -F"\t" "BEGIN{OFS=\"\t\"} {svlen=0;svtype=\"\"; for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){svlen=a[1]+0;if(svlen<0)svlen=-svlen} if(match(\$i,/SVTYPE=([^;]+)/,b))svtype=b[1] } if(svtype!=\"\" && svtype!=\"SNV\") print svlen, svtype, \$1, \$2}" \
  | sort -rn | head -20

echo "--- svim-asm top-20 (h1+h2) ---"
sort -k5,5rn "$TMP/svim_gw.tsv" | head -20

echo
echo "=== Q3: PAV SVs > 500kb (svim-asm cap) ==="
echo "--- PAV variants > 500,000 bp (where svim-asm caps) ---"
$PIX bcftools view -H "$PAV" 2>/dev/null \
  | awk -F"\t" "BEGIN{OFS=\"\t\"} {svlen=0;svtype=\"\"; for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){svlen=a[1]+0;if(svlen<0)svlen=-svlen} if(match(\$i,/SVTYPE=([^;]+)/,b))svtype=b[1] } if(svtype!=\"\" && svtype!=\"SNV\" && svlen>500000) print svlen, svtype, \$1, \$2}" \
  | sort -rn | head -30

echo
echo "=== Q4: Inversions (PAV emphasis) ==="
echo "--- PAV inversions genome-wide ---"
$PIX bcftools view -H "$PAV" 2>/dev/null \
  | grep "SVTYPE=INV" | awk -F"\t" "{for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){sv=a[1]+0;if(sv<0)sv=-sv; print sv, \$1, \$2; break}}}" | sort -rn | head -15
echo "--- svim-asm inversions h1+h2 ---"
for h in h1 h2; do
  if [ "$h" = "h1" ]; then F=$SVIM_H1; else F=$SVIM_H2; fi
  $PIX bcftools view -H "$F" 2>/dev/null | grep "SVTYPE=INV" | awk -F"\t" "{for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){sv=a[1]+0;if(sv<0)sv=-sv; print sv, \$1, \$2; break}}}"
done | sort -rn | head -15

echo
echo "=== DONE $(date) ==="
rm -rf "$TMP"
