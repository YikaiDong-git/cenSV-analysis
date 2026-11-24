#!/bin/bash
#SBATCH --job-name=pav_vs_svim_deep
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/compare-deep-%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
ROOT=${CENSV_ROOT}
cd "$ROOT"
PAV=analyses/sv_call_pav_comparison/HG002.vcf.gz
SVIM_H1=$(ls variant_calling/svim_asm/HG002_hap1_mapq20_*.vcf.gz | head -1)
SVIM_H2=$(ls variant_calling/svim_asm/HG002_hap2_mapq20_*.vcf.gz | head -1)
CEN_BED=resources/chm13v2.0_censat_v2.1.active.bed
WORK=analyses/sv_call_pav_comparison/cmp_deep
mkdir -p "$WORK"
PIX="$HOME/.pixi/bin/pixi run --manifest-path $ROOT/pixi.toml"

$PIX bcftools view "$PAV" -R "$CEN_BED" 2>/dev/null | $PIX bcftools view -H 2>/dev/null \
  | awk -F"\t" "BEGIN{OFS=\"\t\"} {svlen=0;svtype=\"\"; for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=([0-9-]+)/,a)){svlen=a[1]+0;if(svlen<0)svlen=-svlen} if(match(\$i,/SVTYPE=([^;]+)/,b))svtype=b[1] } if(svtype!=\"\" && svtype!=\"SNV\" && svlen>=50){end=\$2+(svtype==\"DEL\"?svlen:1); print \$1,\$2,end,svtype,svlen}}" \
  > "$WORK/pav_cen_50.bed"
echo "PAV ≥50bp cen: $(wc -l < $WORK/pav_cen_50.bed)"

for h in h1 h2; do
  if [ "$h" = "h1" ]; then F=$SVIM_H1; else F=$SVIM_H2; fi
  $PIX bcftools view "$F" -R "$CEN_BED" 2>/dev/null | $PIX bcftools view -H 2>/dev/null \
    | awk -F"\t" -v hap="$h" "BEGIN{OFS=\"\t\"} {svlen=0;svtype=\"\"; for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=(-?[0-9]+)/,a)){svlen=a[1]+0;if(svlen<0)svlen=-svlen} if(match(\$i,/SVTYPE=([^;]+)/,b))svtype=b[1] } if(svtype!=\"\" && svlen>=50){end=\$2+(svtype==\"DEL\"?svlen:1); print \$1,\$2,end,svtype,svlen,hap}}" \
    > "$WORK/svim_${h}_cen_50.bed"
done
cat "$WORK/svim_h1_cen_50.bed" "$WORK/svim_h2_cen_50.bed" | sort -k1,1 -k2,2n -u > "$WORK/svim_union_cen_50.bed"
echo "svim h1 ≥50bp cen: $(wc -l < $WORK/svim_h1_cen_50.bed)"
echo "svim h2 ≥50bp cen: $(wc -l < $WORK/svim_h2_cen_50.bed)"
echo "svim union ≥50bp cen: $(wc -l < $WORK/svim_union_cen_50.bed)"

$PIX bedtools window -a "$WORK/pav_cen_50.bed" -b "$WORK/svim_union_cen_50.bed" -w 100 \
  | awk -F"\t" "\$4 == \$9" | sort -u -k1,1 -k2,2n > "$WORK/concordant_full.bed"
$PIX bedtools window -a "$WORK/pav_cen_50.bed" -b "$WORK/svim_union_cen_50.bed" -w 100 \
  | awk -F"\t" "\$4 == \$9 {print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}" | sort -u > "$WORK/pav_concordant_keys.bed"
$PIX bedtools intersect -v -a "$WORK/pav_cen_50.bed" -b "$WORK/pav_concordant_keys.bed" > "$WORK/pav_only.bed"
$PIX bedtools window -a "$WORK/svim_union_cen_50.bed" -b "$WORK/pav_cen_50.bed" -w 100 \
  | awk -F"\t" "\$4 == \$10 {print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}" | sort -u > "$WORK/svim_concordant_keys.bed"
$PIX bedtools intersect -v -a "$WORK/svim_union_cen_50.bed" -b "$WORK/svim_concordant_keys.bed" > "$WORK/svim_only.bed"

echo
echo "=== HEADLINE: HG002 centromeric SV recovery (≥50bp, same SVTYPE, ±100bp position) ==="
PAV_TOTAL=$(wc -l < $WORK/pav_cen_50.bed)
SVIM_TOTAL=$(wc -l < $WORK/svim_union_cen_50.bed)
CONCORDANT=$(wc -l < $WORK/pav_concordant_keys.bed)
PAV_ONLY=$(wc -l < $WORK/pav_only.bed)
SVIM_ONLY=$(wc -l < $WORK/svim_only.bed)
echo "  PAV total: $PAV_TOTAL | svim total: $SVIM_TOTAL"
echo "  Concordant (PAV ∩ svim, ±100bp, same type): $CONCORDANT"
echo "  PAV-only: $PAV_ONLY  ($((100*PAV_ONLY/PAV_TOTAL))% of PAV)"
echo "  svim-only: $SVIM_ONLY  ($((100*SVIM_ONLY/SVIM_TOTAL))% of svim)"
echo
echo "=== By SVTYPE ==="
for t in DEL INS INV; do
  pn=$(awk -v t="$t" "\$4==t" $WORK/pav_cen_50.bed | wc -l)
  sn=$(awk -v t="$t" "\$4==t" $WORK/svim_union_cen_50.bed | wc -l)
  printf "  %-4s PAV: %5d | svim: %5d\n" "$t" "$pn" "$sn"
done
echo
echo "=== By size bucket (PAV non-SNV / svim) ==="
for cutoff in 50 100 500 1000 5000 10000; do
  pn=$(awk -v c="$cutoff" "\$5>=c" $WORK/pav_cen_50.bed | wc -l)
  sn=$(awk -v c="$cutoff" "\$5>=c" $WORK/svim_union_cen_50.bed | wc -l)
  printf "  ≥%-5dbp  PAV: %5d | svim: %5d\n" "$cutoff" "$pn" "$sn"
done
echo
echo "=== DONE $(date) ==="
