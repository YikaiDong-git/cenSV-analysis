#!/bin/bash
#SBATCH --job-name=pav_vs_svim_hor
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/in-hor-%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -eu
cd ${CENSV_ROOT}
source ~/.bashrc
PIXI=${CENSV_ROOT}/.pixi/envs/default/bin

PAV=analyses/sv_call_pav_comparison/HG002.vcf.gz
SVIM_H1=$(ls variant_calling/svim_asm/HG002_hap1_*.vcf.gz | head -1)
SVIM_H2=$(ls variant_calling/svim_asm/HG002_hap2_*.vcf.gz | head -1)
HOR_BED=resources/chm13v2.0_censat_v2.1.active.bed

echo "=== Inputs ==="
echo "PAV:     $PAV"
echo "svim h1: $SVIM_H1"
echo "svim h2: $SVIM_H2"
echo "HOR bed: $HOR_BED"
echo

$PIXI/bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n" $PAV \
  | awk -v OFS="\t" "{ sz = (\$4 < 0 ? -\$4 : \$4); print \$1, \$2-1, \$2+sz, \$3, sz }" \
  > /tmp/pav_all.bed

$PIXI/bedtools intersect -a /tmp/pav_all.bed -b $HOR_BED -u > /tmp/pav_in_hor.bed
echo "PAV total SVs (genome-wide): $(wc -l < /tmp/pav_all.bed)"
echo "PAV SVs inside active HOR:    $(wc -l < /tmp/pav_in_hor.bed)"
echo

for V in $SVIM_H1 $SVIM_H2; do
  $PIXI/bcftools query -f "%CHROM\t%POS\t%END\t%INFO/SVTYPE\t%INFO/SVLEN\n" $V
done | awk -v OFS="\t" "{ sz = (\$5 < 0 ? -\$5 : \$5); print \$1, \$2-1, \$3, \$4, sz }" \
  | sort -u > /tmp/svim_all.bed
$PIXI/bedtools intersect -a /tmp/svim_all.bed -b $HOR_BED -u > /tmp/svim_in_hor.bed
echo "svim total SVs (h1+h2 union, already cen-restricted at call time): $(wc -l < /tmp/svim_all.bed)"
echo "svim SVs verified in active HOR: $(wc -l < /tmp/svim_in_hor.bed)"
echo

echo "=== APPLES-TO-APPLES: SVs inside active HOR by size bucket ==="
printf "%-12s %10s %10s %10s\n" "Size cutoff" "PAV" "svim-asm" "Ratio"
for cutoff in 50 100 500 1000 5000 10000 50000 100000 500000; do
  pav=$(awk -v c=$cutoff "\$5>=c" /tmp/pav_in_hor.bed | wc -l)
  svm=$(awk -v c=$cutoff "\$5>=c" /tmp/svim_in_hor.bed | wc -l)
  if [ "$svm" -gt 0 ]; then
    ratio=$(awk -v p=$pav -v s=$svm "BEGIN{printf \"%.2fx\", p/s}")
  else
    ratio="(svim=0)"
  fi
  printf "%-12s %10d %10d %10s\n" ">=${cutoff}bp" "$pav" "$svm" "$ratio"
done
echo

echo "=== Any PAV SVs >=500kb inside active HOR? ==="
awk "\$5>=500000" /tmp/pav_in_hor.bed
echo "(empty = none)"
echo

echo "=== Any PAV SVs >=100kb inside active HOR? ==="
awk -v OFS="\t" "\$5>=100000 { print }" /tmp/pav_in_hor.bed | sort -k5,5 -nr
echo "(empty = none)"
echo

echo "=== Any PAV SVs >=50kb inside active HOR? ==="
awk -v OFS="\t" "\$5>=50000 { print }" /tmp/pav_in_hor.bed | sort -k5,5 -nr
echo "(empty = none)"
