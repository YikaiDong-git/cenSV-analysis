#!/bin/bash
#SBATCH --job-name=svim_bnd_hor
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/svim-bnd-hor-%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -eu
cd ${CENSV_ROOT}
source ~/.bashrc
PIXI=${CENSV_ROOT}/.pixi/envs/default/bin
SVIM_H1=$(ls variant_calling/svim_asm/HG002_hap1_*.vcf.gz | head -1)
SVIM_H2=$(ls variant_calling/svim_asm/HG002_hap2_*.vcf.gz | head -1)
HOR_BED=resources/chm13v2.0_censat_v2.1.active.bed
PAV=analyses/sv_call_pav_comparison/HG002.vcf.gz

echo "=== svim-asm SVTYPE breakdown (h1 + h2 union, before HOR filter) ==="
for V in $SVIM_H1 $SVIM_H2; do
  echo "--- $V ---"
  $PIXI/bcftools query -f "%INFO/SVTYPE\n" $V | sort | uniq -c
done

echo
echo "=== Extract all BND records from svim-asm h1+h2, parse partner from ALT ==="
for V in $SVIM_H1 $SVIM_H2; do
  $PIXI/bcftools query -f "%CHROM\t%POS\t%ALT\t%INFO/SVTYPE\n" $V \
    | awk -v OFS="\t" "
      \$4 == \"BND\" {
        alt = \$3
        if (match(alt, /[][]([^][]+)[][]/, m)) {
          split(m[1], parts, \":\")
          pchr = parts[1]
          ppos = parts[2]
          if (pchr == \$1) {
            sz = (ppos > \$2 ? ppos - \$2 : \$2 - ppos)
          } else {
            sz = -1   # interchromosomal — no implied size
          }
          print \$1, \$2, pchr, ppos, sz
        }
      }
    "
done | sort -u > /tmp/svim_bnd_all.tsv
echo "Total BNDs (h1+h2 union, deduplicated): $(wc -l < /tmp/svim_bnd_all.tsv)"
echo "  intra-chrom BNDs (same chrom): $(awk \"\$5 >= 0\" /tmp/svim_bnd_all.tsv | wc -l)"
echo "  inter-chrom BNDs (translocations): $(awk \"\$5 == -1\" /tmp/svim_bnd_all.tsv | wc -l)"

awk -v OFS="\t" "{print \$1, \$2-1, \$2, \$3\":\"\$4, \$5}" /tmp/svim_bnd_all.tsv > /tmp/svim_bnd_all.bed
$PIXI/bedtools intersect -a /tmp/svim_bnd_all.bed -b $HOR_BED -u > /tmp/svim_bnd_in_hor.bed
echo
echo "=== svim-asm BNDs inside active HOR ==="
echo "Total BNDs in HOR: $(wc -l < /tmp/svim_bnd_in_hor.bed)"
echo "  intra-chrom: $(awk \"\$5 >= 0\" /tmp/svim_bnd_in_hor.bed | wc -l)"
echo "  inter-chrom: $(awk \"\$5 == -1\" /tmp/svim_bnd_in_hor.bed | wc -l)"

echo
echo "=== Intra-chrom BND implied-size buckets in HOR ==="
printf "%-12s %10s\n" "Cutoff" "BNDs"
for c in 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do
  n=$(awk -v c=$c "\$5 >= c" /tmp/svim_bnd_in_hor.bed | wc -l)
  printf "%-12s %10d\n" ">=${c}bp" "$n"
done

echo
echo "=== Top-20 largest intra-chrom BNDs in active HOR ==="
sort -k5,5 -nr /tmp/svim_bnd_in_hor.bed | awk "\$5 >= 0" | head -20 | column -t

echo
echo "=== Now build the full apples-to-apples table: PAV (default) vs svim-asm DEL+INV+INS in HOR vs svim-asm DEL+INV+INS+BND in HOR ==="
$PIXI/bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n" $PAV \
  | awk -v OFS="\t" "{ sz = (\$4 < 0 ? -\$4 : \$4); print \$1, \$2-1, \$2+sz, \$3, sz }" \
  > /tmp/pav_all.bed
$PIXI/bedtools intersect -a /tmp/pav_all.bed -b $HOR_BED -u > /tmp/pav_in_hor.bed

for V in $SVIM_H1 $SVIM_H2; do
  $PIXI/bcftools query -f "%CHROM\t%POS\t%END\t%INFO/SVTYPE\t%INFO/SVLEN\n" $V
done | awk -v OFS="\t" "\$4 != \"BND\" { sz = (\$5 < 0 ? -\$5 : \$5); print \$1, \$2-1, \$3, \$4, sz }" \
  | sort -u > /tmp/svim_nonbnd.bed
$PIXI/bedtools intersect -a /tmp/svim_nonbnd.bed -b $HOR_BED -u > /tmp/svim_nonbnd_in_hor.bed

printf "%-12s %10s %10s %10s %10s\n" "Cutoff" "PAV" "svim_noBND" "svim_BND" "svim_total"
for c in 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000; do
  pav=$(awk -v c=$c "\$5 >= c" /tmp/pav_in_hor.bed | wc -l)
  sv1=$(awk -v c=$c "\$5 >= c" /tmp/svim_nonbnd_in_hor.bed | wc -l)
  sv2=$(awk -v c=$c "\$5 >= c" /tmp/svim_bnd_in_hor.bed | wc -l)
  total=$((sv1 + sv2))
  printf "%-12s %10d %10d %10d %10d\n" ">=${c}bp" "$pav" "$sv1" "$sv2" "$total"
done
