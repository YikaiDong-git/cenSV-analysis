#!/bin/bash
#SBATCH --job-name=pav_large_locate
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/large-locate-%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -eu
cd ${CENSV_ROOT}
source ~/.bashrc
PIXI=${CENSV_ROOT}/.pixi/envs/default/bin
PAV=analyses/sv_call_pav_comparison/HG002.vcf.gz
ACTIVE=resources/chm13v2.0_censat_v2.1.active.bed
ALL_CEN=resources/chm13v2.0_censat_v2.1.bed

$PIXI/bcftools view -h $PAV | tail -1 | head -1
$PIXI/bcftools query -f "%CHROM\t%POS\t%END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\n" $PAV \
  | awk -v OFS="\t" "{ sz = (\$5 < 0 ? -\$5 : \$5); if (sz>=500000) print \$0 }" \
  > /tmp/pav_large_500kb.tsv
echo "TOTAL PAV >=500kb: $(wc -l < /tmp/pav_large_500kb.tsv)"
echo

echo "=== ALL 33+ PAV SVs >=500kb (CHROM POS END SVTYPE SVLEN ID) ==="
sort -k5,5 -nr /tmp/pav_large_500kb.tsv | column -t

echo
echo "=== ANNOTATE: which fall in ACTIVE HOR vs other CenSat vs OUTSIDE censat ==="
awk -v OFS="\t" "{ sz = (\$5 < 0 ? -\$5 : \$5); print \$1, \$2-1, \$2+sz, \$4\"_\"sz\"_\"\$6 }" /tmp/pav_large_500kb.tsv > /tmp/pav_large_500kb.bed

echo "--- in ACTIVE HOR (chm13v2.0_censat_v2.1.active.bed) ---"
$PIXI/bedtools intersect -a /tmp/pav_large_500kb.bed -b $ACTIVE -u | column -t || true
echo "--- in ANY CenSat region ---"
$PIXI/bedtools intersect -a /tmp/pav_large_500kb.bed -b $ALL_CEN -u | column -t || true
echo "--- OUTSIDE all CenSat regions (truly non-centromeric) ---"
$PIXI/bedtools intersect -a /tmp/pav_large_500kb.bed -b $ALL_CEN -v | column -t || true
