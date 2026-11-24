#!/bin/bash
#SBATCH --job-name=pav_vs_svim
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/logs/pav_smoke_HG002_HG00099/compare-%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
ROOT=${CENSV_ROOT}
cd "$ROOT"
PAV=analyses/sv_call_pav_comparison/HG002.vcf.gz
SVIM_H1=$(ls variant_calling/svim_asm/HG002_hap1_*.vcf.gz 2>/dev/null | head -1)
SVIM_H2=$(ls variant_calling/svim_asm/HG002_hap2_*.vcf.gz 2>/dev/null | head -1)
CEN_BED=resources/chm13v2.0_censat_v2.1.active.bed
PIXIRUN="$HOME/.pixi/bin/pixi run --manifest-path $ROOT/pixi.toml"

echo "=== PAV HG002.vcf.gz total records ==="
$PIXIRUN bcftools view -H "$PAV" 2>/dev/null | wc -l

echo
echo "=== PAV non-SNV variants by SVTYPE x size ==="
$PIXIRUN bcftools view -H "$PAV" 2>/dev/null | awk -F"\t" "
{
  svlen=0; svtype=\"\";
  for(i=8;i<=NF;i++){
    if(match(\$i,/SVLEN=([0-9-]+)/,a)) { svlen=a[1]+0; if(svlen<0)svlen=-svlen }
    if(match(\$i,/SVTYPE=([^;]+)/,b)) svtype=b[1]
  }
  if(svtype!=\"\" && svtype!=\"SNV\"){
    if(svlen>=50)c50[svtype]++
    if(svlen>=100)c100[svtype]++
    if(svlen>=500)c500[svtype]++
    if(svlen>=1000)c1k[svtype]++
    if(svlen>=10000)c10k[svtype]++
  }
}
END {
  for(t in c50) printf \"  %-10s >=50:%d >=100:%d >=500:%d >=1000:%d >=10kb:%d\n\", t, c50[t]+0, c100[t]+0, c500[t]+0, c1k[t]+0, c10k[t]+0
}
"

echo
echo "=== svim-asm HG002 totals ==="
echo "h1: $SVIM_H1"
echo "h2: $SVIM_H2"
SVIM_H1_N=$($PIXIRUN bcftools view -H "$SVIM_H1" 2>/dev/null | wc -l)
SVIM_H2_N=$($PIXIRUN bcftools view -H "$SVIM_H2" 2>/dev/null | wc -l)
echo "  h1 records: $SVIM_H1_N | h2 records: $SVIM_H2_N | combined: $((SVIM_H1_N + SVIM_H2_N))"

echo
echo "=== centromere BED size (records) ==="
wc -l "$CEN_BED"

echo
echo "=== PAV non-SNV >=50bp inside centromere ==="
PAV_CEN_50=$($PIXIRUN bcftools view -H "$PAV" -R "$CEN_BED" 2>/dev/null | awk -F"\t" "{ for(i=8;i<=NF;i++){ if(match(\$i,/SVLEN=([0-9-]+)/,a)) { sv=a[1]+0; if(sv<0)sv=-sv; if(sv>=50)c++; break } } } END{print c+0}")
echo "  PAV cen >=50bp: $PAV_CEN_50"

echo
echo "=== svim-asm HG002 SVs in centromere (h1+h2) ==="
SVIM_CEN_H1=$($PIXIRUN bcftools view -H "$SVIM_H1" -R "$CEN_BED" 2>/dev/null | wc -l)
SVIM_CEN_H2=$($PIXIRUN bcftools view -H "$SVIM_H2" -R "$CEN_BED" 2>/dev/null | wc -l)
echo "  svim h1 cen: $SVIM_CEN_H1 | h2 cen: $SVIM_CEN_H2 | combined: $((SVIM_CEN_H1 + SVIM_CEN_H2))"

echo
echo "=== DONE ==="; date
