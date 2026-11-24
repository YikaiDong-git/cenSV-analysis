#!/bin/bash
#SBATCH -J a1_truvari
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -t 4:00:00
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
D=${CENSV_ROOT}/analyses/sv_call_reliability
B=${CENSV_COLLAB_ROOT}/exploration/simulate_svim_asm_accuracy/results/multichrom_sv_v1/bundles
REF=${CENSV_ROOT}/resources/chm13v2.0.fa.gz
ENV=${CENSV_ROOT}/.pixi/envs/default/bin
export PATH=$ENV:$PATH
O=$D/output/a1_truvari
mkdir -p $O/truth $O/calls $O/bench
$ENV/python $D/code/a1_truvari/build_truth_vcfs.py $O/truth
for t in $O/truth/*.truth.vcf; do
  bid=$(basename $t .truth.vcf)
  bcftools sort -Oz -o $O/truth/$bid.truth.vcf.gz $t 2>/dev/null && tabix -f -p vcf $O/truth/$bid.truth.vcf.gz
  echo SAMPLE > $O/calls/_s.txt
  bcftools reheader -s $O/calls/_s.txt $B/$bid/alpha_mixed/active.vcf.gz -o $O/calls/$bid.calls.vcf.gz
  tabix -f -p vcf $O/calls/$bid.calls.vcf.gz
  for cfg in hg002 seq0 seq0_big; do
    case $cfg in
      hg002)    ARGS="--refdist 500 --pctseq 0.7 --pctsize 0.7 --pick single" ;;
      seq0)     ARGS="--refdist 500 --pctseq 0 --pctsize 0.7 --pick single" ;;
      seq0_big) ARGS="--refdist 1000 --pctseq 0 --pctsize 0.5 --pick single --sizemax 1000000" ;;
    esac
    rm -rf $O/bench/$cfg/$bid
    mkdir -p $O/bench/$cfg
    truvari bench -b $O/truth/$bid.truth.vcf.gz -c $O/calls/$bid.calls.vcf.gz \
      -o $O/bench/$cfg/$bid $ARGS --reference $REF --includebed $O/truth/$bid.active.bed \
      > $O/bench/$cfg/$bid.log 2>&1
  done
done
echo -e "bundle\tcalls_ge50\tcalls_ge50_le50k" > $O/control_calls.tsv
for d in $B/*/control; do
  bid=$(basename $(dirname $d))
  n=$(bcftools view -H $d/active.vcf.gz | awk -F'\t' '{l=length($4)-length($5); if(l<0)l=-l; m=0; if(match($8,/SVLEN=-?[0-9]+/)){m=substr($8,RSTART+6,RLENGTH-6); if(m<0)m=-m}; L=(m>0?m:l); if(L>=50)a++; if(L>=50&&L<=50000)b++} END{printf "%d\t%d", a+0, b+0}')
  echo -e "$bid\t$n" >> $O/control_calls.tsv
done
$ENV/python $D/code/a1_truvari/summarize.py $O
echo A1_DONE
