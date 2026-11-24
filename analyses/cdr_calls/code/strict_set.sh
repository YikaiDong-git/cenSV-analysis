#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -uo pipefail
D=${CENSV_ROOT}/analyses/cdr_calls
P=${CENSV_ROOT}/subproject_methylation
B=${CENSV_WORKSPACE}/tools/CDR-Finder
PY=$B/env/.pixi/envs/default/bin/python; CF=$B/upstream/workflow/scripts/cdr_finder.py
awk -F'\t' 'NR>1{print $1"_"$2"\t"$6}' $P/manifest/methylation_manifest.tsv | sort -u > /tmp/chem2.tsv
echo -e 'sample\thap\tchemistry\tn_cdrs_strict\tn_cen_strict\tas_hor_flag' > $D/output/cohort_strict_0p5.tsv
for d in $D/output/calls/*/; do
  t=$(basename $d); BAA=$d/binned_aa.bed; [ -s "$BAA" ] || continue
  $PY $CF -i "$BAA" --thr_height_perc_valley 0.5 --thr_prom_perc_valley 0.5 --bp_merge 1 --bp_edge 500000 --baseline_avg_methyl 0.4 --edge_height_heuristic min > $d/cdrs_strict.bed 2>/dev/null
  NC=$(wc -l < $d/cdrs_strict.bed); NCH=$(cut -f1 $d/cdrs_strict.bed|sort -u|wc -l)
  SA=${t%_hap*}; HA=hap${t##*_hap}
  chem=$(awk -F'\t' -v k=$t '$1==k{print $2}' /tmp/chem2.tsv); chem=${chem:-NA}
  case $SA in HG00733|HG02818|NA19036|NA19240) flag=as_hor;; *) flag=globe;; esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$SA" "$HA" "$chem" "$NC" "$NCH" "$flag" >> $D/output/cohort_strict_0p5.tsv
done
echo STRICT_DONE; wc -l $D/output/cohort_strict_0p5.tsv
