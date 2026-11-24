#!/usr/bin/env python3
import os, sys, csv
from google.cloud import bigquery

CDR = os.environ["WORKSPACE_CDR"]; PROJ = os.environ["GOOGLE_CLOUD_PROJECT"]
BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
CATL = f"{BASE}/aou_replication/results_vds/step21_catalog/loinc_catalog.tsv"
OUT = ("/tmp/claude-1000/-home-dataproc-workspace-working-ydon268/"
       "86d47bc2-43bf-4157-b8ed-79f770cbd7e0/scratchpad/phewide_evidence"); os.makedirs(OUT, exist_ok=True)
NMIN = 1000


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "dryrun"
    loincs = [r["loinc"] for r in csv.DictReader(open(CATL), delimiter="\t") if int(r["n_people"]) >= NMIN]
    print(f"quant panel: {len(loincs)} LOINCs (>= {NMIN} people)", flush=True)
    sql = f"""
WITH v AS (
  SELECT m.person_id AS pid, c.concept_code AS loinc, m.value_as_number AS val
  FROM `{CDR}.measurement` m JOIN `{CDR}.concept` c ON m.measurement_concept_id = c.concept_id
  WHERE c.vocabulary_id = 'LOINC' AND c.concept_code IN UNNEST(@l) AND m.value_as_number IS NOT NULL
)
SELECT pid, loinc, APPROX_QUANTILES(val, 2)[OFFSET(1)] AS median_v FROM v GROUP BY pid, loinc
"""
    cfg = bigquery.QueryJobConfig(query_parameters=[bigquery.ArrayQueryParameter("l", "STRING", loincs)])
    cli = bigquery.Client(project=PROJ)
    dry = cli.query(sql, job_config=bigquery.QueryJobConfig(
        query_parameters=cfg.query_parameters, dry_run=True, use_query_cache=False))
    gb = dry.total_bytes_processed / 1e9
    print(f"[quant_medians] scan {gb:.1f} GB -> ${gb/1000*6.25:.3f}", flush=True)
    if mode == "dryrun":
        print("DRYRUN_DONE (dry_run only, $0)", flush=True); return
    p = f"{OUT}/quant_medians.tsv"; n = 0
    with open(p, "w") as f:
        f.write("IID\tloinc\tmedian_value\n")
        for r in cli.query(sql, job_config=cfg).result(page_size=100000):
            f.write(f"{r['pid']}\t{r['loinc']}\t{r['median_v']}\n"); n += 1
            if n % 10_000_000 == 0:
                print(f"  wrote {n:,} rows", flush=True)
    print(f"wrote {p} ({n:,} rows). IN-PERIMETER.", flush=True)
    print("QUANT_PULL_DONE", flush=True)


if __name__ == "__main__":
    main()
