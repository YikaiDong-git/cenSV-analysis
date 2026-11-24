#!/usr/bin/env python3
"""Enlarged phenotype pull: per-person median for the 44 AoU-available LOINCs (49 quant traits)."""
import os
from google.cloud import bigquery
BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
R = f"{BASE}/aou_replication/results/step6_222"
import csv
loincs = sorted({r['loinc'] for r in csv.DictReader(open(f"{R}/testable_trait_map.tsv"), delimiter='\t') if r['type'] == 'quant'})
cdr = os.environ["WORKSPACE_CDR"]; proj = os.environ.get("GOOGLE_CLOUD_PROJECT")
client = bigquery.Client(project=proj)
sql = f"""
WITH v AS (
  SELECT m.person_id AS pid, c.concept_code AS loinc, m.value_as_number AS val
  FROM `{cdr}.measurement` m JOIN `{cdr}.concept` c ON m.measurement_concept_id=c.concept_id
  WHERE c.vocabulary_id='LOINC' AND c.concept_code IN UNNEST(@l) AND m.value_as_number IS NOT NULL
)
SELECT pid, loinc, APPROX_QUANTILES(val,2)[OFFSET(1)] AS median_v FROM v GROUP BY pid, loinc
"""
cfg = bigquery.QueryJobConfig(query_parameters=[bigquery.ArrayQueryParameter("l", "STRING", loincs)])
print(f"pulling {len(loincs)} LOINCs (streaming)...", flush=True)
from collections import Counter
npers = Counter(); ntot = 0
with open(f"{R}/measurement_medians_222.tsv", "w") as f:
    f.write("IID\tloinc\tmedian_value\n")
    for r in client.query(sql, job_config=cfg).result(page_size=100000):
        f.write(f"{r['pid']}\t{r['loinc']}\t{r['median_v']}\n"); npers[r['loinc']] += 1; ntot += 1
        if ntot % 1000000 == 0:
            print(f"  {ntot} rows...", flush=True)
print(f"DONE: {ntot} person-LOINC rows across {len(npers)} LOINCs")
