#!/usr/bin/env python3
"""Step 3b: pull per-person median measurement values for the AoU-testable quant traits."""
import os
from google.cloud import bigquery

cdr = os.environ["WORKSPACE_CDR"]
proj = os.environ.get("GOOGLE_CLOUD_PROJECT")
OUT = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results/step3")
os.makedirs(OUT, exist_ok=True)
client = bigquery.Client(project=proj)

L2T = {'8302-2':'height','39156-5':'BMI','1751-7':'albumin','718-7':'haemoglobin','777-3':'platelet_count',
 '776-5':'platelet_dist_width','742-7':'monocyte_count','736-9':'lymphocyte_pct','5905-5':'monocyte_pct',
 '17849-1':'reticulocyte_pct','60474-4':'retic_mean_vol','1869-7':'ApoA','1988-5':'CRP','33863-2':'cystatinC',
 '2484-4':'IGF1','2571-8':'triglycerides','3084-1':'urate'}
case = " ".join(f"WHEN '{k}' THEN '{v}'" for k, v in L2T.items())
sql = f"""
WITH vals AS (
  SELECT m.person_id AS pid, CASE c.concept_code {case} END AS trait, m.value_as_number AS v
  FROM `{cdr}.measurement` m JOIN `{cdr}.concept` c ON m.measurement_concept_id=c.concept_id
  WHERE c.vocabulary_id='LOINC' AND c.concept_code IN UNNEST(@loincs) AND m.value_as_number IS NOT NULL
)
SELECT pid, trait, APPROX_QUANTILES(v,2)[OFFSET(1)] AS median_v, COUNT(*) AS n_meas
FROM vals GROUP BY pid, trait
"""
cfg = bigquery.QueryJobConfig(query_parameters=[bigquery.ArrayQueryParameter("loincs", "STRING", list(L2T))])
print("running measurement median query (streaming, low-mem)...", flush=True)
from collections import Counter
npers = Counter()
job = client.query(sql, job_config=cfg)
ntot = 0
with open(f"{OUT}/measurement_medians.tsv", "w") as f:
    f.write("IID\ttrait\tmedian_value\tn_meas\n")
    for r in job.result(page_size=50000):
        f.write(f"{r['pid']}\t{r['trait']}\t{r['median_v']}\t{r['n_meas']}\n")
        npers[r['trait']] += 1
        ntot += 1
        if ntot % 500000 == 0:
            print(f"  streamed {ntot} rows...", flush=True)
rows = [None] * ntot
with open(f"{OUT}/measurement_trait_N.tsv", "w") as f:
    f.write("trait\tn_persons\n")
    for t, n in sorted(npers.items(), key=lambda x: -x[1]):
        f.write(f"{t}\t{n}\n")
print(f"DONE: {len(rows)} person-trait rows across {len(npers)} traits")
for t, n in sorted(npers.items(), key=lambda x: -x[1]):
    print(f"  {t:20} {n:>10}")
