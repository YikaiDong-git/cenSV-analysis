#!/usr/bin/env python3
import os
from google.cloud import bigquery

CDR = os.environ["WORKSPACE_CDR"]; PROJ = os.environ["GOOGLE_CLOUD_PROJECT"]
OUT = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/"
                         "disease_replication/aou_replication/results_vds/step21_catalog")
os.makedirs(OUT, exist_ok=True)

Q_MEASURE = f"""
SELECT c.concept_code AS loinc, ANY_VALUE(c.concept_name) AS name,
       COUNT(DISTINCT m.person_id) AS n_people, COUNT(*) AS n_rows
FROM `{CDR}.measurement` m JOIN `{CDR}.concept` c ON m.measurement_concept_id = c.concept_id
WHERE c.vocabulary_id = 'LOINC' AND m.value_as_number IS NOT NULL
GROUP BY loinc
HAVING n_people >= 100
ORDER BY n_people DESC
"""

Q_ICD = f"""
WITH ev AS (
  SELECT co.person_id AS pid, c.vocabulary_id AS vocab, c.concept_code AS code
  FROM `{CDR}.condition_occurrence` co JOIN `{CDR}.concept` c
    ON co.condition_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
  UNION ALL
  SELECT o.person_id, c.vocabulary_id, c.concept_code
  FROM `{CDR}.observation` o JOIN `{CDR}.concept` c ON o.observation_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
  UNION ALL
  SELECT p.person_id, c.vocabulary_id, c.concept_code
  FROM `{CDR}.procedure_occurrence` p JOIN `{CDR}.concept` c ON p.procedure_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
  UNION ALL
  SELECT m.person_id, c.vocabulary_id, c.concept_code
  FROM `{CDR}.measurement` m JOIN `{CDR}.concept` c ON m.measurement_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
)
SELECT vocab, code, COUNT(DISTINCT pid) AS n_people
FROM ev GROUP BY vocab, code HAVING n_people >= 20 ORDER BY n_people DESC
"""


def run(cli, sql, label):
    dry = cli.query(sql, job_config=bigquery.QueryJobConfig(dry_run=True, use_query_cache=False))
    gb = dry.total_bytes_processed / 1e9
    print(f"[{label}] scan {gb:.1f} GB -> ${gb/1000*6.25:.3f}; running...", flush=True)
    return cli.query(sql).result()


def main():
    cli = bigquery.Client(project=PROJ)
    rows = run(cli, Q_MEASURE, "measurement")
    n = 0
    with open(f"{OUT}/loinc_catalog.tsv", "w") as f:
        f.write("loinc\tname\tn_people\tn_rows\n")
        for r in rows:
            f.write(f"{r.loinc}\t{r.name}\t{r.n_people}\t{r.n_rows}\n"); n += 1
    print(f"  wrote loinc_catalog.tsv: {n} LOINCs (>=100 people)", flush=True)
    rows = run(cli, Q_ICD, "icd_4domain")
    n = 0
    with open(f"{OUT}/icd_code_catalog.tsv", "w") as f:
        f.write("vocabulary\tcode\tn_people\n")
        for r in rows:
            f.write(f"{r.vocab}\t{r.code}\t{r.n_people}\n"); n += 1
    print(f"  wrote icd_code_catalog.tsv: {n} ICD codes (>=20 people)", flush=True)
    print("CATALOG_DONE", flush=True)


if __name__ == "__main__":
    main()
