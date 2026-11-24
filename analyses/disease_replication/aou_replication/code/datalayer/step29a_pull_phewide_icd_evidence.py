#!/usr/bin/env python3
import os, sys, csv
from google.cloud import bigquery

CDR = os.environ["WORKSPACE_CDR"]; PROJ = os.environ["GOOGLE_CLOUD_PROJECT"]
BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
XMAP = f"{BASE}/aou_replication/results_vds/step21_catalog/phecodeX/phecodeX_map_expanded.tsv"
CAT = f"{BASE}/aou_replication/results_vds/step21_catalog/icd_code_catalog.tsv"
OUT = ("/tmp/claude-1000/-home-dataproc-workspace-working-ydon268/"
       "86d47bc2-43bf-4157-b8ed-79f770cbd7e0/scratchpad/phewide_evidence"); os.makedirs(OUT, exist_ok=True)


def code_set():
    xm = set()
    for r in csv.DictReader(open(XMAP), delimiter="\t"):
        xm.add((r["vocabulary_id"], r["concept_code"]))
    aou = set()
    for r in csv.DictReader(open(CAT), delimiter="\t"):
        aou.add((r["vocabulary"], r["code"]))
    keep = xm & aou
    return keep


def build_sql(keep):
    lits = ",".join("'" + v + "|" + c.replace("'", "") + "'" for (v, c) in sorted(keep))
    return f"""
WITH keep AS (SELECT k FROM UNNEST([{lits}]) AS k),
ev AS (
  SELECT co.person_id AS pid, c.vocabulary_id AS vocab, c.concept_code AS code,
         CAST(co.condition_start_date AS STRING) AS dt
  FROM `{CDR}.condition_occurrence` co JOIN `{CDR}.concept` c
    ON co.condition_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
  UNION ALL
  SELECT o.person_id, c.vocabulary_id, c.concept_code, CAST(o.observation_date AS STRING)
  FROM `{CDR}.observation` o JOIN `{CDR}.concept` c ON o.observation_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
  UNION ALL
  SELECT p.person_id, c.vocabulary_id, c.concept_code, CAST(p.procedure_date AS STRING)
  FROM `{CDR}.procedure_occurrence` p JOIN `{CDR}.concept` c ON p.procedure_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
  UNION ALL
  SELECT m.person_id, c.vocabulary_id, c.concept_code, CAST(m.measurement_date AS STRING)
  FROM `{CDR}.measurement` m JOIN `{CDR}.concept` c ON m.measurement_source_concept_id = c.concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
)
SELECT DISTINCT pid, vocab, code, dt
FROM ev
WHERE CONCAT(vocab,'|',code) IN (SELECT k FROM keep) AND dt IS NOT NULL
"""


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "dryrun"
    keep = code_set()
    print(f"phecodeX∩AoU codes kept: {len(keep)}", flush=True)
    sql = build_sql(keep)
    cli = bigquery.Client(project=PROJ)
    dry = cli.query(sql, job_config=bigquery.QueryJobConfig(dry_run=True, use_query_cache=False))
    gb = dry.total_bytes_processed / 1e9
    print(f"[icd_evidence] scan {gb:.1f} GB -> ${gb/1000*6.25:.3f}", flush=True)
    if mode == "dryrun":
        print("DRYRUN_DONE (dry_run only, no data pulled, $0)", flush=True)
        return
    res = cli.query(sql).result(page_size=100000)
    p = f"{OUT}/phewide_icd_evidence.tsv"; n = 0
    with open(p, "w") as f:
        f.write("IID\tvocabulary_id\tconcept_code\tcode_date\n")
        for r in res:
            f.write(f"{r.pid}\t{r.vocab}\t{r.code}\t{r.dt}\n"); n += 1
            if n % 20_000_000 == 0:
                print(f"  wrote {n:,} rows", flush=True)
    print(f"wrote {p} ({n:,} rows). IN-PERIMETER.", flush=True)
    print("ICD_PULL_DONE", flush=True)


if __name__ == "__main__":
    main()
