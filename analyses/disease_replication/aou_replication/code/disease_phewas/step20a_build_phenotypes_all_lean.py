#!/usr/bin/env python3
"""Step 20a (Task D): memory-lean, STREAMING rebuild of the all-ancestry disease phenotype matrix."""
import csv, os
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
OUT = f"{BASE}/aou_replication/results_vds/step20_disease_all"
EV = f"{OUT}/aou_target_icd_evidence.tsv"
PRS = f"{HOME}/workspace/working/PRS"
COV = f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv"
DENOM = f"{PRS}/results/phewas/aou_ehr_icd_denominator.tsv"
TARGETS = f"{PRS}/resources/phewas/target_phenotypes.tsv"
TMAP = f"{PRS}/resources/phewas/target_phecode_icd_map.tsv"
MIN_CODE = 1


def norm_sex(s):
    s = (s or "").strip().lower()
    return {"male": "male", "female": "female", "m": "male", "f": "female"}.get(s, "")


def main():
    code2phen = defaultdict(set); phen_has_map = set()
    for r in csv.DictReader(open(TMAP), delimiter="\t"):
        if r.get("map_status") == "ok":
            code2phen[(r["vocabulary_id"], r["concept_code"])].add(r["phenotype_id"])
            phen_has_map.add(r["phenotype_id"])
    print(f"{len(code2phen)} (vocab,code) keys; {len(phen_has_map)} phenotypes with ICD map", flush=True)

    targets = list(csv.DictReader(open(TARGETS), delimiter="\t"))

    denom = set()
    with open(DENOM) as f:
        next(f)
        for line in f:
            denom.add(line.strip().split("\t")[0])
    base_sex = {}
    with open(COV) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ii, si = hdr.index("IID"), hdr.index("sex")
        for line in f:
            p = line.rstrip("\n").split("\t")
            iid = p[ii]; sx = norm_sex(p[si])
            if sx and iid in denom:
                base_sex[iid] = sx
    print(f"base cohort (all-anc cov ∩ denominator, sex-known) = {len(base_sex)}", flush=True)

    cases = defaultdict(set); nrow = 0
    with open(EV) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[3]:
                continue
            phens = code2phen.get((p[1], p[2]))
            if phens and p[0] in base_sex:
                for ph in phens:
                    cases[ph].add(p[0])
            nrow += 1
            if nrow % 5_000_000 == 0:
                print(f"  streamed {nrow:,} evidence rows", flush=True)
    print(f"streamed {nrow:,} evidence rows total", flush=True)

    base_iids = sorted(base_sex)
    incl = [t for t in targets if t["phenotype_id"] in phen_has_map]
    cols = [t["phenotype_id"] for t in incl]
    with open(f"{OUT}/phenotypes_all.tsv", "w") as fo:
        fo.write("IID\t" + "\t".join(cols) + "\n")
        for iid in base_iids:
            sx = base_sex[iid]; out = [iid]
            for t in incl:
                sr = str(t.get("sex_restriction", "both")).strip().lower()
                eligible = sr not in ("male", "female") or sx == sr
                if not eligible:
                    out.append("")
                else:
                    out.append("1" if iid in cases[t["phenotype_id"]] else "0")
            fo.write("\t".join(out) + "\n")
    with open(f"{OUT}/phenotype_summary_all.tsv", "w") as fs:
        fs.write("phenotype_id\tlabel\tphecodes\tsex_restriction\tcategory\tn_eligible\tcases\tcontrols\tinclude\tstatus\n")
        for t in incl:
            pid = t["phenotype_id"]; sr = str(t.get("sex_restriction", "both")).strip().lower()
            elig = [i for i in base_iids if sr not in ("male", "female") or base_sex[i] == sr]
            ncase = sum(1 for i in elig if i in cases[pid])
            fs.write(f"{pid}\t{t['label']}\t{t.get('phecodes','')}\t{sr}\t{t.get('category','')}\t"
                     f"{len(elig)}\t{ncase}\t{len(elig)-ncase}\tTrue\tincluded_all_ancestry\n")
    print(f"wrote phenotypes_all.tsv ({len(base_iids)} IIDs x {len(cols)} diseases) + summary", flush=True)
    for t in incl[:6]:
        pid = t["phenotype_id"]
        print(f"    {pid:24} cases={len(cases[pid])}", flush=True)


if __name__ == "__main__":
    main()
