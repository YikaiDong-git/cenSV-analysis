#!/usr/bin/env python3
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
MIN_CODE = 2


def norm_sex(s):
    s = (s or "").strip().lower()
    return {"male": "male", "female": "female", "m": "male", "f": "female"}.get(s, "")


def main():
    targets = list(csv.DictReader(open(TARGETS), delimiter="\t"))
    phen_has_map = set()
    code2phen = defaultdict(set)
    for r in csv.DictReader(open(TMAP), delimiter="\t"):
        if r.get("map_status") == "ok":
            code2phen[(r["vocabulary_id"], r["concept_code"])].add(r["phenotype_id"])
            phen_has_map.add(r["phenotype_id"])
    incl = [t for t in targets if t["phenotype_id"] in phen_has_map]
    phidx = {t["phenotype_id"]: i for i, t in enumerate(incl)}
    code2idx = {k: [phidx[p] for p in ps if p in phidx] for k, ps in code2phen.items()}
    code2idx = {k: v for k, v in code2idx.items() if v}
    print(f"{len(incl)} phenotypes with map; {len(code2idx)} (vocab,code) keys", flush=True)

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

    state = {}; nrow = 0
    with open(EV) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[3]:
                continue
            idxs = code2idx.get((p[1], p[2]))
            if not idxs or p[0] not in base_sex:
                continue
            try:
                iid = int(p[0]); di = int(p[3].replace("-", ""))
            except (ValueError, IndexError):
                continue
            for pi in idxs:
                k = iid * 100 + pi
                v = state.get(k, 0)
                if v == 0:
                    state[k] = di
                elif v != -1 and v != di:
                    state[k] = -1
            nrow += 1
            if nrow % 5_000_000 == 0:
                print(f"  streamed {nrow:,} rows; {len(state):,} (person,phecode) pairs", flush=True)
    print(f"streamed {nrow:,} matched evidence rows total", flush=True)

    cases = [set() for _ in incl]
    for k, v in state.items():
        if v == -1:
            cases[k % 100].add(k // 100)
    del state

    base_iids = sorted(base_sex)
    cols = [t["phenotype_id"] for t in incl]
    with open(f"{OUT}/phenotypes_all_mcc2.tsv", "w") as fo:
        fo.write("IID\t" + "\t".join(cols) + "\n")
        for iid in base_iids:
            sx = base_sex[iid]; iidi = int(iid); out = [iid]
            for j, t in enumerate(incl):
                sr = str(t.get("sex_restriction", "both")).strip().lower()
                eligible = sr not in ("male", "female") or sx == sr
                out.append("" if not eligible else ("1" if iidi in cases[j] else "0"))
            fo.write("\t".join(out) + "\n")
    with open(f"{OUT}/phenotype_summary_all_mcc2.tsv", "w") as fs:
        fs.write("phenotype_id\tlabel\tphecodes\tsex_restriction\tcategory\tn_eligible\tcases\tcontrols\tinclude\tstatus\n")
        for j, t in enumerate(incl):
            pid = t["phenotype_id"]; sr = str(t.get("sex_restriction", "both")).strip().lower()
            elig = [i for i in base_iids if sr not in ("male", "female") or base_sex[i] == sr]
            ncase = sum(1 for i in elig if int(i) in cases[j])
            fs.write(f"{pid}\t{t['label']}\t{t.get('phecodes','')}\t{sr}\t{t.get('category','')}\t"
                     f"{len(elig)}\t{ncase}\t{len(elig)-ncase}\tTrue\tincluded_all_ancestry_mcc2\n")
    print(f"wrote phenotypes_all_mcc2.tsv ({len(base_iids)} IIDs x {len(cols)} diseases) + summary", flush=True)
    for j, t in enumerate(incl[:6]):
        print(f"    {t['phenotype_id']:24} cases(MCC2)={len(cases[j])}", flush=True)


if __name__ == "__main__":
    main()
