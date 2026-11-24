#!/usr/bin/env python3
import os, csv, subprocess, array
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
XMAP = f"{BASE}/aou_replication/results_vds/step21_catalog/phecodeX/phecodeX_map_expanded.tsv"
XINFO = f"{BASE}/aou_replication/results_vds/step21_catalog/phecodeX/phecodeX_info.tsv"
COV = f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv"
SC = "/tmp/claude-1000/-home-dataproc-workspace-working-ydon268/86d47bc2-43bf-4157-b8ed-79f770cbd7e0/scratchpad"
EVID = f"{SC}/phewide_evidence/phewide_icd_evidence.tsv"
FAM = f"{SC}/phewide_pgen/sv.fam"
INP = f"{SC}/inputs"; os.makedirs(INP, exist_ok=True)
OUTA = f"{BASE}/aou_replication/results_vds/step30_phewide"; os.makedirs(OUTA, exist_ok=True)
MIN_CASE = 50
BASECOV = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
           "sex", *[f"gPC{i}" for i in range(1, 17)]]


def read_fam():
    iids, sex = [], {}
    with open(FAM) as f:
        for line in f:
            p = line.split()
            iids.append(p[1]); sex[p[1]] = p[4]
    return iids, sex


def build_covar(iids):
    rows = {}
    ancs = {}
    with open(COV) as f:
        h = f.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            iid = p[idx["IID"]]
            rows[iid] = [p[idx[c]] for c in BASECOV]
            ancs[iid] = p[idx["ancestry_pred"]]
    anc_levels = sorted({ancs[i] for i in iids if ancs.get(i)})
    drop = anc_levels[0] if anc_levels else None
    dums = [a for a in anc_levels if a != drop]
    covcols = list(BASECOV) + [f"anc_{a}" for a in dums]
    with open(f"{INP}/covar.txt", "w") as f:
        f.write("FID\tIID\t" + "\t".join(covcols) + "\n")
        for iid in iids:
            base = rows.get(iid)
            if base is None:
                continue
            sexnum = "1" if base[3] == "Male" else "0"
            vals = list(base); vals[3] = sexnum
            dv = ["1" if ancs.get(iid) == a else "0" for a in dums]
            f.write(f"{iid}\t{iid}\t" + "\t".join(vals + dv) + "\n")
    with open(f"{INP}/covar_names.txt", "w") as f:
        f.write(",".join(covcols))
    print(f"wrote covar.txt ({len(covcols)} covars: {len(BASECOV)} base + {len(dums)} ancestry dummies)", flush=True)


def main():
    iids, sex = read_fam()
    print(f"sample set (sv.fam) N={len(iids)}", flush=True)
    build_covar(iids)

    code2ph = defaultdict(list); phlist = []
    for r in csv.DictReader(open(XMAP), delimiter="\t"):
        ph = r["phecode"]
        if ph not in phlist:
            phlist.append(ph)
        code2ph[(r["vocabulary_id"], r["concept_code"])].append(ph)
    phidx = {p: i for i, p in enumerate(phlist)}
    code2idx = {k: [phidx[p] for p in v] for k, v in code2ph.items()}
    info = {r["phecode"]: (r["label"], r.get("sex", "Both")) for r in csv.DictReader(open(XINFO), delimiter="\t")}
    nph = len(phlist)
    print(f"{nph} phecodeX phecodes in map", flush=True)

    ev_sorted = f"{SC}/phewide_evidence/icd_evidence_sorted.tsv"
    if not os.path.exists(ev_sorted):
        print("sorting evidence by pid ...", flush=True)
        with open(ev_sorted, "w") as fo:
            subprocess.run(["sort", "-t", "\t", "-k1,1n", "-S", "128M", "--parallel=2", EVID],
                           stdout=fo, check=True)
    case_map = {}
    ccount = [0] * nph
    cur = None; ph_dates = defaultdict(set); nseen = 0
    def flush(pid):
        cases = array.array('i', (pi for pi, ds in ph_dates.items() if len(ds) >= 2))
        if cases:
            case_map[pid] = cases
            for pi in cases:
                ccount[pi] += 1
    with open(ev_sorted) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            pid = p[0]
            if pid != cur:
                if cur is not None:
                    flush(cur)
                cur = pid; ph_dates = defaultdict(set); nseen += 1
            for pi in code2idx.get((p[1], p[2]), ()):
                ph_dates[pi].add(p[3])
        if cur is not None:
            flush(cur)
    print(f"streamed evidence for {nseen:,} persons; {len(case_map):,} have >=1 MCC2 case-phecode", flush=True)

    kept = [i for i in range(nph) if ccount[i] >= MIN_CASE]
    print(f"phecodes with >= {MIN_CASE} cases: {len(kept)} / {nph}", flush=True)
    kept_ph = [phlist[i] for i in kept]

    def sr(ph):
        return (info.get(ph, ("", "Both"))[1] or "Both").strip().lower()
    with open(f"{INP}/disease.pheno", "w") as f:
        f.write("FID\tIID\t" + "\t".join(kept_ph) + "\n")
        for iid in iids:
            cs = set(case_map.get(iid, ()))
            sxn = sex.get(iid, "0")
            out = []
            for i, ph in zip(kept, kept_ph):
                s = sr(ph)
                if (s == "male" and sxn != "1") or (s == "female" and sxn != "2"):
                    out.append("NA")
                else:
                    out.append("1" if i in cs else "0")
            f.write(f"{iid}\t{iid}\t" + "\t".join(out) + "\n")
    with open(f"{OUTA}/phecodeX_pheno_summary.tsv", "w") as f:
        f.write("phecode\tlabel\tsex\tcases_mcc2\n")
        for i in kept:
            ph = phlist[i]; lab, sx = info.get(ph, ("", "Both"))
            f.write(f"{ph}\t{lab}\t{sx}\t{ccount[i]}\n")
    print(f"wrote disease.pheno ({len(iids)} samples x {len(kept)} phecodes) + summary -> {OUTA}", flush=True)
    print("DISEASE_PHENO_DONE", flush=True)


if __name__ == "__main__":
    main()
