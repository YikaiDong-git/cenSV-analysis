#!/usr/bin/env python3
import os, csv

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
D = f"{BASE}/aou_replication/results_vds/step20_disease_all"
MAP = f"{HOME}/workspace/working/PRS/resources/phewas/target_phecode_icd_map.tsv"
OUT = D


def main():
    code2ph = {}
    phlist = []
    for r in csv.DictReader(open(MAP), delimiter="\t"):
        ph = r["phenotype_id"]
        if ph not in phlist:
            phlist.append(ph)
        code2ph.setdefault((r["vocabulary_id"], r["concept_code"]), []).append(phlist.index(ph))
    phidx = {p: i for i, p in enumerate(phlist)}
    label = {}
    for r in csv.DictReader(open(MAP), delimiter="\t"):
        label[r["phenotype_id"]] = r.get("label", r.get("phecode", ""))
    print(f"{len(phlist)} phenotypes, {len(code2ph)} (vocab,code) keys", flush=True)

    state = {}
    nrow = 0
    with open(f"{D}/aou_target_icd_evidence.tsv") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            phs = code2ph.get((p[1], p[2]))
            if not phs:
                continue
            try:
                iid = int(p[0]); di = int(p[3].replace("-", ""))
            except (ValueError, IndexError):
                continue
            for pi in phs:
                k = iid * 100 + pi
                v = state.get(k, 0)
                if v == 0:
                    state[k] = di
                elif v != -1 and v != di:
                    state[k] = -1
            nrow += 1
            if nrow % 5_000_000 == 0:
                print(f"  {nrow:,} rows; {len(state):,} (person,phecode) pairs", flush=True)

    mcc1 = [0] * len(phlist); mcc2 = [0] * len(phlist)
    for k, v in state.items():
        pi = k % 100
        mcc1[pi] += 1
        if v == -1:
            mcc2[pi] += 1

    rows = sorted(range(len(phlist)), key=lambda i: -mcc1[i])
    print(f"\n{'phenotype':26} {'MCC1':>9} {'MCC2':>9} {'kept%':>6}", flush=True)
    with open(f"{OUT}/mcc_compare_31disease.tsv", "w") as fo:
        fo.write("phenotype_id\tlabel\tcases_mcc1\tcases_mcc2\tpct_kept_mcc2\n")
        for i in rows:
            ph = phlist[i]; pct = 100 * mcc2[i] / mcc1[i] if mcc1[i] else 0
            fo.write(f"{ph}\t{label.get(ph,'')}\t{mcc1[i]}\t{mcc2[i]}\t{pct:.1f}\n")
            print(f"{ph[:26]:26} {mcc1[i]:9d} {mcc2[i]:9d} {pct:5.1f}", flush=True)
    tot1, tot2 = sum(mcc1), sum(mcc2)
    print(f"\nTOTAL case-pairs: MCC1={tot1:,}  MCC2={tot2:,}  ({100*tot2/tot1:.1f}% retained)", flush=True)
    print(f"wrote {OUT}/mcc_compare_31disease.tsv", flush=True)


if __name__ == "__main__":
    main()
