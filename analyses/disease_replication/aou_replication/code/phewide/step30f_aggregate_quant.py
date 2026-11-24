#!/usr/bin/env python3
import os, glob, csv
import numpy as np

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
SC = "/tmp/phewide"; INP = f"{SC}/inputs"; GLM = f"{SC}/glm/quant"
OUTA = f"{BASE}/aou_replication/results_vds/step30_phewide"
QSUM = f"{OUTA}/quant_pheno_summary.tsv"


def bh(p):
    ps = np.asarray(p, float); ok = np.isfinite(ps); n = int(ok.sum())
    q = np.full(len(ps), np.nan)
    if n == 0:
        return q
    order = np.argsort(np.where(ok, ps, 2.0)); prev = 1.0
    for rank, i in enumerate(order[::-1]):
        if not ok[i]:
            continue
        prev = min(prev, ps[i] * n / (n - rank)); q[i] = prev
    return q


def main():
    tier = {r["sv_locus"]: r["tier"] for r in csv.DictReader(open(f"{INP}/sv_tier.tsv"), delimiter="\t")}
    name = {r["loinc"]: r["name"] for r in csv.DictReader(open(QSUM), delimiter="\t")}
    files = sorted(glob.glob(f"{GLM}/q.*.glm.linear"))
    print(f"parsing {len(files)} LOINC files ...", flush=True)
    SV, PH = [], []
    BE, SE, PV, NN, AF = [], [], [], [], []
    for fi, fp in enumerate(files):
        ph = os.path.basename(fp)[len("q."):-len(".glm.linear")]
        with open(fp) as f:
            h = f.readline().rstrip("\n").split("\t")
            iid, it, ib, ip = h.index("ID"), h.index("TEST"), h.index("BETA"), h.index("P")
            ise, inn, iaf = h.index("SE"), h.index("OBS_CT"), h.index("A1_FREQ")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p[it] != "ADD":
                    continue
                pv = p[ip]
                SV.append(p[iid]); PH.append(ph)
                BE.append(p[ib]); SE.append(p[ise]); NN.append(p[inn]); AF.append(p[iaf])
                PV.append(pv if pv not in ("", "NA") else "nan")
        if (fi + 1) % 200 == 0:
            print(f"  {fi+1}/{len(files)} files, {len(SV)} rows", flush=True)
    T = len(SV)
    PVf = np.array([float(x) for x in PV], np.float64)
    tiers = np.array([tier.get(s, "") for s in SV])
    print(f"{T} tests parsed; tiered BH ...", flush=True)
    q = np.full(T, np.nan)
    for tn in ("signal", "null_ext"):
        idx = np.nonzero((tiers == tn) & np.isfinite(PVf))[0]
        if idx.size:
            q[idx] = bh(PVf[idx])
    n_q05 = int(np.nansum(q < 0.05))
    keep = np.nonzero((np.nan_to_num(q, nan=1.0) < 0.05) | (PVf < 1e-5))[0]
    keep = keep[np.argsort(PVf[keep])]
    out_c = f"{OUTA}/phewide_quant_full.tsv"
    with open(out_c, "w") as f:
        f.write("sv_locus\ttier\tloinc\tname\tn\ta1_freq\tbeta\tse\tp\tq_bh_tier\n")
        for i in keep:
            f.write(f"{SV[i]}\t{tiers[i]}\t{PH[i]}\t{name.get(PH[i],'')}\t{NN[i]}\t{AF[i]}\t{BE[i]}\t{SE[i]}\t"
                    f"{'' if not np.isfinite(PVf[i]) else f'{PVf[i]:.3g}'}\t"
                    f"{'' if not np.isfinite(q[i]) else f'{q[i]:.4g}'}\n")
    dump = f"{SC}/phewide_quant_full_all.tsv"
    with open(dump, "w") as f:
        f.write("sv_locus\ttier\tloinc\tn\tbeta\tp\n")
        for i in range(T):
            f.write(f"{SV[i]}\t{tiers[i]}\t{PH[i]}\t{NN[i]}\t{BE[i]}\t"
                    f"{'' if not np.isfinite(PVf[i]) else f'{PVf[i]:.3g}'}\n")
    print(f"QUANT_AGG_DONE: {T} tests, {n_q05} q<0.05 (tiered), {len(keep)} written -> {out_c}; "
          f"full dump -> {dump}", flush=True)


if __name__ == "__main__":
    main()
