"""Step 32 (phewide P3): replicate discovery's near-credible pairs in AoU by DIRECTION CONCORDANCE."""
import os, csv, math
import numpy as np

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
PH = f"{BASE}/phewide_all_handoff"
RV = f"{BASE}/aou_replication/results_vds"
OUT = f"{RV}/step31_phewide_ext2"; os.makedirs(OUT, exist_ok=True)
XW = f"{RV}/step30_phewide_ext2/reference/trait_crosswalk.tsv"
MATDIR = f"{RV}/step30_phewide_ext2/disease_matrices"


SCRATCH = f"{RV}/step30_phewide_ext2/quant"


def load_disease_from_matrices():
    if not os.path.exists(f"{MATDIR}/disease_Pv.npy"):
        print(f"  WARN missing disease matrices at {MATDIR}", flush=True); return {}
    needed = {r["aou_id"] for r in csv.DictReader(open(XW), delimiter="\t")
              if r.get("aou_type") == "phecode" and r.get("aou_id")}
    Pv = np.load(f"{MATDIR}/disease_Pv.npy"); Be = np.load(f"{MATDIR}/disease_Be.npy")
    svs = [ln.split("\t")[0] for ln in open(f"{MATDIR}/disease_svorder.txt") if ln.strip()]
    phcol = {r["phecode"]: j for j, r in enumerate(csv.DictReader(open(f"{MATDIR}/disease_meta.tsv"), delimiter="\t"))}
    d = {}
    for ph in needed:
        j = phcol.get(ph)
        if j is None:
            continue
        for i, sv in enumerate(svs):
            pv = Pv[i, j]
            if np.isfinite(pv):
                d[(sv, ph)] = (float(Be[i, j]), float(pv), "", "")
    print(f"  disease matrices: {len(needed)} needed phecodes -> {len(d)} (sv,phecode) cells", flush=True)
    return d


def load_p2():
    d = {}
    qp = f"{SCRATCH}/phewide_quant_allpanel_ext2.tsv"
    if os.path.exists(qp):
        for r in csv.DictReader(open(qp), delimiter="\t"):
            key = "loinc" if "loinc" in r else "phecode" if "phecode" in r else "pheno"
            try:
                b = float(r["beta"])
            except (ValueError, KeyError):
                b = float("nan")
            try:
                pv = float(r["p"])
            except (ValueError, KeyError):
                pv = float("nan")
            d[(r["sv_locus"], r[key])] = (b, pv, r.get("q_bh_tier", ""), r.get("n", ""))
    else:
        print(f"  WARN missing {qp}", flush=True)
    d.update(load_disease_from_matrices())
    return d


def main():
    xw = {(r["study"], r["trait_code"]): r for r in csv.DictReader(open(XW), delimiter="\t")}
    pairs = list(csv.DictReader(open(f"{PH}/sv_trait_pairs.tsv"), delimiter="\t"))
    p2 = load_p2()

    rows = []
    for pr in pairs:
        m = xw.get((pr["study"], pr["trait"]))
        aou_id = m["aou_id"] if m else ""
        dsign = pr.get("discovery_sign", "")
        if not aou_id or dsign not in ("+", "-"):
            continue
        hit = p2.get((pr["sv_locus"], aou_id))
        if hit is None:
            rows.append([pr["sv_locus"], pr["study"], pr["trait"], m["aou_type"], aou_id, m["aou_name"],
                         dsign, "", "", "", "", "not_tested_in_aou"]); continue
        b, pv, q, n = hit
        asign = "+" if b > 0 else "-" if b < 0 else ""
        conc = "concordant" if (asign and asign == dsign) else "discordant" if asign else "na"
        rows.append([pr["sv_locus"], pr["study"], pr["trait"], m["aou_type"], aou_id, m["aou_name"],
                     dsign, asign, f"{b:.4g}" if math.isfinite(b) else "",
                     f"{pv:.3g}" if math.isfinite(pv) else "", q, conc])

    with open(f"{OUT}/replication_concordance.tsv", "w") as f:
        f.write("sv_locus\tstudy\ttrait\taou_type\taou_id\taou_name\tdiscovery_sign\taou_sign\t"
                "aou_beta\taou_p\taou_q_bh_tier\tconcordance\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    tested = [r for r in rows if r[11] in ("concordant", "discordant")]
    conc = [r for r in tested if r[11] == "concordant"]
    sig = [r for r in tested if r[9] and float(r[9]) < 0.05]
    sig_conc = [r for r in sig if r[11] == "concordant"]
    print("=== P3 direction-concordance replication ===", flush=True)
    print(f"  mapped pairs with an AoU test: {len(tested)}", flush=True)
    print(f"  concordant sign: {len(conc)}/{len(tested)} ({100*len(conc)/max(1,len(tested)):.0f}%) "
          f"[null=50%]", flush=True)
    print(f"  nominally significant in AoU (p<0.05): {len(sig)}; of those concordant: "
          f"{len(sig_conc)}/{len(sig)} ({100*len(sig_conc)/max(1,len(sig)):.0f}%)", flush=True)
    print(f"\n  top concordant significant pairs:")
    for r in sorted(sig_conc, key=lambda r: float(r[9]))[:15]:
        print(f"    {r[0][:24]:24} {r[5][:26]:26} d{r[6]}/a{r[7]} p={r[9]} q={r[10]}", flush=True)
    print(f"\nwrote {OUT}/replication_concordance.tsv", flush=True)


if __name__ == "__main__":
    main()
