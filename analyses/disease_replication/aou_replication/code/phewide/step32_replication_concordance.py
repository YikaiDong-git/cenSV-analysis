"""Step 32 (phewide P3): replicate discovery's near-credible pairs in AoU by DIRECTION CONCORDANCE."""
import os, csv, math

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
PH = f"{BASE}/phewide_all_handoff"
RV = f"{BASE}/aou_replication/results_vds"
XW = f"{RV}/step31_phewide/trait_crosswalk.tsv"
OUT = f"{RV}/step31_phewide"


SCRATCH = "/tmp/claude-1000/-home-dataproc-workspace-working-ydon268/86d47bc2-43bf-4157-b8ed-79f770cbd7e0/scratchpad"


def load_p2():
    d = {}
    srcs = [f"{SCRATCH}/phewide_quant_allpanel.tsv", f"{SCRATCH}/phewide_disease_allpanel.tsv"]
    for p in srcs:
        if not os.path.exists(p):
            print(f"  WARN missing {p}", flush=True); continue
        for r in csv.DictReader(open(p), delimiter="\t"):
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
