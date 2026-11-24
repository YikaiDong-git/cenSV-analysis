#!/usr/bin/env python3
"""Step 58 (ext2 arm) — cross-ancestry META + replication verdict for the ext2 pooled disease hits."""
import os, csv, sys
from collections import Counter, defaultdict
import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S
import step56_ancestry_disease_ext2 as P

MAT = P.MATDIR
ANCS = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
WELL = ["EUR", "AFR", "AMR"]
CHI2_SIG = 3.841459
HITS = f"{S.OUTA}/disease_matrices_firth/disease_hits_final.tsv"
CMP = f"{S.OUTA}/comparison/perancestry"; os.makedirs(CMP, exist_ok=True)
OUT = f"{CMP}/ancestry_meta_pooledhits_ext2.tsv"


def load_anc(anc):
    tag = f"_{anc}"
    if not os.path.exists(f"{MAT}/disease_U{tag}.npy"):
        return None
    U = np.load(f"{MAT}/disease_U{tag}.npy"); V = np.load(f"{MAT}/disease_V{tag}.npy")
    svs = [ln.split("\t")[0] for ln in open(f"{MAT}/disease_svorder{tag}.txt") if ln.strip()]
    meta = list(csv.DictReader(open(f"{MAT}/disease_meta{tag}.tsv"), delimiter="\t"))
    return dict(U=U, V=V, phcol={r["phecode"]: j for j, r in enumerate(meta)},
                cases={r["phecode"]: int(r["cases"]) for r in meta}, svrow={s: i for i, s in enumerate(svs)})


def cell_uv(d, sv, ph):
    i = d["svrow"].get(sv); j = d["phcol"].get(ph)
    if i is None or j is None:
        return np.nan, np.nan, 0
    return float(d["U"][i, j]), float(d["V"][i, j]), d["cases"].get(ph, 0)


def analyze(cells):
    well = {a: cells[a] for a in WELL if a in cells and np.isfinite(cells[a][1]) and cells[a][1] > 0}
    if not well:
        return dict(beta=np.nan, chi2=np.nan, p=np.nan, I2=np.nan, k=0, verdict="UNTESTED", strata={})
    sU = sum(U for U, V, c in well.values()); sV = sum(V for U, V, c in well.values())
    bM = sU / sV; chi2M = sU * sU / sV; pM = float(chi2dist.sf(chi2M, 1))
    betas = {a: U / V for a, (U, V, c) in well.items()}
    Q = sum(V * (betas[a] - bM) ** 2 for a, (U, V, c) in well.items()); k = len(well)
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 and k > 1 else 0.0
    strata = {}
    for a, (U, V, c) in well.items():
        obs = U * U / V
        strata[a] = dict(beta=betas[a], chi2=obs, p=float(chi2dist.sf(obs, 1)),
                         sig=obs >= CHI2_SIG, cases=c, V=V)
    _sig = [a for a in well if strata[a]["sig"]]
    ref_beta = max((abs(strata[a]["beta"]) for a in _sig), default=abs(bM))
    for a in well:
        strata[a]["exp_chi2"] = ref_beta * ref_beta * strata[a]["V"]
        strata[a]["powered"] = strata[a]["exp_chi2"] >= CHI2_SIG
    sig = [a for a in well if strata[a]["sig"]]
    sig_signs = {np.sign(strata[a]["beta"]) for a in sig}
    powered_null = [a for a in well if strata[a]["exp_chi2"] >= 8.0 and strata[a]["chi2"] < 3.84]
    if len(sig) == 0:
        v = "UNDERPOWERED"
    elif len(sig) >= 2 and len(sig_signs) > 1:
        v = "DISCORDANT"
    elif powered_null:
        v = "STRAT_SUSPECT"
    elif len(sig) >= 2:
        v = "REPLICATES" if I2 < 50 else "REPLICATES_HET"
    else:
        s = np.sign(strata[sig[0]]["beta"])
        opp = [a for a in well if a not in sig and strata[a]["chi2"] >= 2.0
               and np.sign(strata[a]["beta"]) != s]
        v = "STRAT_SUSPECT" if opp else "REPLICATES_PLIM"
    return dict(beta=bM, chi2=chi2M, p=pM, I2=I2, k=k, verdict=v, strata=strata)


def main():
    hits_file = sys.argv[1] if len(sys.argv) > 1 else HITS
    out_file = sys.argv[2] if len(sys.argv) > 2 else OUT
    anc = {a: load_anc(a) for a in ANCS}
    have = [a for a in ANCS if anc[a] is not None]
    print(f"ancestries with matrices: {have}", flush=True)
    print(f"hits file: {hits_file}\nout file:  {out_file}", flush=True)
    allhits = list(csv.DictReader(open(hits_file), delimiter="\t"))
    hits = [h for h in allhits if h.get("tier_label", "significant") == "significant"] or allhits
    print(f"pooled hits in file: {len(allhits)}; adjudicating {len(hits)} significant\n", flush=True)

    rows = []
    for h in hits:
        sv = h["sv_locus"]; ph = h["phecode"]
        cells = {a: cell_uv(anc[a], sv, ph) for a in have}
        rows.append((h, cells, analyze(cells)))

    with open(out_file, "w") as f:
        hdr = ["sv_locus", "tier", "phecode", "disease", "pooled_p_final", "pooled_beta", "pooled_q"]
        for a in ANCS:
            hdr += [f"{a}_beta", f"{a}_p", f"{a}_cases"]
        hdr += ["meta_beta_well", "meta_p_well", "I2", "k_well", "n_sig_well", "sig_strata", "verdict"]
        f.write("\t".join(hdr) + "\n")
        for h, cells, res in rows:
            sigs = [a for a in WELL if res["strata"].get(a, {}).get("sig")]
            out = [h["sv_locus"], h.get("tier", ""), h["phecode"], h.get("label", ""),
                   h.get("p_final", ""), h.get("beta_final_logOR", ""), h.get("q_bh_tier", "")]
            for a in ANCS:
                U, V, c = cells.get(a, (np.nan, np.nan, 0))
                b = (U / V) if (np.isfinite(V) and V > 0) else np.nan
                p = float(chi2dist.sf(U * U / V, 1)) if (np.isfinite(V) and V > 0) else np.nan
                out += [f"{b:.4g}" if np.isfinite(b) else "", f"{p:.3g}" if np.isfinite(p) else "", str(c)]
            out += [f"{res['beta']:.4g}" if np.isfinite(res['beta']) else "",
                    f"{res['p']:.3g}" if np.isfinite(res['p']) else "",
                    f"{res['I2']:.0f}" if np.isfinite(res['I2']) else "", str(res["k"]),
                    str(len(sigs)), ",".join(sigs), res["verdict"]]
            f.write("\t".join(out) + "\n")
    print(f"wrote {out_file}\n", flush=True)

    tally = Counter(res["verdict"] for _, _, res in rows)
    print("=== replication verdicts over the ext2 pooled significant cells (WELL = EUR/AFR/AMR) ===")
    for v, n in tally.most_common():
        print(f"  {v:18s} {n}")

    by_tier = defaultdict(Counter)
    for h, _c, res in rows:
        by_tier[h.get("tier", "?")][res["verdict"]] += 1
    print("\n=== verdicts BY TIER (does the new-predictor set behave like the old?) ===")
    verds = [v for v, _ in tally.most_common()]
    print(f"  {'tier':<12}{'n':>5}  " + "".join(f"{v[:13]:>15}" for v in verds))
    for tn in sorted(by_tier):
        tot = sum(by_tier[tn].values())
        print(f"  {tn:<12}{tot:>5}  " + "".join(f"{by_tier[tn].get(v,0):>15}" for v in verds))

    amr_only = [(h, res) for h, _c, res in rows
                if [a for a in WELL if res["strata"].get(a, {}).get("sig")] == ["AMR"]]
    eur_pow_null = [(h, res) for h, res in amr_only
                    if res["strata"].get("EUR", {}).get("exp_chi2", 0) >= 8.0
                    and res["strata"].get("EUR", {}).get("chi2", 99) < 3.84]
    print(f"\n=== AMR-restriction check ===")
    print(f"  cells where AMR is the ONLY significant WELL stratum : {len(amr_only)}/{len(rows)}")
    print(f"    ...of which EUR is WELL-POWERED yet NULL           : {len(eur_pow_null)}")
    ch16 = [h for h, _r in amr_only if h["sv_locus"].startswith("chr16_")]
    print(f"    ...on chr16 predictors                             : {len(ch16)}")

    print("\n=== strongest pooled cells: per-ancestry beta / p ===")
    print(f"  {'SV':<20}{'phecode':<10}" + "".join(f"{a:>17}" for a in WELL) + f"{'metaP':>9}{'I2':>4}  verdict")
    top = sorted(rows, key=lambda r: float(r[0].get("p_final") or 1.0))[:25]
    for h, cells, res in top:
        s = f"  {h['sv_locus'][:19]:<20}{h['phecode']:<10}"
        for a in WELL:
            st = res["strata"].get(a)
            cell = "na" if st is None else "{:+.2f}/{:.0e}".format(st["beta"], st["p"])
            s += f"{cell:>17}"
        s += f"{res['p']:>9.0e}" if np.isfinite(res['p']) else f"{'na':>9}"
        s += f"{res['I2']:>4.0f}" if np.isfinite(res['I2']) else f"{'na':>4}"
        s += f"  {res['verdict']}"
        print(s, flush=True)
    print("\nEXT2_ANCESTRY_META_DONE", flush=True)


if __name__ == "__main__":
    main()
