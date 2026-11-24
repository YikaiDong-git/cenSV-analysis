#!/usr/bin/env python
import numpy as np, os, csv
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
AGG = f"{SSIMP}/validation/ssimp_impg/work/phewide_agg"
BUND = f"{SSIMP}/output/bundles"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(f"{HERE}/../aou_handoff"); os.makedirs(OUT, exist_ok=True)
N_TOTAL = 50
CONTROLS = {
    "chr3_sv40_DEL_non_best6_alt_edge500": "POSITIVE_control_PROS1_VTE",
    "chr19_sv1949_INS_mon_best8_midpoint500": "NEGATIVE_control_cvr2~1_SNP_equivalent",
}


def istrue(x):
    return str(x).strip().upper() == "TRUE"


def eligible_pip(r):
    p = []
    try:
        if float(r["impute_pip"]) >= 0.5 and istrue(r["converged_imp"]) and float(r["region_n_sv_ge_robust_imp"]) <= 3:
            p.append((float(r["impute_pip"]), "impute"))
        if float(r["matched_pip"]) >= 0.5 and istrue(r["converged_mat"]) and float(r["region_n_sv_ge_robust_mat"]) <= 3:
            p.append((float(r["matched_pip"]), "matched"))
    except Exception:
        return None
    return max(p) if p else None


robust = {}
with open(f"{AGG}/candidate_audit_v2_owntag.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        robust[r["sv_locus"]] = dict(region=r["region"], study=r["study"], trait=r["rep_trait"],
                                     cv_r2=float(r["cv_r2"]), maxpip=float("nan"))

best = {}
with open(f"{AGG}/phewide_pip_ge_0.01.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ep = eligible_pip(r)
        if ep is None:
            continue
        sv = r["sv_locus"]
        if sv not in best or ep[0] > best[sv]["maxpip"]:
            best[sv] = dict(region=r["region_id"], study=r["study"], trait=r["pheno"],
                            cv_r2=float(r["cv_r2"]), maxpip=ep[0], handling=ep[1])

sel = list(robust.keys())
for c in CONTROLS:
    if c not in sel:
        sel.append(c)
for sv in sorted(best, key=lambda s: -best[s]["maxpip"]):
    if len(sel) >= N_TOTAL:
        break
    if sv not in sel:
        sel.append(sv)
print(f"selected {len(sel)} SVs = {len(robust)} robust-24 + extras to {N_TOTAL} (non-degenerate eligible)")

score_rows, meta, traitcount, cache = [], [], {}, {}
for sv in sel:
    info = best.get(sv) or robust.get(sv)
    region = info["region"]; chrom = region.split("_")[0].replace("chr", "")
    bpath = f"{BUND}/chr{chrom}/{region}_bundle.npz"
    if bpath not in cache:
        cache[bpath] = np.load(bpath, allow_pickle=True)
    d = cache[bpath]; svn = np.array([str(x) for x in d["sv_locus"]])
    note = CONTROLS.get(sv, "")
    if sv not in svn:
        meta.append(dict(sv_locus=sv, region=region, chrom=chrom, cv_r2=info["cv_r2"], n_weight_snps=0,
                         rep_study=info["study"], rep_trait=info["trait"], max_pip=info.get("maxpip"),
                         in_robust24=(sv in robust), note=(note + ";SV_NOT_IN_BUNDLE").strip(";")))
        continue
    ti = int(np.where(svn == sv)[0][0]); w = d["B"][:, ti]; sup = np.where(w != 0)[0]
    chr_, pos, ref, alt = d["snp_chrom"], d["snp_pos"], d["snp_ref"], d["snp_alt"]
    for j in sup:
        cc = str(chr_[j]); cc = cc if cc.startswith("chr") else f"chr{cc}"
        score_rows.append((sv, cc, int(pos[j]), str(ref[j]).upper(), str(alt[j]).upper(),
                           str(alt[j]).upper(), float(w[j])))
    meta.append(dict(sv_locus=sv, region=region, chrom=chrom, cv_r2=info["cv_r2"], n_weight_snps=int(len(sup)),
                     rep_study=info["study"], rep_trait=info["trait"], max_pip=info.get("maxpip"),
                     in_robust24=(sv in robust), note=note))
    traitcount[info["trait"]] = traitcount.get(info["trait"], 0) + 1


def atomic_write(path, header, rows):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if x is None else str(x) for x in r) + "\n")
    os.replace(tmp, path)


score_rows.sort(key=lambda r: (r[0], r[1], r[2]))
atomic_write(f"{OUT}/sv_weights_long.tsv", ["sv_locus", "chrom", "pos", "ref", "alt", "a1_effect", "beta_std"], score_rows)
mcols = ["sv_locus", "region", "chrom", "cv_r2", "n_weight_snps", "rep_study", "rep_trait", "max_pip", "in_robust24", "note"]
atomic_write(f"{OUT}/sv_metadata.tsv", mcols,
             [[m[c] for c in mcols] for m in sorted(meta, key=lambda m: (not m["in_robust24"], -(m["max_pip"] or 0)))])
atomic_write(f"{OUT}/trait_codes.tsv", ["rep_trait", "n_svs"], sorted(traitcount.items(), key=lambda kv: -kv[1]))

print(f"SVs written: {len(meta)} | in_robust24: {sum(m['in_robust24'] for m in meta)} | "
      f"weight-SNP rows: {len(score_rows)} | unique weight SNPs: {len({(r[1], r[2], r[3], r[4]) for r in score_rows})} | "
      f"distinct rep_traits: {len(traitcount)}")
print("rep_traits:", ", ".join(f"{t}({n})" for t, n in sorted(traitcount.items(), key=lambda kv: -kv[1])))
print("DONE_AOU_EXTRACT")
