import csv, os
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]
P     = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
AGG   = f"{SSIMP}/validation/ssimp_impg/work/phewide_agg"
GD    = f"{SSIMP}/../sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix"

rows = [r for r in csv.DictReader(open(f"{AGG}/phewide_pip_ge_0.01.tsv"), delimiter="\t")
        if float(r["impute_pip"]) >= 0.5 and float(r["matched_pip"]) >= 0.5]
svs = {r["sv_locus"] for r in rows}
print(f"222-check: pairs={len(rows)} unique_SVs={len(svs)} unique_regions={len({r['region_id'] for r in rows})}")

traits = {}
for r in rows:
    k = (r["study"], r["pheno"]); traits[k] = traits.get(k, 0) + 1
have = sum(1 for (s, t) in traits if os.path.exists(f"{GD}/{s}/{t}.tsv.bgz"))
print(f"unique (study,trait) pairs among the 222 = {len(traits)} ; with a GWAS .bgz in GD = {have}")
print("study split:", {s: sum(1 for (ss, t) in traits if ss == s) for s in {k[0] for k in traits}})
print("\n-- unique traits (study  trait  n_SVs  gwas_bgz?) --")
for (s, t), n in sorted(traits.items(), key=lambda kv: (-kv[1], kv[0])):
    print(f"  {s:8s} {t:22s} nSV={n:2d} bgz={'Y' if os.path.exists(f'{GD}/{s}/{t}.tsv.bgz') else 'NO'}")

persv = {}
for r in rows:
    persv.setdefault(r["sv_locus"], set()).add(r["pheno"])
print("\nper-SV n_traits_hit (sorted):", sorted((len(v) for v in persv.values()), reverse=True))
print("SVs hitting >=3 traits:", sum(1 for v in persv.values() if len(v) >= 3))

print("\n-- 33 SVs: cv_r2 + region degeneracy --")
seen = set()
for r in sorted(rows, key=lambda r: -float(r["cv_r2"])):
    sv = r["sv_locus"]
    if sv in seen: continue
    seen.add(sv)
    print(f"  {sv:42s} cv_r2={float(r['cv_r2']):.3f} region_n_sv_imp={r['region_n_sv_ge_robust_imp']} mat={r['region_n_sv_ge_robust_mat']}")
print("DONE_RECON")
