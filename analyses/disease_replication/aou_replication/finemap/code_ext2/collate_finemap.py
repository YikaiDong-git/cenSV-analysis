#!/usr/bin/env python3
"""Collate per-fit SuSiE fine-mapping results into master + localization summary tables."""
import json
import os
import sys
import csv

FM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FITS = os.path.join(FM, "results_ext2", "fits")
OUT = os.path.join(FM, "results_ext2")
ENGINES = ("indiv", "rss_lpm", "rss_logit")


def _num(x):
    if x is None or x == "" or x == "NA" or x == "NaN":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _bool(x):
    if x is None:
        return None
    return str(x).strip().upper() in ("TRUE", "1", "T", "YES")


def read_kv_tsv(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if len(rows) < 2:
        return {}
    return dict(zip(rows[0], rows[1]))


def read_cs_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def collate():
    master, local = [], []
    n_skip_no_success = 0
    if not os.path.isdir(FITS):
        sys.exit(f"no fits dir: {FITS}")
    for stratum in sorted(os.listdir(FITS)):
        sdir = os.path.join(FITS, stratum)
        if not os.path.isdir(sdir):
            continue
        for fit_id in sorted(os.listdir(sdir)):
            fdir = os.path.join(sdir, fit_id)
            if not os.path.isdir(fdir):
                continue
            for mode in sorted(os.listdir(fdir)):
                mdir = os.path.join(fdir, mode)
                if not os.path.isdir(mdir):
                    continue
                if not os.path.exists(os.path.join(mdir, "_SUCCESS")):
                    if os.path.exists(os.path.join(mdir, "_SKIP")):
                        continue
                    n_skip_no_success += 1
                    continue
                meta = {}
                mp = os.path.join(mdir, "fit_meta.json")
                if os.path.exists(mp):
                    with open(mp) as fh:
                        meta = json.load(fh)
                summ = read_kv_tsv(os.path.join(mdir, "susie_fit_summary.tsv"))
                cs = read_cs_tsv(os.path.join(mdir, "susie_credible_sets.tsv"))

                sv_cs_by_eng = {e: [r for r in cs if r.get("engine") == e and _num(r.get("n_sv")) and _num(r.get("n_sv")) > 0] for e in ENGINES}
                sv_in_any_cs = any(sv_cs_by_eng[e] for e in ENGINES)
                sv_top_cs = [r for r in cs if r.get("top_type") == "sv"]
                best_sv_cs_pip = max((_num(r.get("top_pip")) or 0.0 for r in sv_top_cs), default=0.0)
                best_sv_cs_variant = ""
                if sv_top_cs:
                    best = max(sv_top_cs, key=lambda r: _num(r.get("top_pip")) or 0.0)
                    best_sv_cs_variant = best.get("top_variant", "")

                sum_sv_pips = [_num(summ.get(f"sum_sv_pip_{e}")) for e in ENGINES]
                top_sv_pips = [_num(summ.get(f"top_sv_pip_{e}")) for e in ("indiv", "rss_logit")]
                best_sum_sv_pip = max((v for v in sum_sv_pips if v is not None), default=None)
                best_top_sv_pip = max((v for v in top_sv_pips if v is not None), default=None)

                n_cs_any = max((_num(summ.get(f"n_cs_{e}")) or 0 for e in ENGINES), default=0)
                if sv_in_any_cs:
                    verdict = "SV_in_CS"
                elif n_cs_any > 0:
                    verdict = "tag_SNP_only"
                else:
                    verdict = "unresolved"

                row = {
                    "stratum": stratum,
                    "chrom": meta.get("chrom", ""),
                    "phecode": meta.get("phecode", ""),
                    "mcc": meta.get("mcc", ""),
                    "mode": mode,
                    "fit_id": fit_id,
                    "n_variants": meta.get("n_variants"),
                    "n_sv": meta.get("n_sv"),
                    "n_snp": meta.get("n_snp"),
                    "n_present": meta.get("n_present"),
                    "ncase": meta.get("ncase"),
                    "ncontrol": meta.get("ncontrol"),
                    "n_eff": meta.get("n_eff"),
                    "null_converged": meta.get("null_converged"),
                    "n_snp_scanned": meta.get("n_snp_scanned"),
                    "n_snp_below_callrate": meta.get("n_snp_below_callrate"),
                    "hit_sv_ids": "|".join(meta.get("hit_sv_ids", []) or []),
                    "max_abs_sv_snp_r": _num(summ.get("max_abs_sv_snp_r")),
                    "L": summ.get("L"),
                }
                for e in ENGINES:
                    row[f"n_cs_{e}"] = summ.get(f"n_cs_{e}")
                    row[f"converged_{e}"] = _bool(summ.get(f"converged_{e}"))
                    row[f"sum_sv_pip_{e}"] = _num(summ.get(f"sum_sv_pip_{e}"))
                    row[f"value_{e}"] = _num(summ.get(f"value_{e}"))
                row["jaccard_indiv_vs_rsslpm"] = _num(summ.get("jaccard_indiv_vs_rsslpm"))
                row["pipcor_indiv_vs_rsslpm"] = _num(summ.get("pipcor_indiv_vs_rsslpm"))
                row["jaccard_indiv_vs_rsslogit"] = _num(summ.get("jaccard_indiv_vs_rsslogit"))
                row["pipcor_indiv_vs_rsslogit"] = _num(summ.get("pipcor_indiv_vs_rsslogit"))
                row["sv_in_any_cs"] = sv_in_any_cs
                row["best_sum_sv_pip"] = best_sum_sv_pip
                row["best_top_sv_pip"] = best_top_sv_pip
                row["best_sv_cs_pip"] = best_sv_cs_pip
                row["best_sv_cs_variant"] = best_sv_cs_variant
                row["verdict"] = verdict
                master.append(row)

                local.append({
                    "stratum": stratum, "chrom": meta.get("chrom", ""),
                    "phecode": meta.get("phecode", ""), "mode": mode,
                    "ncase": meta.get("ncase"), "n_present": meta.get("n_present"),
                    "n_sv": meta.get("n_sv"), "hit_sv_ids": "|".join(meta.get("hit_sv_ids", []) or []),
                    "best_sum_sv_pip": best_sum_sv_pip, "best_top_sv_pip": best_top_sv_pip,
                    "sv_in_any_cs": sv_in_any_cs, "best_sv_cs_pip": best_sv_cs_pip,
                    "best_sv_cs_variant": best_sv_cs_variant,
                    "max_abs_sv_snp_r": _num(summ.get("max_abs_sv_snp_r")),
                    "n_cs_indiv": summ.get("n_cs_indiv"),
                    "pipcor_indiv_vs_rsslogit": _num(summ.get("pipcor_indiv_vs_rsslogit")),
                    "verdict": verdict,
                })
    return master, local, n_skip_no_success


def write_tsv(rows, path, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    master, local, n_missing = collate()
    if not master:
        sys.exit("no completed fits found")
    mcols = list(master[0].keys())
    lcols = list(local[0].keys())
    write_tsv(master, os.path.join(OUT, "finemap_master.tsv"), mcols)
    write_tsv(local, os.path.join(OUT, "finemap_localization.tsv"), lcols)

    n = len(master)
    by_verdict = {}
    for r in master:
        by_verdict[r["verdict"]] = by_verdict.get(r["verdict"], 0) + 1
    sv_cs = sorted([r for r in master if r["sv_in_any_cs"]],
                   key=lambda r: -((r["best_top_sv_pip"] or 0) or (r["best_sv_cs_pip"] or 0)))
    print(f"collated {n} fit-modes ({n_missing} with neither _SUCCESS nor _SKIP)")
    print("  by stratum: " + ", ".join(f"{s}={sum(1 for r in master if r['stratum']==s)}"
                                       for s in sorted({r['stratum'] for r in master})))
    print("  by mode: " + ", ".join(f"{m}={sum(1 for r in master if r['mode']==m)}"
                                     for m in sorted({r['mode'] for r in master})))
    print(f"  verdicts: {by_verdict}")
    print(f"  *** {len(sv_cs)} fit-modes with an SV IN a credible set (clean SV-causal signal) ***")
    for r in sv_cs[:30]:
        pip = (r["best_top_sv_pip"] or 0) or (r["best_sv_cs_pip"] or 0)
        print(f"    {r['stratum']:8s} {r['chrom']:6s} {r['phecode']:10s} {r['mode']:10s} "
              f"ncase={r['ncase']} top_sv_pip={pip:.3f} sv_cs_var={r['best_sv_cs_variant']}")
    print(f"\nwrote {OUT}/finemap_master.tsv ({n} rows), {OUT}/finemap_localization.tsv")


if __name__ == "__main__":
    main()
