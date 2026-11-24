#!/usr/bin/env python
import os, sys, gzip
import numpy as np, pandas as pd
from collections import defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
AUD = f"{P1}/analyses/cdr_calls"
CHAIN_DIR = f"{SUB}/raw/chains_chm13"
OUT_DIR = f"{WU}/output/mqtl/per_hap"
WINDOWS = [10_000, 50_000, 100_000]

TAG = sys.argv[1]
SAMPLE, HAP = TAG.rsplit("_", 1)
os.makedirs(OUT_DIR, exist_ok=True)


def parse_chain_chm13_to_asm(chain_path):
    chains = defaultdict(list)
    opener = gzip.open if chain_path.endswith(".gz") else open
    with opener(chain_path, "rt") as f:
        cur = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                cur = None
                continue
            if line.startswith("chain"):
                parts = line.split()
                cur = dict(t_name=parts[2], t_pos=int(parts[5]),
                           q_name=parts[7], q_size=int(parts[8]), q_strand=parts[9],
                           q_pos=int(parts[10]))
                continue
            if cur is None:
                continue
            parts = line.split()
            try:
                size = int(parts[0])
            except (ValueError, IndexError):
                continue
            t0 = cur["t_pos"]; t1 = t0 + size
            q0 = cur["q_pos"]; q1 = q0 + size
            if cur["q_strand"] == "+":
                cq0, cq1 = q0, q1
            else:
                cq0, cq1 = cur["q_size"] - q1, cur["q_size"] - q0
            chains[cur["q_name"]].append((cq0, cq1, cur["t_name"], t0, t1, cur["q_strand"]))
            if len(parts) >= 3:
                cur["t_pos"] = t1 + int(parts[1])
                cur["q_pos"] = q1 + int(parts[2])
            else:
                cur = None
    for k in chains:
        chains[k].sort(key=lambda x: x[0])
    return chains


def lift_chm13_to_asm(chains, chm13_chrom, lo, hi):
    if chm13_chrom not in chains:
        return []
    out = []
    for cq0, cq1, t_name, t0, t1, q_strand in chains[chm13_chrom]:
        if cq0 >= hi:
            break
        if cq1 <= lo:
            continue
        ovl0, ovl1 = max(cq0, lo), min(cq1, hi)
        d0, d1 = ovl0 - cq0, ovl1 - cq0
        size = cq1 - cq0
        if q_strand == "+":
            ta, tb = t0 + d0, t0 + d1
        else:
            ta, tb = t0 + (size - d1), t0 + (size - d0)
        if ta < tb:
            out.append((t_name, ta, tb))
    return out


MAN = pd.read_csv(f"{SUB}/manifest/methylation_manifest.tsv", sep="\t")
row = MAN[(MAN["sample"] == SAMPLE) & (MAN["hap_local"] == HAP)]
HAP_REMOTE = str(row["hap_remote"].iloc[0]) if len(row) else HAP
CANDS = [HAP_REMOTE, "mat" if HAP == "hap1" else "pat", HAP]
chain_path = None
for cand in dict.fromkeys(CANDS):
    p = f"{CHAIN_DIR}/{SAMPLE}_{cand}_vs_CHM13.chain.gz"
    if os.path.exists(p):
        chain_path, HAP_REMOTE = p, cand
        break
if chain_path is None:
    print(f"[SKIP] no chain for {TAG} (tried {list(dict.fromkeys(CANDS))})")
    sys.exit(0)
chains = parse_chain_chm13_to_asm(chain_path)
CHAIN_CTGS = {t for blocks in chains.values() for (_, _, t, _, _, _) in blocks}

bed = f"{AUD}/output/calls/{TAG}/binned_freq.bed"
if not os.path.exists(bed):
    sys.exit(f"[SKIP] no binned_freq for {TAG}")
M = pd.read_csv(bed, sep="\t", header=None,
                names=["ctg_full", "start", "end", "frac", "covg"])
M["ctg"] = M["ctg_full"].str.rsplit("#", n=1).str[-1]
BY = {c: g.sort_values("start").reset_index(drop=True) for c, g in M.groupby("ctg")}
ARR = {c: (g["start"].values, g["end"].values, g["frac"].values, g["covg"].values)
       for c, g in BY.items()}

vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
vt = vt[["sv_locus", "chrom", "pos", "end", "svtype", "geno_group"]].copy()

bed_ctgs = set(BY)
shared = bed_ctgs & CHAIN_CTGS
frac = len(shared) / max(1, len(bed_ctgs))
print(f"{TAG}: hap_remote={HAP_REMOTE} | chain CHM13 seqs {len(chains)} | methylation bins {len(M)} "
      f"on {len(BY)} contigs | SV loci {len(vt)}")
print(f"  contig-overlap guard: {len(shared)}/{len(bed_ctgs)} bed contigs present in the chain "
      f"({100*frac:.1f}%)")
assert frac >= 0.5, (f"HAPLOTYPE MISPAIRING: only {100*frac:.1f}% of {TAG}'s methylation contigs "
                     f"appear in chain {os.path.basename(chain_path)}")


def window_stats(chrom, lo, hi):
    ivs = lift_chm13_to_asm(chains, chrom, lo, hi)
    if not ivs:
        return np.nan, np.nan, 0, 0
    wsum = fsum = csum = 0.0
    nbin = 0
    for ctg, s, e in ivs:
        if ctg not in ARR:
            continue
        bs, be, bf, bc = ARR[ctg]
        i0, i1 = np.searchsorted(be, s, "right"), np.searchsorted(bs, e, "left")
        if i1 <= i0:
            continue
        ov = np.minimum(be[i0:i1], e) - np.maximum(bs[i0:i1], s)
        ok = ov > 0
        if not ok.any():
            continue
        w = ov[ok].astype(float)
        wsum += w.sum()
        fsum += float(np.sum(bf[i0:i1][ok] * w))
        csum += float(np.sum(bc[i0:i1][ok] * w))
        nbin += int(ok.sum())
    if wsum <= 0:
        return np.nan, np.nan, 0, len(ivs)
    return fsum / wsum, csum / wsum, nbin, len(ivs)


rows = []
for r in vt.itertuples():
    for W in WINDOWS:
        lo, hi = max(0, int(r.pos) - W), int(r.end) + W
        f, c, nbin, nblk = window_stats(r.chrom, lo, hi)
        rows.append((r.sv_locus, r.chrom, W, f, c, nbin, nblk))

D = pd.DataFrame(rows, columns=["sv_locus", "chrom", "window_bp", "meth", "covg", "n_bins", "n_blocks"])
D.insert(0, "tag", TAG)
out = f"{OUT_DIR}/{TAG}.tsv.gz"
D.to_csv(out, sep="\t", index=False, float_format="%.6g", compression="gzip")
cover = D.groupby("window_bp")["meth"].apply(lambda s: 100 * s.notna().mean())
print(f"  wrote {out}   rows {len(D)}")
print("  % SVs with methylation data, by window: " +
      "  ".join(f"{int(w)//1000}kb {v:.1f}%" for w, v in cover.items()))
