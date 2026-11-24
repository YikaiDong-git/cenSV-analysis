#!/usr/bin/env python3
import os, json, time, urllib.request
import pandas as pd
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

SV_PARAM = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
CHROMS = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
BD = CENSV_ROOT + ""
BASE = f"{BD}/results"
OUT_DIR = f"{BD}/results/analysis_qtl/final_tables"
os.makedirs(OUT_DIR, exist_ok=True)
MAGE_PERICEN_FILE = f"{BD}/claude_t/results_v9_combined/mage_pericen_eqtl.tsv"
STRATS = ['edge500','edge1000','edge2000','fksize100k','fksize200k','fksize500k','midpoint500','midpoint1000']

cen = {}
with open(f"{BD}/resources/grch38_centromeres.bed") as f:
    for line in f:
        fs = line.strip().split('\t')
        cen[fs[0]] = (int(fs[1]), int(fs[2]))

mage_lookup = {}
try:
    if os.path.exists(MAGE_PERICEN_FILE) and os.path.getsize(MAGE_PERICEN_FILE) > 0:
        mage = pd.read_csv(MAGE_PERICEN_FILE, sep='\t')
        if 'ensemblID' in mage.columns and len(mage) > 0:
            for gid, grp in mage.groupby("ensemblID"):
                base = str(gid).split('.')[0]
                mage_lookup[base] = (len(grp), grp["pval_nominal"].min(),
                                     str(grp["geneSymbol"].iloc[0]) if "geneSymbol" in grp.columns else "")
except Exception as e:
    print(f"  MAGE pericen unavailable ({e}); skipping")
print(f"MAGE: {len(mage_lookup)} genes")

sv_meta = {}
for strat in STRATS:
    for chrom in CHROMS:
        mf = f"{BASE}/elastic_net/decompose/cosine_nj/{strat}_{SV_PARAM}/{chrom}/sv_metrics.tsv"
        if not os.path.exists(mf):
            continue
        try:
            lines = open(mf).readlines()
            hdr = lines[0].strip().split('\t')
            cols = {c: i for i, c in enumerate(hdr)}
            for l in lines[1:]:
                fs = l.strip().split('\t')
                try:
                    si = int(fs[cols['sv_index']])
                    sv_meta[(strat, chrom, si)] = {
                        'svtype': fs[cols.get('svtype', -1)] if 'svtype' in cols else '',
                        'svlen': fs[cols.get('svlen', -1)] if 'svlen' in cols else '',
                        'svpos': fs[cols.get('svpos', -1)] if 'svpos' in cols else '',
                        'svend': fs[cols.get('svend', -1)] if 'svend' in cols else '',
                        'svfreq': fs[cols.get('svfreq', -1)] if 'svfreq' in cols else '',
                    }
                except:
                    pass
        except:
            pass
print(f"SV metadata: {len(sv_meta)} entries")

sel_data = {}
for ens in ['best6', 'best8']:
    for chrom in CHROMS:
        sf = f"{BASE}/imputation/decompose/cosine_nj/{ens}_{SV_PARAM}/{chrom}/{chrom}_best_strategy_selection.tsv"
        if not os.path.exists(sf):
            continue
        df = pd.read_csv(sf, sep='\t')
        for _, row in df.iterrows():
            dcr2 = row.get('dosage_cv_r2', '')
            if pd.notna(dcr2) and dcr2 != '':
                try:
                    si = int(row['sv_index'])
                    sel_data[(ens, chrom, si)] = {
                        'cv_r2': float(dcr2),
                        'strategy': str(row.get('dosage_strategy', '')),
                    }
                except:
                    pass
print(f"Selection data: {len(sel_data)} entries")

all_signals = []

for ens in ['best6', 'best8']:
    for qtl_type in ['eqtl', 'sqtl']:
        peer = 'peer60' if qtl_type == 'eqtl' else 'peer15'

        for window, wlabel in [('10000000', 'w10M'), ('0', 'w0')]:
            for chrom in CHROMS:
                pf = f"{BASE}/{qtl_type}_overlap/decompose/cosine_nj/{ens}_{SV_PARAM}/{chrom}/sexTrue_gpc5_{peer}/window_{window}/{chrom}_cis_permutation.tsv"
                if not os.path.exists(pf):
                    continue
                df = pd.read_csv(pf, sep='\t')
                if 'qval' not in df.columns:
                    continue
                sig = df[df['qval'] < 0.10]
                for _, row in sig.iterrows():
                    gene = row['phenotype_id']
                    if ':' in gene and 'ENSG' in gene:
                        gene_base = gene.split(':')[-1].split('.')[0]
                    else:
                        gene_base = gene.split('.')[0]
                    try:
                        si = int(row['variant_id'].split('_')[1])
                    except:
                        si = -1
                    sinfo = sel_data.get((ens, chrom, si), {})
                    strat = sinfo.get('strategy', '')
                    smeta = sv_meta.get((strat, chrom, si), {}) if strat else {}

                    all_signals.append({
                        'ens': ens, 'qtl': qtl_type, 'window': wlabel, 'cis_trans': 'cis',
                        'sv_chrom': chrom, 'gene_base': gene_base, 'gene_full': gene,
                        'sv_id': row['variant_id'], 'sv_idx': si,
                        'cv_r2': sinfo.get('cv_r2', 0), 'qval': row['qval'],
                        'slope': row.get('slope', 0),
                        'strategy': strat,
                        'svtype': smeta.get('svtype', ''), 'svlen': smeta.get('svlen', ''),
                        'svpos': smeta.get('svpos', ''), 'svend': smeta.get('svend', ''),
                        'svfreq': smeta.get('svfreq', ''),
                    })

        for chrom in CHROMS:
            tf = f"{BASE}/{qtl_type}_overlap/decompose/cosine_nj/{ens}_{SV_PARAM}/{chrom}/sexTrue_gpc5_{peer}/trans/{chrom}_trans_results.tsv"
            if not os.path.exists(tf):
                continue
            df = pd.read_csv(tf, sep='\t')
            if 'pval' not in df.columns:
                continue
            sig = df[df['pval'] < 5e-8]
            for _, row in sig.iterrows():
                gene = row['phenotype_id']
                if ':' in gene and 'ENSG' in gene:
                    gene_base = gene.split(':')[-1].split('.')[0]
                else:
                    gene_base = gene.split('.')[0]
                try:
                    si = int(row['variant_id'].split('_')[1])
                except:
                    si = -1
                sinfo = sel_data.get((ens, chrom, si), {})
                strat = sinfo.get('strategy', '')
                smeta = sv_meta.get((strat, chrom, si), {}) if strat else {}

                all_signals.append({
                    'ens': ens, 'qtl': qtl_type, 'window': 'trans', 'cis_trans': 'trans',
                    'sv_chrom': chrom, 'gene_base': gene_base, 'gene_full': gene,
                    'sv_id': row['variant_id'], 'sv_idx': si,
                    'cv_r2': sinfo.get('cv_r2', 0), 'qval': row['pval'],
                    'slope': row.get('b', 0),
                    'strategy': strat,
                    'svtype': smeta.get('svtype', ''), 'svlen': smeta.get('svlen', ''),
                    'svpos': smeta.get('svpos', ''), 'svend': smeta.get('svend', ''),
                    'svfreq': smeta.get('svfreq', ''),
                })

sdf = pd.DataFrame(all_signals)
print(f"\nTotal signals: {len(sdf)} (cis: {(sdf['cis_trans']=='cis').sum()}, trans: {(sdf['cis_trans']=='trans').sum()})")
unique_genes = sorted(sdf['gene_base'].unique())
print(f"Unique genes: {len(unique_genes)}")

def query_ensembl(gene_id):
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=0;content-type=application/json"
    try:
        req = urllib.request.Request(url)
        req.add_header('Content-Type', 'application/json')
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            return {
                'symbol': data.get('display_name', ''),
                'biotype': data.get('biotype', ''),
                'description': (data.get('description', '') or '').split(' [Source:')[0][:50],
                'gene_chr': str(data.get('seq_region_name', '')),
                'gene_start': data.get('start', ''),
                'gene_end': data.get('end', ''),
            }
    except:
        return {}

CACHE = f"{OUT_DIR}/ensembl_cache.json"
gene_info = {}
if os.path.exists(CACHE):
    try: gene_info = json.load(open(CACHE))
    except: pass

print(f"\nQuerying Ensembl for {len(unique_genes)} genes (cached={sum(1 for g in unique_genes if g in gene_info)})...")
for i, g in enumerate(unique_genes):
    if i > 0 and i % 50 == 0:
        print(f"  {i}/{len(unique_genes)}...")
        json.dump(gene_info, open(CACHE, 'w'))
    if g not in gene_info:
        gene_info[g] = query_ensembl(g)
        time.sleep(0.1)
json.dump(gene_info, open(CACHE, 'w'))
print(f"  Done. Got symbols for {sum(1 for v in gene_info.values() if v.get('symbol'))} genes")

ARTIFACT_FAMILIES = ['ANKRD20', 'DUX4', 'FAM95', 'AGAP', 'GTF2I', 'LSP1P', 'FRG1',
    'HERC2P', 'IGKV', 'IGHV', 'ZNF33', 'ZNF37', 'ZNF514', 'HSD17B7P',
    'FAM182', 'FAM88', 'FAM27', 'ANKRD36', 'CDRT', 'GPAT2P', 'FAHD2C',
    'PLD5P', 'WHAMM', 'ROCK1P']

def classify_risk(info):
    biotype = info.get('biotype', '')
    sym = info.get('symbol', '')
    if 'pseudogene' in biotype:
        return 'high'
    if 'IG_' in biotype or 'TR_' in biotype:
        return 'high'
    for fam in ARTIFACT_FAMILIES:
        if fam in sym:
            return 'high'
    if biotype == 'lncRNA':
        return 'moderate'
    if biotype == 'protein_coding':
        return 'low'
    if not sym or not biotype:
        return '?'
    return 'moderate'

rows = []
for g in unique_genes:
    gsigs = sdf[sdf['gene_base'] == g]
    info = gene_info.get(g, {})

    cis_sigs = gsigs[gsigs['cis_trans'] == 'cis']
    trans_sigs = gsigs[gsigs['cis_trans'] == 'trans']

    if len(cis_sigs) > 0:
        best_idx = cis_sigs['cv_r2'].idxmax()
        best = cis_sigs.loc[best_idx]
    elif len(trans_sigs) > 0:
        best_idx = trans_sigs['cv_r2'].idxmax()
        best = trans_sigs.loc[best_idx]
    else:
        continue

    n_cis = len(cis_sigs)
    n_trans = len(trans_sigs)

    gene_chr = info.get('gene_chr', '')
    gene_start = info.get('gene_start', '')

    is_trans_only = n_cis == 0
    if is_trans_only:
        dist_str = f"trans_{best['sv_chrom']}"
        wing = 'T'
    elif gene_chr and gene_start and str(gene_chr) != 'nan':
        c = 'chr' + gene_chr if not gene_chr.startswith('chr') else gene_chr
        if c in cen and gene_start:
            cs, ce = cen[c]
            tss = int(gene_start)
            if tss < cs:
                dist = cs - tss; wing = 'L'
            elif tss > ce:
                dist = tss - ce; wing = 'R'
            else:
                dist = 0; wing = 'IN'
            dist_str = f"{dist/1e6:.2f}"
        else:
            dist_str = '-'; wing = '-'
    else:
        dist_str = '-'; wing = '-'

    sym = info.get('symbol', '') or g
    risk = classify_risk(info)

    qtl_types = sorted(gsigs['qtl'].unique())
    ensembles = sorted(gsigs['ens'].unique())
    windows = sorted(gsigs['window'].unique())

    if n_cis > 0 and n_trans > 0:
        ct = 'cis,trans'
    elif n_trans > 0:
        ct = 'trans'
    else:
        ct = 'cis'

    minfo = mage_lookup.get(g)

    rows.append({
        'ENSG': g, 'Symbol': sym, 'gene_chr': gene_chr if gene_chr else '-',
        'dist_Mb': dist_str, 'Wing': wing, 'risk': risk,
        'MAGE': 'Y' if minfo else 'N',
        'MAGE_n': int(minfo[0]) if minfo else 0,
        'MAGE_p': minfo[1] if minfo else 1.0,
        'QTL': ','.join(qtl_types),
        'cis_trans': ct,
        'ens': ','.join(ensembles),
        'windows': ','.join(windows),
        'n_cis': n_cis, 'n_trans': n_trans,
        'cv_r2': best['cv_r2'], 'qval': best['qval'],
        'sv_chrom': best['sv_chrom'], 'sv_id': best['sv_id'],
        'svtype': best.get('svtype', ''), 'svlen': best.get('svlen', ''),
        'svpos': best.get('svpos', ''), 'svend': best.get('svend', ''),
        'svfreq': best.get('svfreq', ''),
        'strategy': best['strategy'],
        'desc': info.get('description', ''),
    })

rdf = pd.DataFrame(rows).sort_values('cv_r2', ascending=False).reset_index(drop=True)

outf = f"{OUT_DIR}/all_genes_v3.tsv"
rdf.to_csv(outf, sep='\t', index=False)
print(f"\nSaved {len(rdf)} genes to {outf}")

print()
print("=" * 80); print("SUMMARY"); print("=" * 80)
print(f"Total genes: {len(rdf)}")
for risk in ['low', 'moderate', 'high', '?']:
    sub = rdf[rdf['risk'] == risk]
    if len(sub) == 0: continue
    print(f"  {risk:<10s}: {len(sub):>4d} | MAGE: {(sub['MAGE']=='Y').sum()}/{len(sub)}")
print()
print(f"  cis only:    {(rdf['cis_trans']=='cis').sum()}")
print(f"  trans only:  {(rdf['cis_trans']=='trans').sum()}")
print(f"  cis+trans:   {(rdf['cis_trans']=='cis,trans').sum()}")
print()
print(f"  eQTL only:   {(rdf['QTL']=='eqtl').sum()}")
print(f"  sQTL only:   {(rdf['QTL']=='sqtl').sum()}")
print(f"  both:        {rdf['QTL'].str.contains(',').sum()}")
