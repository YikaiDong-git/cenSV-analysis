#!/usr/bin/env python3
"""Merge gene annotations from multiple sources into final table."""
import os, json, time, urllib.request, gzip
import pandas as pd
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BD = CENSV_ROOT + ""
SV_PARAM = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
CHROMS = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
BASE = f"{BD}/results/analysis_qtl/final_tables"
STRATS = ['edge500','edge1000','edge2000','fksize100k','fksize200k','fksize500k','midpoint500','midpoint1000']

v3 = pd.read_csv(f"{BASE}/all_genes_v3.tsv", sep='\t')
v3.columns = v3.columns.str.strip()
for col in v3.select_dtypes(include='object').columns:
    v3[col] = v3[col].str.strip()
print(f"Loaded v3 table: {len(v3)} genes")

def query_ensembl(gene_id):
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=0;content-type=application/json"
    try:
        req = urllib.request.Request(url)
        req.add_header('Content-Type', 'application/json')
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            return {
                'api_symbol': data.get('display_name', ''),
                'api_biotype': data.get('biotype', ''),
                'api_desc': (data.get('description', '') or '').split(' [Source:')[0][:60],
                'api_chr': str(data.get('seq_region_name', '')),
                'api_start': data.get('start', ''),
                'api_end': data.get('end', ''),
            }
    except Exception as e:
        return {'api_symbol': '', 'api_biotype': '', 'api_desc': str(e)[:40], 'api_chr': '', 'api_start': '', 'api_end': ''}

unique_genes = v3['ENSG'].unique()
print(f"\nLoading Ensembl info for {len(unique_genes)} genes (cache-first)...")
api_data = {}
_cache_path = f"{BASE}/ensembl_cache.json"
_cache = {}
if os.path.exists(_cache_path):
    try: _cache = json.load(open(_cache_path))
    except: pass
for i, g in enumerate(unique_genes):
    if g in _cache and _cache[g].get('symbol'):
        c = _cache[g]
        api_data[g] = {
            'api_symbol': c.get('symbol',''),
            'api_biotype': c.get('biotype',''),
            'api_desc': c.get('description', c.get('desc','')),
            'api_chr': c.get('gene_chr', c.get('chr','')),
            'api_start': c.get('gene_start', c.get('start','')),
            'api_end': c.get('gene_end', c.get('end','')),
        }
    else:
        if i > 0 and i % 50 == 0: print(f"  {i}/{len(unique_genes)}...")
        api_data[g] = query_ensembl(g)
        time.sleep(0.1)
print(f"  Got symbols for {sum(1 for v in api_data.values() if v.get('api_symbol'))} genes")

print("\nScanning full MAGE eQTL...")
gene_set = set(unique_genes)
full_mage = {}

_full_mage_path = f"{BD}/resources/MAGE/cis_eqtl_nominal.tsv.gz"
if not os.path.exists(_full_mage_path):
    print(f"  Full MAGE scan skipped (file not present)")
else:
  with gzip.open(_full_mage_path, 'rt') as f:
    header = f.readline().strip().split('\t')
    cols = {c: i for i, c in enumerate(header)}
    gene_col = cols.get('ensemblID', cols.get('phenotype_id', 0))
    pval_col = cols.get('pval_nominal', 3)
    sym_col = cols.get('geneSymbol', None)
    n = 0
    for line in f:
        n += 1
        if n % 5000000 == 0:
            print(f"  {n/1e6:.0f}M lines...")
        fs = line.strip().split('\t')
        if len(fs) <= gene_col:
            continue
        gid = fs[gene_col].split('.')[0]
        if gid not in gene_set:
            continue
        try:
            p = float(fs[pval_col])
        except:
            p = 1.0
        sym = fs[sym_col] if sym_col is not None and len(fs) > sym_col else ''
        if gid not in full_mage:
            full_mage[gid] = [0, 1.0, sym]
        full_mage[gid][0] += 1
        if p < full_mage[gid][1]:
            full_mage[gid][1] = p
            if sym:
                full_mage[gid][2] = sym

print(f"  Full MAGE: {len(full_mage)} of our genes found ({locals().get('n',0)} lines scanned)")

pericen_mage = {}
_mp_path = f"{BASE}/mage_pericen_eqtl.tsv"
import os as _os
if _os.path.exists(_mp_path) and _os.path.getsize(_mp_path) > 100:
    try:
        mage_peri = pd.read_csv(_mp_path, sep='\t')
        if 'ensemblID' in mage_peri.columns:
            for gid, grp in mage_peri.groupby("ensemblID"):
                base = str(gid).split('.')[0]
                pericen_mage[base] = (len(grp), grp["pval_nominal"].min(),
                                      str(grp["geneSymbol"].iloc[0]) if "geneSymbol" in grp.columns else "")
    except Exception as _e:
        print(f"  MAGE pericen unavailable: {_e}")

WEB_CORRECTIONS = {
"ENSG00000281741": {"symbol": "ENSG00000281741", "biotype": "TEC", "risk": "high", "note": "TEC - unannotated locus"},    "ENSG00000280043": {"symbol": "ENSG00000280043", "biotype": "TEC", "risk": "high", "note": "TEC - unannotated locus"},    "ENSG00000279425": {"symbol": "ENSG00000279425", "biotype": "TEC", "risk": "high", "note": "TEC - unannotated locus"},    "ENSG00000281571": {"symbol": "ENSG00000281571", "biotype": "TEC", "risk": "high", "note": "TEC - unannotated locus"},    "ENSG00000279352": {"symbol": "ENSG00000279352", "biotype": "TEC", "risk": "high", "note": "TEC - unannotated locus"},    "ENSG00000283352": {"symbol": "ENSG00000283352", "biotype": "TEC", "risk": "high", "note": "TEC/retired - unannotated"},    "ENSG00000280296": {"symbol": "ENSG00000280296", "biotype": "TEC", "risk": "high", "note": "TEC - unannotated locus"},    "ENSG00000272316": {"symbol": "ENSG00000272316", "biotype": "TEC", "risk": "high", "note": "TEC/retired - unannotated"},
    'ENSG00000276805': {'web_symbol': 'GTF2IP', 'web_biotype': 'transcribed_unprocessed_pseudogene', 'web_risk': 'high', 'web_note': 'GTF2I pseudogene on chr10 segdup'},
    'ENSG00000285454': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'moderate', 'web_note': 'Novel lncRNA, no HGNC symbol'},
    'ENSG00000250090': {'web_symbol': '(novel)', 'web_biotype': 'transcribed_unprocessed_pseudogene', 'web_risk': 'high', 'web_note': 'Novel ZNF pseudogene'},
    'ENSG00000214265': {'web_symbol': '(novel)', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Novel protein in 15q11 segdup'},
    'ENSG00000273154': {'web_symbol': 'ZGPAT-LIME1', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Readthrough transcript'},
    'ENSG00000283930': {'web_symbol': 'PLD5P1', 'web_biotype': 'protein_coding', 'web_risk': 'high', 'web_note': 'HGNC says pseudogene despite Ensembl protein_coding'},
    'ENSG00000215190': {'web_symbol': 'LINC00680', 'web_biotype': 'transcribed_unprocessed_pseudogene', 'web_risk': 'high', 'web_note': 'Biotype is pseudogene not lncRNA'},
    'ENSG00000283352': {'web_symbol': '(retired)', 'web_biotype': 'NOT_FOUND', 'web_risk': 'high', 'web_note': 'Retired/merged Ensembl ID'},
    'ENSG00000279352': {'web_symbol': '(TEC)', 'web_biotype': 'TEC', 'web_risk': 'high', 'web_note': 'To be experimentally confirmed'},
    'ENSG00000272316': {'web_symbol': 'LINC00680', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'Retired ID, replaced by ENSG00000290555'},
    'ENSG00000261247': {'web_symbol': 'GOLGA8T', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'chr15 segdup region'},
    'ENSG00000276085': {'web_symbol': 'CCL3L3', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Copy-number variable'},
    'ENSG00000162825': {'web_symbol': 'NBPF20', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'NBPF segdup family'},
    'ENSG00000135976': {'web_symbol': 'ANKRD36', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Segdup-prone repeat family'},
    'ENSG00000196912': {'web_symbol': 'ANKRD36B', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Segdup-prone repeat family'},
    'ENSG00000174501': {'web_symbol': 'ANKRD36C', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Segdup-prone repeat family'},
    'ENSG00000178700': {'web_symbol': 'DHFR2', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'Paralog of DHFR'},
    'ENSG00000272092': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000258924': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'Readthrough between TMEM199 and SARM1'},
    'ENSG00000287908': {'web_symbol': '(novel)', 'web_biotype': 'protein_coding', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000264668': {'web_symbol': '(novel)', 'web_biotype': 'protein_coding', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000267257': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'Antisense to ALPK2'},
    'ENSG00000284116': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000225411': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000287265': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'Antisense to DUSP22'},
    'ENSG00000287300': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000230010': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000227486': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'moderate', 'web_note': 'Novel transcript on chrX'},
    'ENSG00000235893': {'web_symbol': '(novel)', 'web_biotype': 'processed_pseudogene', 'web_risk': 'high', 'web_note': 'No HGNC symbol'},
    'ENSG00000286837': {'web_symbol': '(novel)', 'web_biotype': 'lncRNA', 'web_risk': 'moderate', 'web_note': 'Antisense to THAP1'},
    'ENSG00000273047': {'web_symbol': '(novel)', 'web_biotype': 'protein_coding', 'web_risk': 'moderate', 'web_note': 'LIME1-SLC2A4RG readthrough'},
    'ENSG00000204556': {'web_symbol': '(novel)', 'web_biotype': 'transcribed_unprocessed_pseudogene', 'web_risk': 'high', 'web_note': 'ZNF pseudogene'},
    'ENSG00000283196': {'web_symbol': '(novel)', 'web_biotype': 'transcribed_unprocessed_pseudogene', 'web_risk': 'high', 'web_note': 'LSP1 pseudogene'},
    'ENSG00000281179': {'web_symbol': '(novel)', 'web_biotype': 'IG_V_gene', 'web_risk': 'high', 'web_note': 'Identical to IGHV1OR15-1'},
    'ENSG00000281550': {'web_symbol': '(novel)', 'web_biotype': 'IG_V_pseudogene', 'web_risk': 'high', 'web_note': 'Identical to IGHV1OR15-2'},
}

ARTIFACT_FAMILIES = ['ANKRD20', 'DUX4', 'FAM95', 'AGAP', 'GTF2I', 'LSP1P', 'FRG1',
    'HERC2P', 'IGKV', 'IGHV', 'ZNF33', 'ZNF37', 'HSD17B7P',
    'FAM182', 'FAM88', 'FAM27', 'CDRT', 'GPAT2P', 'FAHD2C',
    'PLD5P', 'WHAMM', 'ROCK1P', 'GUSBP', 'FGF7P', 'LGALS9D']

def auto_risk(symbol, biotype):
    if 'pseudogene' in biotype:
        return 'high'
    if 'IG_' in biotype or 'TR_' in biotype:
        return 'high'
    if biotype == 'TEC' or biotype == 'NOT_FOUND':
        return 'high'
    for fam in ARTIFACT_FAMILIES:
        if fam in symbol:
            return 'high'
    if biotype == 'lncRNA':
        return 'moderate'
    if biotype == 'protein_coding':
        return 'low'
    return 'moderate'

rows = []
for _, r in v3.iterrows():
    g = r['ENSG']
    api = api_data.get(g, {})
    web = WEB_CORRECTIONS.get(g, {})
    fm = full_mage.get(g)
    pm = pericen_mage.get(g)

    final_sym = web.get('web_symbol', '') or api.get('api_symbol', '') or (fm[2] if fm else '') or r.get('Symbol', g)
    if not final_sym or final_sym == g:
        final_sym = g

    # Final biotype: web > api
    final_biotype = web.get('web_biotype', '') or api.get('api_biotype', '')

    if web.get('web_risk'):
        final_risk = web['web_risk']
    elif final_biotype:
        final_risk = auto_risk(final_sym, final_biotype)
    else:
        final_risk = '?'

    final_chr = api.get('api_chr', '') or str(r.get('gene_chr', ''))

    mage_pericen_yn = 'Y' if pm else 'N'
    mage_full_yn = 'Y' if fm else 'N'
    mage_full_n = fm[0] if fm else 0
    mage_full_p = fm[1] if fm else 1.0
    mage_pericen_n = pm[0] if pm else 0
    mage_pericen_p = pm[1] if pm else 1.0

    rows.append({
        'ENSG': g,
        'final_symbol': final_sym,
        'final_biotype': final_biotype,
        'final_risk': final_risk,
        'gene_chr': final_chr,
        'dist_Mb': r.get('dist_Mb', ''),
        'Wing': r.get('Wing', ''),
        'api_symbol': api.get('api_symbol', ''),
        'api_biotype': api.get('api_biotype', ''),
        'api_desc': api.get('api_desc', ''),
        'web_symbol': web.get('web_symbol', ''),
        'web_biotype': web.get('web_biotype', ''),
        'web_risk': web.get('web_risk', ''),
        'web_note': web.get('web_note', ''),
        'mage_pericen': mage_pericen_yn,
        'mage_pericen_n': mage_pericen_n,
        'mage_pericen_p': mage_pericen_p,
        'mage_full': mage_full_yn,
        'mage_full_n': mage_full_n,
        'mage_full_p': mage_full_p,
        'QTL': r.get('QTL', ''),
        'cis_trans': r.get('cis_trans', ''),
        'ens': r.get('ens', ''),
        'windows': r.get('windows', ''),
        'n_cis': r.get('n_cis', 0),
        'n_trans': r.get('n_trans', 0),
        'cv_r2': r.get('cv_r2', 0),
        'qval': r.get('qval', 1),
        'sv_chrom': r.get('sv_chrom', ''),
        'sv_id': r.get('sv_id', ''),
        'svtype': r.get('svtype', ''),
        'svlen': r.get('svlen', ''),
        'svpos': r.get('svpos', ''),
        'svend': r.get('svend', ''),
        'svfreq': r.get('svfreq', ''),
        'strategy': r.get('strategy', ''),
        'desc': api.get('api_desc', ''),
    })

rdf = pd.DataFrame(rows).sort_values('cv_r2', ascending=False).reset_index(drop=True)

outf = f"{BASE}/all_genes_final.tsv"
rdf.to_csv(outf, sep='\t', index=False)
print(f"\nSaved {len(rdf)} genes to {outf}")

print()
print("=" * 200)
print("COMPREHENSIVE GENE TABLE — ALL SOURCES")
print("=" * 200)
hdr = (f"{'#':>3s} {'ENSG':<18s} {'Symbol':<14s} {'Chr':>3s} {'dist':>10s} {'W':>2s} "
       f"{'risk':<6s} {'biotype':<20s} "
       f"{'M_peri':>5s} {'M_full':>5s} {'M_n':>6s} {'M_p':>10s} "
       f"{'QTL':>10s} {'c/t':>10s} {'windows':>14s} {'nc':>2s} {'nt':>2s} "
       f"{'cv_r2':>6s} {'qval':>10s} {'sv_chr':>6s} {'sv_id':>18s} {'svtype':>6s} {'svlen':>7s} {'svfreq':>6s} "
       f"{'strat':<10s} {'desc':<45s}")
print(hdr)
print("-" * 200)

for i, (_, r) in enumerate(rdf.iterrows(), 1):
    mn = str(r['mage_full_n']) if r['mage_full'] == 'Y' else '-'
    mp = f"{r['mage_full_p']:.1e}" if r['mage_full'] == 'Y' and r['mage_full_p'] < 1 else '-'
    svf = f"{float(r['svfreq']):.3f}" if r['svfreq'] and str(r['svfreq']) not in ['', 'nan'] else '-'
    svl = str(r['svlen'])[:7] if r['svlen'] and str(r['svlen']) not in ['', 'nan'] else '-'
    svt = str(r['svtype'])[:6] if r['svtype'] and str(r['svtype']) not in ['', 'nan'] else '-'
    bt = str(r['final_biotype'])[:20] if r['final_biotype'] else '-'
    print(f"{i:>3d} {r['ENSG']:<18s} {str(r['final_symbol']):<14s} {str(r['gene_chr']):>3s} {str(r['dist_Mb']):>10s} {str(r['Wing']):>2s} "
          f"{r['final_risk']:<6s} {bt:<20s} "
          f"{r['mage_pericen']:>5s} {r['mage_full']:>5s} {mn:>6s} {mp:>10s} "
          f"{r['QTL']:>10s} {r['cis_trans']:>10s} {r['windows']:>14s} {r['n_cis']:>2d} {r['n_trans']:>2d} "
          f"{r['cv_r2']:>6.3f} {r['qval']:>10.2e} {str(r['sv_chrom']):>6s} {str(r['sv_id']):>18s} {svt:>6s} {svl:>7s} {svf:>6s} "
          f"{str(r['strategy']):<10s} {str(r['desc'])[:45]:<45s}")

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)
print(f"Total genes: {len(rdf)}")
for risk in ['low', 'moderate', 'high', '?']:
    sub = rdf[rdf['final_risk'] == risk]
    if len(sub) == 0:
        continue
    mf = (sub['mage_full'] == 'Y').sum()
    mp = (sub['mage_pericen'] == 'Y').sum()
    print(f"  {risk:<10s}: {len(sub):>4d} | MAGE_full: {mf}/{len(sub)} | MAGE_pericen: {mp}/{len(sub)}")

print()
print(f"MAGE comparison:")
print(f"  Pericentromeric MAGE: {(rdf['mage_pericen']=='Y').sum()} YES, {(rdf['mage_pericen']=='N').sum()} NO")
print(f"  Full MAGE:            {(rdf['mage_full']=='Y').sum()} YES, {(rdf['mage_full']=='N').sum()} NO")
changed = ((rdf['mage_pericen'] == 'N') & (rdf['mage_full'] == 'Y')).sum()
print(f"  Changed N→Y:          {changed}")
truly_novel = ((rdf['mage_full'] == 'N')).sum()
print(f"  Truly novel:          {truly_novel}")
