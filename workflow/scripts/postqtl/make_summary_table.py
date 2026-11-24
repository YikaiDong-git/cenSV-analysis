#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv("results/analysis_qtl/final_tables/all_genes_final.tsv", sep='\t')
df.columns = df.columns.str.strip()
for col in df.select_dtypes(include='object').columns:
    df[col] = df[col].str.strip()

cols = ["ENSG", "final_symbol", "final_biotype", "cv_r2", "qval", "final_risk",
        "gene_chr", "mage_pericen", "mage_full", "QTL", "ens", "windows", "dist_Mb",
        "sv_id", "svlen", "svfreq", "desc"]
sub = df[cols].copy()
sub.to_csv("results/analysis_qtl/final_tables/all_genes_summary.tsv", sep='\t', index=False)
print(f"Saved {len(sub)} genes to all_genes_summary.tsv")
