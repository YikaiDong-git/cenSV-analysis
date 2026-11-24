#!/usr/bin/env python3
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
CDR = ROOT / "analyses/cdr_qtl"
OUT = CDR / "results/snp_conditioning_chr8"
RSCRIPT = CDR / "code/conditioning/run_susie.R"
SEED = "104729"

res = pd.read_csv(OUT / "conditioning_by_panel.tsv", sep="\t")
out_rows = []
for _, r in res.iterrows():
    panel = r["panel"]
    f = OUT / f"susie_input_{panel}.tsv"
    dat = pd.read_csv(f, sep="\t")
    sv = dat["SV"].to_numpy()
    snps = dat.drop(columns=["phenotype", "SV"]).to_numpy()

    sd = snps.std(axis=0); sd[sd == 0] = 1.0
    z = (snps - snps.mean(axis=0)) / sd
    svz = (sv - sv.mean()) / (sv.std() if sv.std() > 0 else 1.0)
    rs = np.abs(z.T @ svz) / len(sv)
    max_r = float(np.nanmax(rs)); n_r08 = int(np.sum(rs >= 0.8)); n_r09 = int(np.sum(rs >= 0.9))

    sus_out = OUT / f"susie_variables_{panel}.tsv"
    subprocess.run(["Rscript", str(RSCRIPT), str(f), str(sus_out), SEED], check=True)
    sus = pd.read_csv(sus_out, sep="\t")
    sv_row = sus[sus["variable"] == "SV"].iloc[0]
    top_snp = sus[sus["variable"] != "SV"].sort_values("pip", ascending=False).iloc[0]

    out_rows.append({
        "panel": panel, "n_susie_snps": int(r["n_susie_snps"]),
        "n_condition_snps": int(r["n_condition_snps"]),
        "max_abs_r_sv_snp": max_r, "n_snps_r_ge_0.8": n_r08, "n_snps_r_ge_0.9": n_r09,
        "sv_susie_pip": float(sv_row["pip"]),
        "sv_credible_sets": str(sv_row["credible_sets"]),
        "top_snp_pip": float(top_snp["pip"]),
        "p_post": float(r["p_post"]), "attenuation_ratio": float(r["attenuation_ratio"]),
        "ridge_resid_p": float(r["ridge_resid_p"]),
        "ridge_resid_beta": float(r["ridge_resid_beta"]),
        "ridge_oof_r2_sv": float(r["ridge_oof_r2"]),
    })
    print(f"  {panel:15s} max|r|={max_r:.3f} (n>=0.8: {n_r08}, >=0.9: {n_r09})  "
          f"SV_PIP={float(sv_row['pip']):.4g} CS='{sv_row['credible_sets']}'  "
          f"top_snp_PIP={float(top_snp['pip']):.3g}")

df = pd.DataFrame(out_rows)
df.to_csv(OUT / "conditioning_verdict.tsv", sep="\t", index=False)
print()
print(df[["panel", "max_abs_r_sv_snp", "n_condition_snps", "sv_susie_pip", "p_post",
          "attenuation_ratio", "ridge_resid_p"]].to_string(index=False))

summary = {
    "n_panels": len(df),
    "max_abs_r_overall": float(df["max_abs_r_sv_snp"].max()),
    "panels_sv_pip_lt_0.05": int((df["sv_susie_pip"] < 0.05).sum()),
    "panels_sv_pip_ge_0.9": int((df["sv_susie_pip"] >= 0.9).sum()),
    "panels_p_post_ge_0.05": int((df["p_post"] >= 0.05).sum()),
    "panels_ridge_p_lt_0.05": int((df["ridge_resid_p"] < 0.05).sum()),
    "midpoint_panels": df[df["panel"].str.startswith("midpoint")][
        ["panel", "max_abs_r_sv_snp", "sv_susie_pip", "p_post", "attenuation_ratio"]
    ].to_dict("records"),
}
(OUT / "verdict_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print()
print(json.dumps(summary, indent=2))
