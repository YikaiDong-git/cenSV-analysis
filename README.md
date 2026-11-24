# cenSV-analysis

Analysis code for *Common centromeric variation associates with molecular and complex traits in humans*.
Methods packages: [cenSVhap](https://github.com/YikaiDong-git/cenSVhap), [cenSVimp](https://github.com/YikaiDong-git/cenSVimp).

| Folder | Contents |
|---|---|
| `workflow/`, `config/` | Snakemake workflow: SV calling, haplogroup resolution, imputation, QTL mapping |
| `analyses/` | Downstream analyses |
| `external/` | Collaborator code (methylation QTL, simulation) |
| `resources/` | Reference files |

```bash
pixi install
export CENSV_ROOT=/path/to/project   # also CENSV_COLLAB_ROOT, CENSV_AOU_ROOT, CENSV_WORKSPACE, CENSV_LAB, CENSV_SCRATCH
pixi run snakemake --snakefile workflow/Snakefile --profile workflow/profiles/slurm
```

License: MIT. Contact: yikai.dong@emory.edu, josh.weinstock@emory.edu

## Citation

Dong Y, Cutler D, Weinstock JS. Common centromeric variation associates with molecular and complex
traits in humans. *medRxiv* 2026. [10.64898/2026.10.09.26365154](https://doi.org/10.64898/2026.10.09.26365154)

```bibtex
@article{dong2026centromeric,
  title   = {Common centromeric variation associates with molecular and complex traits in humans},
  author  = {Dong, Yikai and Cutler, David and Weinstock, Joshua S.},
  journal = {medRxiv},
  year    = {2026},
  doi     = {10.64898/2026.10.09.26365154}
}
```
