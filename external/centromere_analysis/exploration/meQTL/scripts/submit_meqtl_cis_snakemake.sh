#!/usr/bin/env bash
set -euo pipefail

config_file="config/meqtl_cis_autosomes.yaml"
profile_dir="config/profiles/slurm_meqtl"
snakefile="workflow/meqtl_cis/Snakefile"
max_jobs="120"
retries="3"
dry_run="false"
local_run="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) config_file="$2"; shift 2 ;;
        --profile) profile_dir="$2"; shift 2 ;;
        --snakefile) snakefile="$2"; shift 2 ;;
        --jobs) max_jobs="$2"; shift 2 ;;
        --retries) retries="$2"; shift 2 ;;
        --dry-run) dry_run="true"; shift ;;
        --local) local_run="true"; shift ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if command -v snakemake >/dev/null 2>&1; then
    snakemake_cmd=(snakemake)
elif [[ -x ".venv/bin/snakemake" ]]; then
    snakemake_cmd=(".venv/bin/snakemake")
elif command -v uv >/dev/null 2>&1; then
    export UV_CACHE_DIR="${UV_CACHE_DIR:-/tmp/uv-cache}"
    snakemake_cmd=(uv run snakemake)
else
    echo "snakemake is not available on PATH, in .venv/bin, or via uv" >&2
    exit 1
fi

if [[ -x ".venv/bin/python" ]]; then
    python_cmd=(".venv/bin/python")
    python_cmd_string=".venv/bin/python"
elif command -v uv >/dev/null 2>&1; then
    export UV_CACHE_DIR="${UV_CACHE_DIR:-/tmp/uv-cache}"
    python_cmd=(uv run python)
    python_cmd_string="uv run python"
elif command -v python3 >/dev/null 2>&1; then
    python_cmd=(python3)
    python_cmd_string="python3"
else
    echo "python is not available in .venv/bin, via uv, or on PATH" >&2
    exit 1
fi

workflow_output="$("${python_cmd[@]}" - "$config_file" <<'PY'
import sys
import yaml
from pathlib import Path
cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
print(cfg["workflow_output"])
PY
)"
methylation_root="$("${python_cmd[@]}" - "$config_file" <<'PY'
import sys
import yaml
from pathlib import Path
cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
print(cfg["methylation_root"])
PY
)"
methylation_suffix="$("${python_cmd[@]}" - "$config_file" <<'PY'
import sys
import yaml
from pathlib import Path
cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
print(cfg.get("methylation_file_suffix", "_combined.split5hmc.bed.gz"))
PY
)"
methylation_manifest="$("${python_cmd[@]}" - "$config_file" <<'PY'
import sys
import yaml
from pathlib import Path
cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
print(cfg["methylation_manifest"])
PY
)"

mkdir -p "$workflow_output/manifests" "$workflow_output/slurm_logs"
export MEQTL_SLURM_LOG_DIR="$workflow_output/slurm_logs"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/snakemake-cache}"
export XDG_CONFIG_HOME="${XDG_CONFIG_HOME:-/tmp/snakemake-config}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"

"${python_cmd[@]}" scripts/build_methylation_manifest.py \
    --methylation-root "$methylation_root" \
    --file-suffix "$methylation_suffix" \
    --output "$methylation_manifest"

"${python_cmd[@]}" scripts/prepare_meqtl_manifest.py --config "$config_file"

cmd=(
    "${snakemake_cmd[@]}"
    --snakefile "$snakefile"
    --configfile "$config_file"
    --jobs "$max_jobs"
    --retries "$retries"
    --drop-metadata
    --config python_cmd="$python_cmd_string"
)

if [[ "$local_run" != "true" ]]; then
    cmd+=(--profile "$profile_dir")
fi

if [[ "$dry_run" == "true" ]]; then
    cmd+=(-n)
fi

"${cmd[@]}"
