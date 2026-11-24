#!/usr/bin/env python3
"""Materialize one local methylation BED link/copy and tabix index."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-bed", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--tabix", default=None)
    parser.add_argument("--bgzip", default=None)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--backend", choices=["auto", "duckhts", "tabix"], default="duckhts")
    parser.add_argument("--copy-non-bgzf", action="store_true")
    return parser.parse_args()


def resolve_tool(name: str, override: str | None) -> str:
    if override:
        path = shutil.which(override) or override
        if Path(path).exists() or shutil.which(path):
            return path
    path = shutil.which(name)
    if path:
        return path
    for candidate in [Path.home() / "downloads" / "bin" / name, Path("/home/jwein22/downloads/bin") / name]:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError(f"Could not find {name}; pass --{name} or add it to PATH")


def run_tabix(tabix: str, bed: Path, threads: int) -> subprocess.CompletedProcess[str]:
    cmd = [tabix, "-f", "-p", "bed"]
    if threads > 1:
        cmd.extend(["-@", str(max(0, threads - 1))])
    cmd.append(str(bed))
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def load_duckhts_connection():
    try:
        import duckdb
    except ImportError:
        return None
    con = duckdb.connect()
    try:
        con.execute("LOAD duckhts")
    except Exception:
        con.execute("INSTALL duckhts FROM community")
        con.execute("LOAD duckhts")
    return con


def run_duckhts_tabix_index(con, bed: Path, threads: int) -> None:
    index_path = Path(f"{bed}.tbi")
    con.execute(
        "SELECT * FROM tabix_index(?, preset := 'bed', index_path := ?, threads := ?)",
        [str(bed), str(index_path), int(threads)],
    ).fetchall()


def bgzip_copy(bgzip: str, source: Path, dest: Path, threads: int) -> None:
    cmd = [bgzip, "-c"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(source))
    with dest.open("wb") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.decode(errors="replace").strip() or f"bgzip failed for {source}")


def replace_symlink(link: Path, target: Path) -> None:
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(target)


def choose_index_backend(requested: str, tabix_override: str | None):
    if requested in {"auto", "duckhts"}:
        try:
            con = load_duckhts_connection()
        except Exception as exc:
            if requested == "duckhts":
                raise RuntimeError("DuckHTS index backend requested, but duckdb/duckhts is not available") from exc
            con = None
        if con is not None:
            return "duckhts", con, None
        if requested == "duckhts":
            raise RuntimeError("DuckHTS index backend requested, but duckdb/duckhts is not available")
    return "tabix", None, resolve_tool("tabix", tabix_override)


def build_index(backend: str, duckhts_con, tabix: str | None, bed: Path, threads: int) -> None:
    if backend == "duckhts":
        run_duckhts_tabix_index(duckhts_con, bed, threads)
        return
    proc = run_tabix(tabix, bed, threads)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or f"tabix failed for {bed}")


def main() -> int:
    args = parse_args()
    source = args.source.resolve()
    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    bgzip = resolve_tool("bgzip", args.bgzip) if args.copy_non_bgzf else None
    backend, duckhts_con, tabix = choose_index_backend(args.backend, args.tabix)

    replace_symlink(args.output_bed, source)
    indexed_mode = "symlink"
    try:
        build_index(backend, duckhts_con, tabix, args.output_bed, args.threads)
    except Exception as exc:
        if not args.copy_non_bgzf or bgzip is None:
            raise RuntimeError(f"{backend} indexing failed for {args.output_bed} linked to {source}: {exc}") from exc
        args.output_bed.unlink(missing_ok=True)
        bgzip_copy(bgzip, source, args.output_bed, args.threads)
        indexed_mode = "bgzip_copy"
        build_index(backend, duckhts_con, tabix, args.output_bed, args.threads)

    index_path = Path(f"{args.output_bed}.tbi")
    if not index_path.exists():
        raise FileNotFoundError(f"tabix did not create expected index: {index_path}")

    summary = {
        "source": str(source),
        "output_bed": str(args.output_bed),
        "index": str(index_path),
        "indexed_mode": indexed_mode,
        "index_backend": backend,
        "tabix": tabix,
        "bgzip": bgzip,
        "index_file_size_bytes": index_path.stat().st_size,
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
