#!/usr/bin/env python3
"""Single-edit assembly simulation using the upstream minimap2/svim-asm path."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import time

import pysam


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def digest(data):
    return hashlib.sha256(data.encode()).hexdigest()


def mutate(sequence, event):
    pos, size, kind = event["position"], event["size"], event["kind"]
    if kind == "CONTROL":
        return sequence, ""
    if not 0 < pos < len(sequence) or size <= 0:
        raise ValueError("Invalid edit coordinates or size")
    if kind == "DEL":
        if pos + size > len(sequence):
            raise ValueError("Deletion extends beyond sequence")
        allele = sequence[pos:pos + size]
        return sequence[:pos] + sequence[pos + size:], allele
    if kind == "TANDEM":
        if pos < size:
            raise ValueError("Tandem source extends before sequence")
        allele = sequence[pos - size:pos]
    elif kind == "RANDOM_INS":
        rng = random.Random(event["seed"])
        allele = "".join(rng.choices("ACGT", k=size))
    else:
        raise ValueError(f"Unknown edit kind: {kind}")
    return sequence[:pos] + allele + sequence[pos:], allele


def run(command, log, **kwargs):
    log.write("$ " + " ".join(map(str, command)) + "\n")
    log.flush()
    return subprocess.run(list(map(str, command)), check=True, stderr=log, **kwargs)


def prepare(config_path):
    config = json.loads(Path(config_path).read_text())
    root = Path(config["output_dir"])
    root.mkdir(parents=True, exist_ok=True)
    manifest = root / "manifest.json"
    if manifest.exists():
        raise FileExistsError(f"Refusing to replace {manifest}")
    fa = pysam.FastaFile(config["reference"])
    sequence = fa.fetch(config["chrom"]).upper()
    events = [dict(id="control", kind="CONTROL", context="control", position=0,
                   size=0, seed=config["seed"])]
    for context, pos in config["positions"].items():
        for kind in config["kinds"]:
            for size in config["sizes"]:
                events.append(dict(id=f"{context}_{kind.lower()}_{size}", kind=kind,
                                   context=context, position=pos, size=size,
                                   seed=config["seed"] + len(events)))
    for event in events:
        changed, allele = mutate(sequence, event)
        event.update(chrom=config["chrom"], allele=allele,
                     allele_sha256=digest(allele), query_length=len(changed),
                     query_sha256=digest(changed))
    config.update(events=events, reference_chr_sha256=digest(sequence),
                  reference_chr_length=len(sequence),
                  coordinate_convention="zero-based half-open; insertion before position",
                  created=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    write_json(manifest, config)
    with (root / "truth.tsv").open("w") as handle:
        keys = [k for k in events[0] if k != "allele"]
        writer = csv.DictWriter(handle, keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(events)
    print(f"Prepared {len(events)} queries: {manifest}")


def index(manifest):
    config = json.loads(Path(manifest).read_text())
    root = Path(config["output_dir"])
    target = root / "reference.asm5.mmi"
    if target.exists():
        raise FileExistsError(target)
    with (root / "index.log").open("w") as log:
        run([config["bin_dir"] + "/minimap2", "-x", "asm5", "-t", "3", "-d",
             str(target) + ".tmp", config["reference"]], log)
    Path(str(target) + ".tmp").replace(target)


def simulate(manifest, task):
    config = json.loads(Path(manifest).read_text())
    event = config["events"][task]
    root = Path(config["output_dir"])
    out = root / "cases" / event["id"]
    out.mkdir(parents=True, exist_ok=True)
    if (out / "complete.json").exists():
        print(f"Already complete: {out}")
        return
    start = time.time()
    binary = Path(config["bin_dir"])
    provenance = dict(event_id=event["id"], slurm_job=os.getenv("SLURM_JOB_ID"),
                      node=os.uname().nodename, started=time.time())
    write_json(out / "started.json", provenance)
    with (out / "commands.log").open("w") as log:
        for tool in ("minimap2", "samtools", "bcftools"):
            run([binary / tool, "--version"], log, stdout=log)
        fa = pysam.FastaFile(config["reference"])
        seq = fa.fetch(event["chrom"]).upper()
        if digest(seq) != config["reference_chr_sha256"]:
            raise ValueError("Reference chromosome changed since preparation")
        query, allele = mutate(seq, event)
        if digest(query) != event["query_sha256"] or allele != event["allele"]:
            raise ValueError("Generated query disagrees with truth")
        with tempfile.TemporaryDirectory(prefix="sv_recovery_", dir=os.getenv("TMPDIR", "/tmp")) as scratch:
            query_path = Path(scratch) / "query.fa"
            with query_path.open("w") as handle:
                handle.write(f">{event['id']}\n")
                for offset in range(0, len(query), 80):
                    handle.write(query[offset:offset + 80] + "\n")
            del query, seq
            bam = out / "alignment.bam"
            sam = Path(scratch) / "alignment.sam"
            with sam.open("wb") as handle:
                run([binary / "minimap2", "-a", "-x", "asm5", "--cs", "-r2k",
                     "-t", "3", root / "reference.asm5.mmi", query_path], log, stdout=handle)
            run([binary / "samtools", "sort", "-@", "3", "-m", "1G", "-T",
                 Path(scratch) / "sort", "-o", bam, sam], log)
            run([binary / "samtools", "index", "-@", "3", bam], log)
            env = os.environ.copy()
            env["PYTHONPATH"] = config["svim_source"]
            caller = [sys.executable, str(Path(config["svim_source"]) / "svim_asm/svim-asm"),
                      "haploid", out / "svim", bam, config["reference"], "--sample", event["id"]]
            for key, value in config["caller_parameters"].items():
                caller.extend(["--" + key, str(value)])
            run(caller, log, env=env, stdout=log)
            run([binary / "bcftools", "view", "-T", config["active_bed"], "-Oz", "-o",
                 out / "active.vcf.gz", out / "svim/variants.vcf"], log)
            run([binary / "bcftools", "index", "-t", out / "active.vcf.gz"], log)
    provenance.update(elapsed_seconds=time.time() - start, completed=time.time())
    write_json(out / "complete.json", provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "index", "run"])
    parser.add_argument("path", help="Config for prepare; manifest for index/run")
    parser.add_argument("--task", type=int)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.path)
    elif args.action == "index":
        index(args.path)
    else:
        if args.task is None:
            parser.error("run requires --task")
        simulate(args.path, args.task)


if __name__ == "__main__":
    main()
