#!/usr/bin/env python3
"""Compare individual SVIM calls and local call sets with injected sequence truth."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import pysam


def normalize_alleles(start, ref, alt):
    prefix=0
    while prefix<min(len(ref),len(alt)) and ref[prefix]==alt[prefix]:
        prefix+=1
    ref_end,alt_end=len(ref),len(alt)
    while ref_end>prefix and alt_end>prefix and ref[ref_end-1]==alt[alt_end-1]:
        ref_end-=1;alt_end-=1
    return start+prefix,start+ref_end,alt[prefix:alt_end]


def apply_edits(sequence, origin, edits):
    pieces, cursor = [], origin
    for start, end, alt in sorted(edits):
        if start < cursor or start < origin or end > origin + len(sequence):
            raise ValueError("Edits overlap or fall outside the reference slice")
        pieces.extend([sequence[cursor - origin:start - origin], alt])
        cursor = end
    pieces.append(sequence[cursor - origin:])
    return "".join(pieces)


def equivalent(reference, chrom, truth, calls):
    edits = [truth] + calls
    left = max(0, min(e[0] for e in edits) - 1)
    right = min(reference.get_reference_length(chrom), max(e[1] for e in edits) + 1)
    seq = reference.fetch(chrom, left, right).upper()
    try:
        return apply_edits(seq, left, [truth]) == apply_edits(seq, left, calls)
    except ValueError:
        return False


def parse_calls(path, reference, region=None, kinds=None):
    calls = []
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            if region is not None:
                chrom, left, right = region
                if rec.chrom != chrom or rec.start > right or rec.stop < left:
                    continue
            kind = rec.info.get("SVTYPE", "UNKNOWN")
            if kinds is not None and kind not in kinds:
                continue
            edit = None
            if kind == "DUP:TANDEM":
                cn = next(iter(rec.samples.values())).get("CN", 2)
                if isinstance(cn, tuple):
                    cn = cn[0]
                cn = 2 if cn is None else int(cn)
                allele = reference.fetch(rec.chrom, rec.start, rec.stop).upper() * (cn - 1)
                edit = (rec.stop, rec.stop, allele)
            elif rec.alts and not rec.alts[0].startswith("<") and "[" not in rec.alts[0] and "]" not in rec.alts[0]:
                edit = normalize_alleles(rec.start, rec.ref.upper(), rec.alts[0].upper())
            elif kind == "DEL":
                edit = (rec.pos, rec.stop, "")
            length = rec.info.get("SVLEN", 0)
            if isinstance(length, tuple):
                length = length[0]
            if edit:
                length = len(edit[2]) - (edit[1] - edit[0])
            calls.append(dict(id=rec.id, chrom=rec.chrom, pos=rec.pos, end=rec.stop,
                              kind=kind, length=int(length), filters=list(rec.filter), edit=edit))
    return calls


def signature(call):
    return (call["chrom"], call["pos"], call["end"], call["kind"], call["length"],
            tuple(call["edit"]) if call["edit"] else None)


def evaluate(event, calls, baseline, reference, active_intervals):
    pos, size = event["position"], event["size"]
    kind, chrom = event["kind"], event["chrom"]
    truth = (pos, pos + size, "") if kind == "DEL" else (pos, pos, event["allele"])
    expected_length = -size if kind == "DEL" else size
    base_keys = {signature(c) for c in baseline}
    new = [c for c in calls if signature(c) not in base_keys]
    local = [c for c in new if c["chrom"] == chrom and
             min(abs(c["pos"] - pos), abs(c["end"] - pos)) <= 10000]
    candidates = []
    if kind != "CONTROL":
        for call in new:
            edit = call["edit"]
            if call["chrom"] != chrom or edit is None:
                continue
            if (call["length"] > 0) != (expected_length > 0):
                continue
            displacement = edit[0] - pos
            if abs(displacement) > 1000000:
                continue
            seq_match = equivalent(reference, chrom, truth, [edit])
            length_error = abs(call["length"]) - size
            relaxed = abs(displacement) <= 100 and abs(length_error) <= max(10, size * .1)
            exact = edit == truth
            candidates.append((seq_match, exact, relaxed, -abs(length_error),
                               -abs(displacement), call))
    candidates.sort(key=lambda c: c[:5], reverse=True)
    best = candidates[0] if candidates else None
    sequence_match = bool(best and best[0])
    exact = bool(best and best[1])
    relaxed = any(c[2] for c in candidates)
    combined = False
    compatible_local = [c for c in local if c["edit"] and
                        (c["length"] > 0) == (expected_length > 0)]
    if kind != "CONTROL" and len(compatible_local) > 1:
        combined = equivalent(reference, chrom, truth, [c["edit"] for c in compatible_local])
    if kind == "CONTROL":
        status = "control"
    elif exact:
        status = "exact"
    elif sequence_match:
        status = "repeat_equivalent"
    elif combined:
        status = "fragmented_equivalent"
    elif relaxed:
        status = "approximate"
    else:
        status = "not_recovered"
    best_call = best[-1] if best else None
    explained = 1 if sequence_match else len(compatible_local) if combined else 0
    active = [c for c in calls if any(c["chrom"] == ch and a <= c["pos"] - 1 < b
                                     for ch, a, b in active_intervals)]
    return dict(id=event["id"], context=event["context"], kind=kind, size=size,
                position_0based=pos, status=status, exact=exact,
                single_call_sequence_equivalent=sequence_match,
                local_call_set_sequence_equivalent=combined,
                approximate_100bp_10pct=relaxed, total_calls=len(calls),
                active_calls=len(active), new_calls_vs_control=len(new),
                unexplained_new_calls=len(new) - explained,
                new_local_calls_10kb=len(local), best_call_id=best_call["id"] if best_call else "",
                best_call_type=best_call["kind"] if best_call else "",
                best_call_filters=",".join(best_call["filters"]) if best_call else "",
                breakpoint_shift_bp=best_call["edit"][0] - pos if best_call else "",
                length_error_bp=abs(best_call["length"]) - size if best_call else ""), local


def write_tsv(path, rows):
    if not rows:
        return
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_results(rows, root, config):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap, BoundaryNorm
    statuses = ["not_recovered", "approximate", "fragmented_equivalent", "repeat_equivalent", "exact"]
    colors = ["#d55e00", "#e6ab02", "#56b4e9", "#009e73", "#0072b2"]
    contexts = list(config["positions"])
    fig, axes = plt.subplots(1, len(contexts), figsize=(11, 3.5), sharey=True)
    for ax, context in zip(axes, contexts):
        matrix = np.full((len(config["kinds"]), len(config["sizes"])), np.nan)
        for row in rows:
            if row["context"] == context:
                y, x = config["kinds"].index(row["kind"]), config["sizes"].index(row["size"])
                matrix[y, x] = statuses.index(row["status"])
                if row["single_call_sequence_equivalent"] and not row["production_sequence_equivalent"]:
                    ax.text(x, y, "×", ha="center", va="center", color="black", fontsize=20)
        ax.imshow(matrix, cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-.5, 5), 5), aspect="auto")
        ax.set_title(context.replace("_", " "))
        ax.set_xticks(range(len(config["sizes"])), config["sizes"])
        ax.set_yticks(range(len(config["kinds"])), ["Deletion", "Tandem duplication", "Random insertion"])
        ax.set_xlabel("Injected size (bp)")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c, label=s.replace("_", " ")) for c, s in zip(colors, statuses)],
               loc="lower center", ncol=3, frameon=False)
    fig.suptitle("CHM13 chr1: full caller recovery (× = excluded by active-array filter)")
    fig.tight_layout(rect=(0, .18, 1, .92))
    fig.savefig(root / "recovery.png", dpi=180)
    fig.savefig(root / "recovery.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    config = json.loads(Path(args.manifest).read_text())
    root = Path(config["output_dir"])
    reference = pysam.FastaFile(config["reference"])
    control = root / "cases/control"
    if not (control / "complete.json").exists():
        raise RuntimeError("Unedited control has not completed")
    baseline = parse_calls(control / "svim/variants.vcf", reference)
    filtered_baseline = parse_calls(control / "active.vcf.gz", reference)
    intervals = []
    for line in Path(config["active_bed"]).read_text().splitlines():
        fields = line.split()
        intervals.append((fields[0], int(fields[1]), int(fields[2])))
    rows, missing = [], []
    for event in config["events"]:
        case = root / "cases" / event["id"]
        if not (case / "complete.json").exists():
            missing.append(event["id"])
            continue
        calls = parse_calls(case / "svim/variants.vcf", reference)
        row, local = evaluate(event, calls, baseline, reference, intervals)
        filtered = parse_calls(case / "active.vcf.gz", reference)
        filtered_result, _ = evaluate(event, filtered, filtered_baseline, reference, intervals)
        row["active_calls"] = len(filtered)
        row["truth_in_active_array"] = any(ch == event["chrom"] and a <= event["position"] < b
                                           for ch, a, b in intervals)
        row["production_recovery_status"] = filtered_result["status"]
        row["production_sequence_equivalent"] = (filtered_result["single_call_sequence_equivalent"] or
                                                  filtered_result["local_call_set_sequence_equivalent"])
        with pysam.AlignmentFile(str(case / "alignment.bam")) as bam:
            hits = list(bam.fetch(event["chrom"], max(0, event["position"] - 1), event["position"] + 1))
            row["alignments_overlapping_site"] = len(hits)
            row["alignments_overlapping_site_mapq20"] = sum(r.mapping_quality >= 20 for r in hits)
        rows.append(row)
        (case / "local_calls.json").write_text(json.dumps(local, indent=2) + "\n")
    if missing and not args.allow_partial:
        raise RuntimeError(f"Incomplete cases: {missing}")
    suffix = "_partial" if missing else ""
    write_tsv(root / f"recovery{suffix}.tsv", rows)
    groups = defaultdict(Counter)
    for row in rows:
        if row["kind"] != "CONTROL":
            groups[(row["context"], row["kind"])][row["status"]] += 1
    summary = []
    for (context, kind), counts in groups.items():
        subset = [r for r in rows if r["context"] == context and r["kind"] == kind]
        summary.append(dict(context=context, kind=kind, n=sum(counts.values()),
                            expected_active=sum(r["truth_in_active_array"] for r in subset),
                            retained_sequence_equivalent=sum(r["production_sequence_equivalent"] for r in subset),
                            **{s: counts[s] for s in ["exact", "repeat_equivalent", "fragmented_equivalent", "approximate", "not_recovered"]}))
    write_tsv(root / f"summary{suffix}.tsv", summary)
    if not missing:
        plot_results(rows, root, config)
    edits = [r for r in rows if r["kind"] != "CONTROL"]
    report = dict(completed_queries=len(rows), missing=missing,
                  injected_events=len(edits), outcomes=dict(Counter(r["status"] for r in rows)),
                  sequence_equivalent_recoveries=sum(r["single_call_sequence_equivalent"] or
                                                    r["local_call_set_sequence_equivalent"] for r in edits),
                  expected_active_events=sum(r["truth_in_active_array"] for r in edits),
                  retained_active_sequence_equivalent=sum(r["production_sequence_equivalent"] for r in edits),
                  unexplained_new_calls=sum(r["unexplained_new_calls"] for r in edits),
                  control_calls=next(r["total_calls"] for r in rows if r["kind"] == "CONTROL"))
    (root / f"run_summary{suffix}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
