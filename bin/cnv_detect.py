#!/usr/bin/env python3
"""
cnv_detect.py — per-sample CNV event calling from mosdepth output.

Consumes mosdepth's regions.bed.gz (per-panel-gene depth) and
mosdepth.summary.txt (per-contig depth), and emits:

  <prefix>.cnv_events.tsv   one row per panel gene
                            + optional Chr5_full row if aneuploidy detected
  <prefix>.contig_depth.tsv one row per contig (scan-by-contig)

Classifier (5 categories, applied to normalised depth):
   norm >= 5.0   -> amplification     (10-15x ERG11, Burrack 2022)
   norm >= 3.0   -> dup_high
   norm >= 1.5   -> dup_low           (Carolus 2021 signature)
   norm <= 0.1   -> deletion_total
   norm <= 0.5   -> deletion
   else          -> normal
   no depth      -> no_coverage

Confidence based on genome median depth:
   >=30x -> high     ;   >=10x -> medium    ;    <10x -> low

Aneuploidy detection (per contig):
   ratio_vs_genome >= 2.5 -> aneuploidy_3x
   ratio_vs_genome >= 1.7 -> aneuploidy_2x  (Chr5x2 in C. auris, Li 2024)
   ratio_vs_genome >= 1.3 -> partial_gain
   ratio_vs_genome <= 0.3 -> loss
   ratio_vs_genome <= 0.7 -> partial_loss
   else                   -> normal

If the contig hosting ERG11 (per coords) is flagged aneuploidy_2x,
an extra row gene=Chr5_full is appended to cnv_events.tsv.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path


def classify_locus(norm: float | None) -> str:
    if norm is None:
        return "no_coverage"
    if norm >= 5.0:
        return "amplification"
    if norm >= 3.0:
        return "dup_high"
    if norm >= 1.5:
        return "dup_low"
    if norm <= 0.1:
        return "deletion_total"
    if norm <= 0.5:
        return "deletion"
    return "normal"


def classify_contig(ratio: float | None) -> str:
    if ratio is None:
        return "no_coverage"
    if ratio >= 2.5:
        return "aneuploidy_3x"
    if ratio >= 1.7:
        return "aneuploidy_2x"
    if ratio >= 1.3:
        return "partial_gain"
    if ratio <= 0.3:
        return "loss"
    if ratio <= 0.7:
        return "partial_loss"
    return "normal"


def confidence_from_depth(genome_median: float) -> str:
    if genome_median >= 30:
        return "high"
    if genome_median >= 10:
        return "medium"
    return "low"


def safe_log2(x: float) -> float:
    if x is None or x <= 0:
        return float("-inf")
    return math.log2(x)


def load_panel(path: str) -> dict[str, dict]:
    """Return {gene -> {tier, expected_event, ...}}."""
    out: dict[str, dict] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            g = row.get("gene", "").strip()
            if not g or g.startswith("#"):
                continue
            out[g] = row
    return out


def load_coords(path: str) -> dict[str, tuple[str, int, int]]:
    """Return {gene -> (chrom, start, end)}."""
    out: dict[str, tuple[str, int, int]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            g = row.get("gene", "").strip()
            if not g:
                continue
            out[g] = (row["chrom"], int(row["start"]), int(row["end"]))
    return out


def load_regions_bedgz(path: str) -> dict[str, dict]:
    """Mosdepth regions.bed.gz with 5 columns: chrom start end name mean_depth.

    Returns {gene_name -> {chrom, start, end, mean_depth}}.
    If a gene appears in multiple rows (split panel BED) it averages
    depth weighted by length.
    """
    raw: dict[str, list[tuple[str, int, int, float]]] = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end, name, depth = parts[0], int(parts[1]), int(parts[2]), parts[3], float(parts[4])
            raw[name].append((chrom, start, end, depth))
    merged: dict[str, dict] = {}
    for name, rows in raw.items():
        if len(rows) == 1:
            c, s, e, d = rows[0]
            merged[name] = {"chrom": c, "start": s, "end": e, "mean_depth": d}
        else:
            # weighted average
            total_len = sum(e - s for _, s, e, _ in rows)
            wmean = sum((e - s) * d for _, s, e, d in rows) / total_len if total_len else 0.0
            chrom = rows[0][0]
            start = min(s for _, s, _, _ in rows)
            end = max(e for _, _, e, _ in rows)
            merged[name] = {"chrom": chrom, "start": start, "end": end, "mean_depth": wmean}
    return merged


def load_mosdepth_summary(path: str) -> tuple[float, dict[str, dict]]:
    """Read mosdepth.summary.txt.

    Returns (genome_median_proxy, {contig -> {length, mean, min, max}}).
    Mosdepth's "total" row mean is used as the genome-wide proxy
    (median-vs-mean is a small concession for cleaner contig ratios).
    """
    contigs: dict[str, dict] = {}
    genome_mean = 0.0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            chrom = parts[idx["chrom"]]
            row = {
                "length": int(parts[idx["length"]]),
                "mean": float(parts[idx["mean"]]),
                "min": float(parts[idx.get("min", 4)]),
                "max": float(parts[idx.get("max", 5)]),
            }
            # Mosdepth emits one row per contig PLUS a duplicate '<chrom>_region'
            # row that aggregates over the --by BED regions. Keep only the
            # raw per-contig stat and the global 'total'; skip everything else.
            if chrom == "total":
                genome_mean = row["mean"]
            elif chrom.endswith("_region"):
                continue
            else:
                contigs[chrom] = row
    return genome_mean, contigs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions_bed_gz", required=True,
                    help="mosdepth *.regions.bed.gz (5 cols: chrom start end name mean_depth)")
    ap.add_argument("--summary", required=True,
                    help="mosdepth *.mosdepth.summary.txt (per-contig stats)")
    ap.add_argument("--coords", required=True,
                    help="coords_<reference_tag>.tsv from assets/cnv_loci/coords_pre/")
    ap.add_argument("--panel", required=True,
                    help="assets/cnv_loci/panel.tsv")
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--output_prefix", default=None,
                    help="Output prefix; defaults to sample_id")
    args = ap.parse_args()

    prefix = args.output_prefix or args.sample_id
    panel = load_panel(args.panel)
    coords = load_coords(args.coords)
    regions = load_regions_bedgz(args.regions_bed_gz)
    genome_mean, contig_stats = load_mosdepth_summary(args.summary)

    if genome_mean <= 0:
        print(f"[cnv_detect] WARNING: genome mean depth is 0; outputs will be empty",
              file=sys.stderr)
    conf = confidence_from_depth(genome_mean)
    if conf == "low":
        print(f"[cnv_detect] WARNING: genome mean {genome_mean:.1f}x is below 10x; confidence=low",
              file=sys.stderr)

    # ── contig_depth.tsv ──
    contig_rows = []
    for chrom, stats in sorted(contig_stats.items()):
        ratio = stats["mean"] / genome_mean if genome_mean > 0 else None
        log2 = safe_log2(ratio) if ratio is not None else float("-inf")
        contig_rows.append({
            "sample_id": args.sample_id,
            "contig": chrom,
            "length_bp": stats["length"],
            "mean_depth": round(stats["mean"], 3),
            "ratio_vs_genome": round(ratio, 4) if ratio is not None else "",
            "log2_ratio": round(log2, 4) if log2 != float("-inf") else "",
            "event": classify_contig(ratio),
        })

    contig_path = f"{prefix}.contig_depth.tsv"
    with open(contig_path, "w") as fh:
        cols = ["sample_id", "contig", "length_bp", "mean_depth",
                "ratio_vs_genome", "log2_ratio", "event"]
        fh.write("\t".join(cols) + "\n")
        for r in contig_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[cnv_detect] wrote {contig_path} ({len(contig_rows)} contigs)", file=sys.stderr)

    # ── cnv_events.tsv ──
    event_rows = []
    for gene, panel_row in panel.items():
        if gene == "Chr5_full":
            continue  # handled separately via aneuploidy detection
        tier = panel_row.get("tier", "")
        expected = panel_row.get("expected_event", "")
        coord = coords.get(gene)
        region = regions.get(gene)

        if region is None or coord is None:
            event_rows.append({
                "sample_id": args.sample_id,
                "gene": gene,
                "tier": tier,
                "chrom": coord[0] if coord else "",
                "start": coord[1] if coord else "",
                "end": coord[2] if coord else "",
                "mean_depth": "",
                "genome_median": round(genome_mean, 3),
                "norm": "",
                "log2_ratio": "",
                "event": "no_coverage",
                "confidence": conf,
                "expected_event": expected,
                "agreement_with_expected": "NA",
                "notes": "no mosdepth bin or no coord",
            })
            continue

        depth = region["mean_depth"]
        norm = depth / genome_mean if genome_mean > 0 else None
        log2 = safe_log2(norm) if norm is not None else float("-inf")
        event = classify_locus(norm)
        agree = "true" if expected and event in expected.replace(" ", "").split("_or_") else (
            "false" if expected and event != "normal" else "NA")
        event_rows.append({
            "sample_id": args.sample_id,
            "gene": gene,
            "tier": tier,
            "chrom": coord[0],
            "start": coord[1],
            "end": coord[2],
            "mean_depth": round(depth, 3),
            "genome_median": round(genome_mean, 3),
            "norm": round(norm, 4) if norm is not None else "",
            "log2_ratio": round(log2, 4) if log2 != float("-inf") else "",
            "event": event,
            "confidence": conf,
            "expected_event": expected,
            "agreement_with_expected": agree,
            "notes": "",
        })

    # ── Chr5_full row via aneuploidy heuristic ──
    erg11_coord = coords.get("ERG11")
    if erg11_coord:
        erg11_chrom = erg11_coord[0]
        contig_row = next((r for r in contig_rows if r["contig"] == erg11_chrom), None)
        if contig_row and contig_row["event"] in ("aneuploidy_2x", "aneuploidy_3x"):
            event_rows.append({
                "sample_id": args.sample_id,
                "gene": "Chr5_full",
                "tier": panel.get("Chr5_full", {}).get("tier", "3"),
                "chrom": erg11_chrom,
                "start": 0,
                "end": contig_stats.get(erg11_chrom, {}).get("length", ""),
                "mean_depth": contig_row["mean_depth"],
                "genome_median": round(genome_mean, 3),
                "norm": contig_row["ratio_vs_genome"],
                "log2_ratio": contig_row["log2_ratio"],
                "event": contig_row["event"],
                "confidence": conf,
                "expected_event": "aneuploidy_2x",
                "agreement_with_expected": "true",
                "notes": f"contig hosting ERG11 ({erg11_chrom}) shows whole-chromosome gain",
            })

    events_path = f"{prefix}.cnv_events.tsv"
    with open(events_path, "w") as fh:
        cols = ["sample_id", "gene", "tier", "chrom", "start", "end",
                "mean_depth", "genome_median", "norm", "log2_ratio",
                "event", "confidence", "expected_event",
                "agreement_with_expected", "notes"]
        fh.write("\t".join(cols) + "\n")
        for r in event_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[cnv_detect] wrote {events_path} ({len(event_rows)} rows)", file=sys.stderr)

    # human-readable stderr summary
    by_event: dict[str, int] = defaultdict(int)
    for r in event_rows:
        by_event[r["event"]] += 1
    print(f"[cnv_detect] event distribution: " +
          ", ".join(f"{k}={v}" for k, v in sorted(by_event.items())),
          file=sys.stderr)


if __name__ == "__main__":
    main()
