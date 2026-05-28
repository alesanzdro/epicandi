#!/usr/bin/env python3
"""
qc_gate.py — QC gate v2 (post-IDENTIFICATION).

Five AND criteria:
  1. Sylph Adjusted_ANI for the assigned species   >= qc_min_ani
  2. Sylph taxonomic abundance for that species    >= qc_min_sylph_cov
  3. Read-screen alignment rate to the winning ref >= qc_min_align_rate
  4. Estimated mean depth on the winning ref       >= qc_min_depth
  5. Contamination fraction removed                <= qc_max_contam (legacy)

Depth is extrapolated from total cleaned bases on disk and the alignment rate
reported by READ_SCREEN — matches the user's spec
(`profundidad media = total_bases x align_rate / genome_size`).

Inputs:
  --species_call  TSV from SPECIES_CALL (cols: sylph_ani, sylph_tax_pct,
                  readscreen_top_pct, readscreen_total_reads, classification)
  --genome_size_bp integer (from references_manifest.tsv for meta.reference_tag)
  --clean_r1/r2/nano post-clean fastq.gz paths (used for total base count and
                  contam fraction; pass NO_FILE_* placeholders when absent).

Any single criterion failing  ->  qc_flag = FAIL.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


def count_fastq_bases(path: str | None) -> tuple[int, int]:
    """Return (n_reads, total_bases) for a gzipped fastq, or (0, 0) if missing."""
    if not path:
        return 0, 0
    p = Path(path)
    if not p.exists() or p.name.startswith("NO_FILE"):
        return 0, 0
    n_reads = 0
    n_bases = 0
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:   # sequence line
                n_bases += len(line.rstrip("\n"))
                n_reads += 1
    return n_reads, n_bases


def read_species_call(path: str) -> dict[str, str]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return {}
    return rows[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id",       required=True)
    ap.add_argument("--seq_platform",    required=True, choices=["illumina", "nanopore", "hybrid"])
    ap.add_argument("--species_call",    required=True, help="species_call.tsv from SPECIES_CALL")
    ap.add_argument("--genome_size_bp",  required=True, type=int)
    ap.add_argument("--clean_r1",        default=None)
    ap.add_argument("--clean_r2",        default=None)
    ap.add_argument("--clean_nano",      default=None)
    ap.add_argument("--decont_r1",       default=None, help="post-decont R1 (for contam_fraction; defaults to clean_r1)")
    ap.add_argument("--decont_r2",       default=None)
    ap.add_argument("--decont_nano",     default=None)
    ap.add_argument("--qc_min_ani",        type=float, default=97.0)
    ap.add_argument("--qc_min_sylph_cov",  type=float, default=85.0)
    ap.add_argument("--qc_min_align_rate", type=float, default=90.0)
    ap.add_argument("--qc_min_depth",      type=float, default=30.0)
    ap.add_argument("--qc_max_contam",     type=float, default=30.0)
    ap.add_argument("--output",          required=True)
    args = ap.parse_args()

    sc = read_species_call(args.species_call)
    if not sc:
        print(f"[qc_gate] ERROR: empty species_call for {args.sample_id}", file=sys.stderr)
        return 1

    sylph_ani         = float(sc.get("sylph_ani") or 0)
    sylph_cov         = float(sc.get("sylph_tax_pct") or 0)
    classification    = sc.get("classification", "")
    species_assigned  = sc.get("species_assigned", "")

    # When READ_SCREEN was skipped (params.skip_readscreen=true), the species_call
    # has readscreen_total_reads=0 and readscreen_top_pct=0. In that case we
    # skip the align-rate criterion (default to 100%) so the gate falls back to
    # the other 4 criteria only. The depth estimate then assumes all reads map.
    rs_total = float(sc.get("readscreen_total_reads") or 0)
    rs_top   = float(sc.get("readscreen_top_pct") or 0)
    readscreen_ran = rs_total > 0
    align_rate = rs_top if readscreen_ran else 100.0

    # Total bases (used for depth + contam_fraction)
    n_illu_r1, b_illu_r1 = count_fastq_bases(args.clean_r1)
    _,         b_illu_r2 = count_fastq_bases(args.clean_r2)
    n_nano,    b_nano    = count_fastq_bases(args.clean_nano)

    # Depth: total_bases on the chosen platform × (align_rate/100) / genome_size_bp
    # For hybrid we use max(illu_depth, nano_depth) since the assembly uses both.
    illu_bases = b_illu_r1 + b_illu_r2
    illu_depth = (illu_bases * align_rate / 100.0) / args.genome_size_bp if args.genome_size_bp else 0
    nano_depth = (b_nano   * align_rate / 100.0) / args.genome_size_bp if args.genome_size_bp else 0
    if args.seq_platform == "illumina":
        est_depth = illu_depth
    elif args.seq_platform == "nanopore":
        est_depth = nano_depth
    else:  # hybrid
        est_depth = max(illu_depth, nano_depth)

    # Legacy contam_fraction (clean_total - decont_total) / clean_total × 100
    decont_b = 0
    clean_b  = 0
    if args.seq_platform in ("illumina", "hybrid"):
        _, db = count_fastq_bases(args.decont_r1 or args.clean_r1)
        _, db2 = count_fastq_bases(args.decont_r2 or args.clean_r2)
        decont_b += db + db2
        clean_b  += illu_bases
    if args.seq_platform in ("nanopore", "hybrid"):
        _, db = count_fastq_bases(args.decont_nano or args.clean_nano)
        decont_b += db
        clean_b  += b_nano
    contam_b = max(clean_b - decont_b, 0)
    contam_pct = round(100.0 * contam_b / max(clean_b, 1), 3)

    # AND-gate on 5 criteria
    reasons: list[str] = []
    if sylph_ani < args.qc_min_ani:
        reasons.append(f"sylph_ani={sylph_ani:.2f}<{args.qc_min_ani}")
    if sylph_cov < args.qc_min_sylph_cov:
        reasons.append(f"sylph_cov={sylph_cov:.2f}<{args.qc_min_sylph_cov}")
    if readscreen_ran and align_rate < args.qc_min_align_rate:
        reasons.append(f"align_rate={align_rate:.2f}<{args.qc_min_align_rate}")
    if est_depth < args.qc_min_depth:
        reasons.append(f"est_depth={est_depth:.1f}<{args.qc_min_depth}")
    if contam_pct > args.qc_max_contam:
        reasons.append(f"contam={contam_pct}>{args.qc_max_contam}")

    flag = "PASS" if not reasons else "FAIL"
    notes = "; ".join(reasons) if reasons else "ok"

    header = [
        "sample_id", "seq_platform", "species_assigned", "classification",
        "sylph_ani", "sylph_cov_pct", "align_rate_pct", "est_depth_x",
        "contam_fraction_pct", "genome_size_bp",
        "clean_bases_illumina", "clean_bases_nanopore",
        "qc_flag", "qc_notes",
    ]
    row = [
        args.sample_id, args.seq_platform, species_assigned, classification,
        f"{sylph_ani:.3f}", f"{sylph_cov:.3f}", f"{align_rate:.3f}", f"{est_depth:.2f}",
        f"{contam_pct}", str(args.genome_size_bp),
        str(illu_bases), str(b_nano),
        flag, notes,
    ]
    with open(args.output, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")

    print(f"[qc_gate] {args.sample_id}: {flag} ({notes})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
