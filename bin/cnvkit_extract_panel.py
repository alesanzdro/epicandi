#!/usr/bin/env python3
"""
cnvkit_extract_panel.py
=======================
Extrae las llamadas CNVkit (.call.cns segmentado y con --ploidy 1)
SOLO en los loci del panel, en formato compatible con cnv_detect.py.

Esto permite que aggregate_cnv.py compare cell-to-cell:
   mosdepth_call vs cnvkit_call → si difieren, marca LOW_CONFIDENCE

Input (.call.cns columnas estándar CNVkit):
   chromosome  start  end  gene  log2  cn  depth  probes  weight

Output TSV:
   sample_id  gene  chrom  start  end  cnvkit_cn  cnvkit_log2
   cnvkit_event(dup/del/normal/no_cov)  cnvkit_depth

Mapping cn → event para haploides (ploidy 1):
   cn = 0     → del
   cn = 1     → normal
   cn >= 2    → dup
   cn = NA    → no_cov
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def cn_to_event(cn, ploidy=1):
    if pd.isna(cn):
        return "no_cov"
    cn = int(round(float(cn)))
    if cn < ploidy:
        return "del"
    if cn > ploidy:
        return "dup"
    return "normal"


def find_overlap(cns_df, chrom, start, end):
    """Find CNVkit segments overlapping the gene locus. Return best (longest)."""
    hits = cns_df[(cns_df["chromosome"] == chrom)
                  & (cns_df["start"] <= end)
                  & (cns_df["end"]   >= start)]
    if hits.empty:
        return None
    # If multiple segments overlap → pick the one with max overlap length
    def ov(row):
        return min(row["end"], end) - max(row["start"], start)
    hits = hits.copy()
    hits["_ov"] = hits.apply(ov, axis=1)
    return hits.loc[hits["_ov"].idxmax()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls",     required=True, help="<sample>.call.cns from cnvkit call --ploidy 1")
    ap.add_argument("--coords",    required=True, help="coords_<ref>.tsv from derive_cnv_coords or GFF")
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--ploidy",    type=int, default=1)
    ap.add_argument("--output",    required=True)
    args = ap.parse_args()

    if not Path(args.calls).exists():
        print(f"[cnvkit_extract] {args.calls} not found — emitting empty",
              file=sys.stderr)
        with open(args.output, "w") as fh:
            fh.write("sample_id\tgene\tchrom\tstart\tend\tcnvkit_cn\t"
                     "cnvkit_log2\tcnvkit_event\tcnvkit_depth\n")
        return

    # CNVkit .call.cns columns: chromosome, start, end, gene, log2, cn, depth, probes, weight
    cns = pd.read_csv(args.calls, sep="\t")
    if "cn" not in cns.columns:
        print(f"[cnvkit_extract] WARN: 'cn' column missing — was --ploidy 1 used in 'cnvkit call'?",
              file=sys.stderr)
        cns["cn"] = pd.NA

    coords = pd.read_csv(args.coords, sep="\t")

    rows = []
    for _, c in coords.iterrows():
        gene = c["gene"]
        seg = find_overlap(cns, c["chrom"], int(c["start"]), int(c["end"]))
        if seg is None:
            rows.append([args.sample_id, gene, c["chrom"], c["start"], c["end"],
                         "NA", "NA", "no_cov", "NA"])
            continue
        cn = seg.get("cn", pd.NA)
        log2 = seg.get("log2", pd.NA)
        depth = seg.get("depth", pd.NA)
        event = cn_to_event(cn, args.ploidy)
        rows.append([args.sample_id, gene, c["chrom"], c["start"], c["end"],
                     cn if not pd.isna(cn) else "NA",
                     f"{log2:.3f}" if not pd.isna(log2) else "NA",
                     event,
                     f"{depth:.2f}" if not pd.isna(depth) else "NA"])

    with open(args.output, "w") as fh:
        fh.write("sample_id\tgene\tchrom\tstart\tend\tcnvkit_cn\t"
                 "cnvkit_log2\tcnvkit_event\tcnvkit_depth\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"[cnvkit_extract] {args.sample_id}: {len(rows)} loci, "
          f"{sum(r[7] == 'dup' for r in rows)} dup, "
          f"{sum(r[7] == 'del' for r in rows)} del", file=sys.stderr)


if __name__ == "__main__":
    main()
