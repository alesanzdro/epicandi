#!/usr/bin/env python3
"""
cnv_amr_atlas.py — high-resolution clinical atlas of CNV + AMR + mutations
for a single sample.

Three modes:

  --mode genome   one figure with all contigs concatenated on X axis, log2
                  ratio scatter (CNVkit .cnr) + segments (.cns) on Y, panel
                  gene markers above, all detected mismatches (assembly vs
                  reference) plotted as vertical ticks coloured by FungAMR
                  classification. Header carries QC + AMR call summary.

  --mode chrom    same layout, zoomed to one contig (--chrom).

  --mode gene     three vertical panels for one (sample, gene):
                    (a) mosdepth depth profile across the gene ± padding
                    (b) CNVkit log2 in bins covering the gene
                    (c) protein alignment ref vs sample with mismatch
                        annotation and FungAMR-known mutations highlighted

Outputs PNG (300 dpi) + SVG (vector) so the user can iterate in Inkscape.

Inputs (any optional unless listed in the mode docstrings):
  --sample_id      str
  --species        clinical name string (header)
  --qc_flag        PASS/WARN/FAIL string (header)
  --cnr            CNVkit *.cnr (chromosome start end gene log2 depth weight)
  --cns            CNVkit *.call.cns segments (incl. cn column)
  --cnv_events     cnv_events.tsv from cnv_detect.py
  --contig_depth   contig_depth.tsv from cnv_detect.py
  --gene_proteins  <sample>.gene_proteins.tsv from extract_gene_proteins.py
  --aln_dir        directory of <sample>__<gene>.aa_alignment.fasta
  --resistance     <sample>.resistance_report.tsv from ChroQueTaS join
  --coords         coords_<reference_tag>.tsv (gene chrom start end strand)
  --panel          assets/cnv_loci/panel.tsv
  --mutation_catalog  assets/cnv_loci/mutation_catalog_auris.tsv
  --output         output prefix (we write <prefix>.png and <prefix>.svg)
  --mode           {genome,chrom,gene}
  --chrom          required for mode=chrom
  --gene           required for mode=gene
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict, OrderedDict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Disable mathtext globally so '_' in sample/gene names does not become a
# subscript renderer. We keep \bf{...} explicit where wanted.
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams["text.parse_math"]  = False


# ── colour palette ────────────────────────────────────────────────────────
TIER_COLOR = {"1": "#C0392B", "2": "#E67E22", "3": "#3498DB", "": "#7F8C8D"}
MISMATCH_COLOR = {
    "known_resistance":   "#C0392B",   # red
    "known_sensitivity":  "#2980B9",   # blue
    "unknown_missense":   "#F1C40F",   # yellow
    "frameshift_or_LoF":  "#8E44AD",   # purple
}
CNV_EVENT_COLOR = {
    "amplification": "#7B241C",
    "dup_high":      "#C0392B",
    "dup_low":       "#E67E22",
    "normal":        "#7F8C8D",
    "deletion":      "#2980B9",
    "deletion_total":"#1F3A93",
    "aneuploidy_2x": "#9B59B6",
    "aneuploidy_3x": "#5B2C6F",
    "no_coverage":   "#BDC3C7",
}


# ── loaders ───────────────────────────────────────────────────────────────
def read_tsv(path: str | None) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_cns(path: str | None) -> list[dict]:
    """CNVkit segments file (.cns or .call.cns, tab-delimited with header)."""
    return read_tsv(path)


def load_panel_meta(panel_path: str | None) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not panel_path or not Path(panel_path).exists():
        return out
    with open(panel_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("gene") or "").strip()
            if g:
                out[g] = row
    return out


def load_coords(path: str | None) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not path or not Path(path).exists():
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("gene") or "").strip()
            if g:
                out[g] = {
                    "chrom":  row["chrom"],
                    "start":  int(row["start"]),
                    "end":    int(row["end"]),
                    "strand": (row.get("strand") or "+").strip(),
                }
    return out


def build_contig_offsets(cns_rows: list[dict],
                         contig_depth_rows: list[dict],
                         coords: dict[str, dict]) -> "OrderedDict[str, tuple[int, int]]":
    """Establish per-contig (length, x_offset) from whichever source has data.
    Returns OrderedDict preserving natural order."""
    lengths: dict[str, int] = {}

    # CNVkit segments give max(end) per contig
    for r in cns_rows:
        c = r["chromosome"]
        e = int(r["end"])
        lengths[c] = max(lengths.get(c, 0), e)
    # contig_depth.tsv: authoritative length_bp
    for r in contig_depth_rows:
        c = r["contig"]
        try:
            l = int(r["length_bp"])
        except (KeyError, ValueError):
            continue
        lengths[c] = max(lengths.get(c, 0), l)
    # coords (just to make sure all gene-bearing contigs exist)
    for g, c in coords.items():
        lengths[c["chrom"]] = max(lengths.get(c["chrom"], 0), c["end"])

    # Order by natural contig name (CP043531.1 < CP043532.1, etc.)
    sorted_contigs = sorted(lengths.keys())
    offsets: OrderedDict[str, tuple[int, int]] = OrderedDict()
    cursor = 0
    for c in sorted_contigs:
        offsets[c] = (lengths[c], cursor)
        cursor += lengths[c]
    return offsets


def x_coord(contig: str, pos: int,
            offsets: "OrderedDict[str, tuple[int, int]]") -> float | None:
    if contig not in offsets:
        return None
    return offsets[contig][1] + pos


# ── plotting modes ────────────────────────────────────────────────────────
def plot_genome(args, panel_meta, coords, cnr, cns, contig_depth,
                cnv_events, gene_proteins, resistance_report) -> None:
    offsets = build_contig_offsets(cns, contig_depth, coords)
    if not offsets:
        print("[atlas] no contigs to plot", file=sys.stderr)
        return
    total_bp = sum(l for l, _ in offsets.values())

    # Chrom-mode figures are narrower because they show only one contig
    fig_w = 12 if (args.mode == "chrom" or len(offsets) == 1) else 18
    fig = plt.figure(figsize=(fig_w, 9), dpi=300)
    gs = fig.add_gridspec(4, 1, height_ratios=[0.7, 0.6, 4, 1.2], hspace=0.05)
    ax_hdr  = fig.add_subplot(gs[0]); ax_hdr.axis("off")
    ax_pan  = fig.add_subplot(gs[1])
    ax_log2 = fig.add_subplot(gs[2], sharex=ax_pan)
    ax_mut  = fig.add_subplot(gs[3], sharex=ax_pan)

    # ── header ──
    drugs_called = []
    if resistance_report:
        # one row per (sample, mutation); pick whatever drug columns ≠ ''
        for r in resistance_report:
            for d in ["Fluconazole", "Voriconazole", "Itraconazole",
                      "Posaconazole", "Amphotericin", "Echinocandins"]:
                v = (r.get(d) or "").strip()
                if v and v != "":
                    # Tuple format pos/neg
                    pos_s = v.split("/")[0]
                    try:
                        if int(pos_s) <= 4:
                            drugs_called.append(f"{d}=R")
                    except ValueError:
                        pass
        drugs_called = sorted(set(drugs_called))
    n_known_R = sum(1 for r in gene_proteins if r.get("classification") == "known_resistance")
    n_cnv_events = sum(1 for r in cnv_events if r.get("event") not in ("normal", "no_coverage", ""))

    mode_lbl = args.mode.upper()
    if args.mode == "chrom" and args.chrom:
        mode_lbl = f"CHROM ZOOM: {args.chrom}"
    ax_hdr.text(0.01, 0.95, args.sample_id, ha="left", va="top",
                fontsize=13, fontweight="bold", family="monospace",
                transform=ax_hdr.transAxes)
    ax_hdr.text(0.01, 0.55,
                f"{args.species or 'C. auris'}   ·   QC = {args.qc_flag or '?'}   ·   {mode_lbl}",
                ha="left", va="top", fontsize=10, family="monospace",
                transform=ax_hdr.transAxes)
    ax_hdr.text(0.01, 0.18,
                f"contigs = {len(offsets)}    panel_genes = {len(coords)}    "
                f"FungAMR-known R = {n_known_R}    Non-normal CNV = {n_cnv_events}    "
                f"AMR = {', '.join(drugs_called) if drugs_called else 'none'}",
                ha="left", va="top", fontsize=9, family="monospace",
                transform=ax_hdr.transAxes)

    # ── panel gene markers ──
    ax_pan.set_xlim(0, total_bp)
    ax_pan.set_ylim(-0.5, 0.5)
    ax_pan.set_yticks([])
    ax_pan.spines["top"].set_visible(False)
    ax_pan.spines["right"].set_visible(False)
    ax_pan.spines["left"].set_visible(False)
    for gene, c in coords.items():
        x0 = x_coord(c["chrom"], c["start"], offsets)
        x1 = x_coord(c["chrom"], c["end"], offsets)
        if x0 is None:
            continue
        tier = panel_meta.get(gene, {}).get("tier", "")
        col = TIER_COLOR.get(tier, "#7F8C8D")
        # Make marker visible even for small genes
        width = max(x1 - x0, total_bp * 0.0015)
        rect = mpatches.Rectangle((x0, -0.2), width, 0.4,
                                   facecolor=col, edgecolor="black", lw=0.3, alpha=0.85)
        ax_pan.add_patch(rect)
        ax_pan.text(x0 + width / 2, 0.35, gene, ha="center", va="bottom",
                    fontsize=6.5, rotation=45)
    ax_pan.set_ylabel("Panel", fontsize=8, rotation=0, ha="right", va="center")

    # ── log2 ratio scatter + segments ──
    cnr_x, cnr_y = [], []
    for r in cnr:
        c = r["chromosome"]
        if c not in offsets:
            continue
        x = x_coord(c, (int(r["start"]) + int(r["end"])) // 2, offsets)
        try:
            y = float(r["log2"])
        except (ValueError, KeyError):
            continue
        cnr_x.append(x)
        cnr_y.append(y)
    if cnr_x:
        ax_log2.scatter(cnr_x, cnr_y, s=1.2, color="#34495E", alpha=0.3,
                        rasterized=True, edgecolor="none")
    for s in cns:
        c = s["chromosome"]
        if c not in offsets:
            continue
        x0 = x_coord(c, int(s["start"]), offsets)
        x1 = x_coord(c, int(s["end"]), offsets)
        try:
            y = float(s["log2"])
            cn = int(float(s.get("cn", "1")))
        except (ValueError, KeyError):
            continue
        col = "#7F8C8D" if cn == 1 else ("#C0392B" if cn > 1 else "#2980B9")
        ax_log2.plot([x0, x1], [y, y], color=col, lw=2.5, solid_capstyle="butt")
    ax_log2.axhline(0, color="black", lw=0.4, alpha=0.6)
    ax_log2.axhline(1.0, color="#C0392B", lw=0.4, ls="--", alpha=0.5)
    ax_log2.axhline(-1.0, color="#2980B9", lw=0.4, ls="--", alpha=0.5)
    ax_log2.set_ylim(-3, 3)
    ax_log2.set_ylabel("log2 ratio\n(CNVkit)", fontsize=9)
    ax_log2.spines["top"].set_visible(False)
    ax_log2.spines["right"].set_visible(False)

    # ── mismatches / mutations row ──
    ax_mut.set_xlim(0, total_bp)
    ax_mut.set_ylim(0, 1)
    ax_mut.set_yticks([])
    ax_mut.spines["top"].set_visible(False)
    ax_mut.spines["right"].set_visible(False)
    ax_mut.spines["left"].set_visible(False)
    legend_seen: set[str] = set()

    # gene_proteins rows: ref_pos is AA-coord relative to the gene; convert
    # to genomic position via coords (approx: gene_start + (ref_pos-1)*3 for '+'
    # strand, gene_end - (ref_pos-1)*3 for '-').
    for r in gene_proteins:
        g = r.get("gene", "")
        c = coords.get(g)
        if not c:
            continue
        try:
            aa_pos = int(r.get("ref_pos") or 0)
        except ValueError:
            continue
        if aa_pos <= 0:
            continue
        if c["strand"] == "+":
            genome_pos = c["start"] + (aa_pos - 1) * 3
        else:
            genome_pos = c["end"] - (aa_pos - 1) * 3
        x = x_coord(c["chrom"], genome_pos, offsets)
        if x is None:
            continue
        cls = r.get("classification", "unknown_missense")
        col = MISMATCH_COLOR.get(cls, "#7F8C8D")
        h = 0.85 if cls == "known_resistance" else (0.65 if cls == "known_sensitivity"
                                                     else 0.45)
        ax_mut.plot([x, x], [0, h], color=col, lw=1.2 if cls.startswith("known") else 0.6,
                    alpha=0.95 if cls.startswith("known") else 0.5,
                    label=cls if cls not in legend_seen else None)
        legend_seen.add(cls)
        if cls == "known_resistance":
            ax_mut.text(x, h + 0.02, r.get("mutation_id", ""), fontsize=7,
                        rotation=90, ha="center", va="bottom",
                        color=col, fontweight="bold")
    ax_mut.set_ylabel("Mutations", fontsize=9, rotation=0, ha="right", va="center")
    # Figure-level legend below the panel (avoids overlapping data area)
    if legend_seen:
        legend_handles = [
            mpatches.Patch(color=MISMATCH_COLOR["known_resistance"], label="FungAMR known R"),
            mpatches.Patch(color=MISMATCH_COLOR["known_sensitivity"], label="FungAMR known S"),
            mpatches.Patch(color=MISMATCH_COLOR["unknown_missense"], label="unknown missense"),
            mpatches.Patch(color=MISMATCH_COLOR["frameshift_or_LoF"], label="frameshift / LoF"),
        ]
        fig.legend(handles=legend_handles, loc="lower center", ncol=4,
                   fontsize=8, bbox_to_anchor=(0.5, -0.02), frameon=False)

    # ── contig separators + labels (on bottom axis) ──
    cursor = 0
    for c, (l, off) in offsets.items():
        # vertical separator at start of each contig
        for ax in (ax_pan, ax_log2, ax_mut):
            ax.axvline(off, color="black", lw=0.4, alpha=0.6)
        ax_mut.text(off + l / 2, -0.18, c, ha="center", va="top",
                    fontsize=7, transform=ax_mut.get_xaxis_transform(), rotation=0)
    ax_mut.set_xlabel("")  # contig names act as labels
    ax_mut.set_xticks([])

    out_png = f"{args.output}.png"
    out_svg = f"{args.output}.svg"
    plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(out_svg, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"[atlas] wrote {out_png} + {out_svg}", file=sys.stderr)


def plot_chrom(args, panel_meta, coords, cnr, cns, contig_depth,
               cnv_events, gene_proteins, resistance_report) -> None:
    """Same layout as genome, restricted to args.chrom."""
    if not args.chrom:
        sys.exit("--chrom required for mode=chrom")
    # Filter inputs to this contig
    coords_f = {g: c for g, c in coords.items() if c["chrom"] == args.chrom}
    cnr_f = [r for r in cnr if r.get("chromosome") == args.chrom]
    cns_f = [r for r in cns if r.get("chromosome") == args.chrom]
    contig_depth_f = [r for r in contig_depth if r.get("contig") == args.chrom]
    cnv_events_f = [r for r in cnv_events if r.get("chrom") == args.chrom]
    panel_proteins_f = [r for r in gene_proteins
                        if coords.get(r.get("gene", ""), {}).get("chrom") == args.chrom]
    # Re-use the genome plotter
    args.output = f"{args.output}_{args.chrom}"
    plot_genome(args, panel_meta, coords_f, cnr_f, cns_f, contig_depth_f,
                cnv_events_f, panel_proteins_f, resistance_report)


def plot_gene(args, panel_meta, coords, cnr, cns, mosdepth_regions,
              gene_proteins, alignment_fasta) -> None:
    """Per-gene panel: depth + log2 + protein alignment.

    Layout (tall figure so AA alignment is readable):
       row 0: mosdepth depth profile across gene ± padding
       row 1: CNVkit log2 ratio in same window
       row 2+: protein alignment in 60-AA blocks (one block = 4 axes rows:
                ref letters / match bar / sample letters / position ticks)
    """
    if not args.gene or args.gene not in coords:
        sys.exit(f"--gene must be one of: {sorted(coords.keys())}")
    c = coords[args.gene]
    pad = args.padding_bp
    start = max(0, c["start"] - pad)
    end = c["end"] + pad

    # Load alignment fasta first so we know how many blocks of 60 AA
    ref_aa = sample_aa = ""
    if alignment_fasta and Path(alignment_fasta).exists():
        seqs: dict[str, list[str]] = {}
        with open(alignment_fasta) as fh:
            cur = None
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    seqs[cur] = []
                elif cur:
                    seqs[cur].append(line)
        ks = list(seqs.keys())
        if len(ks) >= 2:
            ref_aa = "".join(seqs[ks[0]])
            sample_aa = "".join(seqs[ks[1]])

    # Build per-position class map
    mm_pos_to_class: dict[int, str] = {}
    mm_pos_to_meta: dict[int, dict] = {}
    for r in gene_proteins:
        if r.get("gene") != args.gene:
            continue
        try:
            p = int(r.get("ref_pos") or 0)
        except ValueError:
            continue
        if p <= 0:
            continue
        mm_pos_to_class[p] = r.get("classification", "unknown_missense")
        mm_pos_to_meta[p] = r

    block_w = 60
    n_blocks = (len(ref_aa) + block_w - 1) // block_w if ref_aa else 0
    # Figure height grows with number of blocks
    fig_h = 4.5 + n_blocks * 0.9
    fig = plt.figure(figsize=(15, fig_h), dpi=300)
    gs = fig.add_gridspec(2 + n_blocks, 1,
                          height_ratios=[1.5, 1.5] + [1.0] * n_blocks,
                          hspace=0.55 if n_blocks else 0.3)
    ax_cov  = fig.add_subplot(gs[0])
    ax_log2 = fig.add_subplot(gs[1], sharex=ax_cov)

    # ── depth (mosdepth regions overlapping the window) ──
    mosdepth_xs, mosdepth_ys = [], []
    for r in mosdepth_regions:
        if r.get("chrom") != c["chrom"]:
            continue
        try:
            s, e = int(r["start"]), int(r["end"])
            d = float(r["depth"])
        except (ValueError, KeyError):
            continue
        if e < start or s > end:
            continue
        mosdepth_xs.append((s + e) / 2)
        mosdepth_ys.append(d)
    if mosdepth_xs:
        ax_cov.bar(mosdepth_xs, mosdepth_ys, width=200, color="#34495E", edgecolor="none")
    ax_cov.axvspan(c["start"], c["end"], alpha=0.15, color="green", zorder=-1)
    ax_cov.set_ylabel("mean depth", fontsize=9)
    ax_cov.set_xlim(start, end)
    tier = panel_meta.get(args.gene, {}).get("tier", "?")
    mech = panel_meta.get(args.gene, {}).get("mechanism", "")
    drug = panel_meta.get(args.gene, {}).get("drug_class", "")
    n_known_R = sum(1 for p, cls in mm_pos_to_class.items() if cls == "known_resistance")
    n_unknown = sum(1 for p, cls in mm_pos_to_class.items() if cls == "unknown_missense")
    title_l1 = (f"{args.sample_id}    {args.gene}    Tier {tier}    "
                f"{c['chrom']}:{c['start']:,}-{c['end']:,} ({c['strand']})")
    title_l2 = (f"mechanism: {mech}   |   drug: {drug}   |   "
                f"FungAMR-known R: {n_known_R}   |   unknown missense: {n_unknown}")
    # Title rendered in 2 lines via separate text artists so we can bold the gene name
    ax_cov.set_title("", loc="left")
    fig.text(0.02, 0.985, title_l1, fontsize=11, fontweight="bold",
             family="monospace", ha="left", va="top")
    fig.text(0.02, 0.965, title_l2, fontsize=9, family="monospace",
             ha="left", va="top", color="#444")

    # ── log2 from CNVkit .cnr ──
    cnr_xs, cnr_ys = [], []
    for r in cnr:
        if r.get("chromosome") != c["chrom"]:
            continue
        try:
            s, e = int(r["start"]), int(r["end"])
            y = float(r["log2"])
        except (ValueError, KeyError):
            continue
        if e < start or s > end:
            continue
        cnr_xs.append((s + e) / 2)
        cnr_ys.append(y)
    if cnr_xs:
        # Bars when bins are scarce (per-gene zoom) so they read as real data
        # rather than a few floating dots.
        if len(cnr_xs) <= 8:
            bin_w = (end - start) / max(len(cnr_xs) * 2, 4)
            ax_log2.bar(cnr_xs, cnr_ys, width=bin_w, color="#34495E",
                        edgecolor="black", lw=0.4, alpha=0.7,
                        bottom=0, zorder=2)
        else:
            ax_log2.scatter(cnr_xs, cnr_ys, s=12, color="#34495E", alpha=0.7)
    # Overlay any CNVkit segment that crosses the window
    for s in cns:
        if s.get("chromosome") != c["chrom"]:
            continue
        try:
            x0, x1 = int(s["start"]), int(s["end"])
            y = float(s["log2"])
        except (ValueError, KeyError):
            continue
        if x1 < start or x0 > end:
            continue
        ax_log2.plot([max(x0, start), min(x1, end)], [y, y],
                     color="#C0392B", lw=2.5)
    ax_log2.axhline(0, color="black", lw=0.5)
    ax_log2.axvspan(c["start"], c["end"], alpha=0.15, color="green", zorder=-1)
    ax_log2.set_ylim(-3, 3)
    ax_log2.set_xlim(start, end)
    ax_log2.set_ylabel("log2 (CNVkit)", fontsize=9)
    ax_log2.set_xlabel(f"Position on {c['chrom']}")

    # ── AA alignment, one block of 60 AA per row, monospace via axis coords ──
    for bi in range(n_blocks):
        ax = fig.add_subplot(gs[2 + bi])
        s, e = bi * block_w, min((bi + 1) * block_w, len(ref_aa))
        block_len = e - s
        ax.set_xlim(0, block_w)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([])
        for spine in ("top", "right", "bottom", "left"):
            ax.spines[spine].set_visible(False)

        # Highlight rectangles for each mismatch in this block
        for j in range(block_len):
            pos = s + j + 1
            ra = ref_aa[s + j]
            sa = sample_aa[s + j] if s + j < len(sample_aa) else "-"
            if ra != sa and pos in mm_pos_to_class:
                cls = mm_pos_to_class[pos]
                col = MISMATCH_COLOR.get(cls, "#F1C40F")
                # Highlight the whole column (ref + sample rows)
                ax.add_patch(mpatches.Rectangle(
                    (j, 0.15), 1.0, 0.55,
                    facecolor=col, alpha=0.45, edgecolor="none", zorder=0))
                # Annotate the mutation id above (e.g. Y132F)
                if cls in ("known_resistance", "known_sensitivity"):
                    mid = mm_pos_to_meta[pos].get("mutation_id", "")
                    ax.text(j + 0.5, 0.94, mid, ha="center", va="bottom",
                            fontsize=7.5, color=col, fontweight="bold")

        # Row labels with position numbers
        ax.text(-1.4, 0.62, f"REF {s+1:>4}", ha="right", va="center",
                fontsize=8, family="monospace", color="#555")
        ax.text(-1.4, 0.32, f"SMP {s+1:>4}", ha="right", va="center",
                fontsize=8, family="monospace", color="#555")
        ax.text(block_w + 1.4, 0.62, f"{e:>4}", ha="left", va="center",
                fontsize=8, family="monospace", color="#555")
        ax.text(block_w + 1.4, 0.32, f"{e:>4}", ha="left", va="center",
                fontsize=8, family="monospace", color="#555")

        # AA letters — monospaced via axis coords (1 unit per AA on x)
        for j in range(block_len):
            ra = ref_aa[s + j]
            sa = sample_aa[s + j] if s + j < len(sample_aa) else "-"
            is_mm = (ra != sa)
            ax.text(j + 0.5, 0.62, ra, ha="center", va="center",
                    fontsize=8, family="monospace",
                    color="#000" if is_mm else "#444")
            ax.text(j + 0.5, 0.32, sa, ha="center", va="center",
                    fontsize=8, family="monospace",
                    color="#000" if is_mm else "#444",
                    fontweight="bold" if is_mm else "normal")

        # Tick marks every 10 AA at the bottom
        for tick in range(10, block_len + 1, 10):
            ax.plot([tick, tick], [0.05, 0.13], color="#999", lw=0.5)
            ax.text(tick, 0.02, str(s + tick), ha="center", va="bottom",
                    fontsize=6, color="#888")

    # Legend at the bottom
    legend_handles = [
        mpatches.Patch(color=MISMATCH_COLOR["known_resistance"], label="FungAMR known R"),
        mpatches.Patch(color=MISMATCH_COLOR["known_sensitivity"], label="FungAMR known S"),
        mpatches.Patch(color=MISMATCH_COLOR["unknown_missense"], label="unknown missense"),
        mpatches.Patch(color=MISMATCH_COLOR["frameshift_or_LoF"], label="frameshift / LoF"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=4, fontsize=8, bbox_to_anchor=(0.5, 0.005))

    out_png = f"{args.output}.png"
    out_svg = f"{args.output}.svg"
    plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(out_svg, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"[atlas] wrote {out_png} + {out_svg}", file=sys.stderr)


# ── main ─────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["genome", "chrom", "gene"])
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--species", default="")
    ap.add_argument("--qc_flag", default="")
    ap.add_argument("--cnr", default=None)
    ap.add_argument("--cns", default=None)
    ap.add_argument("--cnv_events", default=None)
    ap.add_argument("--contig_depth", default=None)
    ap.add_argument("--gene_proteins", default=None)
    ap.add_argument("--aln_dir", default=None)
    ap.add_argument("--resistance", default=None)
    ap.add_argument("--coords", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--mutation_catalog", default=None)
    ap.add_argument("--mosdepth_regions_bed", default=None,
                    help="<sample>.regions.bed.gz from mosdepth (mode=gene only)")
    ap.add_argument("--chrom", default=None)
    ap.add_argument("--gene", default=None)
    ap.add_argument("--padding_bp", type=int, default=20000)
    ap.add_argument("--output", required=True,
                    help="output prefix (PNG and SVG will be written)")
    args = ap.parse_args()

    panel_meta = load_panel_meta(args.panel)
    coords = load_coords(args.coords)
    cnr = read_tsv(args.cnr)
    cns = read_cns(args.cns)
    contig_depth = read_tsv(args.contig_depth)
    cnv_events = read_tsv(args.cnv_events)
    gene_proteins = read_tsv(args.gene_proteins)
    resistance = read_tsv(args.resistance)

    if args.mode == "genome":
        plot_genome(args, panel_meta, coords, cnr, cns, contig_depth,
                    cnv_events, gene_proteins, resistance)
    elif args.mode == "chrom":
        plot_chrom(args, panel_meta, coords, cnr, cns, contig_depth,
                   cnv_events, gene_proteins, resistance)
    elif args.mode == "gene":
        # mosdepth regions for the per-gene depth panel
        mosdepth_regions = []
        if args.mosdepth_regions_bed and Path(args.mosdepth_regions_bed).exists():
            import gzip
            op = gzip.open if args.mosdepth_regions_bed.endswith(".gz") else open
            with op(args.mosdepth_regions_bed, "rt") as fh:
                for line in fh:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 5:
                        continue
                    mosdepth_regions.append({"chrom": parts[0], "start": parts[1],
                                             "end": parts[2], "name": parts[3],
                                             "depth": parts[4]})
        aln_fasta = None
        if args.aln_dir:
            aln_fasta = str(Path(args.aln_dir) / f"{args.sample_id}__{args.gene}.aa_alignment.fasta")
        args.output = f"{args.output}_{args.gene}"
        plot_gene(args, panel_meta, coords, cnr, cns, mosdepth_regions,
                  gene_proteins, aln_fasta)


if __name__ == "__main__":
    main()
