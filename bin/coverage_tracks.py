#!/usr/bin/env python3
"""
04_coverage_tracks.py
======================
Genera coverage tracks por locus ± padding para visualización clínica.

El gráfico más útil para validar un evento CNV: una figura por (muestra ×
locus_panel) con un padding de ±15 kb donde se ve la depth normalizada por
ventana. Una dup limpia es un escalón cuadrado; Chr5x2 es una elevación
plana de todo el contig; una amplificación 10× es un pico vertical.

Inputs:
  --mosdepth_regions  globs a *.regions.bed.gz (uno por muestra)
  --coords            coords_<reference_slug>.tsv
  --panel             assets/cnv_loci/panel.tsv
  --samples           (opcional) lista de sample_id a procesar (default: todos)
  --genes             (opcional) lista de genes (default: todos los del panel)
  --padding_bp        default 15000
  --outdir            output dir

Output:
  <outdir>/<sample_id>__<gene>.coverage_track.png
  <outdir>/_index.tsv        (sample, gene, n_windows, has_event, evidence)

Patrones visuales que vas a buscar:

  Caso 1 (Carolus 2021, dup segmental ERG11+TAC1B):
    - Escalón limpio dentro del locus, depth ~2× vs flancos
    - Breakpoints dentro de ±15 kb del padding visibles

  Caso 2 (Li 2024, Chr5x2):
    - TODA la ventana visible (locus + padding) elevada ~2×
    - SIN escalón. Es una pista de aneuploidía completa, no segmental.
    - Para confirmar Chr5x2 necesitas plot_chr5_aneuploidy (plot 3 del visualizador).

  Caso 3 (Burrack 2022, ERG11 10-15 copias):
    - Pico estrecho con depth 10-15× sobre el locus
    - A veces multiples picos internos (heterogeneidad clonal)

  Caso 4 (knockout ERG3/ERG6/NCP1):
    - Valle de depth 0 dentro de la CDS
    - Si el valle se extiende fuera del CDS → posible deleción mayor

  Caso 5 (CDR1 deletion → S itraconazol, Rybak 2019):
    - Valle profundo, breakpoints claros, padding intacto
"""
from __future__ import annotations
import argparse
import glob
import gzip
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


def read_mosdepth_regions(path):
    """Reads mosdepth .regions.bed.gz → DataFrame(chrom, start, end, depth)."""
    if not Path(path).exists():
        return pd.DataFrame()
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            try:
                rows.append((parts[0], int(parts[1]), int(parts[2]), float(parts[3])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "depth"])


def sample_id_from_path(path):
    """e.g. /results/cnv/SAMPLE_A/SAMPLE_A.mosdepth.regions.bed.gz → SAMPLE_A"""
    base = Path(path).name
    return re.sub(r"\.(mosdepth\.)?regions\.bed\.gz$", "", base)


def compute_genome_median(df):
    """Median over non-zero windows (matches cnv_detect.py)."""
    nz = df[df["depth"] > 0]["depth"]
    return float(nz.median()) if len(nz) else 0.0


def plot_locus(sample_id, gene, coord_row, df_regions, padding, outdir, panel_info):
    chrom = coord_row["chrom"]
    g_start = int(coord_row["start"])
    g_end   = int(coord_row["end"])
    pad_start = max(0, g_start - padding)
    pad_end   = g_end + padding

    median = compute_genome_median(df_regions)
    if median == 0:
        return {"sample_id": sample_id, "gene": gene, "n_windows": 0,
                "has_event": False, "evidence": "no_coverage"}

    sub = df_regions[(df_regions["chrom"] == chrom)
                      & (df_regions["end"] > pad_start)
                      & (df_regions["start"] < pad_end)].copy()
    if sub.empty:
        return {"sample_id": sample_id, "gene": gene, "n_windows": 0,
                "has_event": False, "evidence": "no_windows_in_region"}

    sub["mid"] = (sub["start"] + sub["end"]) // 2
    sub["norm"] = sub["depth"] / median

    fig, ax = plt.subplots(figsize=(10, 3.2))

    # Color por categoría depth
    colors = []
    for v in sub["norm"]:
        if v >= 5.0:    colors.append("#7B241C")   # amplification massive
        elif v >= 3.0:  colors.append("#C0392B")   # high dup
        elif v >= 1.5:  colors.append("#E67E22")   # low dup
        elif v <= 0.05: colors.append("#1F3A93")   # deletion (~0)
        elif v <= 0.5:  colors.append("#2980B9")   # del
        else:            colors.append("#7F8C8D")   # normal

    ax.bar(sub["mid"], sub["norm"], width=(sub["end"] - sub["start"]).iloc[0],
           color=colors, edgecolor="none", align="center")

    # Líneas de referencia
    ax.axhline(1.0, color="black", lw=0.6, ls="-",   label="1× (normal)")
    ax.axhline(1.5, color="orange", lw=0.6, ls="--", label="1.5× (dup low)")
    ax.axhline(3.0, color="red",    lw=0.6, ls="--", label="3× (dup high)")
    ax.axhline(0.5, color="blue",   lw=0.6, ls="--", label="0.5× (del)")

    # Sombrear el CDS
    ax.axvspan(g_start, g_end, alpha=0.15, color="green", zorder=-1)
    ax.text((g_start + g_end) / 2, ax.get_ylim()[1] * 0.93,
            f"{gene} CDS\n{(g_end - g_start)/1000:.1f} kb",
            ha="center", va="top", fontsize=8, fontweight="bold",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="green", alpha=0.8))

    # Llamada simple
    locus_depth = sub[(sub["start"] >= g_start) & (sub["end"] <= g_end)]["norm"]
    locus_mean = float(locus_depth.mean()) if not locus_depth.empty else 0.0

    if locus_mean >= 5.0:
        evidence = f"amplification_{locus_mean:.1f}x"; has_event = True
    elif locus_mean >= 3.0:
        evidence = f"dup_high_{locus_mean:.1f}x"; has_event = True
    elif locus_mean >= 1.5:
        evidence = f"dup_low_{locus_mean:.1f}x"; has_event = True
    elif locus_mean <= 0.1:
        evidence = "deletion_total"; has_event = True
    elif locus_mean <= 0.5:
        evidence = f"deletion_{locus_mean:.2f}x"; has_event = True
    else:
        evidence = f"normal_{locus_mean:.2f}x"; has_event = False

    # Título con contexto
    pinfo = panel_info.get(gene, {})
    mech = pinfo.get("mechanism", "?")
    drug = pinfo.get("drug_class", "?")
    title = (f"{sample_id} · {gene} ({chrom}:{g_start:,}-{g_end:,}) · "
             f"median locus depth = {locus_mean:.2f}×\n"
             f"mechanism: {mech} · drug class: {drug}")
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_xlabel(f"Position on {chrom} (bp)")
    ax.set_ylabel("Normalized depth (vs genome median)")
    ax.set_xlim(pad_start, pad_end)

    # Anota auto-escalado con call
    badge_color = {"dup_low": "#E67E22", "dup_high": "#C0392B",
                   "amplification": "#7B241C", "deletion": "#1F3A93",
                   "normal": "#7F8C8D"}
    bc = next((c for prefix, c in badge_color.items() if evidence.startswith(prefix)),
              "#7F8C8D")
    ax.text(0.99, 0.97, evidence.upper(),
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            fontweight="bold", color="white",
            bbox=dict(boxstyle="round", facecolor=bc, edgecolor="black", lw=0.5))

    ax.legend(loc="upper left", fontsize=7, framealpha=0.85)
    plt.tight_layout()
    out = Path(outdir) / f"{sample_id}__{gene}.coverage_track.png"
    plt.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close()

    return {"sample_id": sample_id, "gene": gene,
            "n_windows": len(sub),
            "has_event": has_event,
            "evidence": evidence,
            "locus_mean_depth": locus_mean,
            "plot_path": str(out)}


def load_panel_info(panel_tsv):
    """Returns dict gene → {mechanism, drug_class, tier, expected_event}."""
    out = {}
    if not Path(panel_tsv).exists():
        return out
    df = pd.read_csv(panel_tsv, sep="\t")
    for _, r in df.iterrows():
        out[r["gene"]] = {
            "tier": r.get("tier", "?"),
            "mechanism": r.get("mechanism", "?"),
            "drug_class": r.get("drug_class", "?"),
            "expected_event": r.get("expected_event", "?"),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mosdepth_regions", nargs="+", required=True,
                    help="Globs to mosdepth *.regions.bed.gz files")
    ap.add_argument("--coords", required=True,
                    help="coords_<reference_slug>.tsv (or coords_pre/*.tsv)")
    ap.add_argument("--panel",  required=True)
    ap.add_argument("--samples", nargs="*", default=None,
                    help="Filter to these sample_ids")
    ap.add_argument("--genes",   nargs="*", default=None,
                    help="Filter to these gene names")
    ap.add_argument("--padding_bp", type=int, default=15000)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--only_events", action="store_true",
                    help="Only plot tracks where the locus shows a CNV event")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    coords = pd.read_csv(args.coords, sep="\t")
    panel_info = load_panel_info(args.panel)

    # Resolve glob list of bed.gz files → unique paths
    bed_files = []
    for g in args.mosdepth_regions:
        bed_files.extend(glob.glob(g))
    bed_files = sorted(set(bed_files))

    if not bed_files:
        print(f"[cov_tracks] No bed.gz files matched. Globs: {args.mosdepth_regions}",
              file=sys.stderr); sys.exit(1)
    print(f"[cov_tracks] Found {len(bed_files)} mosdepth regions files", file=sys.stderr)

    # Filtros
    genes_to_plot = (set(args.genes)
                      if args.genes
                      else set(coords["gene"]))
    samples_filter = set(args.samples) if args.samples else None

    index_rows = []
    for bed in bed_files:
        sid = sample_id_from_path(bed)
        if samples_filter and sid not in samples_filter:
            continue
        df = read_mosdepth_regions(bed)
        if df.empty:
            print(f"  [skip] {sid}: empty regions", file=sys.stderr)
            continue
        for _, c in coords.iterrows():
            if c["gene"] not in genes_to_plot:
                continue
            res = plot_locus(sid, c["gene"], c, df,
                              args.padding_bp, args.outdir, panel_info)
            if args.only_events and not res["has_event"]:
                # remove the plot that was just written
                try: Path(res["plot_path"]).unlink()
                except (KeyError, FileNotFoundError): pass
                continue
            index_rows.append(res)

    if index_rows:
        idx = pd.DataFrame(index_rows)
        idx.to_csv(Path(args.outdir) / "_index.tsv", sep="\t", index=False)
        n_ev = idx["has_event"].sum()
        print(f"[cov_tracks] {len(idx)} tracks plotted, {n_ev} with event flag",
              file=sys.stderr)
        print(f"[cov_tracks] index: {args.outdir}/_index.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
