#!/usr/bin/env python3
"""
03_cnv_visualization.py
========================
Genera gráficos combinando llamadas AMR (de ChroQueTas v2 / call_matrix.tsv)
con eventos CNV (de cnv_events.tsv del módulo CNV_DETECT).

Outputs en results/figures/:

  1. plot_cnv_heatmap.png
       muestras × genes_panel; color = log2(norm)
       diverging colormap (rojo dup, azul del, blanco normal)

  2. plot_amr_cnv_combined.png
       muestras × (drugs|genes), faceted: izquierda = AMR call (R/S/MIXED/NE),
       derecha = CNV event (dup/normal/del/no_cov). Tira de especie/clado a la izq.

  3. plot_chr5_aneuploidy.png
       barras: median_depth_per_contig vs genome_median, log2-scale.
       Resalta automáticamente contigs con ratio > 1.7× (aneuploidía sospechosa).
       UNA figura por muestra que tenga al menos un contig con ratio>1.5.

  4. plot_carolus_signature.png
       Heatmap específico para el genotipo Carolus 2021:
       ERG11 dup + TAC1B dup + FKS1 F635del + CIS2 A27T → ¿cuántas muestras
       reproducen el patrón completo o parcial?

Uso:
    python3 03_cnv_visualization.py \\
        --cnv_events results/cnv/*.cnv_events.tsv \\
        --amr_call results/resistance/aggregated/call_matrix.tsv \\
        --species_id results/species_id_FIXED.tsv \\
        --depth_per_contig results/cnv/*.contig_depth.tsv \\
        --outdir results/figures/

Requisitos: pandas, numpy, matplotlib (Agg), seaborn opcional.
"""
from __future__ import annotations
import argparse
import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


# ============================================================================
# Helpers — coloreado de especies (reusa el mapeo del script de heatmaps AMR)
# ============================================================================
SP_PRETTY = {
    "cladeI_B8441":    "C. auris I",   "cladeII_B11220":  "C. auris II",
    "cladeIII_B11221": "C. auris III", "cladeIV_B11243":  "C. auris IV",
    "cladeIV_B11245":  "C. auris IV",  "cladeV_B18474":   "C. auris V",
    "cladeVI_F3485":   "C. auris VI",
    "Calbicans_SC5314":   "C. albicans",
    "Cparapsilosis_CDC317": "C. parapsilosis",
    "Cglabrata_CBS138":   "C. glabrata",
}
SP_COLOR = {
    "C. albicans":     "#1f77b4", "C. parapsilosis": "#ff7f0e",
    "C. glabrata":     "#8c564b",
    "C. auris I":      "#e377c2", "C. auris II":     "#7f7f7f",
    "C. auris III":    "#bcbd22", "C. auris IV":     "#17becf",
    "C. auris V":      "#ff9896", "C. auris VI":     "#000000",
}
CLADE_ORDER = ["C. albicans", "C. parapsilosis", "C. glabrata",
               "C. auris I", "C. auris II", "C. auris III",
               "C. auris IV", "C. auris V", "C. auris VI"]


# ============================================================================
# Loaders
# ============================================================================
def load_cnv_events(globs):
    """Lee múltiples cnv_events.tsv → DF largo."""
    dfs = []
    for g in globs:
        for f in glob.glob(g):
            d = pd.read_csv(f, sep="\t")
            if not d.empty:
                dfs.append(d)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def load_call_matrix(path):
    if not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t").set_index("sample_id")


def load_species(path):
    if not Path(path).exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["sample_id"], df.get("top_genome_slug", df.iloc[:, 1])))


def encode_amr_call(v):
    if pd.isna(v) or v == "" or v == "NE": return np.nan
    if v == "S":     return 0
    if v == "MIXED": return 1
    if v == "R":     return 2
    if v == "NS":    return -1
    return np.nan


# ============================================================================
# Plot 1 — CNV depth heatmap (muestras × genes)
# ============================================================================
def _fmt_cn_delta(cn: float, ploidy: int) -> str:
    """Cell annotation: '<cn>cn (+/-<delta>)' relative to baseline ploidy.
    Empty string for NaN; pure '<cn>cn' when the delta is zero."""
    if pd.isna(cn):
        return ""
    cn_int = int(round(cn))
    delta = cn_int - ploidy
    if delta == 0:
        return f"{cn_int}cn"
    sign = "+" if delta > 0 else "-"
    return f"{cn_int}cn ({sign}{abs(delta)})"


EMPTY_EVENT_LABELS = {"normal", "", "nan", "no_cov", "none"}


def _genes_with_real_events(events_df):
    """Return the set of genes that have at least one non-`normal` event
    across all samples in the events DataFrame.  Empty / nan / no_cov are
    treated as 'no event detected'.  Returns None if no `event` column
    exists (so caller falls back to keeping all genes)."""
    if events_df is None or events_df.empty or "event" not in events_df.columns:
        return None
    ev = events_df["event"].astype(str).str.strip().str.lower()
    real = events_df[~ev.isin(EMPTY_EVENT_LABELS)]
    return set(real["gene"].dropna().astype(str).unique())


def plot_cnv_heatmap(events, species_of, outdir, ploidy: int = 1,
                     hide_empty_genes: bool = True):
    if events.empty:
        print("[plot_cnv_heatmap] no events, skipping")
        return

    # Pivot: filas=sample_id, cols=gene, valor=log2(norm)
    df = events.copy()
    df["log2_ratio"] = np.log2(df["norm"].replace(0, np.nan))
    df["cn_estimate"] = df["norm"] * ploidy
    pivot = df.pivot_table(index="sample_id", columns="gene",
                            values="log2_ratio", aggfunc="mean")
    cn_pivot = df.pivot_table(index="sample_id", columns="gene",
                              values="cn_estimate", aggfunc="mean")

    # Drop gene columns where no sample has a real (≠ normal) CNV event.
    if hide_empty_genes:
        genes_keep = _genes_with_real_events(events)
        if genes_keep is not None:
            keep = [c for c in pivot.columns if c in genes_keep]
            if not keep:
                print("[plot_cnv_heatmap] no CNV events ≠ normal across cohort — skipping plot")
                return
            n_drop = pivot.shape[1] - len(keep)
            if n_drop:
                print(f"[plot_cnv_heatmap] hiding {n_drop} gene column(s) with no events")
            pivot = pivot[keep]
            cn_pivot = cn_pivot.reindex(columns=keep)

    # Ordenar por especie
    pivot["_sp"] = pivot.index.map(lambda s: SP_PRETTY.get(species_of.get(s, ""), "?"))
    pivot["_sort"] = pivot["_sp"].apply(lambda x: CLADE_ORDER.index(x) if x in CLADE_ORDER else 99)
    pivot = pivot.sort_values(["_sort", "_sp"])
    sp_strip = list(pivot["_sp"])
    pivot = pivot.drop(columns=["_sp", "_sort"])

    n, m = pivot.shape
    fig = plt.figure(figsize=(max(8, m*0.5), max(6, n*0.18)))
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.025, 1.0], wspace=0.02)

    # Tira de especie
    ax_sp = fig.add_subplot(gs[0, 0])
    sp_rgb = np.array([[mcolors.to_rgb(SP_COLOR.get(s, "#cccccc"))
                        for s in sp_strip]]).transpose(1, 0, 2)
    ax_sp.imshow(sp_rgb, aspect="auto", interpolation="nearest")
    ax_sp.set_xticks([]); ax_sp.set_yticks([])

    # Heatmap principal
    ax = fig.add_subplot(gs[0, 1])
    # Diverging: rojo (dup, log2>0) ↔ blanco (normal) ↔ azul (del, log2<0)
    vmax = max(abs(np.nanmin(pivot.values)), abs(np.nanmax(pivot.values)), 1.5)
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest")

    # Cell annotation: estimated copy number and delta vs baseline ploidy.
    cn_aligned = cn_pivot.reindex(index=pivot.index, columns=pivot.columns)
    for i in range(n):
        for j in range(m):
            cell = cn_aligned.iat[i, j] if (i < cn_aligned.shape[0]
                                            and j < cn_aligned.shape[1]) else np.nan
            txt = _fmt_cn_delta(cell, ploidy)
            if not txt:
                continue
            log2v = pivot.iat[i, j]
            # Switch to white on saturated cells for legibility.
            color = "white" if (not np.isnan(log2v) and abs(log2v) > vmax * 0.55) else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=6, color=color)

    ax.set_xticks(range(m))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(pivot.index, fontsize=6)
    # Colorea labels muestra según especie
    for tick, sp in zip(ax.get_yticklabels(), sp_strip):
        tick.set_color(SP_COLOR.get(sp, "#000000"))
    ax.set_title(f"CNV depth heatmap — {n} samples × {m} genes "
                 f"(log2 ratio vs genome median)\n"
                 f"Rojo=duplicación · Azul=deleción · Blanco=normal",
                 fontsize=11, pad=10)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2(norm)")
    cbar.ax.axhline(np.log2(1.5), color="black", lw=0.7, linestyle="--")
    cbar.ax.axhline(np.log2(0.5), color="black", lw=0.7, linestyle="--")
    cbar.ax.text(2.5, np.log2(1.5), " dup thr", fontsize=7, va="center")
    cbar.ax.text(2.5, np.log2(0.5), " del thr", fontsize=7, va="center")

    out = Path(outdir) / "plot_cnv_heatmap.png"
    plt.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.savefig(out.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"→ {out}")


# ============================================================================
# Plot 2 — AMR call × CNV event, side-by-side
# ============================================================================
def plot_amr_cnv_combined(events, call, species_of, outdir, ploidy: int = 1,
                          hide_empty_genes: bool = True):
    if call.empty or events.empty:
        print("[plot_amr_cnv_combined] missing inputs, skipping")
        return

    # Encode AMR
    samples = sorted(set(call.index) & set(events["sample_id"].unique()))
    if not samples:
        print("[plot_amr_cnv_combined] no overlap samples")
        return

    sp_list = [SP_PRETTY.get(species_of.get(s, ""), "?") for s in samples]
    order = sorted(range(len(samples)),
                   key=lambda i: (CLADE_ORDER.index(sp_list[i])
                                  if sp_list[i] in CLADE_ORDER else 99, samples[i]))
    samples = [samples[i] for i in order]
    sp_list = [sp_list[i] for i in order]

    # AMR matrix
    drugs = list(call.columns)
    amr_mat = np.full((len(samples), len(drugs)), np.nan)
    for i, s in enumerate(samples):
        for j, d in enumerate(drugs):
            amr_mat[i, j] = encode_amr_call(call.loc[s, d]) if s in call.index else np.nan

    # CNV matrix — log2 ratio + estimated cn for cell annotation.
    cnv_p = events[events["sample_id"].isin(samples)].copy()
    cnv_p["log2"] = np.log2(cnv_p["norm"].replace(0, np.nan))
    cnv_p["cn_estimate"] = cnv_p["norm"] * ploidy
    cnv_mat = cnv_p.pivot_table(index="sample_id", columns="gene",
                                values="log2", aggfunc="mean")
    cnv_mat = cnv_mat.reindex(index=samples)
    cn_mat = cnv_p.pivot_table(index="sample_id", columns="gene",
                               values="cn_estimate", aggfunc="mean")
    cn_mat = cn_mat.reindex(index=samples, columns=cnv_mat.columns)

    # Drop gene columns with zero non-`normal` events across the overlap cohort.
    # If the cohort has zero real events the whole combined plot is skipped.
    if hide_empty_genes:
        genes_keep = _genes_with_real_events(cnv_p)
        if genes_keep is not None:
            keep = [c for c in cnv_mat.columns if c in genes_keep]
            n_drop = cnv_mat.shape[1] - len(keep)
            if not keep:
                print("[plot_amr_cnv_combined] no CNV events ≠ normal across "
                      "cohort — skipping combined heatmap")
                return
            if n_drop:
                print(f"[plot_amr_cnv_combined] hiding {n_drop} gene column(s) "
                      f"with no events")
            cnv_mat = cnv_mat[keep]
            cn_mat = cn_mat.reindex(columns=keep)

    genes = list(cnv_mat.columns)

    fig = plt.figure(figsize=(max(14, (len(drugs) + len(genes)) * 0.5),
                              max(6, len(samples) * 0.18)))
    gs = gridspec.GridSpec(1, 4, width_ratios=[0.025, 1.0, 0.1, 1.4], wspace=0.05)

    # 1. Tira especie
    ax_sp = fig.add_subplot(gs[0, 0])
    sp_rgb = np.array([[mcolors.to_rgb(SP_COLOR.get(s, "#ccc")) for s in sp_list]]).transpose(1, 0, 2)
    ax_sp.imshow(sp_rgb, aspect="auto", interpolation="nearest")
    ax_sp.set_xticks([]); ax_sp.set_yticks([])

    # 2. AMR call panel
    ax_amr = fig.add_subplot(gs[0, 1])
    amr_cmap = mcolors.ListedColormap(["#999999", "#2ECC71", "#F39C12", "#C0392B"])
    bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, amr_cmap.N)
    ax_amr.imshow(np.ma.masked_invalid(amr_mat), aspect="auto",
                  cmap=amr_cmap, norm=norm, interpolation="nearest")
    for i in range(amr_mat.shape[0]):
        for j in range(amr_mat.shape[1]):
            if not np.isnan(amr_mat[i, j]):
                lab = {-1: "?", 0: "S", 1: "M", 2: "R"}.get(int(amr_mat[i, j]), "")
                ax_amr.text(j, i, lab, ha="center", va="center", color="white",
                            fontsize=7, fontweight="bold")
    ax_amr.set_xticks(range(len(drugs)))
    ax_amr.set_xticklabels(drugs, rotation=45, ha="right", fontsize=8)
    ax_amr.set_yticks(range(len(samples)))
    ax_amr.set_yticklabels(samples, fontsize=6)
    for tick, sp in zip(ax_amr.get_yticklabels(), sp_list):
        tick.set_color(SP_COLOR.get(sp, "#000000"))
    ax_amr.set_title("AMR (ChroQueTas)", fontsize=10)

    # 3. spacer (gs[0,2])
    # 4. CNV panel
    ax_cnv = fig.add_subplot(gs[0, 3])
    vmax = max(1.5, np.nanmax(np.abs(cnv_mat.values)) if cnv_mat.size else 1.5)
    im = ax_cnv.imshow(cnv_mat.values, aspect="auto", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, interpolation="nearest")
    # Cell annotation: estimated copy number and delta vs baseline ploidy.
    for i in range(len(samples)):
        for j in range(len(genes)):
            cell = cn_mat.iat[i, j] if (i < cn_mat.shape[0]
                                        and j < cn_mat.shape[1]) else np.nan
            txt = _fmt_cn_delta(cell, ploidy)
            if not txt:
                continue
            log2v = cnv_mat.iat[i, j]
            color = "white" if (not np.isnan(log2v) and abs(log2v) > vmax * 0.55) else "black"
            ax_cnv.text(j, i, txt, ha="center", va="center",
                        fontsize=6, color=color)
    ax_cnv.set_xticks(range(len(genes)))
    ax_cnv.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax_cnv.set_yticks([])
    ax_cnv.set_title("CNV (mosdepth log2 ratio + estimated cn)", fontsize=10)

    plt.colorbar(im, ax=ax_cnv, fraction=0.03, pad=0.02, label="log2 ratio")

    fig.suptitle("AMR × CNV combined view",
                 fontsize=13, y=0.98)
    out = Path(outdir) / "plot_amr_cnv_combined.png"
    plt.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.savefig(out.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"→ {out}")


# ============================================================================
# Plot 3 — Aneuploidía cromosómica (Chr5x2 et al.)
# ============================================================================
def plot_chr5_aneuploidy(contig_depth_files, species_of, outdir):
    """Espera que cnv_detect.py emita ALSO un 'contig_depth.tsv' con
       sample_id, contig, length_bp, median_depth, ratio_vs_genome.
       Esta funcionalidad hay que añadirla al python (ver prompt)."""
    rows = []
    for f in contig_depth_files:
        if not Path(f).exists():
            continue
        d = pd.read_csv(f, sep="\t")
        rows.append(d)
    if not rows:
        print("[plot_chr5_aneuploidy] no contig_depth files yet "
              "(extend cnv_detect.py first)")
        return
    df = pd.concat(rows, ignore_index=True)
    df["log2_ratio"] = np.log2(df["ratio_vs_genome"].replace(0, np.nan))

    # Filtrar muestras con al menos un contig con ratio > 1.5x (sospecha aneuploidía)
    suspect_samples = df[df["ratio_vs_genome"] > 1.5]["sample_id"].unique()
    if len(suspect_samples) == 0:
        print("[plot_chr5_aneuploidy] no samples with any contig dup > 1.5x")
        return

    sub = df[df["sample_id"].isin(suspect_samples)]
    n = len(suspect_samples)
    fig, axes = plt.subplots(n, 1, figsize=(10, max(3, n*1.5)),
                             sharex=False, squeeze=False)
    for i, sid in enumerate(suspect_samples):
        ax = axes[i, 0]
        s = sub[sub["sample_id"] == sid].sort_values("contig")
        colors = ["#C0392B" if r > 1.5 else "#7F8C8D" for r in s["ratio_vs_genome"]]
        ax.bar(range(len(s)), s["log2_ratio"], color=colors, edgecolor="black", lw=0.4)
        ax.axhline(0, color="black", lw=0.5)
        ax.axhline(np.log2(1.5), color="red", ls="--", lw=0.7, label="1.5×")
        ax.axhline(np.log2(2.0), color="darkred", ls="--", lw=0.7, label="2×")
        ax.set_xticks(range(len(s)))
        ax.set_xticklabels(s["contig"], rotation=45, ha="right", fontsize=7)
        sp = SP_PRETTY.get(species_of.get(sid, ""), "?")
        ax.set_title(f"{sid} ({sp})", fontsize=9, loc="left",
                     color=SP_COLOR.get(sp, "#000"))
        ax.set_ylabel("log2 ratio", fontsize=8)
        if i == 0:
            ax.legend(loc="upper right", fontsize=7)

    fig.suptitle(f"Aneuploidía sospechosa — {n} muestras con ≥1 contig > 1.5×",
                 fontsize=12)
    plt.tight_layout()
    out = Path(outdir) / "plot_chr5_aneuploidy.png"
    plt.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"→ {out}")


# ============================================================================
# Plot 4 — Firma Carolus 2021
# ============================================================================
def plot_carolus_signature(events, call, species_of, outdir):
    """Detecta cuántas muestras reproducen total/parcialmente el genotipo
       de Carolus 2021: ERG11 dup + TAC1B dup + FKS1 F635del (mut) +
       CIS2 A27T (mut). Como la mutación puntual no la mide CNV, aquí solo
       validamos las DUP de ERG11/TAC1B + posible DEL/dup en FKS1."""
    if events.empty:
        return

    pivot = events.pivot_table(index="sample_id", columns="gene",
                                values="norm", aggfunc="mean")
    target_genes = ["ERG11", "TAC1B", "FKS1"]
    pivot = pivot[[g for g in target_genes if g in pivot.columns]]
    if pivot.empty:
        print("[plot_carolus_signature] target genes ausentes en eventos")
        return

    # Clasifica cada celda en Dup (>=1.5), Normal, Del (<=0.5), No-cov
    def cat(v):
        if pd.isna(v): return "no_cov"
        if v >= 1.5:   return "dup"
        if v <= 0.5:   return "del"
        return "normal"

    cat_df = pivot.applymap(cat)
    # Score Carolus = nº de genes con DUP (max 2: ERG11, TAC1B)
    cat_df["carolus_dup_score"] = cat_df[["ERG11", "TAC1B"]].apply(
        lambda r: sum(v == "dup" for v in r), axis=1)
    cat_df = cat_df.sort_values("carolus_dup_score", ascending=False)

    n = len(cat_df)
    fig, ax = plt.subplots(figsize=(6, max(3, n*0.2)))
    cmap = {"dup": "#C0392B", "normal": "#ECF0F1",
            "del": "#2980B9", "no_cov": "#7F8C8D"}
    mat_rgb = np.array([[mcolors.to_rgb(cmap[v]) for v in cat_df.iloc[i, :3]]
                         for i in range(n)])
    ax.imshow(mat_rgb, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["ERG11", "TAC1B", "FKS1"], fontsize=10)
    ax.set_yticks(range(n))
    ax.set_yticklabels(cat_df.index, fontsize=6)
    for tick, sid in zip(ax.get_yticklabels(), cat_df.index):
        sp = SP_PRETTY.get(species_of.get(sid, ""), "?")
        tick.set_color(SP_COLOR.get(sp, "#000"))

    # Score badge a la derecha
    for i, score in enumerate(cat_df["carolus_dup_score"]):
        ax.text(2.6, i, f"·{int(score)}", fontsize=7, va="center",
                color="#C0392B" if score >= 1 else "#999")

    ax.set_title("Firma Carolus 2021 — ERG11+TAC1B dup co-amplicón\n"
                 "(score 2 = ambos dup; 1 = uno; 0 = ninguno)",
                 fontsize=10, pad=8)
    out = Path(outdir) / "plot_carolus_signature.png"
    plt.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"→ {out}")


# ============================================================================
# main
# ============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnv_events", nargs="+",
                    help="Globs a cnv_events.tsv (uno por muestra)")
    ap.add_argument("--amr_call", default=None,
                    help="call_matrix.tsv del módulo AMR")
    ap.add_argument("--species_id", default=None,
                    help="TSV con sample_id, top_genome_slug")
    ap.add_argument("--contig_depth", nargs="*", default=[],
                    help="Globs a *.contig_depth.tsv (extensión cnv_detect v2.0.1)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ploidy", type=int, default=1,
                    help="Baseline ploidy to compute the cn delta annotation "
                         "(default: 1 for C. auris).")
    # --panel is accepted for forward compatibility with the Nextflow module that
    # passes the AMR panel TSV; the script does not consume it yet but argparse
    # would otherwise abort with 'unrecognized arguments'.
    ap.add_argument("--panel", default=None,
                    help="Optional AMR panel TSV (currently unused; kept to "
                         "stay compatible with the Nextflow wrapper).")
    ap.add_argument("--show-all-genes", dest="show_all_genes", action="store_true",
                    help="Show every gene column even when no sample has a "
                         "non-`normal` CNV event.  Default is to hide such "
                         "empty columns so the figure stays compact.")
    args = ap.parse_args()
    hide_empty = not args.show_all_genes

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    events    = load_cnv_events(args.cnv_events or [])
    call      = load_call_matrix(args.amr_call) if args.amr_call else pd.DataFrame()
    species_of = load_species(args.species_id) if args.species_id else {}
    contigs    = [f for g in args.contig_depth for f in glob.glob(g)]

    print(f"[main] {len(events)} eventos CNV (de {events['sample_id'].nunique() if not events.empty else 0} muestras)")
    print(f"[main] {len(call)} muestras en call_matrix" if not call.empty else "[main] no call_matrix")
    print(f"[main] {len(species_of)} mapeos sample→clado")

    plot_cnv_heatmap(events, species_of, args.outdir, ploidy=args.ploidy,
                     hide_empty_genes=hide_empty)
    plot_amr_cnv_combined(events, call, species_of, args.outdir, ploidy=args.ploidy,
                          hide_empty_genes=hide_empty)
    plot_chr5_aneuploidy(contigs, species_of, args.outdir)
    plot_carolus_signature(events, call, species_of, args.outdir)


if __name__ == "__main__":
    main()
