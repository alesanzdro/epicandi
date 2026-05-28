#!/usr/bin/env python3
"""
plot_amr_heatmap.py — clinical heatmap of AMR calls across samples × drugs.

Input:  call_matrix.tsv from build_master_table.py (rows = samples,
        columns = drugs, values = R/I/S).

Output: amr_heatmap.png (300 dpi) + amr_heatmap.svg (vector).

Color map:
   R (resistant)   → red
   I (intermediate)→ amber
   S (susceptible) → green
   empty / NA      → grey
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Disable mathtext: sample IDs may contain underscores
plt.rcParams["text.parse_math"] = False

COLORS = {
    "R": "#C0392B",   # red
    "I": "#E67E22",   # amber
    "S": "#27AE60",   # green
    "":  "#BDC3C7",   # grey
    "NA":"#BDC3C7",
}
NUMERIC = {"R": 2, "I": 1, "S": 0, "": -1, "NA": -1}


def load_call_matrix(path: str):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = [r for r in reader if r]
    drugs = header[1:]
    samples = [r[0] for r in rows]
    matrix = [[r[i + 1] if i + 1 < len(r) else "" for i in range(len(drugs))] for r in rows]
    return samples, drugs, matrix


def plot(samples, drugs, matrix, output_prefix: str):
    n_s, n_d = len(samples), len(drugs)
    # Figure size: grow with N samples (rows) but cap drug-side width
    fig_w = max(6, 0.5 * n_d + 4)
    fig_h = max(3, 0.25 * n_s + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)

    # Build color array
    col_array = np.zeros((n_s, n_d, 3))
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            hexcol = COLORS.get(val.strip().upper() if val else "", COLORS[""])
            r, g, b = int(hexcol[1:3], 16), int(hexcol[3:5], 16), int(hexcol[5:7], 16)
            col_array[i, j] = [r / 255, g / 255, b / 255]

    ax.imshow(col_array, aspect="auto", interpolation="nearest")

    # Cell labels (R/I/S in white text, grey blank)
    for i in range(n_s):
        for j in range(n_d):
            val = matrix[i][j].strip().upper() if matrix[i][j] else ""
            if val in ("R", "I", "S"):
                ax.text(j, i, val, ha="center", va="center",
                        color="white", fontsize=8, fontweight="bold")

    ax.set_xticks(range(n_d))
    ax.set_xticklabels(drugs, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_s))
    ax.set_yticklabels(samples, fontsize=8)
    ax.set_xlabel("Antifungal drug")
    ax.set_title(f"AMR call matrix · {n_s} samples × {n_d} drugs",
                 loc="left", fontsize=11, fontweight="bold")

    # Grid between cells
    ax.set_xticks(np.arange(-0.5, n_d, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_s, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)

    # Legend
    legend_handles = [
        mpatches.Patch(color=COLORS["R"], label="R (resistant)"),
        mpatches.Patch(color=COLORS["I"], label="I (intermediate)"),
        mpatches.Patch(color=COLORS["S"], label="S (susceptible)"),
        mpatches.Patch(color=COLORS[""],  label="no data"),
    ]
    ax.legend(handles=legend_handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.15 - 0.001 * n_s),
              ncol=4, fontsize=8, frameon=False)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{output_prefix}.svg",          bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"[amr_heatmap] wrote {output_prefix}.png + .svg "
          f"({n_s} samples × {n_d} drugs)", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--call_matrix", required=True,
                    help="call_matrix.tsv (sample_id\\tdrug1\\tdrug2\\t...)")
    ap.add_argument("--output_prefix", required=True,
                    help="output path prefix (writes .png and .svg)")
    args = ap.parse_args()

    samples, drugs, matrix = load_call_matrix(args.call_matrix)
    if not samples:
        sys.exit(f"[amr_heatmap] empty call_matrix: {args.call_matrix}")
    plot(samples, drugs, matrix, args.output_prefix)


if __name__ == "__main__":
    main()
