#!/usr/bin/env python3
"""
Genera 4 JSON de MultiQC con scatter RAW vs TRIM mejorados
"""

import json
import math
from pathlib import Path
import pandas as pd


def build_scatter_json(df: pd.DataFrame, metric_key: str, raw_col: str, trim_col: str,
                       section_name: str, title: str, xlab: str, ylab: str, out_path: Path):
    
    # Calcular límites para el gráfico
    all_values = list(df[raw_col]) + list(df[trim_col])
    all_values = [v for v in all_values if pd.notnull(v) and v > 0]
    
    if all_values:
        min_val = min(all_values) * 0.9
        max_val = max(all_values) * 1.1
    else:
        min_val = 0
        max_val = 100

    points = []
    for _, row in df.iterrows():
        sample = str(row["sample"]) if "sample" in row else str(row.get("Sample", ""))
        raw_v = float(row[raw_col]) if pd.notnull(row[raw_col]) else 0.0
        trim_v = float(row[trim_col]) if pd.notnull(row[trim_col]) else 0.0

        points.append({
            "x": trim_v,
            "y": raw_v,
            "name": sample
        })

    # Configuración del gráfico mejorada
    payload = {
        "id": metric_key,
        "section_name": section_name,
        "description": f"Diagonal line (y=x) shows no change. Points above: lost in filtering. Points below: improved.",
        "plot_type": "scatter",
        "pconfig": {
            "id": f"{metric_key}_plot",
            "title": title,
            "xlab": xlab,
            "ylab": ylab,
            "xmin": min_val,
            "xmax": max_val,
            "ymin": min_val,
            "ymax": max_val,
            "square": True,  # Fuerza aspecto cuadrado
            "xlog": False,
            "ylog": False,
            "marker_size": 8,
            "marker_line_width": 1,
            "tt_label": "{point.name}<br/>Raw: {point.y:.0f}<br/>Filtered: {point.x:.0f}",
            # Línea diagonal y=x
            "extra_series": [{
                "name": "y=x line",
                "data": [
                    {"x": min_val, "y": min_val},
                    {"x": max_val, "y": max_val}
                ],
                "dash": "dash",
                "width": 2,
                "color": "#FF0000",
                "marker": {"enabled": False},
                "showInLegend": False
            }]
        },
        "data": {
            "Samples": points
        }
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def main():
    # Entradas / salidas vía snakemake
    tsv_path = Path(snakemake.input.tsv)
    out_reads = Path(snakemake.output.reads)
    out_mean_len = Path(snakemake.output.mean_len)
    out_mean_q = Path(snakemake.output.mean_q)
    out_n50 = Path(snakemake.output.n50)

    df = pd.read_csv(tsv_path, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    # 1) Número de lecturas
    build_scatter_json(
        df,
        metric_key="nanostats_trim_reads",
        raw_col="raw_number_of_reads",
        trim_col="trim_number_of_reads",
        section_name="Read Count Comparison",
        title="Read Count: Raw vs Filtered",
        xlab="Filtered Reads",
        ylab="Raw Reads",
        out_path=out_reads,
    )

    # 2) Longitud media
    build_scatter_json(
        df,
        metric_key="nanostats_trim_mean_len",
        raw_col="raw_mean_read_length",
        trim_col="trim_mean_read_length",
        section_name="Mean Length Comparison",
        title="Mean Read Length: Raw vs Filtered",
        xlab="Filtered Mean Length (bp)",
        ylab="Raw Mean Length (bp)",
        out_path=out_mean_len,
    )

    # 3) Calidad media
    build_scatter_json(
        df,
        metric_key="nanostats_trim_mean_q",
        raw_col="raw_mean_read_quality",
        trim_col="trim_mean_read_quality",
        section_name="Mean Quality Comparison",
        title="Mean Quality: Raw vs Filtered",
        xlab="Filtered Mean Quality (Q)",
        ylab="Raw Mean Quality (Q)",
        out_path=out_mean_q,
    )

    # 4) N50
    build_scatter_json(
        df,
        metric_key="nanostats_trim_n50",
        raw_col="raw_n50",
        trim_col="trim_n50",
        section_name="N50 Comparison",
        title="N50: Raw vs Filtered",
        xlab="Filtered N50 (bp)",
        ylab="Raw N50 (bp)",
        out_path=out_n50,
    )


if __name__ == "__main__":
    main()