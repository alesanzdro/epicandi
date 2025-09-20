#!/usr/bin/env python3
"""Genera métricas de BUSCO para MultiQC.

Lee los ficheros normalizados en la colección:
    output/04_report/04.1_collection_renamed/busco/{sample}.busco_summary.txt
que contienen una línea tipo:
    C:96.3%[S:96.2%,D:0.1%],F:0.1%,M:3.6%,n:758

Modos de salida:
    - Por defecto (sin flags): genera JSON general_stats en
                output/.../busco/busco_multiqc.json
    - Con --tsv: genera un TSV con cabecera comentada (custom content 'table') en
                output/.../busco/busco_summary_mqc.tsv
"""
import re
import json
import os
from pathlib import Path

import sys
import argparse

COLLECTION_DIR = Path('output/04_report/04.1_collection_renamed/busco')
OUTPUT_JSON = Path('output/04_report/04.1_collection_renamed/busco/busco_multiqc.json')
OUTPUT_TSV = Path('output/04_report/04.1_collection_renamed/busco/busco_summary_mqc.tsv')

summary_pattern = re.compile(r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>\d+)")

metrics = {}

if not COLLECTION_DIR.exists():
    print(f"BUSCO collection dir not found: {COLLECTION_DIR}", file=sys.stderr)
    sys.exit(0)

for file in COLLECTION_DIR.glob('*.busco_summary.txt'):
    sample = file.name.replace('.busco_summary.txt','')
    try:
        text = file.read_text()
        match = summary_pattern.search(text)
        if not match:
            print(f"WARN: No summary pattern found in {file}")
            continue
        data = {k: float(v) for k,v in match.groupdict().items() if k in ['C','S','D','F','M']}
        metrics[sample] = {
            'BUSCO_C': data['C'],
            'BUSCO_S': data['S'],
            'BUSCO_D': data['D'],
            'BUSCO_F': data['F'],
            'BUSCO_M': data['M']
        }
    except Exception as e:
        print(f"ERROR parsing {file}: {e}")

def write_json(metrics: dict):
    multiqc_json = {
        "id": "busco_general_stats",
        "section_name": "BUSCO General Stats",
        "description": "BUSCO completeness metrics (percentages)",
        "plot_type": "general_stats",
        "data": {s: {k+"(%)": v for k, v in vals.items()} for s, vals in metrics.items()}
    }
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_JSON, 'w') as fh:
        json.dump(multiqc_json, fh, indent=2)
    print(f"Written BUSCO MultiQC JSON: {OUTPUT_JSON}")


def write_tsv(metrics: dict):
    """Escribe TSV con cabecera comentada para MultiQC custom content (table)."""
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_TSV, 'w') as fh:
        # Cabecera comentada exactamente como se solicitó
        fh.write("# plot_type: 'table'\n")
        fh.write("# section_name: 'BUSCO Summary'\n")
        fh.write("# description: 'Completeness assessment with BUSCO'\n")
        fh.write("# headers:\n")
        fh.write("#     BUSCO_C:\n")
        fh.write("#         title: 'Complete (%)'\n")
        fh.write("#         format: '{:.1f}'\n")
        fh.write("#         suffix: '%'\n")
        fh.write("#     BUSCO_S:\n")
        fh.write("#         title: 'Single-copy (%)'\n")
        fh.write("#         format: '{:.1f}'\n")
        fh.write("#         suffix: '%'\n")
        fh.write("#     BUSCO_D:\n")
        fh.write("#         title: 'Duplicated (%)'\n")
        fh.write("#         format: '{:.1f}'\n")
        fh.write("#         suffix: '%'\n")
        fh.write("#     BUSCO_F:\n")
        fh.write("#         title: 'Fragmented (%)'\n")
        fh.write("#         format: '{:.1f}'\n")
        fh.write("#         suffix: '%'\n")
        fh.write("#         max: 5\n")
        fh.write("#     BUSCO_M:\n")
        fh.write("#         title: 'Missing (%)'\n")
        fh.write("#         format: '{:.1f}'\n")
        fh.write("#         suffix: '%'\n")
        fh.write("#         max: 5\n")
        # Cabecera de tabla
        fh.write("Sample\tBUSCO_C\tBUSCO_S\tBUSCO_D\tBUSCO_F\tBUSCO_M\n")
        for sample in sorted(metrics.keys()):
            vals = metrics[sample]
            fh.write(
                f"{sample}\t{vals['BUSCO_C']:.1f}\t{vals['BUSCO_S']:.1f}\t{vals['BUSCO_D']:.1f}\t{vals['BUSCO_F']:.1f}\t{vals['BUSCO_M']:.1f}\n"
            )
    print(f"Written BUSCO MultiQC TSV: {OUTPUT_TSV}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv', action='store_true', help='Escribir TSV en lugar de JSON')
    args = parser.parse_args()

    if args.tsv:
        write_tsv(metrics)
    else:
        write_json(metrics)


if __name__ == '__main__':
    main()
