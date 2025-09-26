import os
import re
import csv
from pathlib import Path

# Campos a extraer de NanoStats
FIELDS = [
    ("Mean read length", "mean_read_length"),
    ("Mean read quality", "mean_read_quality"),
    ("Median read length", "median_read_length"),
    ("Median read quality", "median_read_quality"),
    ("Number of reads", "number_of_reads"),
    ("Read length N50", "n50"),
]

number_re = re.compile(r"([-+]?[0-9]*\.?[0-9]+(?:,[0-9]{3})*)")


def parse_nanostats(path: str) -> dict:
    res = {}
    if not os.path.exists(path):
        return res
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            for label, key in FIELDS:
                if line.startswith(label + ":"):
                    # Extraer primer número en la línea
                    m = number_re.search(line)
                    if m:
                        val = m.group(1).replace(",", "")
                        try:
                            res[key] = float(val)
                        except ValueError:
                            pass
    return res


def safe_ratio(trim: float | int | None, raw: float | int | None):
    if raw in (None, 0):
        return None
    if trim is None:
        return None
    return trim / raw


base_collection = snakemake.params.collection_dir  # output/04_report/04.1_collection_renamed
output_tsv = snakemake.output[0]
logfile = snakemake.log[0]
samples = snakemake.params.samples

# Rutas de entrada (raw/trim) producidas por NanoPlot
raw_fmt = os.path.join(snakemake.params.output_dir, "01_data/01.2_nanoplot_raw/{sample}/NanoStats.txt")
trim_fmt = os.path.join(snakemake.params.output_dir, "01_data/01.4_nanoplot_filtered/{sample}/NanoStats.txt")

os.makedirs(os.path.dirname(output_tsv), exist_ok=True)

headers = [
    "sample",
    # raw
    "raw_mean_read_length", "raw_mean_read_quality", "raw_median_read_length",
    "raw_median_read_quality", "raw_number_of_reads", "raw_n50",
    # trim
    "trim_mean_read_length", "trim_mean_read_quality", "trim_median_read_length",
    "trim_median_read_quality", "trim_number_of_reads", "trim_n50",
    # ratios trim/raw
    "ratio_mean_read_length", "ratio_mean_read_quality", "ratio_median_read_length",
    "ratio_median_read_quality", "ratio_number_of_reads", "ratio_n50",
]

with open(output_tsv, 'w', newline='') as out, open(logfile, 'w') as log:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(headers)

    for s in samples:
        raw_stats = parse_nanostats(raw_fmt.format(sample=s))
        trim_stats = parse_nanostats(trim_fmt.format(sample=s))

        row = [s]
        # raw values
        row.extend([
            raw_stats.get("mean_read_length"),
            raw_stats.get("mean_read_quality"),
            raw_stats.get("median_read_length"),
            raw_stats.get("median_read_quality"),
            raw_stats.get("number_of_reads"),
            raw_stats.get("n50"),
        ])
        # trim values
        row.extend([
            trim_stats.get("mean_read_length"),
            trim_stats.get("mean_read_quality"),
            trim_stats.get("median_read_length"),
            trim_stats.get("median_read_quality"),
            trim_stats.get("number_of_reads"),
            trim_stats.get("n50"),
        ])
        # ratios trim/raw
        row.extend([
            safe_ratio(trim_stats.get("mean_read_length"), raw_stats.get("mean_read_length")),
            safe_ratio(trim_stats.get("mean_read_quality"), raw_stats.get("mean_read_quality")),
            safe_ratio(trim_stats.get("median_read_length"), raw_stats.get("median_read_length")),
            safe_ratio(trim_stats.get("median_read_quality"), raw_stats.get("median_read_quality")),
            safe_ratio(trim_stats.get("number_of_reads"), raw_stats.get("number_of_reads")),
            safe_ratio(trim_stats.get("n50"), raw_stats.get("n50")),
        ])

        writer.writerow(row)
        log.write(f"Processed {s}: raw={raw_stats} trim={trim_stats}\n")
