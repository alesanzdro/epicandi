#!/usr/bin/env python3
"""
Generate sample metrics table for EpiCandi pipeline.
Extracts and processes QC metrics from fastp and NanoPlot outputs.
"""

import json
import re
import pandas as pd
from datetime import datetime
import os
import sys
from pathlib import Path
from typing import Dict, Any, List


def parse_fastp_metrics(fastp_file: str) -> Dict[str, Any]:
    """
    Parse fastp JSON output and extract key metrics.

    Args:
        fastp_file: Path to fastp JSON output file

    Returns:
        Dictionary with processed Illumina metrics
    """
    try:
        with open(fastp_file, 'r') as f:
            fastp_data = json.load(f)

            # Data before filtering
            raw_reads = fastp_data['summary']['before_filtering']['total_reads'] // 2

            # Data after filtering
            clean_reads = fastp_data['summary']['after_filtering']['total_reads'] // 2
            total_bases = fastp_data['summary']['after_filtering']['total_bases']
            q30_bases = fastp_data['summary']['after_filtering']['q30_bases']
            q30_rate = (q30_bases / total_bases) if total_bases > 0 else 0
            gc_content = fastp_data['summary']['after_filtering']['gc_content']
            keep_rate = clean_reads / raw_reads if raw_reads > 0 else 0

            # Mean length
            len_mean = total_bases / clean_reads / 2 if clean_reads > 0 else 0

            return {
                'Illumina_total_reads': clean_reads,
                'Illumina_len_mean': round(len_mean, 1),
                'Illumina_q30_bases': q30_bases,
                'Illumina_q30_rate': round(q30_rate, 6),
                'Illumina_gc_content': round(gc_content, 6),
                'Illumina_keep_rate': round(keep_rate, 6)
            }
    except Exception as e:
        print(f"Error parsing fastp file {fastp_file}: {e}")
        return {}


def parse_nanostats(nano_stats_file: str) -> Dict[str, Any]:
    """
    Parse NanoStats.txt file and extract key metrics.

    Args:
        nano_stats_file: Path to NanoStats.txt file

    Returns:
        Dictionary with processed Nanopore metrics
    """
    try:
        with open(nano_stats_file, 'r') as f:
            content = f.read()

            # Extract general metrics with regex
            reads_match = re.search(r'Number of reads:\s*([0-9,]+)', content)
            len_mean_match = re.search(r'Mean read length:\s*([0-9,.]+)', content)
            q15_match = re.search(r'>Q15:\s*([0-9,]+)', content)
            gc_match = re.search(r'Average GC:\s*([0-9.]+)', content)

            # Parse values
            total_reads = int(reads_match.group(1).replace(',', '')) if reads_match else 0
            len_mean = float(len_mean_match.group(1).replace(',', '')) if len_mean_match else 0
            q15_bases = int(q15_match.group(1).replace(',', '')) if q15_match else 0
            gc_content = float(gc_match.group(1)) / 100 if gc_match else 0

            return {
                'total_reads': total_reads,
                'len_mean': len_mean,
                'q15_bases': q15_bases,
                'gc_content': gc_content
            }
    except Exception as e:
        print(f"Error parsing NanoStats file {nano_stats_file}: {e}")
        return {}


def calculate_nanopore_metrics(filtered_stats: Dict, raw_stats: Dict) -> Dict[str, Any]:
    """
    Calculate Nanopore metrics from filtered and raw statistics.

    Args:
        filtered_stats: Metrics from filtered reads
        raw_stats: Metrics from raw reads

    Returns:
        Dictionary with calculated Nanopore metrics
    """
    if not filtered_stats:
        return {}

    # Calculate keep rate
    keep_rate = 0
    if raw_stats and raw_stats.get('total_reads', 0) > 0:
        keep_rate = filtered_stats['total_reads'] / raw_stats['total_reads']

    # Calculate Q15 rate
    q15_rate = 0
    total_bases = filtered_stats['total_reads'] * filtered_stats['len_mean']
    if total_bases > 0:
        q15_rate = filtered_stats['q15_bases'] / total_bases

    return {
        'Nanopore_total_reads': filtered_stats['total_reads'],
        'Nanopore_len_mean': round(filtered_stats['len_mean'], 1),
        'Nanopore_q15_bases': filtered_stats['q15_bases'],
        'Nanopore_q15_rate': round(q15_rate, 6),
        'Nanopore_gc_content': round(filtered_stats['gc_content'], 6),
        'Nanopore_keep_rate': round(keep_rate, 6)
    }


def determine_assembly_strategy(illumina_pass: bool, nanopore_pass: bool) -> str:
    """
    Determine the optimal assembly strategy based on QC results.

    Args:
        illumina_pass: Whether Illumina data passes QC
        nanopore_pass: Whether Nanopore data passes QC

    Returns:
        Assembly strategy string
    """
    if illumina_pass and nanopore_pass:
        return "hybrid"
    elif illumina_pass:
        return "illumina"
    elif nanopore_pass:
        return "nanopore"
    else:
        return "fail"


def generate_sample_metrics(samples_dict: Dict, qc_thresholds: Dict,
                          output_file: str, log_file: str) -> None:
    """
    Generate comprehensive sample metrics table.

    Args:
        samples_dict: Dictionary with sample information
        qc_thresholds: QC thresholds configuration
        output_file: Output TSV file path
        log_file: Log file path
    """
    # Get thresholds
    illumina_min_reads = qc_thresholds.get("illumina_min_reads", 50000)
    nanopore_min_reads = qc_thresholds.get("nanopore_min_reads", 5000)

    # Create log directory
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    with open(log_file, "w") as log:
        log.write("=== GENERATING SAMPLE METRICS TABLE ===\n")
        log.write(f"Date: {datetime.now()}\n\n")

        # List to store data for each sample
        metrics_data = []

        # Process each sample
        for sample_id, data in samples_dict.items():
            sample_type = data['type']

            # Default values
            metrics = {
                'Sample': sample_id,
                'target_assembly': sample_type,  # Default to sample type
                'final_assembly': None,  # Will be determined later
                'Nanopore_total_reads': 'NA',
                'Nanopore_len_mean': 'NA',
                'Nanopore_q15_bases': 'NA',
                'Nanopore_q15_rate': 'NA',
                'Nanopore_gc_content': 'NA',
                'Nanopore_keep_rate': 'NA',
                'Illumina_total_reads': 'NA',
                'Illumina_len_mean': 'NA',
                'Illumina_q30_bases': 'NA',
                'Illumina_q30_rate': 'NA',
                'Illumina_gc_content': 'NA',
                'Illumina_keep_rate': 'NA',
                'Nanopore_sample': 'NA',
                'Illumina_sample': 'NA'
            }

            # Variables to track QC pass/fail
            illumina_pass = False
            nanopore_pass = False

            # Process Illumina metrics
            if sample_type in ['illumina', 'hybrid']:
                fastp_file = f"output/01_data/01.2_filtered/fastp/{sample_id}_fastp.json"

                if os.path.exists(fastp_file):
                    illumina_metrics = parse_fastp_metrics(fastp_file)
                    if illumina_metrics:
                        metrics.update(illumina_metrics)

                        # Determine if passes QC
                        illumina_pass = metrics['Illumina_total_reads'] >= illumina_min_reads
                        metrics['Illumina_sample'] = 'pass' if illumina_pass else 'fail'

                        log.write(f"{sample_id} (Illumina): {metrics['Illumina_total_reads']:,} reads - "
                                f"{'PASS' if illumina_pass else 'FAIL'}\n")
                else:
                    log.write(f"Warning: fastp file not found for {sample_id}: {fastp_file}\n")

            # Process Nanopore metrics
            if sample_type in ['nanopore', 'hybrid']:
                nano_stats_file = f"output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample_id}/NanoStats.txt"
                nano_raw_file = f"output/01_data/01.1_fastq_raw_qc/nanoplot/{sample_id}/NanoStats.txt"

                if os.path.exists(nano_stats_file):
                    filtered_stats = parse_nanostats(nano_stats_file)
                    raw_stats = parse_nanostats(nano_raw_file) if os.path.exists(nano_raw_file) else {}

                    if filtered_stats:
                        nanopore_metrics = calculate_nanopore_metrics(filtered_stats, raw_stats)
                        metrics.update(nanopore_metrics)

                        # Determine if passes QC
                        nanopore_pass = metrics['Nanopore_total_reads'] >= nanopore_min_reads
                        metrics['Nanopore_sample'] = 'pass' if nanopore_pass else 'fail'

                        log.write(f"{sample_id} (Nanopore): {metrics['Nanopore_total_reads']:,} reads - "
                                f"{'PASS' if nanopore_pass else 'FAIL'}\n")
                else:
                    log.write(f"Warning: NanoStats file not found for {sample_id}: {nano_stats_file}\n")

            # Determine final assembly strategy
            final_assembly = determine_assembly_strategy(illumina_pass, nanopore_pass)
            metrics['final_assembly'] = final_assembly

            log.write(f"{sample_id}: Final assembly strategy = {final_assembly}\n\n")

            metrics_data.append(metrics)

        # Create DataFrame and save
        df = pd.DataFrame(metrics_data)

        # Create output directory
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Save as TSV
        df.to_csv(output_file, sep='\t', index=False)

        log.write(f"Metrics table saved to: {output_file}\n")
        log.write(f"Total samples processed: {len(metrics_data)}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate sample metrics table")
    parser.add_argument("--samples-dict", required=True, help="Path to samples dictionary JSON")
    parser.add_argument("--qc-thresholds", required=True, help="Path to QC thresholds JSON")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--log", required=True, help="Log file")

    args = parser.parse_args()

    # Load samples dictionary
    with open(args.samples_dict, 'r') as f:
        samples_dict = json.load(f)

    # Load QC thresholds
    with open(args.qc_thresholds, 'r') as f:
        qc_thresholds = json.load(f)

    generate_sample_metrics(samples_dict, qc_thresholds, args.output, args.log)