#!/usr/bin/env python3
"""
Sample QC validation script for EpiCandi pipeline.
Processes metrics table to generate validation reports and filtered sample lists.
"""

import pandas as pd
from datetime import datetime
import os
import sys
from typing import Dict, List, Tuple


def generate_validation_report(metrics_df: pd.DataFrame, report_file: str) -> None:
    """
    Generate detailed validation report from metrics DataFrame.

    Args:
        metrics_df: DataFrame with sample metrics
        report_file: Output report file path
    """
    os.makedirs(os.path.dirname(report_file), exist_ok=True)

    with open(report_file, "w") as f:
        f.write("=== SAMPLE VALIDATION REPORT ===\n")
        f.write(f"Date: {datetime.now()}\n\n")

        # Assembly strategy summary
        hybrid_count = len(metrics_df[metrics_df['final_assembly'] == 'hybrid'])
        illumina_count = len(metrics_df[metrics_df['final_assembly'] == 'illumina'])
        nanopore_count = len(metrics_df[metrics_df['final_assembly'] == 'nanopore'])
        failed_count = len(metrics_df[metrics_df['final_assembly'] == 'fail'])

        f.write("ASSEMBLY STRATEGY SUMMARY:\n")
        f.write(f"Total samples: {len(metrics_df)}\n")
        f.write(f"Hybrid assembly: {hybrid_count}\n")
        f.write(f"Illumina-only assembly: {illumina_count}\n")
        f.write(f"Nanopore-only assembly: {nanopore_count}\n")
        f.write(f"Failed samples: {failed_count}\n\n")

        f.write("DETAIL BY SAMPLE:\n")
        for _, row in metrics_df.iterrows():
            sample_id = row['Sample']
            target = row['target_assembly']
            final = row['final_assembly']

            f.write(f"\n{sample_id}:\n")
            f.write(f"  Target strategy: {target}\n")
            f.write(f"  Final strategy: {final}\n")

            # Illumina details if available
            if pd.notna(row['Illumina_sample']) and row['Illumina_sample'] != 'NA':
                ill_status = str(row['Illumina_sample'])
                ill_reads = row['Illumina_total_reads']
                if pd.notna(ill_reads) and ill_reads != 'NA':
                    f.write(f"  Illumina: {ill_status.upper()} ({ill_reads:,} reads)\n")
                else:
                    f.write(f"  Illumina: {ill_status.upper()}\n")

            # Nanopore details if available
            if pd.notna(row['Nanopore_sample']) and row['Nanopore_sample'] != 'NA':
                nano_status = str(row['Nanopore_sample'])
                nano_reads = row['Nanopore_total_reads']
                if pd.notna(nano_reads) and nano_reads != 'NA':
                    f.write(f"  Nanopore: {nano_status.upper()} ({nano_reads:,} reads)\n")
                else:
                    f.write(f"  Nanopore: {nano_status.upper()}\n")

            # Indicate if strategy changed
            if target != final and final != 'fail':
                f.write(f"  Note: Strategy changed from {target} to {final} based on QC\n")


def create_valid_samples_csv(metrics_df: pd.DataFrame, original_csv_path: str,
                           output_csv_path: str) -> List[str]:
    """
    Create filtered samplesinfo.csv with only valid samples.

    Args:
        metrics_df: DataFrame with sample metrics
        original_csv_path: Path to original samplesinfo.csv
        output_csv_path: Path for filtered CSV output

    Returns:
        List of valid sample IDs
    """
    # Define valid assemblies
    valid_assemblies = ['hybrid', 'illumina', 'nanopore']
    valid_samples_df = metrics_df[metrics_df['final_assembly'].isin(valid_assemblies)]

    # Load original samplesinfo.csv
    original_csv = pd.read_csv(original_csv_path)

    # Filter by valid samples
    valid_ids = set(valid_samples_df['Sample'].tolist())
    filtered_csv = original_csv[original_csv['id'].isin(valid_ids)]

    # Create output directory
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)

    # Save filtered CSV
    filtered_csv.to_csv(output_csv_path, index=False)

    return list(valid_ids)


def log_validation_summary(metrics_df: pd.DataFrame, log_file: str) -> None:
    """
    Generate validation summary log.

    Args:
        metrics_df: DataFrame with sample metrics
        log_file: Log file path
    """
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    with open(log_file, "w") as log:
        log.write("=== SAMPLE VALIDATION ===\n")
        log.write(f"Date: {datetime.now()}\n\n")

        # Load metrics table
        passed_samples = metrics_df[metrics_df['final_assembly'] != 'fail']
        failed_samples = metrics_df[metrics_df['final_assembly'] == 'fail']

        log.write(f"Total samples: {len(metrics_df)}\n")
        log.write(f"Samples passed: {len(passed_samples)}\n")
        log.write(f"Samples failed: {len(failed_samples)}\n\n")

        log.write("Passed samples:\n")
        for _, row in passed_samples.iterrows():
            sample_id = row['Sample']
            target = row['target_assembly']
            final = row['final_assembly']
            log.write(f"{sample_id}: {target} → {final}\n")

        log.write("\nFailed QC samples:\n")
        for _, row in failed_samples.iterrows():
            sample_id = row['Sample']
            target = row['target_assembly']
            log.write(f"{sample_id}: {target} → FAIL\n")


def validate_samples(metrics_table: str, original_csv: str,
                    output_report: str, output_csv: str, log_file: str) -> Dict[str, any]:
    """
    Main validation function that processes all outputs.

    Args:
        metrics_table: Path to input metrics table (TSV)
        original_csv: Path to original samplesinfo.csv
        output_report: Path for validation report
        output_csv: Path for filtered CSV
        log_file: Path for log file

    Returns:
        Dictionary with validation results
    """
    # Load metrics table
    metrics_df = pd.read_csv(metrics_table, sep='\t')

    # Generate log summary
    log_validation_summary(metrics_df, log_file)

    # Generate detailed report
    generate_validation_report(metrics_df, output_report)

    # Create filtered CSV
    valid_ids = create_valid_samples_csv(metrics_df, original_csv, output_csv)

    # Calculate summary statistics
    total_samples = len(metrics_df)
    valid_samples = len(valid_ids)
    failed_samples = total_samples - valid_samples

    return {
        'total_samples': total_samples,
        'valid_samples': valid_samples,
        'failed_samples': failed_samples,
        'valid_sample_ids': valid_ids,
        'assembly_strategies': {
            'hybrid': len(metrics_df[metrics_df['final_assembly'] == 'hybrid']),
            'illumina': len(metrics_df[metrics_df['final_assembly'] == 'illumina']),
            'nanopore': len(metrics_df[metrics_df['final_assembly'] == 'nanopore']),
            'failed': len(metrics_df[metrics_df['final_assembly'] == 'fail'])
        }
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Validate samples based on QC metrics")
    parser.add_argument("--metrics-table", required=True, help="Input metrics table (TSV)")
    parser.add_argument("--original-csv", required=True, help="Original samplesinfo.csv")
    parser.add_argument("--output-report", required=True, help="Output validation report")
    parser.add_argument("--output-csv", required=True, help="Output filtered CSV")
    parser.add_argument("--log", required=True, help="Log file")

    args = parser.parse_args()

    results = validate_samples(
        args.metrics_table,
        args.original_csv,
        args.output_report,
        args.output_csv,
        args.log
    )

    print(f"Validation completed: {results['valid_samples']}/{results['total_samples']} samples passed QC")