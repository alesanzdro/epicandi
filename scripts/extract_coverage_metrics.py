#!/usr/bin/env python3
"""
Extract coverage metrics from samtools and mosdepth outputs
Generate TSV summary and JSON for MultiQC
"""

import json
import os
import sys
import logging
from pathlib import Path
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_coverage_file(coverage_file):
    """Parse samtools coverage output"""
    genome_length = 0
    weighted_coverage = 0
    
    try:
        with open(coverage_file, 'r') as f:
            # Skip header
            header = next(f)
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    contig_length = int(parts[2])
                    mean_depth = float(parts[6])
                    
                    genome_length += contig_length
                    weighted_coverage += mean_depth * contig_length
        
        mean_coverage = weighted_coverage / genome_length if genome_length > 0 else 0
        return genome_length, round(mean_coverage, 2)
    
    except Exception as e:
        logger.error(f"Error parsing {coverage_file}: {e}")
        return 0, 0

def parse_stats_file(stats_file):
    """Parse samtools stats output"""
    bases_mapped = 0
    reads_mapped = 0
    
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith('SN'):
                    if 'bases mapped:' in line:
                        bases_mapped = int(line.split('\t')[2])
                    elif 'reads mapped:' in line:
                        reads_mapped = int(line.split('\t')[2])
        
        return bases_mapped, reads_mapped
    
    except Exception as e:
        logger.error(f"Error parsing {stats_file}: {e}")
        return 0, 0

def parse_mosdepth_summary(mosdepth_file):
    """Parse mosdepth summary file for median coverage"""
    try:
        with open(mosdepth_file, 'r') as f:
            for line in f:
                if line.startswith('total'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        return float(parts[3])  # median coverage
        return None
    
    except Exception as e:
        logger.warning(f"Could not parse mosdepth file {mosdepth_file}: {e}")
        return None

def main():
    # Get parameters from snakemake
    samples = snakemake.params.samples
    output_dir = snakemake.params.output_dir
    
    # Initialize data structures
    metrics_data = []
    coverage_plot_data = {}
    
    logger.info(f"Processing {len(samples)} samples...")
    
    for sample in samples:
        logger.info(f"Processing sample: {sample}")
        
        # Define file paths
        coverage_file = Path(output_dir) / "02_assembly" / "02.6_coverage" / sample / f"{sample}.coverage.txt"
        stats_file = Path(output_dir) / "02_assembly" / "02.6_coverage" / sample / f"{sample}.stats.txt"
        mosdepth_file = Path(output_dir) / "02_assembly" / "02.7_mosdepth" / sample / f"{sample}.mosdepth.summary.txt"
        
        # Parse coverage data
        genome_length, mean_coverage = parse_coverage_file(coverage_file)
        
        # Parse stats data
        bases_mapped, reads_mapped = parse_stats_file(stats_file)
        
        # Parse mosdepth data
        median_coverage = parse_mosdepth_summary(mosdepth_file)
        if median_coverage is None:
            median_coverage = mean_coverage  # Fallback to mean if median not available
        
        # Add to metrics table
        metrics_data.append({
            'Sample': sample,
            'Genome_Length': genome_length,
            'Mean_Coverage': mean_coverage,
            'Median_Coverage': median_coverage,
            'Bases_Mapped': bases_mapped,
            'Reads_Mapped': reads_mapped
        })
        
        # Add to plot data (only if valid)
        if genome_length > 0 and mean_coverage > 0:
            genome_mb = round(genome_length / 1_000_000, 2)
            coverage_plot_data[sample] = {
                "x": genome_mb,
                "y": mean_coverage
            }
            logger.info(f"  - Genome: {genome_mb} Mb, Coverage: {mean_coverage}X")
    
    # Write TSV file
    logger.info("Writing TSV file...")
    df = pd.DataFrame(metrics_data)
    df.to_csv(snakemake.output.tsv, sep='\t', index=False)
    
    # Create MultiQC JSON structure
    logger.info("Creating MultiQC JSON...")
    multiqc_json = {
        "id": "assembly_coverage",
        "section_name": "Assembly Coverage Analysis",
        "description": "Coverage statistics for Oxford Nanopore assembled genomes",
        "plot_type": "scatter",
        "pconfig": {
            "id": "assembly_coverage_plot",
            "title": "Mean Coverage vs Genome Size",
            "xlab": "Genome Length (Mb)",
            "ylab": "Mean Coverage Depth (X)",
            "xlog": False,
            "ylog": False,
            "marker_size": 10,
            "marker_line_width": 1.5,
            "xmin": 0,
            "ymin": 0
        },
        "data": coverage_plot_data
    }
    
    # Write JSON file
    with open(snakemake.output.json, 'w') as f:
        json.dump(multiqc_json, f, indent=4)
    
    # Summary statistics
    if metrics_data:
        avg_genome = sum(m['Genome_Length'] for m in metrics_data) / len(metrics_data) / 1_000_000
        avg_coverage = sum(m['Mean_Coverage'] for m in metrics_data) / len(metrics_data)
        logger.info(f"\nSummary:")
        logger.info(f"  Average genome size: {avg_genome:.2f} Mb")
        logger.info(f"  Average coverage: {avg_coverage:.1f}X")
        logger.info(f"  Samples processed: {len(metrics_data)}")
    
    logger.info("Coverage metrics extraction completed successfully!")

if __name__ == "__main__":
    main()