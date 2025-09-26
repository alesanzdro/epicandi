#!/usr/bin/env python3
"""
Sample discovery utilities for EpiCandi pipeline.
Extracts sample identification and path resolution logic from Snakefile.
"""

import os
import glob
from typing import List


def get_samples(input_dir: str = "input") -> List[str]:
    """
    Get list of samples from input directory, excluding control samples.

    Args:
        input_dir: Directory to search for input files

    Returns:
        List of sample names found in input directory
    """
    samples = []

    # Tags to exclude (controls and blanks)
    exclude_tags = [
        'negative', 'Negative', 'NEGATIVE',
        'control', 'Control', 'CONTROL',
        'blanco', 'Blanco', 'BLANCO',
        'unclassified', 'Unclassified', 'UNCLASSIFIED',
        'negativo', 'Negativo', 'NEGATIVO'
    ]

    # Search patterns including subdirectories
    patterns = [
        f"{input_dir}/*.fastq.gz",
        f"{input_dir}/*.fq.gz",
        f"{input_dir}/*/*.fastq.gz",
        f"{input_dir}/*/*.fq.gz"
    ]

    for pattern in patterns:
        for file_path in glob.glob(pattern):
            basename = os.path.basename(file_path)
            # Remove extensions
            sample = basename.replace('.fastq.gz', '').replace('.fq.gz', '')

            # Exclude control samples
            if any(tag in sample for tag in exclude_tags):
                print(f"⚠️  Excluding control sample: {sample}")
                continue

            if sample not in samples:
                samples.append(sample)
                print(f"✓ Sample added: {sample}")

    return samples


def get_sample_path(sample_name: str, input_dir: str = "input") -> str:
    """
    Find the complete path of a sample file.

    Searches first in input/ and immediate subdirectories with patterns:
      input/{sample}.fastq.gz
      input/{sample}.fq.gz
      input/*/{sample}.fastq.gz
      input/*/{sample}.fq.gz

    Args:
        sample_name: Name of the sample to find
        input_dir: Directory to search in

    Returns:
        Absolute path to the first matching file

    Raises:
        FileNotFoundError: If no matching file is found
    """
    candidate_patterns = [
        f"{input_dir}/{sample_name}.fastq.gz",
        f"{input_dir}/{sample_name}.fq.gz",
        f"{input_dir}/*/{sample_name}.fastq.gz",
        f"{input_dir}/*/{sample_name}.fq.gz",
    ]

    for pattern in candidate_patterns:
        matches = glob.glob(pattern)
        if matches:
            # Prefer shortest match (less ambiguous) if multiple
            matches.sort(key=len)
            return matches[0]

    raise FileNotFoundError(f"No FASTQ/FQ file found for sample: {sample_name}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Discover samples in input directory")
    parser.add_argument("--input-dir", default="input", help="Input directory to search")
    parser.add_argument("--list-samples", action="store_true",
                       help="List all discovered samples")
    parser.add_argument("--get-path", help="Get path for specific sample")

    args = parser.parse_args()

    if args.list_samples:
        samples = get_samples(args.input_dir)
        print(f"\nFound {len(samples)} samples:")
        for sample in samples:
            print(f"  {sample}")

    if args.get_path:
        try:
            path = get_sample_path(args.get_path, args.input_dir)
            print(f"Path for {args.get_path}: {path}")
        except FileNotFoundError as e:
            print(f"Error: {e}")