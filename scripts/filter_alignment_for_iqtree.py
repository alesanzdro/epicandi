#!/usr/bin/env python3
"""
Filter concatenated alignment to ensure sufficient non-gap characters for IQ-TREE.

Usage:
    python filter_alignment_for_iqtree.py <input_alignment> <partition_file> <output_alignment> <output_partition>
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict
import logging


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('filter_alignment_for_iqtree.log')
        ]
    )
    return logging.getLogger(__name__)


def parse_partition_file(partition_file):
    """
    Parse IQ-TREE partition file.

    Args:
        partition_file (Path): Path to partition file

    Returns:
        list: List of (partition_name, start, end) tuples
    """
    partitions = []

    with open(partition_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                if '=' in line:
                    parts = line.split('=')
                    if len(parts) == 2:
                        left_part = parts[0].strip()
                        range_part = parts[1].strip()

                        # Check if format includes model (e.g., "DNA, partition_name")
                        if ',' in left_part:
                            partition_name = left_part.split(',')[1].strip()
                        else:
                            partition_name = left_part

                        if '-' in range_part:
                            start, end = map(int, range_part.split('-'))
                            partitions.append((partition_name, start, end))

    return partitions


def parse_fasta_alignment(fasta_file):
    """
    Parse FASTA alignment file efficiently.

    Args:
        fasta_file (Path): Path to FASTA file

    Returns:
        dict: Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)

                # Start new sequence
                current_id = line[1:]  # Remove '>' prefix
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def analyze_sequence_characters(sequences, partitions, min_char_fraction=0.3):
    """
    Analyze character content in sequences and identify sequences to exclude.

    Args:
        sequences (dict): Dictionary of sequences
        partitions (list): List of partition tuples
        min_char_fraction (float): Minimum fraction of non-gap characters required

    Returns:
        tuple: (valid_sequences, excluded_sequences, partition_stats)
    """
    logger = logging.getLogger(__name__)

    valid_sequences = set(sequences.keys())
    excluded_sequences = []
    partition_stats = []

    logger.info(f"Analyzing {len(partitions)} partitions...")

    for i, (partition_name, start, end) in enumerate(partitions):
        if i % 100 == 0:
            logger.info(f"Processing partition {i+1}/{len(partitions)}: {partition_name}")

        partition_length = end - start + 1
        min_chars_required = int(partition_length * min_char_fraction)

        partition_char_counts = {}

        for seq_id in list(valid_sequences):
            if seq_id in sequences:
                sequence = sequences[seq_id]
                if len(sequence) >= end:
                    # Extract partition sequence (convert to 0-based indexing)
                    partition_seq = sequence[start-1:end]
                    non_gap_chars = len(partition_seq) - partition_seq.count('-')
                    partition_char_counts[seq_id] = non_gap_chars

                    if non_gap_chars < min_chars_required:
                        if seq_id not in excluded_sequences:
                            excluded_sequences.append(seq_id)
                            logger.warning(f"Excluding {seq_id}: only {non_gap_chars}/{partition_length} "
                                         f"characters in partition {partition_name}")

        # Update valid sequences
        valid_sequences = valid_sequences - set(excluded_sequences)

        # Store partition statistics
        if partition_char_counts:
            avg_chars = sum(partition_char_counts.values()) / len(partition_char_counts)
            min_chars = min(partition_char_counts.values())
            max_chars = max(partition_char_counts.values())

            partition_stats.append({
                'name': partition_name,
                'length': partition_length,
                'min_required': min_chars_required,
                'avg_chars': avg_chars,
                'min_chars': min_chars,
                'max_chars': max_chars,
                'num_sequences': len(partition_char_counts)
            })

    return valid_sequences, excluded_sequences, partition_stats


def filter_sequences_by_overall_quality(sequences, min_total_chars=500000):
    """
    Filter sequences by overall character content.

    Args:
        sequences (dict): Dictionary of sequences
        min_total_chars (int): Minimum total non-gap characters

    Returns:
        tuple: (valid_sequences, excluded_sequences)
    """
    logger = logging.getLogger(__name__)

    valid_sequences = []
    excluded_sequences = []

    for seq_id, sequence in sequences.items():
        total_chars = len(sequence) - sequence.count('-')
        if total_chars >= min_total_chars:
            valid_sequences.append(seq_id)
        else:
            excluded_sequences.append(seq_id)
            logger.warning(f"Excluding {seq_id}: only {total_chars} total characters "
                         f"(minimum {min_total_chars} required)")

    return valid_sequences, excluded_sequences


def write_filtered_alignment(sequences, valid_sequences, output_file):
    """
    Write filtered sequences to FASTA format.

    Args:
        sequences (dict): Dictionary of all sequences
        valid_sequences (set): Set of valid sequence IDs
        output_file (Path): Output file path
    """
    logger = logging.getLogger(__name__)

    with open(output_file, 'w') as f:
        for seq_id in sorted(valid_sequences):
            if seq_id in sequences:
                sequence = sequences[seq_id]
                f.write(f">{seq_id}\n")
                # Write sequence in lines of 80 characters
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')

    logger.info(f"Written {len(valid_sequences)} sequences to {output_file}")


def write_filtered_partitions(partitions, output_file):
    """
    Write partition file with original coordinates.

    Args:
        partitions (list): List of partition tuples
        output_file (Path): Output file path
    """
    with open(output_file, 'w') as f:
        for partition_name, start, end in partitions:
            f.write(f"{partition_name} = {start}-{end}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter concatenated alignment for IQ-TREE compatibility"
    )
    parser.add_argument("input_alignment", help="Input FASTA alignment file")
    parser.add_argument("partition_file", help="IQ-TREE partition file")
    parser.add_argument("output_alignment", help="Output filtered FASTA alignment file")
    parser.add_argument("output_partition", help="Output partition file")
    parser.add_argument("--min-char-fraction", type=float, default=0.3,
                        help="Minimum fraction of non-gap characters per partition (default: 0.3)")
    parser.add_argument("--min-total-chars", type=int, default=500000,
                        help="Minimum total non-gap characters per sequence (default: 500000)")

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging()

    try:
        # Parse input files
        logger.info(f"Reading alignment from: {args.input_alignment}")
        sequences = parse_fasta_alignment(Path(args.input_alignment))
        logger.info(f"Loaded {len(sequences)} sequences")

        logger.info(f"Reading partitions from: {args.partition_file}")
        partitions = parse_partition_file(Path(args.partition_file))
        logger.info(f"Found {len(partitions)} partitions")

        if not partitions:
            logger.error("No partitions found in partition file")
            sys.exit(1)

        # Filter by overall quality first
        logger.info("Filtering sequences by overall character content...")
        valid_by_total, excluded_by_total = filter_sequences_by_overall_quality(
            sequences, args.min_total_chars
        )

        # Keep only sequences that pass total character filter
        filtered_sequences = {seq_id: sequences[seq_id] for seq_id in valid_by_total}

        # Analyze partition-specific character content
        logger.info("Analyzing partition-specific character content...")
        valid_sequences, excluded_by_partition, partition_stats = analyze_sequence_characters(
            filtered_sequences, partitions, args.min_char_fraction
        )

        # Final filtering
        final_sequences = {seq_id: sequences[seq_id] for seq_id in valid_sequences}

        logger.info(f"\n=== FILTERING SUMMARY ===")
        logger.info(f"Original sequences: {len(sequences)}")
        logger.info(f"Excluded by total characters: {len(excluded_by_total)}")
        logger.info(f"Excluded by partition analysis: {len(excluded_by_partition)}")
        logger.info(f"Final valid sequences: {len(final_sequences)}")

        if len(final_sequences) < 4:
            logger.error("Too few sequences remain for phylogenetic analysis")
            sys.exit(1)

        # Write filtered alignment and partition file
        write_filtered_alignment(final_sequences, valid_sequences, Path(args.output_alignment))
        write_filtered_partitions(partitions, Path(args.output_partition))

        # Write summary statistics
        stats_file = Path(args.output_alignment).parent / "filtering_stats.txt"
        with open(stats_file, 'w') as f:
            f.write("=== ALIGNMENT FILTERING SUMMARY ===\n")
            f.write(f"Original sequences: {len(sequences)}\n")
            f.write(f"Excluded by total characters: {len(excluded_by_total)}\n")
            f.write(f"Excluded by partition analysis: {len(excluded_by_partition)}\n")
            f.write(f"Final valid sequences: {len(final_sequences)}\n")
            f.write(f"Partitions: {len(partitions)}\n")
            f.write(f"\nExcluded sequences (total chars):\n")
            for seq_id in excluded_by_total:
                if seq_id in sequences:
                    total_chars = len(sequences[seq_id]) - sequences[seq_id].count('-')
                    f.write(f"  {seq_id}: {total_chars} characters\n")

            f.write(f"\nExcluded sequences (partition analysis):\n")
            for seq_id in excluded_by_partition:
                f.write(f"  {seq_id}\n")

            f.write(f"\nFinal valid sequences:\n")
            for seq_id in sorted(valid_sequences):
                if seq_id in sequences:
                    total_chars = len(sequences[seq_id]) - sequences[seq_id].count('-')
                    f.write(f"  {seq_id}: {total_chars} characters\n")

        logger.info(f"Filtering complete! Statistics written to: {stats_file}")

    except Exception as e:
        logger.error(f"Error filtering alignment: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()