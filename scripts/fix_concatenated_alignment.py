#!/usr/bin/env python3
"""
Fix concatenated alignment sequence length inconsistencies for IQ-TREE partitioned analysis.

Usage:
    python fix_concatenated_alignment.py <input_alignment> <partition_file> <output_alignment>
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
            logging.FileHandler('fix_concatenated_alignment.log')
        ]
    )
    return logging.getLogger(__name__)


def parse_partition_file(partition_file):
    """
    Parse IQ-TREE partition file to get expected partition lengths.

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
                # Handle different partition file formats:
                # Format 1: DNA, partition_name = start-end
                # Format 2: partition_name = start-end
                if '=' in line:
                    parts = line.split('=')
                    if len(parts) == 2:
                        left_part = parts[0].strip()
                        range_part = parts[1].strip()

                        # Check if format includes model (e.g., "DNA, partition_name")
                        if ',' in left_part:
                            partition_name = left_part.split(',')[1].strip()
                        else:
                            # Simple format: partition_name = start-end
                            partition_name = left_part

                        if '-' in range_part:
                            start, end = map(int, range_part.split('-'))
                            partitions.append((partition_name, start, end))

    return partitions


def parse_fasta_alignment(fasta_file):
    """
    Parse FASTA alignment file.

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


def calculate_expected_length(partitions):
    """
    Calculate expected total alignment length from partitions.

    Args:
        partitions (list): List of partition tuples

    Returns:
        int: Expected total length
    """
    if not partitions:
        return 0

    # Find maximum end position
    max_end = max(end for _, _, end in partitions)
    return max_end


def fix_sequence_lengths(sequences, expected_length, gap_char='-'):
    """
    Fix sequence lengths by padding with gaps or truncating.

    Args:
        sequences (dict): Dictionary of sequences
        expected_length (int): Expected sequence length
        gap_char (str): Character to use for padding

    Returns:
        dict: Fixed sequences
    """
    logger = logging.getLogger(__name__)
    fixed_sequences = {}

    for seq_id, sequence in sequences.items():
        current_length = len(sequence)

        if current_length < expected_length:
            # Pad with gaps
            padding_needed = expected_length - current_length
            fixed_sequence = sequence + (gap_char * padding_needed)
            logger.info(f"Padded {seq_id}: {current_length} -> {expected_length} (+{padding_needed})")

        elif current_length > expected_length:
            # Truncate (remove trailing gaps first)
            sequence_trimmed = sequence.rstrip(gap_char)
            if len(sequence_trimmed) <= expected_length:
                # Safe to truncate trailing gaps
                fixed_sequence = sequence_trimmed + (gap_char * (expected_length - len(sequence_trimmed)))
                logger.info(f"Truncated trailing gaps {seq_id}: {current_length} -> {expected_length}")
            else:
                # Need to truncate actual sequence - this might be problematic
                fixed_sequence = sequence[:expected_length]
                logger.warning(f"Truncated sequence data {seq_id}: {current_length} -> {expected_length}")

        else:
            # Length is correct
            fixed_sequence = sequence

        fixed_sequences[seq_id] = fixed_sequence

    return fixed_sequences


def write_fasta_alignment(sequences, output_file):
    """
    Write sequences to FASTA format.

    Args:
        sequences (dict): Dictionary of sequences
        output_file (Path): Output file path
    """
    with open(output_file, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in lines of 80 characters
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Fix concatenated alignment sequence length inconsistencies"
    )
    parser.add_argument("input_alignment", help="Input FASTA alignment file")
    parser.add_argument("partition_file", help="IQ-TREE partition file")
    parser.add_argument("output_alignment", help="Output fixed FASTA alignment file")
    parser.add_argument("--gap-char", default='-', help="Gap character (default: -)")

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

        # Calculate expected length
        expected_length = calculate_expected_length(partitions)
        logger.info(f"Expected alignment length: {expected_length}")

        # Check current sequence lengths
        seq_lengths = [len(seq) for seq in sequences.values()]
        logger.info(f"Current sequence lengths: {min(seq_lengths)} - {max(seq_lengths)}")

        # Fix sequence lengths
        logger.info("Fixing sequence lengths...")
        fixed_sequences = fix_sequence_lengths(sequences, expected_length, args.gap_char)

        # Verify all sequences have correct length
        fixed_lengths = [len(seq) for seq in fixed_sequences.values()]
        if len(set(fixed_lengths)) > 1:
            logger.error(f"Still have inconsistent lengths: {set(fixed_lengths)}")
            sys.exit(1)

        logger.info(f"All sequences now have length: {fixed_lengths[0]}")

        # Write fixed alignment
        logger.info(f"Writing fixed alignment to: {args.output_alignment}")
        write_fasta_alignment(fixed_sequences, Path(args.output_alignment))

        logger.info("Alignment fixed successfully!")

        # Print summary
        print("\n=== ALIGNMENT FIX SUMMARY ===")
        print(f"Input sequences: {len(sequences)}")
        print(f"Expected length: {expected_length}")
        print(f"Fixed alignment: {args.output_alignment}")
        print(f"All sequences now have length: {fixed_lengths[0]}")

    except Exception as e:
        logger.error(f"Error fixing alignment: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()