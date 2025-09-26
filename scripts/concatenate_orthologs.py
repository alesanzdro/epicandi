#!/usr/bin/env python3
"""
Extract and concatenate single-copy orthologous sequences from OrthoFinder results.

Usage:
    python concatenate_orthologs.py <orthofinder_results_dir> <output_dir> [--max-orthologs N]
"""

import os
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
            logging.FileHandler('concatenate_orthologs.log')
        ]
    )
    return logging.getLogger(__name__)


def find_orthofinder_results(results_dir):
    """
    Find OrthoFinder results directory and single-copy ortholog files.

    Args:
        results_dir (str): Path to OrthoFinder results directory

    Returns:
        tuple: (results_path, single_copy_dir)
    """
    results_path = Path(results_dir)

    # Look for Results_* directory
    result_subdirs = list(results_path.glob("Results_*"))
    if not result_subdirs:
        # Maybe results_dir is already the Results_* directory
        if results_path.name.startswith("Results_"):
            ortho_results = results_path
        else:
            raise FileNotFoundError(f"No OrthoFinder Results_* directory found in {results_dir}")
    else:
        ortho_results = result_subdirs[0]

    single_copy_dir = ortho_results / "Single_Copy_Orthologue_Sequences"

    if not single_copy_dir.exists():
        raise FileNotFoundError(f"Single-copy ortholog directory not found: {single_copy_dir}")

    return ortho_results, single_copy_dir


def parse_fasta_file(fasta_file):
    """
    Parse a FASTA file and return sequences by species.

    Args:
        fasta_file (Path): Path to FASTA file

    Returns:
        dict: Dictionary mapping species names to sequences
    """
    sequences = {}
    current_seq = []
    current_id = None

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None and current_seq:
                    sequences[current_id] = ''.join(current_seq)

                # Extract species name from header
                # Format: >SpeciesName|ProteinID or >SpeciesName_additional_info
                header = line[1:]
                if '|' in header:
                    current_id = header.split('|')[0]
                else:
                    current_id = header.split()[0]

                # Clean up species name
                current_id = current_id.split('_')[0]  # Remove additional suffixes
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None and current_seq:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def concatenate_orthologs(single_copy_dir, output_dir, max_orthologs=None, min_taxa=3):
    """
    Concatenate single-copy orthologous sequences.

    Args:
        single_copy_dir (Path): Directory containing single-copy ortholog files
        output_dir (Path): Output directory for concatenated alignment
        max_orthologs (int): Maximum number of orthologs to use (None for all)
        min_taxa (int): Minimum number of taxa required per ortholog

    Returns:
        tuple: (concatenated_file, statistics)
    """
    logger = logging.getLogger(__name__)

    # Find all ortholog files
    og_files = list(single_copy_dir.glob("OG*.fa"))
    if not og_files:
        raise FileNotFoundError(f"No ortholog files found in {single_copy_dir}")

    logger.info(f"Found {len(og_files)} single-copy ortholog files")

    # Limit number of orthologs if specified
    if max_orthologs and max_orthologs < len(og_files):
        og_files = sorted(og_files)[:max_orthologs]
        logger.info(f"Using first {max_orthologs} orthologs for concatenation")

    # Concatenate sequences
    concatenated_seqs = defaultdict(str)
    used_orthologs = []
    alignment_lengths = []

    for og_file in og_files:
        try:
            sequences = parse_fasta_file(og_file)

            # Check if ortholog has minimum required taxa
            if len(sequences) < min_taxa:
                logger.warning(f"Skipping {og_file.name}: only {len(sequences)} taxa "
                             f"(minimum {min_taxa} required)")
                continue

            # Get alignment length and check sequence quality
            seq_lengths = [len(seq.replace('-', '')) for seq in sequences.values()]  # Count only non-gap characters
            min_seq_length = min(seq_lengths)

            # Skip orthologs with very short sequences (less than 100 aa)
            if min_seq_length < 100:
                logger.warning(f"Skipping {og_file.name}: shortest sequence only {min_seq_length} aa")
                continue

            # Get actual alignment length
            alignment_lengths_raw = [len(seq) for seq in sequences.values()]
            if len(set(alignment_lengths_raw)) > 1:
                logger.warning(f"Unequal sequence lengths in {og_file.name}: {alignment_lengths_raw}")

            alignment_length = max(alignment_lengths_raw)
            alignment_lengths.append(alignment_length)

            # Add sequences to concatenated alignment - ensure all taxa are represented
            all_taxa = set(concatenated_seqs.keys()) | set(sequences.keys())

            for taxon in all_taxa:
                if taxon in sequences:
                    sequence = sequences[taxon]
                    # Pad sequence to alignment length if needed
                    if len(sequence) < alignment_length:
                        sequence += '-' * (alignment_length - len(sequence))
                    concatenated_seqs[taxon] += sequence
                else:
                    # Add missing data for taxa not in this ortholog
                    concatenated_seqs[taxon] += '-' * alignment_length

            used_orthologs.append(og_file.name)

        except Exception as e:
            logger.error(f"Error processing {og_file}: {e}")

    if not concatenated_seqs:
        raise ValueError("No valid orthologs found for concatenation")

    # Filter out sequences that are mostly gaps (>70% gaps)
    filtered_seqs = {}
    for species, sequence in concatenated_seqs.items():
        gap_percentage = sequence.count('-') / len(sequence)
        if gap_percentage < 0.7:  # Keep sequences with <70% gaps
            filtered_seqs[species] = sequence
        else:
            logger.warning(f"Removing {species}: {gap_percentage:.1%} gaps")

    if len(filtered_seqs) < 3:
        raise ValueError("Too few sequences remain after filtering")

    # Write concatenated alignment
    output_dir.mkdir(parents=True, exist_ok=True)
    concatenated_file = output_dir / "concatenated_orthologs.fasta"

    with open(concatenated_file, 'w') as f:
        for species in sorted(filtered_seqs.keys()):
            sequence = filtered_seqs[species]
            f.write(f">{species}\n{sequence}\n")

    # Calculate statistics
    total_length = len(list(filtered_seqs.values())[0])
    num_taxa = len(filtered_seqs)
    num_orthologs = len(used_orthologs)

    statistics = {
        'concatenated_file': concatenated_file,
        'num_taxa': num_taxa,
        'num_orthologs': num_orthologs,
        'total_length': total_length,
        'used_orthologs': used_orthologs,
        'alignment_lengths': alignment_lengths
    }

    logger.info(f"Concatenation complete:")
    logger.info(f"  Output file: {concatenated_file}")
    logger.info(f"  Taxa: {num_taxa}")
    logger.info(f"  Orthologs used: {num_orthologs}")
    logger.info(f"  Total alignment length: {total_length} bp")

    return concatenated_file, statistics


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate single-copy orthologs from OrthoFinder results"
    )
    parser.add_argument("orthofinder_results", help="OrthoFinder results directory")
    parser.add_argument("output_dir", help="Output directory for concatenated alignment")
    parser.add_argument("--max-orthologs", type=int, default=None,
                        help="Maximum number of orthologs to use (default: all)")
    parser.add_argument("--min-taxa", type=int, default=3,
                        help="Minimum number of taxa required per ortholog (default: 3)")

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging()

    try:
        # Find OrthoFinder results
        results_dir, single_copy_dir = find_orthofinder_results(args.orthofinder_results)
        logger.info(f"Using OrthoFinder results from: {results_dir}")

        # Concatenate orthologs
        output_dir = Path(args.output_dir)
        concatenated_file, stats = concatenate_orthologs(
            single_copy_dir, output_dir, args.max_orthologs, args.min_taxa
        )

        # Write statistics
        stats_file = output_dir / "concatenation_stats.txt"
        with open(stats_file, 'w') as f:
            f.write(f"Concatenation Statistics\n")
            f.write(f"========================\n")
            f.write(f"Number of taxa: {stats['num_taxa']}\n")
            f.write(f"Number of orthologs: {stats['num_orthologs']}\n")
            f.write(f"Total alignment length: {stats['total_length']} bp\n")
            f.write(f"Average ortholog length: {sum(stats['alignment_lengths'])/len(stats['alignment_lengths']):.1f} bp\n")
            f.write(f"\nUsed orthologs:\n")
            for og in stats['used_orthologs']:
                f.write(f"  {og}\n")

        logger.info(f"Statistics written to: {stats_file}")

    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()