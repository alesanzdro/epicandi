#!/usr/bin/env python3
"""
Select core single-copy orthologs that are present in all or most taxa for phylogenetic analysis.

Usage:
    python select_core_orthologs.py <single_copy_dir> <output_dir> [--min-taxa-fraction 1.0] [--max-orthologs 100]
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
            logging.FileHandler('select_core_orthologs.log')
        ]
    )
    return logging.getLogger(__name__)


def parse_fasta_file(fasta_file):
    """
    Parse a FASTA file and return sequence information.

    Args:
        fasta_file (Path): Path to FASTA file

    Returns:
        dict: Dictionary mapping species names to sequence info (id, length, sequence)
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
                    sequence = ''.join(current_seq)
                    # Extract species name from header
                    header = line[1:]
                    if '|' in header:
                        species = header.split('|')[0]
                    else:
                        species = header.split()[0]

                    # Clean up species name
                    species = species.split('_')[0]  # Remove additional suffixes

                    sequences[species] = {
                        'id': current_id,
                        'length': len(sequence),
                        'non_gap_length': len(sequence.replace('-', '')),
                        'sequence': sequence
                    }

                # Start new sequence
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None and current_seq:
            sequence = ''.join(current_seq)
            header = current_id
            if '|' in header:
                species = header.split('|')[0]
            else:
                species = header.split()[0]

            species = species.split('_')[0]

            sequences[species] = {
                'id': current_id,
                'length': len(sequence),
                'non_gap_length': len(sequence.replace('-', '')),
                'sequence': sequence
            }

    return sequences


def analyze_ortholog_coverage(single_copy_dir, min_seq_length=100):
    """
    Analyze coverage of each ortholog across taxa.

    Args:
        single_copy_dir (Path): Directory containing ortholog files
        min_seq_length (int): Minimum non-gap sequence length required

    Returns:
        tuple: (ortholog_stats, all_taxa)
    """
    logger = logging.getLogger(__name__)

    # Find all ortholog files
    og_files = list(single_copy_dir.glob("OG*.fa"))
    if not og_files:
        raise FileNotFoundError(f"No ortholog files found in {single_copy_dir}")

    logger.info(f"Found {len(og_files)} ortholog files")

    ortholog_stats = []
    all_taxa = set()

    for og_file in og_files:
        try:
            sequences = parse_fasta_file(og_file)

            # Filter sequences by minimum length
            valid_sequences = {}
            for species, seq_info in sequences.items():
                if seq_info['non_gap_length'] >= min_seq_length:
                    valid_sequences[species] = seq_info
                    all_taxa.add(species)

            if valid_sequences:
                # Calculate statistics
                alignment_length = max(seq_info['length'] for seq_info in valid_sequences.values())
                avg_coverage = sum(seq_info['non_gap_length'] for seq_info in valid_sequences.values()) / len(valid_sequences)

                ortholog_stats.append({
                    'file': og_file,
                    'ortholog': og_file.stem,
                    'num_taxa': len(valid_sequences),
                    'taxa': set(valid_sequences.keys()),
                    'alignment_length': alignment_length,
                    'avg_coverage': avg_coverage,
                    'sequences': valid_sequences
                })

        except Exception as e:
            logger.warning(f"Error processing {og_file}: {e}")

    logger.info(f"Processed orthologs with valid sequences: {len(ortholog_stats)}")
    logger.info(f"Total taxa found: {len(all_taxa)}")

    return ortholog_stats, all_taxa


def select_core_orthologs(ortholog_stats, all_taxa, min_taxa_fraction=1.0, max_orthologs=100):
    """
    Select core orthologs based on taxonomic coverage.

    Args:
        ortholog_stats (list): List of ortholog statistics
        all_taxa (set): Set of all taxa names
        min_taxa_fraction (float): Minimum fraction of taxa required
        max_orthologs (int): Maximum number of orthologs to select

    Returns:
        list: Selected ortholog statistics
    """
    logger = logging.getLogger(__name__)

    min_taxa_required = int(len(all_taxa) * min_taxa_fraction)
    logger.info(f"Minimum taxa required per ortholog: {min_taxa_required}/{len(all_taxa)} "
               f"({min_taxa_fraction:.1%})")

    # Filter orthologs by taxonomic coverage
    core_orthologs = []
    for stat in ortholog_stats:
        if stat['num_taxa'] >= min_taxa_required:
            core_orthologs.append(stat)

    logger.info(f"Orthologs meeting coverage requirement: {len(core_orthologs)}")

    if not core_orthologs:
        # Relax requirement if no orthologs meet the strict criteria
        logger.warning("No orthologs meet the strict coverage requirement. Relaxing criteria...")
        min_taxa_required = max(3, int(len(all_taxa) * 0.8))  # At least 80% or 3 taxa
        logger.info(f"Relaxed minimum taxa requirement: {min_taxa_required}/{len(all_taxa)}")

        for stat in ortholog_stats:
            if stat['num_taxa'] >= min_taxa_required:
                core_orthologs.append(stat)

    # Sort by taxonomic coverage (descending) and average coverage (descending)
    core_orthologs.sort(key=lambda x: (x['num_taxa'], x['avg_coverage']), reverse=True)

    # Limit number of orthologs
    if max_orthologs and len(core_orthologs) > max_orthologs:
        core_orthologs = core_orthologs[:max_orthologs]
        logger.info(f"Limited to top {max_orthologs} orthologs")

    return core_orthologs


def write_selected_alignments(selected_orthologs, output_dir):
    """
    Write selected ortholog alignments to output directory.

    Args:
        selected_orthologs (list): List of selected ortholog statistics
        output_dir (Path): Output directory
    """
    logger = logging.getLogger(__name__)

    alignments_dir = output_dir / "selected_alignments"
    alignments_dir.mkdir(parents=True, exist_ok=True)

    for ortholog_stat in selected_orthologs:
        output_file = alignments_dir / f"{ortholog_stat['ortholog']}.fa"

        with open(output_file, 'w') as f:
            for species, seq_info in ortholog_stat['sequences'].items():
                f.write(f">{seq_info['id']}\n{seq_info['sequence']}\n")

        logger.debug(f"Written {ortholog_stat['ortholog']}: {ortholog_stat['num_taxa']} taxa")

    logger.info(f"Written {len(selected_orthologs)} selected alignments to {alignments_dir}")


def create_concatenated_alignment(selected_orthologs, all_taxa, output_dir):
    """
    Create concatenated alignment from selected orthologs.

    Args:
        selected_orthologs (list): List of selected ortholog statistics
        all_taxa (set): Set of all taxa names
        output_dir (Path): Output directory

    Returns:
        tuple: (concatenated_file, statistics)
    """
    logger = logging.getLogger(__name__)

    # Build concatenated sequences
    concatenated_seqs = defaultdict(str)
    used_orthologs = []
    alignment_lengths = []

    # Get taxa that appear in selected orthologs
    represented_taxa = set()
    for ortholog_stat in selected_orthologs:
        represented_taxa.update(ortholog_stat['taxa'])

    logger.info(f"Taxa represented in selected orthologs: {len(represented_taxa)}")

    for ortholog_stat in selected_orthologs:
        sequences = ortholog_stat['sequences']
        alignment_length = ortholog_stat['alignment_length']

        # Add sequences for all represented taxa
        for taxon in represented_taxa:
            if taxon in sequences:
                concatenated_seqs[taxon] += sequences[taxon]['sequence']
            else:
                # Add missing data for taxa not in this ortholog
                concatenated_seqs[taxon] += '-' * alignment_length

        used_orthologs.append(ortholog_stat['ortholog'])
        alignment_lengths.append(alignment_length)

    # Write concatenated alignment
    concatenated_file = output_dir / "concatenated_core_orthologs.fasta"

    with open(concatenated_file, 'w') as f:
        for species in sorted(concatenated_seqs.keys()):
            sequence = concatenated_seqs[species]
            f.write(f">{species}\n{sequence}\n")

    # Calculate statistics
    total_length = len(list(concatenated_seqs.values())[0]) if concatenated_seqs else 0
    num_taxa = len(concatenated_seqs)
    num_orthologs = len(used_orthologs)

    statistics = {
        'concatenated_file': concatenated_file,
        'num_taxa': num_taxa,
        'num_orthologs': num_orthologs,
        'total_length': total_length,
        'used_orthologs': used_orthologs,
        'alignment_lengths': alignment_lengths,
        'avg_ortholog_length': sum(alignment_lengths) / len(alignment_lengths) if alignment_lengths else 0
    }

    logger.info(f"Concatenated alignment created:")
    logger.info(f"  Taxa: {num_taxa}")
    logger.info(f"  Orthologs: {num_orthologs}")
    logger.info(f"  Total length: {total_length} bp")
    logger.info(f"  Average ortholog length: {statistics['avg_ortholog_length']:.1f} bp")

    return concatenated_file, statistics


def main():
    parser = argparse.ArgumentParser(
        description="Select core single-copy orthologs for phylogenetic analysis"
    )
    parser.add_argument("single_copy_dir", help="Directory containing single-copy ortholog files")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--min-taxa-fraction", type=float, default=1.0,
                        help="Minimum fraction of taxa required per ortholog (default: 1.0)")
    parser.add_argument("--max-orthologs", type=int, default=100,
                        help="Maximum number of orthologs to select (default: 100)")
    parser.add_argument("--min-seq-length", type=int, default=100,
                        help="Minimum non-gap sequence length (default: 100)")

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging()

    try:
        single_copy_dir = Path(args.single_copy_dir)
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Analyzing orthologs from: {single_copy_dir}")
        logger.info(f"Output directory: {output_dir}")

        # Analyze ortholog coverage
        ortholog_stats, all_taxa = analyze_ortholog_coverage(single_copy_dir, args.min_seq_length)

        if not ortholog_stats:
            logger.error("No valid orthologs found")
            sys.exit(1)

        # Select core orthologs
        selected_orthologs = select_core_orthologs(
            ortholog_stats, all_taxa, args.min_taxa_fraction, args.max_orthologs
        )

        if not selected_orthologs:
            logger.error("No orthologs selected")
            sys.exit(1)

        # Write selected alignments
        write_selected_alignments(selected_orthologs, output_dir)

        # Create concatenated alignment
        concatenated_file, statistics = create_concatenated_alignment(
            selected_orthologs, all_taxa, output_dir
        )

        # Write summary statistics
        stats_file = output_dir / "core_orthologs_stats.txt"
        with open(stats_file, 'w') as f:
            f.write("=== CORE ORTHOLOGS SELECTION SUMMARY ===\n")
            f.write(f"Total taxa found: {len(all_taxa)}\n")
            f.write(f"Total orthologs analyzed: {len(ortholog_stats)}\n")
            f.write(f"Selected orthologs: {statistics['num_orthologs']}\n")
            f.write(f"Taxa in concatenated alignment: {statistics['num_taxa']}\n")
            f.write(f"Total alignment length: {statistics['total_length']} bp\n")
            f.write(f"Average ortholog length: {statistics['avg_ortholog_length']:.1f} bp\n")
            f.write(f"Minimum taxa fraction required: {args.min_taxa_fraction}\n")
            f.write(f"Maximum orthologs limit: {args.max_orthologs}\n")
            f.write(f"\nSelected orthologs:\n")
            for i, ortholog_stat in enumerate(selected_orthologs, 1):
                f.write(f"  {i:3d}. {ortholog_stat['ortholog']} "
                       f"({ortholog_stat['num_taxa']} taxa, "
                       f"{ortholog_stat['avg_coverage']:.1f} avg coverage)\n")

            f.write(f"\nTaxa in alignment:\n")
            for taxon in sorted(all_taxa):
                f.write(f"  {taxon}\n")

        logger.info(f"Summary written to: {stats_file}")
        logger.info(f"Concatenated alignment: {concatenated_file}")

        print(f"\n=== CORE ORTHOLOGS SELECTION COMPLETED ===")
        print(f"Selected orthologs: {statistics['num_orthologs']}")
        print(f"Taxa: {statistics['num_taxa']}")
        print(f"Total length: {statistics['total_length']} bp")
        print(f"Concatenated file: {concatenated_file}")

    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()