#!/usr/bin/env python3
"""
Comprehensive phylogenetic alignment processor for IQ-TREE compatibility.

This script consolidates functionality for:
1. Fixing concatenated alignment length inconsistencies
2. Quality-based sequence filtering
3. Core ortholog selection and concatenation

Usage:
    python phylogenetic_alignment_processor.py <mode> [options]

Modes:
    fix-lengths     - Fix sequence length inconsistencies
    filter-quality  - Filter sequences by character content quality
    core-orthologs  - Select and concatenate core orthologs (RECOMMENDED)
"""

import sys
import os
import argparse
from pathlib import Path
from collections import defaultdict
import logging


def setup_logging(log_file="phylogenetic_alignment_processor.log"):
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file)
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


def write_fasta_alignment(sequences, output_file):
    """
    Write sequences to FASTA format.

    Args:
        sequences (dict): Dictionary of sequences
        output_file (Path): Output file path
    """
    with open(output_file, 'w') as f:
        for seq_id in sorted(sequences.keys()):
            sequence = sequences[seq_id]
            f.write(f">{seq_id}\n")
            # Write sequence in lines of 80 characters
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')


def fix_sequence_lengths(args):
    """Fix sequence length inconsistencies for IQ-TREE compatibility."""
    logger = logging.getLogger(__name__)

    logger.info(f"Mode: Fix sequence lengths")
    logger.info(f"Input alignment: {args.input_alignment}")
    logger.info(f"Partition file: {args.partition_file}")
    logger.info(f"Output alignment: {args.output_alignment}")

    # Parse input files
    logger.info("Reading alignment...")
    sequences = parse_fasta_alignment(Path(args.input_alignment))
    logger.info(f"Loaded {len(sequences)} sequences")

    logger.info("Reading partitions...")
    partitions = parse_partition_file(Path(args.partition_file))
    logger.info(f"Found {len(partitions)} partitions")

    if not partitions:
        logger.error("No partitions found in partition file")
        return False

    # Calculate expected length
    expected_length = max(end for _, _, end in partitions)
    logger.info(f"Expected alignment length: {expected_length}")

    # Check current sequence lengths
    seq_lengths = [len(seq) for seq in sequences.values()]
    logger.info(f"Current sequence lengths: {min(seq_lengths)} - {max(seq_lengths)}")

    # Fix sequence lengths
    logger.info("Fixing sequence lengths...")
    fixed_sequences = {}
    gap_char = args.gap_char if hasattr(args, 'gap_char') else '-'

    for seq_id, sequence in sequences.items():
        current_length = len(sequence)

        if current_length < expected_length:
            # Pad with gaps
            padding_needed = expected_length - current_length
            fixed_sequence = sequence + (gap_char * padding_needed)
            logger.debug(f"Padded {seq_id}: {current_length} -> {expected_length} (+{padding_needed})")

        elif current_length > expected_length:
            # Truncate (remove trailing gaps first)
            sequence_trimmed = sequence.rstrip(gap_char)
            if len(sequence_trimmed) <= expected_length:
                # Safe to truncate trailing gaps
                fixed_sequence = sequence_trimmed + (gap_char * (expected_length - len(sequence_trimmed)))
                logger.debug(f"Truncated trailing gaps {seq_id}: {current_length} -> {expected_length}")
            else:
                # Need to truncate actual sequence - this might be problematic
                fixed_sequence = sequence[:expected_length]
                logger.warning(f"Truncated sequence data {seq_id}: {current_length} -> {expected_length}")

        else:
            # Length is correct
            fixed_sequence = sequence

        fixed_sequences[seq_id] = fixed_sequence

    # Verify all sequences have correct length
    fixed_lengths = [len(seq) for seq in fixed_sequences.values()]
    if len(set(fixed_lengths)) > 1:
        logger.error(f"Still have inconsistent lengths: {set(fixed_lengths)}")
        return False

    logger.info(f"All sequences now have length: {fixed_lengths[0]}")

    # Write fixed alignment
    logger.info(f"Writing fixed alignment to: {args.output_alignment}")
    write_fasta_alignment(fixed_sequences, Path(args.output_alignment))

    logger.info("Length fixing completed successfully!")
    return True


def filter_by_quality(args):
    """Filter sequences by character content quality."""
    logger = logging.getLogger(__name__)

    logger.info(f"Mode: Filter by quality")
    logger.info(f"Input alignment: {args.input_alignment}")
    logger.info(f"Output alignment: {args.output_alignment}")

    # Parse input files
    logger.info("Reading alignment...")
    sequences = parse_fasta_alignment(Path(args.input_alignment))
    logger.info(f"Loaded {len(sequences)} sequences")

    partitions = []
    if hasattr(args, 'partition_file') and args.partition_file:
        logger.info("Reading partitions...")
        partitions = parse_partition_file(Path(args.partition_file))
        logger.info(f"Found {len(partitions)} partitions")

    # Filter by overall quality first
    logger.info("Filtering sequences by overall character content...")
    min_total_chars = args.min_total_chars if hasattr(args, 'min_total_chars') else 500000

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

    # Keep only sequences that pass total character filter
    filtered_sequences = {seq_id: sequences[seq_id] for seq_id in valid_sequences}

    # Additional partition-based filtering if partitions available
    if partitions and hasattr(args, 'min_char_fraction'):
        logger.info("Analyzing partition-specific character content...")
        min_char_fraction = args.min_char_fraction

        for i, (partition_name, start, end) in enumerate(partitions):
            if i % 100 == 0:
                logger.info(f"Processing partition {i+1}/{len(partitions)}: {partition_name}")

            partition_length = end - start + 1
            min_chars_required = int(partition_length * min_char_fraction)

            for seq_id in list(valid_sequences):
                if seq_id in filtered_sequences:
                    sequence = filtered_sequences[seq_id]
                    if len(sequence) >= end:
                        # Extract partition sequence (convert to 0-based indexing)
                        partition_seq = sequence[start-1:end]
                        non_gap_chars = len(partition_seq) - partition_seq.count('-')

                        if non_gap_chars < min_chars_required:
                            if seq_id not in excluded_sequences:
                                excluded_sequences.append(seq_id)
                                logger.warning(f"Excluding {seq_id}: only {non_gap_chars}/{partition_length} "
                                             f"characters in partition {partition_name}")

            # Update valid sequences
            valid_sequences = [seq_id for seq_id in valid_sequences if seq_id not in excluded_sequences]

    # Final filtering
    final_sequences = {seq_id: sequences[seq_id] for seq_id in valid_sequences}

    logger.info(f"\n=== FILTERING SUMMARY ===")
    logger.info(f"Original sequences: {len(sequences)}")
    logger.info(f"Excluded sequences: {len(excluded_sequences)}")
    logger.info(f"Final valid sequences: {len(final_sequences)}")

    if len(final_sequences) < 3:
        logger.warning("Less than 3 sequences remain - insufficient for phylogenetic analysis")

    # Write filtered alignment
    write_fasta_alignment(final_sequences, Path(args.output_alignment))

    # Write partition file if provided
    if partitions and hasattr(args, 'output_partition'):
        with open(args.output_partition, 'w') as f:
            for partition_name, start, end in partitions:
                f.write(f"{partition_name} = {start}-{end}\n")

    logger.info("Quality filtering completed successfully!")
    return True


def parse_fasta_file_with_species(fasta_file):
    """
    Parse a FASTA file and return sequence information with species extraction.

    Args:
        fasta_file (Path): Path to FASTA file

    Returns:
        dict: Dictionary mapping species names to sequence info
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


def select_core_orthologs(args):
    """Select and concatenate core orthologs for phylogenetic analysis."""
    logger = logging.getLogger(__name__)

    logger.info(f"Mode: Core orthologs selection")
    logger.info(f"Single copy directory: {args.single_copy_dir}")
    logger.info(f"Output directory: {args.output_dir}")

    single_copy_dir = Path(args.single_copy_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parameters
    min_taxa_fraction = args.min_taxa_fraction if hasattr(args, 'min_taxa_fraction') else 0.8
    max_orthologs = args.max_orthologs if hasattr(args, 'max_orthologs') else 100
    min_seq_length = args.min_seq_length if hasattr(args, 'min_seq_length') else 150

    logger.info(f"Parameters:")
    logger.info(f"  Min taxa fraction: {min_taxa_fraction}")
    logger.info(f"  Max orthologs: {max_orthologs}")
    logger.info(f"  Min sequence length: {min_seq_length}")

    # Find all ortholog files
    og_files = list(single_copy_dir.glob("OG*.fa"))
    if not og_files:
        logger.error(f"No ortholog files found in {single_copy_dir}")
        return False

    logger.info(f"Found {len(og_files)} ortholog files")

    # Analyze ortholog coverage
    ortholog_stats = []
    all_taxa = set()

    for og_file in og_files:
        try:
            sequences = parse_fasta_file_with_species(og_file)

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

    if not ortholog_stats:
        logger.error("No valid orthologs found")
        return False

    # Select core orthologs
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

    if not core_orthologs:
        logger.error("No orthologs selected after filtering")
        return False

    # Write selected alignments
    alignments_dir = output_dir / "selected_alignments"
    alignments_dir.mkdir(parents=True, exist_ok=True)

    for ortholog_stat in core_orthologs:
        output_file = alignments_dir / f"{ortholog_stat['ortholog']}.fa"

        with open(output_file, 'w') as f:
            for species, seq_info in ortholog_stat['sequences'].items():
                f.write(f">{seq_info['id']}\n{seq_info['sequence']}\n")

    logger.info(f"Written {len(core_orthologs)} selected alignments to {alignments_dir}")

    # Create concatenated alignment
    logger.info("Creating concatenated alignment...")

    # Build concatenated sequences
    concatenated_seqs = defaultdict(str)
    used_orthologs = []
    alignment_lengths = []

    # Get taxa that appear in selected orthologs
    represented_taxa = set()
    for ortholog_stat in core_orthologs:
        represented_taxa.update(ortholog_stat['taxa'])

    logger.info(f"Taxa represented in selected orthologs: {len(represented_taxa)}")

    for ortholog_stat in core_orthologs:
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
    avg_ortholog_length = sum(alignment_lengths) / len(alignment_lengths) if alignment_lengths else 0

    # Write summary statistics
    stats_file = output_dir / "core_orthologs_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("=== CORE ORTHOLOGS SELECTION SUMMARY ===\n")
        f.write(f"Total taxa found: {len(all_taxa)}\n")
        f.write(f"Total orthologs analyzed: {len(ortholog_stats)}\n")
        f.write(f"Selected orthologs: {num_orthologs}\n")
        f.write(f"Taxa in concatenated alignment: {num_taxa}\n")
        f.write(f"Total alignment length: {total_length} bp\n")
        f.write(f"Average ortholog length: {avg_ortholog_length:.1f} bp\n")
        f.write(f"Minimum taxa fraction required: {min_taxa_fraction}\n")
        f.write(f"Maximum orthologs limit: {max_orthologs}\n")
        f.write(f"\nSelected orthologs:\n")
        for i, ortholog_stat in enumerate(core_orthologs, 1):
            f.write(f"  {i:3d}. {ortholog_stat['ortholog']} "
                   f"({ortholog_stat['num_taxa']} taxa, "
                   f"{ortholog_stat['avg_coverage']:.1f} avg coverage)\n")

        f.write(f"\nTaxa in alignment:\n")
        for taxon in sorted(represented_taxa):
            f.write(f"  {taxon}\n")

    logger.info(f"Concatenated alignment created:")
    logger.info(f"  Taxa: {num_taxa}")
    logger.info(f"  Orthologs: {num_orthologs}")
    logger.info(f"  Total length: {total_length} bp")
    logger.info(f"  Average ortholog length: {avg_ortholog_length:.1f} bp")
    logger.info(f"  Output file: {concatenated_file}")
    logger.info(f"  Statistics: {stats_file}")

    logger.info("Core orthologs selection completed successfully!")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive phylogenetic alignment processor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fix sequence length inconsistencies
  python phylogenetic_alignment_processor.py fix-lengths \\
    --input-alignment concat.aln \\
    --partition-file concat.part \\
    --output-alignment concat_fixed.aln

  # Filter sequences by quality
  python phylogenetic_alignment_processor.py filter-quality \\
    --input-alignment concat_fixed.aln \\
    --output-alignment concat_filtered.aln \\
    --min-total-chars 1800000 \\
    --min-char-fraction 0.8 \\
    --partition-file concat.part \\
    --output-partition concat_filtered.part

  # Select core orthologs (RECOMMENDED)
  python phylogenetic_alignment_processor.py core-orthologs \\
    --single-copy-dir single_copy_alignments \\
    --output-dir core_orthologs \\
    --min-taxa-fraction 0.8 \\
    --max-orthologs 100 \\
    --min-seq-length 150
        """
    )

    subparsers = parser.add_subparsers(dest='mode', help='Processing mode')

    # Fix lengths mode
    fix_parser = subparsers.add_parser('fix-lengths', help='Fix sequence length inconsistencies')
    fix_parser.add_argument('--input-alignment', required=True, help='Input FASTA alignment file')
    fix_parser.add_argument('--partition-file', required=True, help='IQ-TREE partition file')
    fix_parser.add_argument('--output-alignment', required=True, help='Output fixed FASTA alignment file')
    fix_parser.add_argument('--gap-char', default='-', help='Gap character (default: -)')

    # Filter quality mode
    filter_parser = subparsers.add_parser('filter-quality', help='Filter sequences by character content quality')
    filter_parser.add_argument('--input-alignment', required=True, help='Input FASTA alignment file')
    filter_parser.add_argument('--output-alignment', required=True, help='Output filtered FASTA alignment file')
    filter_parser.add_argument('--min-total-chars', type=int, default=500000,
                               help='Minimum total non-gap characters per sequence (default: 500000)')
    filter_parser.add_argument('--min-char-fraction', type=float, default=0.3,
                               help='Minimum fraction of non-gap characters per partition (default: 0.3)')
    filter_parser.add_argument('--partition-file', help='IQ-TREE partition file (optional)')
    filter_parser.add_argument('--output-partition', help='Output partition file (optional)')

    # Core orthologs mode
    core_parser = subparsers.add_parser('core-orthologs', help='Select and concatenate core orthologs (RECOMMENDED)')
    core_parser.add_argument('--single-copy-dir', required=True,
                             help='Directory containing single-copy ortholog files')
    core_parser.add_argument('--output-dir', required=True, help='Output directory')
    core_parser.add_argument('--min-taxa-fraction', type=float, default=0.8,
                             help='Minimum fraction of taxa required per ortholog (default: 0.8)')
    core_parser.add_argument('--max-orthologs', type=int, default=100,
                             help='Maximum number of orthologs to select (default: 100)')
    core_parser.add_argument('--min-seq-length', type=int, default=150,
                             help='Minimum non-gap sequence length (default: 150)')

    args = parser.parse_args()

    if not args.mode:
        parser.print_help()
        sys.exit(1)

    # Set up logging
    logger = setup_logging()

    try:
        success = False
        if args.mode == 'fix-lengths':
            success = fix_sequence_lengths(args)
        elif args.mode == 'filter-quality':
            success = filter_by_quality(args)
        elif args.mode == 'core-orthologs':
            success = select_core_orthologs(args)
        else:
            logger.error(f"Unknown mode: {args.mode}")
            sys.exit(1)

        if success:
            logger.info(f"Mode '{args.mode}' completed successfully!")
        else:
            logger.error(f"Mode '{args.mode}' failed!")
            sys.exit(1)

    except Exception as e:
        logger.error(f"Error in mode '{args.mode}': {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()