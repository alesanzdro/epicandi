# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**nanocheck** is a Snakemake-based bioinformatics pipeline for analyzing Oxford Nanopore long-read sequencing data, specifically optimized for fungal genome assembly and quality control. The pipeline focuses on processing unclassified/problematic nanopore reads and performing comprehensive genome assembly workflows.

## Key Commands

### Environment Setup
```bash
# Create and activate the conda environment
conda env create -f nanocheck.yml
conda activate nanocheck

# Alternative using the resources environment file
conda env create -f resources/envs/nanocheck.yml
```

### Running the Pipeline
```bash
# Basic pipeline execution
snakemake --cores 8

# Dry run to check workflow
snakemake -n

# Generate workflow visualization
snakemake --dag | dot -Tsvg > pipeline_dag.svg

# Clean temporary files
snakemake clean

# Remove all output
snakemake clean_all
```

### Helper Scripts
```bash
# Analyze FASTQ read lengths (multi-threaded)
python3 scripts/fastq_len_check.py -f input_file.fastq.gz -t 8

# Alternative read length checker
python3 scripts/fastq_len_check_2.py input_file.fastq.gz
```

## Architecture

### Core Pipeline (Snakefile)
- **Language**: Python 3 with Snakemake workflow management
- **Target organisms**: Candida species (C. albicans, Candidozyma auris, Nakaseomyces glabratus)
- **Workflow stages**:
  1. Data preparation and copying
  2. Initial quality control (NanoPlot)
  3. Adapter trimming (Porechop)
  4. Quality filtering (Filtlong)
  5. Post-filtering QC
  6. Genome assembly (Flye)
  7. Consensus polishing (Medaka)
  8. Assembly evaluation (Quast)

### Key Configuration Parameters (Snakefile:16-26)
- `GENOME_SIZE`: "13m" (estimated genome size)
- `THREADS`: 8 (default thread count)
- `MIN_LENGTH`: 1000 (minimum read length)
- `MIN_MEAN_Q`: 12 (minimum mean quality score)
- `MEDAKA_MODEL`: "r1041_e82_400bps_sup_v5.2.0"

### Directory Structure
```
input/                  # Raw FASTQ files (.fastq.gz, .fq.gz)
output/
├── 00_data/           # Copied raw data
├── 01_nanoplot/       # Quality control reports
│   ├── raw/          # Initial QC
│   └── filtered/     # Post-filtering QC
├── 01_data_filtered/  # Processed reads
│   └── fastq/        # Filtered FASTQ files
└── 02_assembly/       # Assembly results
    ├── flye/         # Initial assembly
    ├── medaka/       # Polished consensus
    └── quast/        # Quality evaluation
```

### Reference Genomes (references/)
- Candida albicans SC5314 (FASTA + GFF3)
- Candidozyma auris B11220 and B11221 (FASTA + GFF3)
- Nakaseomyces glabratus BG2 (FASTA + GFF3)

### Dependencies
The pipeline uses conda/bioconda packages including:
- Core tools: NanoPlot, Porechop, Filtlong, Flye, Medaka, Quast
- Classification tools: Kraken2, BLAST, Centrifuge
- Assembly tools: Hifiasm, TGS-GapCloser
- Utilities: vsearch, uchime, R packages for visualization

## Development Notes

### Sample Detection Logic (Snakefile:32-48)
The pipeline automatically detects samples from input directory by scanning for `.fastq.gz` and `.fq.gz` files, removing extensions to generate sample names.

### Error Handling
- NanoPlot failures are caught and logged but don't stop the pipeline
- Default retry mechanism: 2 retries with increasing memory allocation
- Memory scaling function: `attempt * 8000 MB`

### Data Analysis Context
Based on `resources/Analysis_long_reads_nanopore.md`, this pipeline addresses the challenge of analyzing "unclassified" nanopore reads (typically 10-20% of total reads, but up to 80-90% for reads >75kb). The methodology includes:
- Structural analysis for detecting chimeras and concatenated reads
- Taxonomic classification for contamination detection  
- Quality-based filtering optimized for long reads
- Assembly-focused workflow for genome reconstruction

### Python Helper Scripts
- `fastq_len_check.py`: Parallel processing of FASTQ files for read length analysis
- Uses multiprocessing to avoid Python GIL limitations
- Chunk-based processing for memory efficiency (400,000 lines per chunk)

## Typical Workflow

1. Place FASTQ files in `input/` directory
2. Activate conda environment: `conda activate nanocheck`
3. Run pipeline: `snakemake --cores 8`
4. Check results in `output/` subdirectories
5. Review HTML reports from NanoPlot and Quast for quality assessment

## Output Files of Interest
- `{sample}_porechop_filtlong.fastq.gz`: Final processed reads
- `assembly.fasta`: Initial Flye assembly
- `consensus.fasta`: Medaka-polished final assembly  
- `NanoPlot-report.html`: Quality control reports
- `report.html`: Quast assembly evaluation