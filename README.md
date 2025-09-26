# EpiCandi – Comprehensive *Candida auris* Genomic Surveillance Pipeline

<div align="center">
  <img src="conf/logo_epicandi.png" alt="EpiCandi Logo" width="300"/>
</div>

**EpiCandi** is a state-of-the-art Snakemake-based bioinformatics pipeline designed for comprehensive whole genome sequencing (WGS) analysis of *Candida auris*. From raw sequencing reads to publication-ready phylogenetic trees, EpiCandi automates the complete genomic surveillance workflow with intelligent sample processing, quality-driven assembly selection, and dual phylogenetic strategies optimized for outbreak investigation and epidemiological studies.

## Key Features

- **Intelligent Assembly Strategy**: Automatic selection between Illumina-only, Nanopore-only, or hybrid approaches based on data availability
- **Quality-Driven Workflow**: Dynamic sample filtering with configurable QC thresholds ensures high-quality results
- **Comprehensive Characterization**: Integrated resistance profiling, clade classification, and phylogenetic analysis
- **Dual Phylogenetic Analysis**: Fast triage trees (5-10 min) for rapid outbreak assessment plus robust ortholog-based trees (2-4 hours) for detailed studies
- **Professional Reporting**: Automated MultiQC reports with customizable metadata and visualizations
- **Epidemiological Focus**: Specifically optimized for *C. auris* genomic surveillance and outbreak investigation

## Quick Start

```bash
# Setup environment and databases (one-time)
mamba env create -f envs/epicandi.yml
snakemake -s Snakefile_robust --cores 8 --use-conda setup_references

# Run complete pipeline
snakemake -s Snakefile_robust --cores 16 --use-conda --rerun-triggers mtime
```

## Pipeline Modules

### Sample Quality Control and Validation
EpiCandi implements intelligent sample filtering with configurable quality thresholds to ensure reliable downstream analysis:

- **Dynamic Sample Detection**: Automatic sample discovery from `samplesinfo.csv` or directory scanning
- **Quality Gating**: Configurable minimum read count thresholds (50K Illumina, 5K Nanopore)
- **Control Exclusion**: Automatic filtering of negative controls and blank samples
- **Contamination Assessment**: Kraken2 taxonomic classification for quality validation

### Intelligent Assembly Strategy Selection
The pipeline automatically determines the optimal assembly approach based on available sequencing data:

#### **A) Illumina-only Strategy**
- **Assembly**: SPAdes with diploid-aware parameters
- **Polishing**: Pilon error correction with multiple iterations
- **Variant Calling**: FreeBayes/BCFtools for heterozygous variants

#### **B) Nanopore-only Strategy**
- **Assembly**: Flye assembly optimized for fungal genomes
- **Polishing**: Medaka consensus with species-specific models
- **Variant Calling**: Clair3 optimized for long-read indel detection

#### **C) Hybrid Strategy (Gold Standard)**
- **Assembly**: Flye (Nanopore) for structural accuracy
- **Polishing**: Medaka + PyPolca (Illumina) for base-level precision
- **Variant Calling**: FreeBayes on Illumina reads for maximum accuracy

### Assembly Integrity and Quality Assessment
Comprehensive evaluation ensuring assembly reliability for downstream analysis:

- **BUSCO Analysis**: Conserved orthologous genes assessment (saccharomycetes_odb12)
- **CheckM2**: Genome completeness and contamination estimation
- **QUAST**: Assembly statistics and structural evaluation
- **Coverage Analysis**: BWA/minimap2 mapping with depth profiling

### Antifungal Resistance Profiling
Multi-layered approach to resistance characterization:

- **ChroQueTas**: Automated resistance gene detection and annotation
- **Custom Analysis**: Deep analysis of 5 key resistance genes (FUR1, FCY1, FKS1, CDR1, ERG11)
  - TBLASTN homology search with *Candida*-specific genetic code
  - Precise mutation calling with clinical resistance correlation
  - Coverage-aware confidence scoring

### Clade Classification
Accurate phylogenetic placement within *C. auris* global diversity:

- **AuriClass**: Machine learning-based clade assignment
- **FastANI**: Average Nucleotide Identity comparison with reference genomes
- **Reference Database**: Comprehensive collection spanning all major clades

### Dual Phylogenetic Analysis

#### **Triage Trees (Fast: ~5-10 minutes)**
Rapid outbreak assessment and initial classification:
- **FastANI**: Average Nucleotide Identity matrix
- **Mash**: K-mer based genomic distances
- **Sourmash**: Scalable k-mer comparison
- **Visualization**: Professional tree rendering with outbreak context

#### **Robust Phylogenetic Trees (Comprehensive: ~2-4 hours)**
Publication-quality phylogenetic analysis:
- **FunAnnotate**: Complete genome annotation with functional domains
- **OrthoFinder**: Orthologous gene family identification
- **Single-Copy Orthologs**: Concatenated alignment of conserved genes
- **IQ-TREE**: Maximum likelihood phylogeny with model selection
- **Bootstrap Support**: Statistical confidence assessment
- **Professional Visualization**: ETE3-based publication-ready figures

### Integrated Reporting
Comprehensive analysis summary with customizable metadata:

- **MultiQC Integration**: Unified HTML reports with sample metrics
- **Run Metadata**: Configurable sequencing platform and run information
- **Quality Metrics**: Assembly statistics, coverage, and contamination assessment
- **Resistance Summary**: Comprehensive mutation and gene presence tables

## Input Requirements

Create a `samplesinfo.csv` file defining your samples:

```csv
id,illumina_r1,illumina_r2,nanopore,dorado_model
EPI003361,data/illumina/sample_R1.fastq.gz,data/illumina/sample_R2.fastq.gz,NA,NA
EPI003102,data/illumina/sample_R1.fastq.gz,data/illumina/sample_R2.fastq.gz,data/nanopore/sample.fastq.gz,r1041_e82_400bps_sup_v5.2.0
EPI003135,NA,NA,data/nanopore/sample.fastq.gz,r1041_e82_400bps_sup_v5.2.0
```

## Output Structure

```
output/
├── 01_data/                    # Quality control and data processing
│   ├── 01.1_fastq_raw_qc/      # Raw data quality assessment
│   ├── 01.2_filtered/          # Processed and filtered reads
│   ├── 01.3_fastq_filtered_qc/ # Post-filtering quality control
│   ├── 01.4_kraken2/           # Taxonomic classification
│   └── 01.5_samples_pass/      # QC validation checkpoint
│
├── 02_assembly/                # Genome assembly and evaluation
│   ├── 02.1_draft_assembly/    # Initial assemblies (Flye/SPAdes)
│   ├── 02.2_polishing/         # Polished assemblies (Medaka/Pilon/PyPolca)
│   ├── 02.3_consensus/         # Final high-quality assemblies
│   ├── 02.4_evaluation/        # Assembly metrics (QUAST/BUSCO/CheckM2)
│   └── 02.5_coverage/          # Coverage and depth analysis
│
├── 03_characterization/        # Genomic characterization
│   ├── 03.1_fastani/           # ANI-based phylogenetic analysis
│   ├── 03.2_auriclass/         # Clade classification
│   ├── 03.3_chroquetas/        # Resistance gene analysis
│   ├── 03.4_resistance/        # Custom resistance profiling
│   ├── 03.5_triage/            # Fast triage phylogenetic trees
│   ├── 03.6_funannotate/       # Genome annotation
│   ├── 03.7_orthofinder/       # Orthology analysis
│   └── 03.8_phylogeny/         # Robust phylogenetic trees
│
└── 04_report/                  # Comprehensive reporting
    └── multiqc_report.html     # Integrated analysis dashboard
```

## Usage Examples

### Basic Analysis
```bash
# Analyze specific samples
snakemake -s Snakefile_robust --cores 16 --use-conda

# Generate triage trees only (fast outbreak assessment)
snakemake -s Snakefile_robust --cores 8 --use-conda triage_phylogeny

# Generate robust phylogenetic analysis only
snakemake -s Snakefile_robust --cores 24 --use-conda robust_phylogeny
```

### High-Performance Computing
```bash
# High-memory workloads (up to 389GB RAM)
snakemake -s Snakefile_robust --cores 64 --use-conda --resources mem_mb=389120

# Background execution with progress monitoring
screen -S epicandi -dm bash -lc 'snakemake -s Snakefile_robust --cores 32 --use-conda --printshellcmds'
```

## Software Dependencies

### **Quality Control and Preprocessing**
- **FastQC** (v0.12.1): Illumina read quality assessment
- **NanoPlot** (v1.42.0): Nanopore read quality visualization
- **fastp** (v0.23.4): Illumina read filtering and adapter trimming
- **Porechop** (v0.2.4): Nanopore adapter trimming
- **Filtlong** (v0.2.1): Nanopore read quality filtering
- **Kraken2** (v2.1.2): Taxonomic classification and contamination detection

### **Genome Assembly**
- **SPAdes** (v4.0.0): Illumina short-read assembly
- **Flye** (v2.9.3): Nanopore long-read assembly
- **Medaka** (v1.11.3): Nanopore assembly consensus polishing
- **Pilon** (v1.24): Illumina-based assembly polishing
- **PyPolca** (v0.3.1): Hybrid assembly polishing

### **Assembly Quality Assessment**
- **QUAST** (v5.2.0): Assembly statistics and evaluation
- **BUSCO** (v5.7.1): Conserved gene completeness assessment
- **CheckM2** (v1.0.2): Genome quality and contamination estimation
- **BWA** (v0.7.18): Read mapping for coverage analysis
- **minimap2** (v2.28): Long-read mapping
- **mosdepth** (v0.3.8): Fast coverage depth calculation

### **Variant Calling**
- **FreeBayes** (v1.3.7): Small variant detection
- **BCFtools** (v1.20): VCF manipulation and filtering
- **Clair3** (v1.0.9): Nanopore-optimized variant calling

### **Resistance and Clade Analysis**
- **AuriClass** (v0.5.4): *C. auris* clade classification
- **ChroQueTas** (v1.0): Antifungal resistance gene analysis
- **BLAST+** (v2.15.0): Homology searching for resistance genes
- **MAFFT** (v7.525): Multiple sequence alignment
- **EMBOSS** (v6.6.0): Pairwise sequence alignment

### **Phylogenetic Analysis**

#### **Triage Phylogeny (Fast)**
- **FastANI** (v1.34): Average Nucleotide Identity calculation
- **Mash** (v2.3): K-mer based genomic distance
- **sourmash** (v4.9.4): Scalable k-mer comparison
- **SciPy** (v1.16.2): Hierarchical clustering and tree construction

#### **Robust Phylogeny (Comprehensive)**
- **FunAnnotate** (v1.8.17): Genome annotation and gene prediction
- **OrthoFinder** (v3.1.0): Orthologous gene family identification
- **IQ-TREE** (v3.0.1): Maximum likelihood phylogenetic inference
- **ETE3** (v3.1.3): Phylogenetic tree visualization

### **Reporting and Visualization**
- **MultiQC** (v1.21): Unified analysis reporting
- **Matplotlib** (v3.8.2): Scientific plotting
- **Seaborn** (v0.13.0): Statistical data visualization

### **Workflow Management**
- **Snakemake** (v9.11.1): Workflow management system
- **conda/mamba**: Environment and dependency management
- **Docker**: Containerized tool execution (InterProScan)

## Documentation

- **[Directory Structure](docs/DIRECTORY_STRUCTURE.md)**: Detailed output organization
- **[Reference Genomes](docs/REFERENCES.md)**: Available *C. auris* reference collection
- **[Installation Guide](docs/INSTALL_CONDA.md)**: Environment setup instructions

## Use Cases

### **Outbreak Investigation**
- **Rapid Assessment**: Use triage phylogeny (5-10 min) for immediate outbreak characterization
- **Contact Tracing**: ANI-based clustering for transmission inference
- **Public Health Response**: Fast turnaround for epidemiological decision-making

### **Genomic Surveillance**
- **Resistance Monitoring**: Comprehensive antifungal resistance profiling
- **Clade Tracking**: Global phylogenetic context and emergence patterns
- **Quality Assurance**: Standardized assembly and annotation workflows

### **Research Applications**
- **Comparative Genomics**: Robust ortholog-based phylogenetic analysis
- **Evolution Studies**: High-resolution phylogenies with bootstrap support
- **Population Genetics**: Variant calling and diversity analysis

## Key Advantages

- **Speed**: Dual phylogenetic approach balances speed with accuracy
- **Precision**: *C. auris*-specific optimizations and reference databases
- **Reproducibility**: Standardized workflows with version-controlled environments
- **Flexibility**: Supports single-platform or hybrid sequencing approaches
- **Scalability**: Efficient resource management from laptop to HPC systems

## Citation

If you use EpiCandi in your research, please cite:

```
EpiCandi: A comprehensive Candida auris genomic surveillance pipeline
[Authors and publication details to be added]
```

## Contributing

We welcome contributions! Please see our contributing guidelines and submit issues or pull requests via GitHub.

## Support

For technical support and questions:
- **Issues**: [GitHub Issues](https://github.com/username/epicandi/issues)
- **Email**: [alejandro.sanz@fisabio.es]

---

**EpiCandi** - Empowering *Candida auris* genomic surveillance and outbreak investigation worldwide.
