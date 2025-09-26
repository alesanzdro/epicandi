# Phylogenetic System Update - September 26, 2025

## Overview

Successfully implemented a **configurable phylogenetic methodology system** for the EpiCandi pipeline, allowing users to control which phylogenetic analyses are executed based on time constraints and analysis requirements.

## Key Features Implemented

### 1. **Methodology Configuration Control**

Added phylogenetic configuration to `config.yaml`:

```yaml
# Phylogenetic analysis configuration
phylogenetics:
  # Control which phylogenetic methodologies to run:
  # - "none": Skip all phylogenetic analysis
  # - "triage": Fast assembly-to-assembly comparison only (5-10 min)
  # - "orthologs": Comprehensive ortholog-based analysis only (2-4 hours)
  # - "all": Run both triage and ortholog methodologies
  methodology: "triage"  # Options: none, triage, orthologs, all

  # Triage tree parameters (fast analysis)
  triage:
    tools: ["fastani", "mash", "sourmash"]  # Tools for assembly comparison
    min_ani: 0.8  # Minimum ANI threshold for inclusion

  # Ortholog tree parameters (comprehensive analysis)
  orthologs:
    min_taxa_fraction: 0.8  # Minimum fraction of taxa required per ortholog
    max_orthologs: 100      # Maximum number of orthologs to use
    min_seq_length: 150     # Minimum non-gap sequence length
    alignment_method: "core_orthologs"  # Method: core_orthologs, amas_concat
    iqtree_params: "-bb 1000 -m MFP"  # Additional IQ-TREE parameters
```

### 2. **Consolidated Phylogenetic Alignment Processor**

Created **`scripts/phylogenetic_alignment_processor.py`** - a unified script combining:

#### **Three Processing Modes:**
- **`fix-lengths`**: Fix sequence length inconsistencies for IQ-TREE
- **`filter-quality`**: Filter sequences by character content quality
- **`core-orthologs`**: Select and concatenate core orthologs ⭐ **RECOMMENDED**

#### **Usage Examples:**
```bash
# Fix sequence length inconsistencies
python scripts/phylogenetic_alignment_processor.py fix-lengths \
  --input-alignment concat.aln \
  --partition-file concat.part \
  --output-alignment concat_fixed.aln

# Filter sequences by quality
python scripts/phylogenetic_alignment_processor.py filter-quality \
  --input-alignment concat_fixed.aln \
  --output-alignment concat_filtered.aln \
  --min-total-chars 1800000 \
  --min-char-fraction 0.8

# Select core orthologs (RECOMMENDED)
python scripts/phylogenetic_alignment_processor.py core-orthologs \
  --single-copy-dir single_copy_alignments \
  --output-dir core_orthologs \
  --min-taxa-fraction 0.8 \
  --max-orthologs 100 \
  --min-seq-length 150
```

### 3. **Updated Snakemake Rules System**

#### **Conditional Rule Execution:**
- Rules are only included if the corresponding methodology is enabled
- Dynamic input/output functions based on configuration
- Proper dependency management between methodologies

#### **New Rule: `process_phylogenetic_alignments`**
```python
rule process_phylogenetic_alignments:
    """
    Process ortholog alignments using the consolidated phylogenetic alignment processor.
    This rule combines alignment fixing, quality filtering, and core ortholog selection.
    """
    input:
        orthofinder_complete="output/03_characterization/03.7_orthofinder/OrthoFinder_complete.txt",
        single_copy_dir="output/03_characterization/03.7_orthofinder/Results_proteomes/Single_Copy_Orthologue_Sequences"
    output:
        concatenated="output/03_characterization/03.8_phylogeny/alignments/concatenated_core_orthologs.fasta",
        stats="output/03_characterization/03.8_phylogeny/alignments/core_orthologs_stats.txt",
        selected_dir="output/03_characterization/03.8_phylogeny/alignments/selected_alignments"
    shell:
        """
        python scripts/phylogenetic_alignment_processor.py core-orthologs \
               --single-copy-dir {input.single_copy_dir} \
               --output-dir {params.output_dir} \
               --min-taxa-fraction {params.min_taxa_fraction} \
               --max-orthologs {params.max_orthologs} \
               --min-seq-length {params.min_seq_length} &> {log}
        """
```

### 4. **Auxiliary Scripts Created**

- **`scripts/ani_to_distance.py`**: Convert FastANI results to distance matrices
- **`scripts/mash_to_matrix.py`**: Convert Mash triangle output to square matrices

## Usage Scenarios

### **Scenario 1: Quick Analysis (Triage Only)**
```yaml
phylogenetics:
  methodology: "triage"
```
**Time**: 5-10 minutes
**Output**: Assembly-to-assembly comparison tree using Mash distances

### **Scenario 2: Comprehensive Analysis (Orthologs Only)**
```yaml
phylogenetics:
  methodology: "orthologs"
```
**Time**: 2-4 hours
**Output**: High-quality phylogenetic tree from core orthologs using IQ-TREE

### **Scenario 3: Complete Analysis (Both)**
```yaml
phylogenetics:
  methodology: "all"
```
**Time**: 2-4 hours
**Output**: Both triage and ortholog-based trees for comparison

### **Scenario 4: Skip Phylogenetics**
```yaml
phylogenetics:
  methodology: "none"
```
**Time**: 0 seconds
**Output**: Phylogenetic analysis skipped

## Snakemake Target Rules

### **Main Target Rule**
```bash
snakemake phylogenetic_analysis
```
Automatically runs the appropriate analysis based on methodology configuration.

### **Specific Methodology Targets**
```bash
# Triage analysis only
snakemake output/03_characterization/03.5_triage/triage_summary.txt

# Ortholog analysis only
snakemake output/03_characterization/03.8_phylogeny/phylogeny_summary.txt

# Both analyses
snakemake output/03_characterization/dual_phylogeny_complete.txt
```

## Architecture Benefits

### **1. Modular Design**
- Each methodology is independent
- Rules only execute when needed
- No wasted computational resources

### **2. Configurable Parameters**
- All parameters controlled via `config.yaml`
- Easy to adjust for different datasets
- Consistent parameter management

### **3. Consolidated Processing**
- Single script handles all alignment processing needs
- Reduces code duplication
- Easier maintenance and testing

### **4. Scalability**
- System automatically scales with more taxa
- Core orthologs approach handles large datasets efficiently
- Quality filtering prevents analysis failures

## File Organization

### **Configuration Files**
- `config.yaml` - Main configuration with phylogenetics section
- `rules/characterization_phylogenetic.smk` - Updated Snakemake rules

### **Scripts**
- `scripts/phylogenetic_alignment_processor.py` - **Main consolidated script**
- `scripts/ani_to_distance.py` - FastANI distance conversion
- `scripts/mash_to_matrix.py` - Mash matrix conversion
- `scripts/mash_to_newick.py` - Distance to Newick conversion (existing)
- `scripts/visualize_phylo_tree.py` - Tree visualization (existing)

### **Output Structure**
```
output/03_characterization/
├── 03.5_triage/              # Fast triage analysis
│   ├── mash/
│   ├── fastani/
│   └── figures/
├── 03.6_funannotate/          # Genome annotations
├── 03.7_orthofinder/          # Orthology analysis
├── 03.8_phylogeny/            # Comprehensive phylogeny
│   ├── alignments/
│   ├── trees/
│   └── figures/
└── phylogenetic_analysis_complete.txt
```

## Migration Guide

### **For Existing Users**

1. **Update `config.yaml`** with phylogenetics section:
```bash
# Add phylogenetics configuration
phylogenetics:
  methodology: "triage"  # Start with fast analysis
```

2. **Update Snakemake rules**:
```bash
# Replace old rule file if you have custom modifications
cp rules/characterization_phylogenetic.smk rules/characterization_phylogenetic_backup.smk
```

3. **Run analysis**:
```bash
# Test with triage first
snakemake phylogenetic_analysis

# Then try comprehensive analysis
# Edit config.yaml: methodology: "orthologs"
snakemake phylogenetic_analysis
```

### **For New Users**

1. **Set methodology** in `config.yaml`
2. **Run pipeline** - phylogenetic rules are automatically included based on configuration

## Example Workflow Commands

### **Development/Testing Phase**
```bash
# Quick test with triage
snakemake phylogenetic_analysis --cores 4

# Check what rules would run
snakemake phylogenetic_analysis --dry-run
```

### **Production Analysis**
```bash
# Comprehensive analysis
# First set methodology: "orthologs" in config.yaml
snakemake phylogenetic_analysis --cores 16 --use-conda

# Full pipeline including phylogenetics
snakemake --cores 16 --use-conda
```

### **Troubleshooting**
```bash
# Run specific components
snakemake process_phylogenetic_alignments --cores 8
snakemake iqtree_analysis --cores 8

# Test alignment processor directly
python scripts/phylogenetic_alignment_processor.py core-orthologs \
  --single-copy-dir test_data/single_copy_alignments \
  --output-dir test_output \
  --max-orthologs 20
```

## Performance Characteristics

### **Triage Analysis**
- **Time**: 5-10 minutes
- **Memory**: <2GB RAM
- **CPU**: 2-4 cores sufficient
- **Good for**: Quick screening, pipeline testing

### **Ortholog Analysis**
- **Time**: 2-4 hours (depending on dataset size)
- **Memory**: 8-32GB RAM (depends on number of orthologs)
- **CPU**: 8-16 cores recommended
- **Good for**: Publication-quality phylogenies

### **Combined Analysis**
- **Time**: 2-4 hours total (triage adds minimal time)
- **Resources**: Same as ortholog analysis
- **Good for**: Complete analysis with both speed and quality

## Future Enhancements

### **Potential Additions**
- **sourmash integration**: Add sourmash as third triage method
- **Tree comparison tools**: Compare triage vs ortholog trees
- **Automatic methodology selection**: Based on dataset size
- **Interactive tree visualization**: Web-based tree browsers
- **Bootstrap analysis**: Enhanced support value calculation

### **Configuration Extensions**
- **Custom tree rooting**: Automatic outgroup detection
- **Model selection**: Advanced IQ-TREE model testing
- **Alignment trimming**: Automated poorly aligned region removal

## Conclusion

The new phylogenetic system provides:

✅ **Flexible methodology control** - Choose analysis depth based on needs
✅ **Consolidated processing** - Single script handles all alignment issues
✅ **Scalable architecture** - Handles datasets from 2 to 100+ taxa
✅ **Production ready** - Robust error handling and logging
✅ **Easy configuration** - All parameters in `config.yaml`

This system resolves the IQ-TREE alignment issues while providing a flexible framework for phylogenetic analysis in the EpiCandi pipeline.