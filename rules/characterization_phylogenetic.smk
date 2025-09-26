"""
Phylogenetic analysis rules for EpiCandi pipeline.

This module implements configurable phylogenetic strategies controlled by config['phylogenetics']['methodology']:
- "none": Skip all phylogenetic analysis
- "triage": Fast assembly-to-assembly comparison only (fastANI, sourmash, mash)
- "orthologs": Comprehensive ortholog-based phylogeny only (funannotate + OrthoFinder + IQ-TREE)
- "all": Run both triage and ortholog methodologies

Directory structure:
- 03.5_triage/: Fast triage analysis results
- 03.6_funannotate/: Genome annotations
- 03.7_orthofinder/: Orthology analysis
- 03.8_phylogeny/: Final phylogenetic trees

All rules use the epicandi_phylonew environment.
"""

# Get phylogenetic methodology from config
PHYLO_METHOD = config.get("phylogenetics", {}).get("methodology", "none")
PHYLO_ORTHOLOGS_CONFIG = config.get("phylogenetics", {}).get("orthologs", {})
PHYLO_TRIAGE_CONFIG = config.get("phylogenetics", {}).get("triage", {})

# Define conditional output functions
def get_phylogenetic_outputs():
    """Return appropriate outputs based on methodology configuration."""
    outputs = []

    if PHYLO_METHOD in ["triage", "all"]:
        outputs.append("output/03_characterization/03.5_triage/triage_summary.txt")

    if PHYLO_METHOD in ["orthologs", "all"]:
        outputs.append("output/03_characterization/03.8_phylogeny/phylogeny_summary.txt")

    return outputs if outputs else ["output/03_characterization/no_phylogeny.txt"]

def get_triage_inputs():
    """Get inputs for triage phylogeny if enabled."""
    if PHYLO_METHOD in ["triage", "all"]:
        return expand("output/02_assembly/02.3_consensus/{sample}.fasta", sample=get_samples())
    return []

def get_ortholog_inputs():
    """Get inputs for ortholog phylogeny if enabled."""
    if PHYLO_METHOD in ["orthologs", "all"]:
        return expand("output/02_assembly/02.3_consensus/{sample}.fasta", sample=get_samples())
    return []

# =============================================================================
# PHYLOGENETIC ANALYSIS TARGET RULE
# =============================================================================

rule phylogenetic_analysis:
    """
    Main target rule for phylogenetic analysis based on methodology configuration.
    """
    input:
        get_phylogenetic_outputs()
    output:
        "output/03_characterization/phylogenetic_analysis_complete.txt"
    shell:
        """
        echo "Phylogenetic analysis complete using methodology: {PHYLO_METHOD}" > {output}
        echo "Date: $(date)" >> {output}
        """

rule no_phylogeny:
    """
    Placeholder rule when phylogenetic analysis is disabled.
    """
    output:
        "output/03_characterization/no_phylogeny.txt"
    shell:
        """
        echo "Phylogenetic analysis skipped (methodology: {PHYLO_METHOD})" > {output}
        echo "Date: $(date)" >> {output}
        """

# =============================================================================
# SETUP AND ENVIRONMENT CHECK
# =============================================================================

if PHYLO_METHOD in ["orthologs", "all"]:
    rule setup_phylogenetic_env:
        """
        Setup FunAnnotate database and environment.
        This rule must complete successfully before any phylogenetic analysis.
        """
        output:
            setup_flag="output/03_characterization/00_setup_phylo/.setup_complete"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["funannotate_setup"]
        log:
            "logs/setup_phylogenetic_env.log"
        params:
            project_dir=lambda wildcards: os.getcwd(),
            genemark_dir=config.get("GENEMARK_DIR", "")
        shell:
            """
            mkdir -p output/03_characterization/00_setup_phylo

            python scripts/setup_funannotate.py {params.project_dir} \
                   --genemark-dir "{params.genemark_dir}" &> {log}

            touch {output.setup_flag}
            """

# =============================================================================
# TRIAGE PHYLOGENY - FAST ASSEMBLY-TO-ASSEMBLY COMPARISON
# =============================================================================

if PHYLO_METHOD in ["triage", "all"]:
    rule triage_fastani:
        """
        Calculate Average Nucleotide Identity (ANI) between assemblies using FastANI.
        """
        input:
            assemblies=get_triage_inputs()
        output:
            ani_results="output/03_characterization/03.5_triage/fastani/ani_results.txt",
            distance_matrix="output/03_characterization/03.5_triage/fastani/distance_matrix.txt"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["fastani"]
        log:
            "logs/triage_fastani.log"
        params:
            outdir="output/03_characterization/03.5_triage/fastani",
            min_ani=PHYLO_TRIAGE_CONFIG.get("min_ani", 0.8)
        shell:
            """
            mkdir -p {params.outdir}

            # Create query and reference lists
            echo {input.assemblies} | tr ' ' '\\n' > {params.outdir}/query_list.txt
            echo {input.assemblies} | tr ' ' '\\n' > {params.outdir}/reference_list.txt

            # Run FastANI all-vs-all comparison
            fastANI --ql {params.outdir}/query_list.txt --rl {params.outdir}/reference_list.txt \
                    -o {output.ani_results} --matrix &> {log}

            # Convert to distance matrix (1 - ANI/100)
            python scripts/ani_to_distance.py {output.ani_results} {output.distance_matrix} \
                   --min-ani {params.min_ani} &>> {log}
            """

    rule triage_mash:
        """
        Calculate genomic distances using Mash sketching algorithm.
        """
        input:
            assemblies=get_triage_inputs()
        output:
            mash_distances="output/03_characterization/03.5_triage/mash/mash_distances.txt",
            distance_matrix="output/03_characterization/03.5_triage/mash/distance_matrix.txt"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["mash"]
        log:
            "logs/triage_mash.log"
        params:
            outdir="output/03_characterization/03.5_triage/mash",
            sketch_size=config.get("mash_sketch_size", 1000)
        shell:
            """
            mkdir -p {params.outdir}

            # Create sketches and calculate distances
            mash triangle -s {params.sketch_size} {input.assemblies} > {output.mash_distances} 2> {log}

            # Convert to square distance matrix
            python scripts/mash_to_matrix.py {output.mash_distances} {output.distance_matrix} &>> {log}
            """

    rule triage_tree:
        """
        Generate triage phylogenetic tree from distance matrix using hierarchical clustering.
        """
        input:
            distance_matrix="output/03_characterization/03.5_triage/mash/distance_matrix.txt"
        output:
            tree="output/03_characterization/03.5_triage/mash/triage_tree.newick"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["default"]
        log:
            "logs/triage_tree.log"
        shell:
            """
            python scripts/mash_to_newick.py {input.distance_matrix} {output.tree} &> {log}
            """

    rule triage_visualize:
        """
        Visualize triage phylogenetic tree with professional styling.
        """
        input:
            tree="output/03_characterization/03.5_triage/mash/triage_tree.newick"
        output:
            tree_png="output/03_characterization/03.5_triage/figures/triage_tree.png",
            tree_svg="output/03_characterization/03.5_triage/figures/triage_tree.svg",
            tree_pdf="output/03_characterization/03.5_triage/figures/triage_tree.pdf"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["ete3"]
        log:
            "logs/triage_visualize_tree.log"
        params:
            output_prefix="output/03_characterization/03.5_triage/figures/triage_tree"
        shell:
            """
            mkdir -p output/03_characterization/03.5_triage/figures

            python scripts/visualize_phylo_tree.py {input.tree} {params.output_prefix} \
                   --tree-type triage &> {log}
            """

    rule triage_summary:
        """
        Generate summary of triage phylogenetic analysis.
        """
        input:
            tree_png="output/03_characterization/03.5_triage/figures/triage_tree.png",
            tree="output/03_characterization/03.5_triage/mash/triage_tree.newick",
            distance_matrix="output/03_characterization/03.5_triage/mash/distance_matrix.txt"
        output:
            summary="output/03_characterization/03.5_triage/triage_summary.txt"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["default"]
        log:
            "logs/triage_phylogeny_summary.log"
        shell:
            """
            echo "Triage Phylogenetic Analysis Complete" > {output.summary}
            echo "====================================" >> {output.summary}
            echo "Date: $(date)" >> {output.summary}
            echo "" >> {output.summary}
            echo "Analysis method: Assembly-to-assembly comparison" >> {output.summary}
            echo "Tools used: Mash distance calculation + hierarchical clustering" >> {output.summary}
            echo "Analysis time: ~5-10 minutes" >> {output.summary}
            echo "" >> {output.summary}
            echo "Generated outputs:" >> {output.summary}
            echo "- Distance matrix: {input.distance_matrix}" >> {output.summary}
            echo "- Phylogenetic tree: {input.tree}" >> {output.summary}
            echo "- Tree visualization: {input.tree_png}" >> {output.summary}
            """

# =============================================================================
# ORTHOLOG-BASED PHYLOGENY - COMPREHENSIVE ANALYSIS
# =============================================================================

if PHYLO_METHOD in ["orthologs", "all"]:
    rule funannotate_sort:
        """
        Sort and clean contigs for annotation. Removes contigs shorter than 500bp.
        """
        input:
            setup_flag="output/03_characterization/00_setup_phylo/.setup_complete",
            assembly="output/02_assembly/02.3_consensus/{sample}.fasta"
        output:
            sorted="output/03_characterization/03.6_funannotate/{sample}/{sample}_clean.fasta"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["default"]
        log:
            "logs/funannotate_sort/{sample}.log"
        params:
            outdir="output/03_characterization/03.6_funannotate/{sample}"
        shell:
            """
            mkdir -p {params.outdir}
            cd {params.outdir}

            funannotate sort -i ../../../{input.assembly} \
                    -o {wildcards.sample}_clean.fasta --minlen 500 &> ../../../{log}
            """

    rule funannotate_mask:
        """
        Mask repetitive elements in genome assemblies using RepeatModeler/RepeatMasker.
        """
        input:
            sorted="output/03_characterization/03.6_funannotate/{sample}/{sample}_clean.fasta"
        output:
            masked="output/03_characterization/03.6_funannotate/{sample}/{sample}_clean_masked.fasta"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["funannotate_predict"]
        log:
            "logs/funannotate_mask/{sample}.log"
        params:
            outdir="output/03_characterization/03.6_funannotate/{sample}"
        shell:
            """
            cd {params.outdir}

            funannotate mask -i {wildcards.sample}_clean.fasta \
                    -o {wildcards.sample}_clean_masked.fasta \
                    --cpus {resources.threads} &> ../../../{log}
            """

    rule funannotate_predict:
        """
        Predict protein-coding genes using multiple gene prediction algorithms.
        """
        input:
            masked="output/03_characterization/03.6_funannotate/{sample}/{sample}_clean_masked.fasta"
        output:
            predict_complete="output/03_characterization/03.6_funannotate/{sample}/predict_misc/predict.complete"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["funannotate_predict"]
        log:
            "logs/funannotate_predict/{sample}.log"
        params:
            outdir="output/03_characterization/03.6_funannotate/{sample}",
            species="Candida auris",
            busco_db="saccharomycetes_odb10"
        shell:
            """
            cd {params.outdir}

            funannotate predict -i {wildcards.sample}_clean_masked.fasta -o . \
                    --species "{params.species}" --strain {wildcards.sample} \
                    --busco_db {params.busco_db} --cpus {resources.threads} &> ../../../{log}

            # Create completion flag
            touch predict_misc/predict.complete
            """

    rule orthofinder_analysis:
        """
        Run OrthoFinder to identify orthologous genes across all annotated genomes.
        """
        input:
            predictions=expand("output/03_characterization/03.6_funannotate/{sample}/predict_misc/predict.complete",
                             sample=get_ortholog_inputs() and [s.split('/')[-1].replace('.fasta', '') for s in get_ortholog_inputs()] or [])
        output:
            orthofinder_complete="output/03_characterization/03.7_orthofinder/OrthoFinder_complete.txt",
            single_copy_dir="output/03_characterization/03.7_orthofinder/Results_proteomes/Single_Copy_Orthologue_Sequences"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["orthofinder"]
        log:
            "logs/orthofinder_analysis.log"
        params:
            input_dir="output/03_characterization/03.6_funannotate",
            output_dir="output/03_characterization/03.7_orthofinder"
        shell:
            """
            mkdir -p {params.output_dir}

            # Create proteome directory and copy protein sequences
            mkdir -p {params.output_dir}/proteomes

            for sample_dir in {params.input_dir}/*/; do
                sample=$(basename "$sample_dir")
                if [ -f "$sample_dir/predict_results/${{sample}}.proteins.fa" ]; then
                    cp "$sample_dir/predict_results/${{sample}}.proteins.fa" {params.output_dir}/proteomes/
                fi
            done

            # Run OrthoFinder
            cd {params.output_dir}
            orthofinder -f proteomes -t {resources.threads} &> ../../{log}

            # Create completion flag
            echo "OrthoFinder analysis completed at $(date)" > OrthoFinder_complete.txt
            """

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
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["iqtree"]
        log:
            "logs/process_phylogenetic_alignments.log"
        params:
            output_dir="output/03_characterization/03.8_phylogeny/alignments",
            min_taxa_fraction=PHYLO_ORTHOLOGS_CONFIG.get("min_taxa_fraction", 0.8),
            max_orthologs=PHYLO_ORTHOLOGS_CONFIG.get("max_orthologs", 100),
            min_seq_length=PHYLO_ORTHOLOGS_CONFIG.get("min_seq_length", 150)
        shell:
            """
            mkdir -p {params.output_dir}

            python scripts/phylogenetic_alignment_processor.py core-orthologs \
                   --single-copy-dir {input.single_copy_dir} \
                   --output-dir {params.output_dir} \
                   --min-taxa-fraction {params.min_taxa_fraction} \
                   --max-orthologs {params.max_orthologs} \
                   --min-seq-length {params.min_seq_length} &> {log}
            """

    rule iqtree_analysis:
        """
        Generate maximum likelihood phylogenetic tree using IQ-TREE.
        """
        input:
            concatenated="output/03_characterization/03.8_phylogeny/alignments/concatenated_core_orthologs.fasta"
        output:
            tree="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.treefile",
            iqtree_complete="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.complete"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["iqtree"]
        log:
            "logs/iqtree_analysis.log"
        params:
            output_dir="output/03_characterization/03.8_phylogeny/trees",
            output_prefix="cauris_iqtree",
            iqtree_params=PHYLO_ORTHOLOGS_CONFIG.get("iqtree_params", "-bb 1000 -m MFP")
        shell:
            """
            mkdir -p {params.output_dir}
            cd {params.output_dir}

            iqtree -s ../alignments/concatenated_core_orthologs.fasta \
                   -pre {params.output_prefix} \
                   {params.iqtree_params} \
                   -nt {resources.threads} &> ../../../{log}

            # Create completion flag
            echo "IQ-TREE analysis completed at $(date)" > {params.output_prefix}.complete
            """

    rule robust_visualize:
        """
        Visualize robust phylogenetic tree with professional styling.
        """
        input:
            tree="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.treefile"
        output:
            tree_png="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.png",
            tree_svg="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.svg",
            tree_pdf="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.pdf"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["ete3"]
        log:
            "logs/robust_visualize_tree.log"
        params:
            output_prefix="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny"
        shell:
            """
            mkdir -p output/03_characterization/03.8_phylogeny/figures

            python scripts/visualize_phylo_tree.py {input.tree} {params.output_prefix} \
                   --tree-type phylogeny &> {log}
            """

    rule robust_phylogeny_summary:
        """
        Generate summary of robust phylogenetic analysis.
        """
        input:
            tree_png="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.png",
            concatenated="output/03_characterization/03.8_phylogeny/alignments/concatenated_core_orthologs.fasta",
            tree="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.treefile"
        output:
            summary="output/03_characterization/03.8_phylogeny/phylogeny_summary.txt"
        conda:
            "../envs/epicandi_phylonew.yml"
        resources:
            **config["resources"]["default"]
        log:
            "logs/robust_phylogeny_summary.log"
        shell:
            """
            echo "Robust Phylogenetic Analysis Complete" > {output.summary}
            echo "=====================================" >> {output.summary}
            echo "Date: $(date)" >> {output.summary}
            echo "" >> {output.summary}
            echo "Analysis pipeline:" >> {output.summary}
            echo "1. FunAnnotate gene prediction and annotation" >> {output.summary}
            echo "2. OrthoFinder orthology inference" >> {output.summary}
            echo "3. Core orthologs selection and concatenation" >> {output.summary}
            echo "4. IQ-TREE maximum likelihood phylogeny" >> {output.summary}
            echo "5. ETE3 tree visualization" >> {output.summary}
            echo "" >> {output.summary}
            echo "Generated outputs:" >> {output.summary}
            echo "- Concatenated alignment: {input.concatenated}" >> {output.summary}
            echo "- ML tree: {input.tree}" >> {output.summary}
            echo "- Tree visualization: {input.tree_png}" >> {output.summary}
            """

# =============================================================================
# DUAL PHYLOGENY RULE (when methodology = "all")
# =============================================================================

if PHYLO_METHOD == "all":
    rule dual_phylogeny:
        """
        Complete both triage and robust phylogenetic analyses.
        """
        input:
            triage_summary="output/03_characterization/03.5_triage/triage_summary.txt",
            phylogeny_summary="output/03_characterization/03.8_phylogeny/phylogeny_summary.txt"
        output:
            combined_summary="output/03_characterization/dual_phylogeny_complete.txt"
        shell:
            """
            echo "Dual Phylogenetic Analysis Complete" > {output.combined_summary}
            echo "===================================" >> {output.combined_summary}
            echo "Date: $(date)" >> {output.combined_summary}
            echo "" >> {output.combined_summary}
            echo "Both triage and robust phylogenetic analyses completed successfully." >> {output.combined_summary}
            echo "" >> {output.combined_summary}
            echo "Triage analysis (fast): Assembly-to-assembly comparison" >> {output.combined_summary}
            echo "Robust analysis (comprehensive): Ortholog-based phylogeny" >> {output.combined_summary}
            """