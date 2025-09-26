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

# =============================================================================
# SETUP AND ENVIRONMENT CHECK
# =============================================================================

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

rule triage_prepare_assemblies:
    """
    Copy clean assemblies to triage directory for analysis.
    """
    input:
        setup_flag="output/03_characterization/00_setup_phylo/.setup_complete",
        assembly="output/02_assembly/02.3_consensus/{sample}.fasta"
    output:
        assembly="output/03_characterization/03.5_triage/assemblies/{sample}.fasta"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["default"]
    log:
        "logs/triage_prepare_assemblies/{sample}.log"
    shell:
        """
        cp {input.assembly} {output.assembly} 2> {log}
        """

rule triage_fastani_matrix:
    """
    Calculate Average Nucleotide Identity (ANI) matrix using fastANI.
    """
    input:
        assemblies=expand("output/03_characterization/03.5_triage/assemblies/{sample}.fasta", sample=SAMPLES)
    output:
        matrix="output/03_characterization/03.5_triage/fastani/ani_matrix.txt",
        file_list="output/03_characterization/03.5_triage/fastani/genome_list.txt"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["fastani"]
    log:
        "logs/triage_fastani_matrix.log"
    shell:
        """
        mkdir -p output/03_characterization/03.5_triage/fastani

        # Create absolute path list of assemblies
        ls -d "$PWD"/output/03_characterization/03.5_triage/assemblies/*.fasta > {output.file_list}

        # Run fastANI with matrix output
        fastANI --ql {output.file_list} --rl {output.file_list} \
                -o {output.matrix} --matrix &> {log}
        """

rule triage_sourmash_sketch:
    """
    Create k-mer sketches for rapid genomic comparison using sourmash.
    """
    input:
        assembly="output/03_characterization/03.5_triage/assemblies/{sample}.fasta"
    output:
        sketch="output/03_characterization/03.5_triage/sourmash/sketches/{sample}.sig"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["sourmash"]
    log:
        "logs/triage_sourmash_sketch/{sample}.log"
    params:
        k=31,
        scale=1000
    shell:
        """
        mkdir -p output/03_characterization/03.5_triage/sourmash/sketches

        sourmash sketch dna -p k={params.k},scaled={params.scale} \
                -o {output.sketch} {input.assembly} &> {log}
        """

rule triage_sourmash_compare:
    """
    Compare k-mer sketches to generate similarity matrix.
    """
    input:
        sketches=expand("output/03_characterization/03.5_triage/sourmash/sketches/{sample}.sig", sample=SAMPLES)
    output:
        matrix="output/03_characterization/03.5_triage/sourmash/sourmash_matrix.csv"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["sourmash"]
    log:
        "logs/triage_sourmash_compare.log"
    shell:
        """
        mkdir -p output/03_characterization/03.5_triage/sourmash

        sourmash compare -o {output.matrix} --csv {input.sketches} \
                -k 31 -p {resources.threads} &> {log}
        """

rule triage_mash_triangle:
    """
    Generate triangular distance matrix using Mash for phylogenetic inference.
    """
    input:
        assemblies=expand("output/03_characterization/03.5_triage/assemblies/{sample}.fasta", sample=SAMPLES)
    output:
        distances="output/03_characterization/03.5_triage/mash/mash_distances.txt"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["mash"]
    log:
        "logs/triage_mash_triangle.log"
    shell:
        """
        mkdir -p output/03_characterization/03.5_triage/mash

        mash triangle -p {resources.threads} {input.assemblies} > {output.distances} 2> {log}
        """

rule triage_mash_to_newick:
    """
    Convert Mash distance matrix to Newick phylogenetic tree format.
    """
    input:
        distances="output/03_characterization/03.5_triage/mash/mash_distances.txt"
    output:
        tree="output/03_characterization/03.5_triage/mash/tree.nwk"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["python_script"]
    log:
        "logs/triage_mash_to_newick.log"
    shell:
        """
        python scripts/mash_to_newick.py {input.distances} {output.tree} &> {log}
        """

rule triage_visualize_tree:
    """
    Create publication-ready visualizations of the triage phylogenetic tree.
    """
    input:
        tree="output/03_characterization/03.5_triage/mash/tree.nwk"
    output:
        png="output/03_characterization/03.5_triage/mash/tree.png",
        svg="output/03_characterization/03.5_triage/mash/tree.svg",
        pdf="output/03_characterization/03.5_triage/mash/tree.pdf"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["ete3"]
    log:
        "logs/triage_visualize_tree.log"
    params:
        output_prefix="output/03_characterization/03.5_triage/mash/tree"
    shell:
        """
        python scripts/visualize_phylo_tree.py {input.tree} {params.output_prefix} \
               --tree-type triage &> {log}
        """

# =============================================================================
# FUNANNOTATE - GENOME ANNOTATION
# =============================================================================

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

rule funannotate_iprscan:
    """
    Perform functional domain annotation using InterProScan.
    This is the most computationally intensive step.
    """
    input:
        predict_complete="output/03_characterization/03.6_funannotate/{sample}/predict_misc/predict.complete"
    output:
        iprscan_complete="output/03_characterization/03.6_funannotate/{sample}/annotate_misc/iprscan.complete"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["funannotate_iprscan"]
    log:
        "logs/funannotate_iprscan/{sample}.log"
    params:
        outdir="output/03_characterization/03.6_funannotate/{sample}"
    shell:
        """
        cd {params.outdir}

        funannotate iprscan -i . -m docker --cpus {resources.threads} &> ../../../{log}

        # Create completion flag
        mkdir -p annotate_misc
        touch annotate_misc/iprscan.complete
        """

rule funannotate_annotate:
    """
    Integrate functional annotations from multiple databases.
    """
    input:
        iprscan_complete="output/03_characterization/03.6_funannotate/{sample}/annotate_misc/iprscan.complete"
    output:
        annotation_complete="output/03_characterization/03.6_funannotate/{sample}/annotate_misc/annotation.complete",
        proteins="output/03_characterization/03.6_funannotate/{sample}/predict_results/{sample}.proteins.fa"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["funannotate_annotate"]
    log:
        "logs/funannotate_annotate/{sample}.log"
    params:
        outdir="output/03_characterization/03.6_funannotate/{sample}",
        busco_db="saccharomycetes_odb10"
    shell:
        """
        cd {params.outdir}

        funannotate annotate -i . --busco_db {params.busco_db} \
                --cpus {resources.threads} &> ../../../{log}

        # Create completion flag
        touch annotate_misc/annotation.complete

        # Ensure proteins file exists (create symlink if needed)
        if [ ! -f predict_results/{wildcards.sample}.proteins.fa ]; then
            find predict_results -name "*.proteins.fa" | head -1 | xargs -I{{}} ln -sf {{}} predict_results/{wildcards.sample}.proteins.fa
        fi
        """

# =============================================================================
# ORTHOFINDER - ORTHOLOGY ANALYSIS
# =============================================================================

rule orthofinder_prepare_proteomes:
    """
    Extract and prepare protein sequences from FunAnnotate results for orthology analysis.
    Also download reference genomes and outgroup species.
    """
    input:
        annotations=expand("output/03_characterization/03.6_funannotate/{sample}/annotate_misc/annotation.complete", sample=SAMPLES)
    output:
        proteome_ready="output/03_characterization/03.7_orthofinder/proteomes/.proteomes_ready"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["default"]
    log:
        "logs/orthofinder_prepare_proteomes.log"
    params:
        samples=SAMPLES,
        proteome_dir="output/03_characterization/03.7_orthofinder/proteomes"
    shell:
        """
        mkdir -p {params.proteome_dir}

        # Extract proteins from FunAnnotate results
        for sample in {params.samples}; do
            echo "Processing $sample..." >> {log}

            protein_file=output/03_characterization/03.6_funannotate/$sample/predict_results/$sample.proteins.fa
            if [ -f "$protein_file" ]; then
                cp "$protein_file" {params.proteome_dir}/$sample.faa
                # Clean headers for OrthoFinder compatibility
                sed -i 's/ .*//' {params.proteome_dir}/$sample.faa
                echo "✓ $sample: $(grep -c '>' {params.proteome_dir}/$sample.faa) proteins" >> {log}
            else
                echo "✗ Warning: No proteins found for $sample" >> {log}
            fi
        done

        # Download reference genomes
        cd {params.proteome_dir}

        # C. auris reference (Clade I - South Asian)
        echo "Downloading C. auris B8441 reference..." >> ../../../{log}
        wget -q -O B8441_ref.faa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/775/015/GCF_002775015.1_Cand_auris_B8441_V2/GCF_002775015.1_Cand_auris_B8441_V2_protein.faa.gz" 2>> ../../../{log}
        gunzip B8441_ref.faa.gz 2>> ../../../{log}

        # C. haemulonii as outgroup
        echo "Downloading C. haemulonii outgroup..." >> ../../../{log}
        wget -q -O Chaemulonii_out.faa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/926/085/GCF_002926085.2_Cand_haemulonii/GCF_002926085.2_Cand_haemulonii_protein.faa.gz" 2>> ../../../{log}
        gunzip Chaemulonii_out.faa.gz 2>> ../../../{log}

        # Create completion flag
        touch .proteomes_ready

        echo "Proteome preparation completed" >> ../../../{log}
        """

rule orthofinder_analysis:
    """
    Identify orthologous gene groups using OrthoFinder with comprehensive analysis.
    """
    input:
        proteome_ready="output/03_characterization/03.7_orthofinder/proteomes/.proteomes_ready"
    output:
        results_flag="output/03_characterization/03.7_orthofinder/results/.orthofinder_complete"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["orthofinder"]
    log:
        "logs/orthofinder_analysis.log"
    params:
        proteome_dir="output/03_characterization/03.7_orthofinder/proteomes",
        results_dir="output/03_characterization/03.7_orthofinder/results"
    shell:
        """
        mkdir -p {params.results_dir}

        orthofinder -f {params.proteome_dir} \
                -t {resources.threads} \
                -M msa \
                -S diamond_ultra_sens \
                -A mafft \
                -T fasttree \
                -o {params.results_dir} \
                &> {log}

        # Create completion flag
        touch {params.results_dir}/.orthofinder_complete
        """

rule phylogeny_concatenate_orthologs:
    """
    Extract and concatenate single-copy orthologous sequences for phylogenetic analysis.
    """
    input:
        orthofinder_complete="output/03_characterization/03.7_orthofinder/results/.orthofinder_complete"
    output:
        concatenated="output/03_characterization/03.8_phylogeny/alignments/concatenated_orthologs.fasta",
        statistics="output/03_characterization/03.8_phylogeny/alignments/concatenation_stats.txt"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["python_script"]
    log:
        "logs/phylogeny_concatenate_orthologs.log"
    params:
        results_dir="output/03_characterization/03.7_orthofinder/results",
        output_dir="output/03_characterization/03.8_phylogeny/alignments",
        max_orthologs=100
    shell:
        """
        mkdir -p {params.output_dir}

        python scripts/concatenate_orthologs.py {params.results_dir} {params.output_dir} \
               --max-orthologs {params.max_orthologs} &> {log}
        """

rule phylogeny_iqtree:
    """
    Construct maximum likelihood phylogenetic tree using IQ-TREE with model selection.
    """
    input:
        concatenated="output/03_characterization/03.8_phylogeny/alignments/concatenated_orthologs.fasta"
    output:
        tree="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.treefile",
        log_file="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.log"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["iqtree"]
    log:
        "logs/phylogeny_iqtree.log"
    params:
        tree_dir="output/03_characterization/03.8_phylogeny/trees",
        prefix="cauris_iqtree"
    shell:
        """
        mkdir -p {params.tree_dir}
        cd {params.tree_dir}

        iqtree -s ../../alignments/concatenated_orthologs.fasta \
               -m MFP \
               -bb 1000 \
               -nt AUTO \
               -ntmax {resources.threads} \
               -pre {params.prefix} \
               -quiet &> ../../../{log}
        """

rule phylogeny_visualize_tree:
    """
    Create publication-ready visualizations of the robust phylogenetic tree.
    """
    input:
        tree="output/03_characterization/03.8_phylogeny/trees/cauris_iqtree.treefile"
    output:
        png="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.png",
        svg="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.svg",
        pdf="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.pdf"
    conda:
        "../envs/epicandi_phylonew.yml"
    resources:
        **config["resources"]["ete3"]
    log:
        "logs/phylogeny_visualize_tree.log"
    params:
        output_prefix="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny"
    shell:
        """
        mkdir -p output/03_characterization/03.8_phylogeny/figures

        python scripts/visualize_phylo_tree.py {input.tree} {params.output_prefix} \
               --tree-type phylogeny &> {log}
        """

# =============================================================================
# AGGREGATE RULES
# =============================================================================

rule triage_phylogeny:
    """
    Complete triage phylogenetic analysis - fast assembly-to-assembly comparison.
    """
    input:
        ani_matrix="output/03_characterization/03.5_triage/fastani/ani_matrix.txt",
        sourmash_matrix="output/03_characterization/03.5_triage/sourmash/sourmash_matrix.csv",
        tree_png="output/03_characterization/03.5_triage/mash/tree.png"
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
        echo "=====================================" >> {output.summary}
        echo "Date: $(date)" >> {output.summary}
        echo "Samples analyzed: $(echo '{input}' | tr ' ' '\\n' | wc -l)" >> {output.summary}
        echo "" >> {output.summary}
        echo "Generated outputs:" >> {output.summary}
        echo "- fastANI matrix: {input.ani_matrix}" >> {output.summary}
        echo "- Sourmash matrix: {input.sourmash_matrix}" >> {output.summary}
        echo "- Triage tree: {input.tree_png}" >> {output.summary}
        """

rule robust_phylogeny:
    """
    Complete robust phylogenetic analysis - comprehensive ortholog-based approach.
    """
    input:
        tree_png="output/03_characterization/03.8_phylogeny/figures/cauris_phylogeny.png",
        concatenated="output/03_characterization/03.8_phylogeny/alignments/concatenated_orthologs.fasta",
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
        echo "3. IQ-TREE maximum likelihood phylogeny" >> {output.summary}
        echo "4. ETE3 tree visualization" >> {output.summary}
        echo "" >> {output.summary}
        echo "Generated outputs:" >> {output.summary}
        echo "- Concatenated alignment: {input.concatenated}" >> {output.summary}
        echo "- ML tree: {input.tree}" >> {output.summary}
        echo "- Tree visualization: {input.tree_png}" >> {output.summary}
        """

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