// EXTRACT_GENE_PROTEINS — assembly-based per-gene protein extraction with
// FungAMR cross-reference.
//
// Why assembly-based and not VCF-aware:
//   The smoke run showed that GATK4's filtered VCF can drop AMR-causing
//   mutations like ERG11 Y132F (stringent priors / mask), even though they
//   are present in the polished assembly. ChroQueTaS reaches them by working
//   on the assembly directly — we do the same here.
//
// Steps per (sample, gene):
//   1. miniprot canonical_<gene>.faa → assembly.fasta  → contig:start-end + strand
//   2. samtools faidx region (+ revcomp on strand -)
//   3. translate with genetic code 12 (Yeast Alternative — Candida CUG → Ser)
//   4. global pairwise align vs canonical protein (BLOSUM62)
//   5. cross every mismatch with mutation_catalog_auris.tsv → classify
//      known_resistance / known_sensitivity / unknown_missense / LoF
//
// Outputs:
//   <sample>.gene_proteins.tsv          long format, one row per mismatch
//   <sample>.gene_proteins_summary.tsv  one row per (sample, gene)
//   alignments/<sample>__<gene>.aa_alignment.fasta
//
process EXTRACT_GENE_PROTEINS {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::miniprot bioconda::samtools conda-forge::biopython'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d65f6e9ce8b5c5d6:abc' :
        'community.wave.seqera.io/library/miniprot_samtools_biopython:abc' }"

    input:
    tuple val(meta), path(assembly)
    path canonical_proteins_dir
    path panel_tsv
    path mutation_catalog

    output:
    tuple val(meta), path("${meta.id}.gene_proteins.tsv"),         emit: long
    tuple val(meta), path("${meta.id}.gene_proteins_summary.tsv"), emit: summary
    tuple val(meta), path("alignments/"),                          emit: alignments
    tuple val("${task.process}"), val('miniprot'),
          eval('miniprot --version 2>&1 | head -1'),
          topic: versions, emit: versions_miniprot

    when:
    task.ext.when == null || task.ext.when

    script:
    def code = params.candida_genetic_code ?: 12
    """
    set -euo pipefail
    extract_gene_proteins.py \\
        --sample_id ${meta.id} \\
        --assembly ${assembly} \\
        --canonical_proteins_dir ${canonical_proteins_dir} \\
        --panel ${panel_tsv} \\
        --mutation_catalog ${mutation_catalog} \\
        --genetic_code ${code} \\
        --threads ${task.cpus} \\
        --outdir .
    """

    stub:
    """
    mkdir -p alignments
    printf "sample_id\\tgene\\tref_pos\\tref_aa\\tsample_aa\\talignment_class\\tmutation_id\\tfungamr_match\\tdrugs\\tevidence_strength\\tconfidence_score\\tcompanion_mutations\\tn_reports_fungamr\\tclassification\\n%s\\tERG11\\t132\\tY\\tF\\tmissense\\tY132F\\tyes\\tfluconazole\\tR-strong\\t2.0\\t\\t226\\tknown_resistance\\n" "${meta.id}" > ${meta.id}.gene_proteins.tsv
    printf "sample_id\\tgene\\tref_aa_length\\tsample_aa_length\\tn_mismatches\\tn_known_resistance\\tn_known_sensitivity\\tn_unknown_missense\\talignment_pct_identity\\tminiprot_status\\n%s\\tERG11\\t524\\t524\\t1\\t1\\t0\\t0\\t99.81\\tstub\\n" "${meta.id}" > ${meta.id}.gene_proteins_summary.tsv
    touch alignments/${meta.id}__ERG11.aa_alignment.fasta
    """
}
