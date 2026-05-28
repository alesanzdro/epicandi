// QC_GATE v2: post-IDENTIFICATION, AND-logic on 5 criteria.
process QC_GATE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(species_call), path(clean_r1), path(clean_r2), path(clean_nano)
    val   genome_size_bp

    output:
    tuple val(meta), path("${meta.id}.qc_flag.tsv"), emit: flag_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_illu = !clean_r1.name.startsWith('NO_FILE')
    def has_nano = !clean_nano.name.startsWith('NO_FILE')
    def illu_args = has_illu ? "--clean_r1 ${clean_r1} --clean_r2 ${clean_r2}" : ""
    def nano_args = has_nano ? "--clean_nano ${clean_nano}" : ""
    """
    qc_gate.py \\
        --sample_id      ${meta.id} \\
        --seq_platform   ${meta.assembly_type} \\
        --species_call   ${species_call} \\
        --genome_size_bp ${genome_size_bp} \\
        ${illu_args} \\
        ${nano_args} \\
        --qc_min_ani        ${params.qc_min_ani} \\
        --qc_min_sylph_cov  ${params.qc_min_sylph_cov} \\
        --qc_min_align_rate ${params.qc_min_align_rate} \\
        --qc_min_depth      ${params.qc_min_depth} \\
        --qc_max_contam     ${params.qc_max_contam} \\
        --output ${meta.id}.qc_flag.tsv
    """

    stub:
    """
    printf "sample_id\\tseq_platform\\tspecies_assigned\\tclassification\\tsylph_ani\\tsylph_cov_pct\\talign_rate_pct\\test_depth_x\\tcontam_fraction_pct\\tgenome_size_bp\\tclean_bases_illumina\\tclean_bases_nanopore\\tqc_flag\\tqc_notes\\n%s\\t%s\\tstub\\tstub\\t99\\t99\\t99\\t999\\t0\\t12000000\\t0\\t0\\tPASS\\tstub\\n" "${meta.id}" "${meta.assembly_type}" > ${meta.id}.qc_flag.tsv
    """
}
