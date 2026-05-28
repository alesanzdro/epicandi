// Wraps bin/species_call_v2.py — combines sylph + read_screen + auriclass into a single species call.
process SPECIES_CALL {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::pandas conda-forge::numpy"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(sylph_profile), path(rs_illu), path(rs_nano), path(auriclass_tsv)
    path manifest

    output:
    tuple val(meta), path("${meta.id}.species_call.tsv"), emit: tsv
    tuple val("${task.process}"), val('species_call'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_species_call

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_rs_illu = !rs_illu.name.startsWith('NO_FILE')
    def has_rs_nano = !rs_nano.name.startsWith('NO_FILE')
    def has_auri    = !auriclass_tsv.name.startsWith('NO_FILE')
    def rs_illu_arg = has_rs_illu ? "--readscreen_illumina ${rs_illu}" : ""
    def rs_nano_arg = has_rs_nano ? "--readscreen_nanopore ${rs_nano}" : ""
    def auri_arg    = has_auri    ? "--auriclass ${auriclass_tsv}"      : ""
    """
    species_call_v2.py \\
        --sample_id ${meta.id} \\
        --sylph ${sylph_profile} \\
        ${rs_illu_arg} ${rs_nano_arg} ${auri_arg} \\
        --manifest ${manifest} \\
        --ani_high ${params.ani_high_conf} \\
        --ani_low  ${params.ani_low_conf} \\
        --readscreen_pct_top ${params.readscreen_pct_top} \\
        --readscreen_pct_unmapped_max ${params.readscreen_pct_unmapped_max} \\
        --output ${meta.id}.species_call.tsv
    """

    stub:
    """
    printf "sample_id\\tspecies_assigned\\tspecies_clinical\\tclassification\\tsylph_ani\\tsylph_tax_pct\\treadscreen_top_pct\\treadscreen_other_pct\\treadscreen_unmapped_pct\\treadscreen_total_reads\\treference_slug\\tploidy\\tploidy_conf\\tclade\\n%s\\tCandidozyma auris\\tC. auris\\thigh_conf\\t99.99\\t99.0\\t99.0\\t1.0\\t0.5\\t100000\\tcladeII_B11220\\t1\\talta\\tII\\n" "${meta.id}" > ${meta.id}.species_call.tsv
    """
}
