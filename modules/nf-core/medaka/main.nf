process MEDAKA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0' :
        'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.fa.gz"), emit: assembly
    tuple val("${task.process}"), val('medaka'), eval('medaka --version 2>&1 | sed "s/medaka //g"'), topic: versions, emit: versions_medaka

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # FLYE (upstream) emits a gzip-compressed fasta. samtools faidx (called by
    # medaka_consensus → mini_align) requires bgzip, not plain gzip. Decompress
    # to plain fasta first.
    if [[ "${assembly}" == *.gz ]]; then
        zcat ${assembly} > input_assembly.fa
        DRAFT=input_assembly.fa
    else
        DRAFT=${assembly}
    fi

    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -d \$DRAFT \\
        -o ./

    mv consensus.fasta ${prefix}.fa
    gzip -n ${prefix}.fa
    rm -f input_assembly.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fa.gz

    """
}
