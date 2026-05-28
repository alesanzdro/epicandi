// AuriClass — C. auris clade typing from Illumina reads. Local because no nf-core module exists.
process AURICLASS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::auriclass=0.5.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/auriclass:0.5.4--pyhdfd78af_0' :
        'quay.io/biocontainers/auriclass:0.5.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads, stageAs: 'reads/*')

    output:
    tuple val(meta), path("${meta.id}.auriclass.tsv"), emit: tsv
    tuple val("${task.process}"), val('auriclass'), eval('auriclass --version 2>&1 | sed "s/.*auriclass //"'), topic: versions, emit: versions_auriclass

    when:
    task.ext.when == null || task.ext.when

    script:
    def r1 = reads instanceof List ? reads[0] : reads
    def r2 = reads instanceof List && reads.size() > 1 ? reads[1] : ''
    """
    auriclass \\
        --name ${meta.id} \\
        --output ${meta.id}.auriclass.tsv \\
        ${r1} ${r2}
    """

    stub:
    """
    printf "Sample\\tClade\\n%s\\tNA\\n" "${meta.id}" > ${meta.id}.auriclass.tsv
    """
}
