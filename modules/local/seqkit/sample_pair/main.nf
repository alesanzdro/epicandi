// SEQKIT_SAMPLE_PAIR — subsample paired R1/R2 fastq.gz with a single shared
// seed so the same read indices are picked in both files (pairing preserved).
// The upstream nf-core seqkit/sample module only accepts a single fastx input
// and breaks on pairs; this local wrapper is the minimal fix.

process SEQKIT_SAMPLE_PAIR {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.13"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'oras://community.wave.seqera.io/library/seqkit:2.13.0--205358a3675c7775' :
        'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751' }"

    input:
    tuple val(meta), path(reads)        // [R1, R2]

    output:
    tuple val(meta), path("${prefix}_R{1,2}.fastq.gz"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def r1 = reads instanceof List ? reads[0] : reads
    def r2 = reads instanceof List ? reads[1] : reads
    """
    seqkit sample --threads ${task.cpus} ${args} ${r1} -o ${prefix}_R1.fastq.gz
    seqkit sample --threads ${task.cpus} ${args} ${r2} -o ${prefix}_R2.fastq.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}_R1.fastq.gz
    echo '' | gzip > ${prefix}_R2.fastq.gz
    """
}
