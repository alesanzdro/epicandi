// FILTLONG_SNP: length + quality filter for Nanopore reads on the SNP-calling path.
// No --target_bases (filter-only). The assembly path will use the same tool with
// `--target_bases` via ext.args2 in a future iteration (Fase F).

process FILTLONG_SNP {
    tag { meta.id }
    label 'process_low'

    conda 'bioconda::filtlong bioconda::nanoq conda-forge::pigz'
    container 'community.wave.seqera.io/library/filtlong_nanoq_pigz:5c0f5717a1aaf80f'

    input:
    tuple val(meta), path(nanopore)

    output:
    tuple val(meta), path("${meta.id}.filt.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.nanoq.json"),    emit: stats, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_len  = params.filtlong_min_length
    def min_q    = params.filtlong_min_mean_q
    def args     = task.ext.args ?: ''
    """
    filtlong --min_length ${min_len} --min_mean_q ${min_q} ${args} ${nanopore} \\
        | pigz -p ${task.cpus} > ${meta.id}.filt.fastq.gz
    nanoq -i ${meta.id}.filt.fastq.gz --json > ${meta.id}.nanoq.json || true
    """

    stub:
    """
    touch ${meta.id}.filt.fastq.gz ${meta.id}.nanoq.json
    """
}
