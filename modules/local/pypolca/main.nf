// PYPOLCA — short-read polishing after Polypolish (Wick et al. 2024 order).
// Bayesian polish on uniquely-mapped Illumina reads; complements Polypolish
// which handles repetitive regions via multi-mapping.

process PYPOLCA {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::pypolca=0.3.1'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/pypolca:0.3.1--pyhdfd78af_0' :
        'quay.io/biocontainers/pypolca:0.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(draft), path(illumina_reads)

    output:
    tuple val(meta), path("${meta.id}.pypolca.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.pypolca.log"),   emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def r1 = illumina_reads instanceof List ? illumina_reads[0] : illumina_reads
    def r2 = illumina_reads instanceof List && illumina_reads.size() > 1 ? illumina_reads[1] : ''
    """
    set -euo pipefail
    pypolca run \\
        --assembly ${draft} \\
        --reads1   ${r1} \\
        --reads2   ${r2} \\
        --threads  ${task.cpus} \\
        --output   pypolca_out \\
        --prefix   ${meta.id}
    cp pypolca_out/${meta.id}_corrected.fasta ${meta.id}.pypolca.fasta
    cp pypolca_out/${meta.id}.report          ${meta.id}.pypolca.log
    """

    stub:
    """
    touch ${meta.id}.pypolca.fasta ${meta.id}.pypolca.log
    """
}
