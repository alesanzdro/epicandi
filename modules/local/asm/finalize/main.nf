// Filter contigs >= params.asm_minlen and rename to ${id}_contig{nr}. Mirrors steps/12_asm_finalize.sh.
process ASM_FINALIZE {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.13"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.13.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(draft_fa)

    output:
    tuple val(meta), path("${meta.id}.fasta"),     emit: fasta
    tuple val(meta), path("${meta.id}.finalize.tsv"), emit: stats
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def minlen = params.asm_minlen
    def cat_cmd = draft_fa.toString().endsWith('.gz') ? 'zcat' : 'cat'
    """
    set -euo pipefail
    PRE=\$(${cat_cmd} ${draft_fa} | seqkit stats -T | awk 'NR==2{print \$4}')
    ${cat_cmd} ${draft_fa} \\
        | seqkit seq -m ${minlen} \\
        | seqkit replace -p '^.+\$' -r '${meta.id}_contig{nr}' \\
        > ${meta.id}.fasta
    POST=\$(seqkit stats -T ${meta.id}.fasta | awk 'NR==2{print \$4}')
    printf "sample_id\\tcontigs_pre\\tcontigs_post\\tminlen\\n%s\\t%s\\t%s\\t%d\\n" "${meta.id}" "\$PRE" "\$POST" ${minlen} > ${meta.id}.finalize.tsv
    """

    stub:
    """
    echo ">${meta.id}_contig1\\nACGT" > ${meta.id}.fasta
    printf "sample_id\\tcontigs_pre\\tcontigs_post\\tminlen\\n%s\\t0\\t0\\t${params.asm_minlen}\\n" "${meta.id}" > ${meta.id}.finalize.tsv
    """
}
