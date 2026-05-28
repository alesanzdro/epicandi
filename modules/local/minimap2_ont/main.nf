// MINIMAP2_ONT: ONT-mode alignment with explicit @RG for downstream Clair3.
// We don't reuse modules/nf-core/minimap2/align/ here because we need precise
// ReadGroup control (Clair3 expects PL:ONT).

process MINIMAP2_ONT {
    tag { "${meta.id} → ${meta.reference_tag}" }
    label 'process_medium'

    conda 'bioconda::minimap2 bioconda::samtools'
    container 'community.wave.seqera.io/library/minimap2_samtools:b3a3da06d6cf2f4b'

    input:
    tuple val(meta), path(reads), path(fasta), path(fai), path(mmi)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT\\tLB:${meta.id}\\tPU:${meta.batch ?: 'na'}"
    """
    minimap2 -ax map-ont -t ${task.cpus} -R "${rg}" ${mmi} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index -@ ${task.cpus} ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    """
}
