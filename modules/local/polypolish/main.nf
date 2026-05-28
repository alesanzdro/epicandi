// All-in-one Polypolish: bwa index → bwa mem -a (per mate) → polypolish filter → polypolish polish.
// Mirrors steps/11c_polypolish.sh. One round only (Wick 2022).
process POLYPOLISH {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::polypolish=0.6.1 bioconda::bwa=0.7.18 bioconda::samtools=1.22"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a3e6e2c0f0d8d8b7d2f6e1c0a4e9c2f3a8f5c5d6:abc' :
        'quay.io/biocontainers/polypolish:0.6.1--h0d85af4_0' }"

    input:
    tuple val(meta), path(draft), path(illumina_reads)

    output:
    tuple val(meta), path("${meta.id}.polypolish.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.polypolish.log"),   emit: log
    tuple val("${task.process}"), val('polypolish'), eval('polypolish --version 2>&1 | sed "s/Polypolish v//"'), topic: versions, emit: versions_polypolish
    tuple val("${task.process}"), val('bwa'),        eval('bwa 2>&1 | sed -n "s/^Version: //p"'),                topic: versions, emit: versions_bwa

    when:
    task.ext.when == null || task.ext.when

    script:
    def r1 = illumina_reads instanceof List ? illumina_reads[0] : illumina_reads
    def r2 = illumina_reads instanceof List && illumina_reads.size() > 1 ? illumina_reads[1] : ''
    def cat_cmd = draft.toString().endsWith('.gz') ? 'zcat' : 'cat'
    """
    set -euo pipefail
    ${cat_cmd} ${draft} > draft.fasta
    bwa index draft.fasta
    bwa mem -t ${task.cpus} -a draft.fasta ${r1} > alignments_R1.sam 2>> ${meta.id}.polypolish.log
    bwa mem -t ${task.cpus} -a draft.fasta ${r2} > alignments_R2.sam 2>> ${meta.id}.polypolish.log
    polypolish filter --in1 alignments_R1.sam --in2 alignments_R2.sam \\
        --out1 filtered_R1.sam --out2 filtered_R2.sam 2>> ${meta.id}.polypolish.log
    # TODO: re-enable --careful once the SNP-only re-run experiment is complete.
    polypolish polish draft.fasta filtered_R1.sam filtered_R2.sam \\
        > ${meta.id}.polypolish.fasta 2>> ${meta.id}.polypolish.log
    rm draft.fasta* alignments_*.sam filtered_*.sam
    """

    stub:
    """
    echo ">${meta.id}_polypolish_contig1\\nACGT" > ${meta.id}.polypolish.fasta
    echo stub > ${meta.id}.polypolish.log
    """
}
