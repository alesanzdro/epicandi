// PREPARE_REFERENCE: index + dict + mask in one shot for a given reference tag.
// Equivalent to the PREPARE_REFERENCE block of docs/gatk4/main.nf:93-117.

process PREPARE_REFERENCE {
    tag { meta.id }
    label 'process_low'

    conda 'bioconda::samtools bioconda::gatk4 bioconda::bwa-mem2 bioconda::minimap2'
    container 'community.wave.seqera.io/library/bwa-mem2_gatk4_minimap2_samtools:7ce93feccf2a3a3a'

    input:
    tuple val(meta), path(fasta, stageAs: 'input.fasta'), path(mask_bed, stageAs: 'input.mask.bed')

    output:
    tuple val(meta),
          path("ref.fasta"),
          path("ref.fasta.fai"),
          path("ref.dict"),
          path("bwamem2"),
          path("ref.fasta.mmi"),
          path("mask.bed"),                                                       emit: prepared

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp ${fasta} ref.fasta
    samtools faidx ref.fasta
    gatk --java-options "-Xmx${task.memory.giga - 1}g" CreateSequenceDictionary \\
        -R ref.fasta -O ref.dict --TMP_DIR .

    mkdir -p bwamem2
    bwa-mem2 index -p bwamem2/ref ref.fasta

    minimap2 -d ref.fasta.mmi -x map-ont ref.fasta

    if [ -s ${mask_bed} ]; then
        # Drop zero-length and negative-length BED intervals (start >= stop) so
        # GATK HaplotypeCaller does not abort with "Badly formed unclippedLoc".
        awk 'BEGIN{FS=OFS="\\t"} \$2+0 < \$3+0 {print}' ${mask_bed} > mask.bed
    else
        : > mask.bed
    fi
    """

    stub:
    """
    cp ${fasta} ref.fasta
    touch ref.fasta.fai ref.dict ref.fasta.mmi mask.bed
    mkdir -p bwamem2 && touch bwamem2/ref.0123 bwamem2/ref.ann bwamem2/ref.bwt.2bit.64
    """
}
