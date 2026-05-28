// Tally mapped reads per reference slug from a screen BAM.
// Slug is taken from RNAME up to '__' (header convention from PREPARE_SCREEN_DB).
process SCREEN_TALLY {
    tag "${meta.id}/${platform}"
    label 'process_low'

    conda "bioconda::samtools=1.22"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), val(platform)

    output:
    tuple val(meta), val(platform), path("readscreen.${platform}.tsv"), emit: tally
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | awk \'{print $2}\''), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -F 4 -F 0x100 -F 0x800 ${bam} \\
      | awk -v sid='${meta.id}' 'BEGIN{OFS="\\t"; print "sample","reference","reads_mapped"} \\
            { split(\$3,a,"__"); cnt[a[1]]++ } \\
            END { for (r in cnt) print sid, r, cnt[r] }' \\
      > readscreen.${platform}.tsv
    """

    stub:
    """
    printf "sample\\treference\\treads_mapped\\n%s\\tstub_slug\\t0\\n" "${meta.id}" > readscreen.${platform}.tsv
    """
}
