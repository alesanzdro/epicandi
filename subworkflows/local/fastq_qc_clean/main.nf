// QC the cleaned reads (post-fastp / post-chopper).
include { FASTQC   } from '../../../modules/nf-core/fastqc/main'
include { NANOPLOT } from '../../../modules/nf-core/nanoplot/main'

workflow FASTQ_QC_CLEAN {

    take:
    ch_clean_illu  // tuple val(meta), path([R1,R2])
    ch_clean_nano  // tuple val(meta), path(clean.fq.gz)

    main:
    FASTQC(ch_clean_illu)
    NANOPLOT(ch_clean_nano)

    emit:
    fastqc_zip   = FASTQC.out.zip
    nanoplot_txt = NANOPLOT.out.txt
}
