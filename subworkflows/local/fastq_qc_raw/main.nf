// QC raw reads: FASTQC on Illumina, NANOPLOT on Nanopore.
// Input shape: ch_samplesheet = [meta, r1, r2, nano] (paths may be empty list)

include { FASTQC   } from '../../../modules/nf-core/fastqc/main'
include { NANOPLOT } from '../../../modules/nf-core/nanoplot/main'

workflow FASTQ_QC_RAW {

    take:
    ch_samplesheet

    main:

    ch_illu = ch_samplesheet
        .filter { meta, r1, r2, nano -> meta.has_illumina }
        .map    { meta, r1, r2, nano -> [meta + [single_end: false], [r1, r2]] }

    ch_nano = ch_samplesheet
        .filter { meta, r1, r2, nano -> meta.has_nanopore }
        .map    { meta, r1, r2, nano -> [meta + [single_end: true],  nano] }

    FASTQC(ch_illu)
    NANOPLOT(ch_nano)

    emit:
    fastqc_zip   = FASTQC.out.zip               // tuple val(meta), path(*.zip)
    nanoplot_txt = NANOPLOT.out.txt             // tuple val(meta), path(*.txt)
}
