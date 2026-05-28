// CLEAN: FASTP (Illumina) + PORECHOP_ABI → CHOPPER (Nanopore).
// Input  ch_samplesheet shape: [meta, r1, r2, nano]
// Emits  ch_clean_illu = [meta(single_end:false), [trim_R1, trim_R2]]
//        ch_clean_nano = [meta(single_end:true),  clean.fastq.gz]

include { FASTP        } from '../../../modules/nf-core/fastp/main'
include { PORECHOP_ABI } from '../../../modules/nf-core/porechop/abi/main'
include { CHOPPER      } from '../../../modules/nf-core/chopper/main'

workflow FASTQ_CLEAN {

    take:
    ch_samplesheet

    main:

    ch_illu = ch_samplesheet
        .filter { meta, r1, r2, nano -> meta.has_illumina }
        .map    { meta, r1, r2, nano -> [meta + [single_end: false], [r1, r2], []] }

    ch_nano_raw = ch_samplesheet
        .filter { meta, r1, r2, nano -> meta.has_nanopore }
        .map    { meta, r1, r2, nano -> [meta + [single_end: true], nano] }

    // Illumina branch — FASTP. discard_trimmed_pass=false, save_trimmed_fail=false, save_merged=false
    FASTP(ch_illu, false, false, false)

    // Nanopore branch — PORECHOP_ABI then CHOPPER (no contam fasta passed).
    PORECHOP_ABI(ch_nano_raw, [])
    CHOPPER(PORECHOP_ABI.out.reads, [])

    emit:
    clean_illu   = FASTP.out.reads             // tuple val(meta), path([R1, R2])
    clean_nano   = CHOPPER.out.fastq           // tuple val(meta), path(clean.fastq.gz)
    fastp_json   = FASTP.out.json
    fastp_log    = FASTP.out.log
    porechop_log = PORECHOP_ABI.out.log
}
