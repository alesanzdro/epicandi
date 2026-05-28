// QC_GATE v2 — post-IDENTIFICATION 5-criteria AND gate.
//
// Inputs:
//   ch_samplesheet   [meta, r1, r2, nano]               (raw, enriched with reference_tag)
//   ch_clean_illu    [meta, [trim_R1, trim_R2]]         (post-FASTQ_CLEAN)
//   ch_clean_nano    [meta, clean.fq.gz]                (post-FASTQ_CLEAN)
//   ch_species_call  [meta, species_call.tsv]           (from IDENTIFICATION)
//   refs_manifest    path to assets/references_manifest.tsv
//
// Emits:
//   samplesheet_pass [meta(+qc_flag), r1, r2, nano]     (FAIL samples removed)
//   flag_tsv         [meta, qc_flag.tsv]

include { QC_GATE } from '../../../modules/local/qc/gate/main'

workflow QC_GATE_SWF {

    take:
    ch_samplesheet
    ch_clean_illu
    ch_clean_nano
    ch_species_call
    refs_manifest

    main:

    def no_r1   = file("${projectDir}/assets/NO_FILE_R1",   checkIfExists: true)
    def no_r2   = file("${projectDir}/assets/NO_FILE_R2",   checkIfExists: true)
    def no_nano = file("${projectDir}/assets/NO_FILE_NANO", checkIfExists: true)

    // Build tag → genome_size_bp map from the manifest, used by QC_GATE for depth estimation.
    ch_genome_size_by_tag = channel
        .fromPath(refs_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.tag, row.genome_size_bp as Integer] }

    // Per-sample join: meta + clean reads + species_call.
    ch_meta_keyed     = ch_samplesheet  .map { meta, r1, r2, nano -> [meta.id, meta] }
    ch_illu_keyed     = ch_clean_illu   .map { meta, files        -> [meta.id, files] }
    ch_nano_keyed     = ch_clean_nano   .map { meta, file_        -> [meta.id, file_] }
    ch_species_keyed  = ch_species_call .map { meta, tsv          -> [meta.id, tsv] }

    ch_gate_input = ch_meta_keyed
        .join(ch_species_keyed)
        .join(ch_illu_keyed, remainder: true)
        .join(ch_nano_keyed, remainder: true)
        .map { id, meta, sp_tsv, illu_files, nano_file ->
            [
                meta.reference_tag,                       // join key for genome_size
                meta,
                sp_tsv,
                illu_files ? illu_files[0] : no_r1,
                illu_files ? illu_files[1] : no_r2,
                nano_file  ?: no_nano
            ]
        }
        .combine(ch_genome_size_by_tag, by: 0)
        .map { _tag, meta, sp_tsv, r1, r2, nano, gsize ->
            [[meta, sp_tsv, r1, r2, nano], gsize]
        }

    QC_GATE(
        ch_gate_input.map { it -> it[0] },
        ch_gate_input.map { it -> it[1] }
    )

    ch_flagged = QC_GATE.out.flag_tsv
        .map { meta, tsv ->
            // qc_flag is the 13th column (index 12) on row 2 of the TSV
            def flag = tsv.readLines()[1].split('\t')[12].trim()
            [meta + [qc_flag: flag], tsv]
        }

    ch_flagged
        .branch { meta, _tsv ->
            pass: meta.qc_flag == 'PASS'
            fail: meta.qc_flag == 'FAIL'
        }
        .set { ch_branched }

    // Re-pair surviving meta with the raw [r1, r2, nano] for downstream subworkflows.
    ch_pass_keyed = ch_branched.pass.map { meta, _tsv -> [meta.id, meta] }

    ch_pass_samplesheet = ch_samplesheet
        .map { meta, r1, r2, nano -> [meta.id, r1, r2, nano] }
        .join(ch_pass_keyed)
        .map { id, r1, r2, nano, enriched_meta -> [enriched_meta, r1, r2, nano] }

    emit:
    samplesheet_pass = ch_pass_samplesheet
    flag_tsv         = QC_GATE.out.flag_tsv
    failed_meta      = ch_branched.fail.map { meta, _tsv -> meta }
}
