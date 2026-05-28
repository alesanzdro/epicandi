// IDENTIFICATION: Sylph (ANI + tax abundance) + READ_SCREEN (per-ref purity) +
// AuriClass (C. auris clade) → SPECIES_CALL (combined classification).
//
// Inputs:
//   ch_samplesheet  : [meta, r1, r2, nano] (PASSED through qc_gate)
//   ch_clean_illu   : [meta(single_end:false), [trim_R1, trim_R2]]
//   ch_clean_nano   : [meta(single_end:true),  clean.fq.gz]
//   sylph_db        : path to epicandi.syldb
//   refs_dir        : per-slug reference directory (only consulted if screen indices need to be built)
//   manifest        : ref_candida_info.tsv

include { SYLPH_PROFILE       } from '../../../modules/nf-core/sylph/profile/main'
include { BOWTIE2_ALIGN       } from '../../../modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN      } from '../../../modules/nf-core/minimap2/align/main'
include { SEQKIT_SAMPLE_PAIR  } from '../../../modules/local/seqkit/sample_pair/main'
include { SEQKIT_SAMPLE as SEQKIT_SAMPLE_NANO } from '../../../modules/nf-core/seqkit/sample/main'
include { PREPARE_SCREEN_DB   } from '../prepare_screen_db'
include { SCREEN_TALLY as TALLY_ILLU } from '../../../modules/local/screen/tally/main'
include { SCREEN_TALLY as TALLY_NANO } from '../../../modules/local/screen/tally/main'
include { AURICLASS                  } from '../../../modules/local/auriclass/main'
include { SPECIES_CALL               } from '../../../modules/local/species/call/main'

workflow IDENTIFICATION {

    take:
    ch_samplesheet
    ch_clean_illu
    ch_clean_nano
    sylph_db          // value channel of file
    refs_dir          // value channel of path
    manifest          // value channel of file

    main:

    def no_rs_illu = file("${projectDir}/assets/NO_FILE_RS_ILLU", checkIfExists: true)
    def no_rs_nano = file("${projectDir}/assets/NO_FILE_RS_NANO", checkIfExists: true)
    def no_auri    = file("${projectDir}/assets/NO_FILE_AURI",    checkIfExists: true)
    def no_file    = file("${projectDir}/assets/NO_FILE",         checkIfExists: true)

    // ── SYLPH: prefer Illumina reads, fall back to Nanopore for nanopore-only samples ──
    ch_sylph_input = ch_samplesheet.flatMap { meta, r1, r2, nano ->
        meta.has_illumina
            ? [[ meta + [single_end: false], [r1, r2] ]]
            : [[ meta + [single_end: true],  [nano]   ]]
    }
    SYLPH_PROFILE(ch_sylph_input, sylph_db)

    // ── READ_SCREEN (gated by params.skip_readscreen) ──
    ch_rs_illu = channel.empty()
    ch_rs_nano = channel.empty()

    if (!params.skip_readscreen) {

        // Resolve / build the screen indices
        if (!params.screen_bt2_index || !params.screen_mmi) {
            PREPARE_SCREEN_DB(refs_dir, channel.value(file(params.refs_manifest, checkIfExists: true)))
            ch_bt2 = PREPARE_SCREEN_DB.out.bt2
            ch_mmi = PREPARE_SCREEN_DB.out.mmi
            ch_screen_fa = PREPARE_SCREEN_DB.out.fasta
        } else {
            ch_bt2 = channel.value([[id: 'preset'], file(params.screen_bt2_index, checkIfExists: true)])
            ch_mmi = channel.value([[id: 'preset'], file(params.screen_mmi,        checkIfExists: true)])
            ch_screen_fa = channel.value([[id: 'preset'], no_file])
        }

        // Subsample then align — keeps screening fast on big libraries.
        ch_seqkit_illu = ch_clean_illu
        ch_seqkit_nano = ch_clean_nano

        SEQKIT_SAMPLE_PAIR(ch_seqkit_illu)
        SEQKIT_SAMPLE_NANO(ch_seqkit_nano)

        BOWTIE2_ALIGN(
            SEQKIT_SAMPLE_PAIR.out.fastx.map { meta, files -> [meta, files instanceof List ? files : [files]] },
            ch_bt2,
            ch_screen_fa,
            false,   // save_unaligned
            true     // sort_bam (we just need a BAM the tally can read)
        )
        TALLY_ILLU(BOWTIE2_ALIGN.out.bam.map { meta, bam -> [meta, bam, 'illumina'] })

        MINIMAP2_ALIGN(
            SEQKIT_SAMPLE_NANO.out.fastx,
            ch_mmi,
            true,    // bam_format
            '',      // bam_index_extension (none)
            false,   // cigar_paf_format
            false    // cigar_bam
        )
        TALLY_NANO(MINIMAP2_ALIGN.out.bam.map { meta, bam -> [meta, bam, 'nanopore'] })

        ch_rs_illu = TALLY_ILLU.out.tally.map { meta, _platform, tsv -> [meta.id, tsv] }
        ch_rs_nano = TALLY_NANO.out.tally.map { meta, _platform, tsv -> [meta.id, tsv] }
    }

    // ── AURICLASS (Illumina only) ──
    ch_auri_input = ch_clean_illu  // [meta, [R1, R2]]
    AURICLASS(ch_auri_input)
    ch_auri_keyed = AURICLASS.out.tsv.map { meta, tsv -> [meta.id, tsv] }

    // ── SPECIES_CALL — join all per-sample inputs by id ──
    // To avoid `.join(remainder:true)` edge cases (Nextflow's remainder pads
    // unmatched RHS entries with nulls in unexpected positions), we pre-fill
    // every per-sample channel with a placeholder file, then mix+groupTuple
    // to take the real file when present and the placeholder otherwise.
    ch_sylph_keyed = SYLPH_PROFILE.out.profile_out.map { meta, tsv -> [meta.id, tsv] }
    ch_meta_keyed  = ch_samplesheet.map { meta, _r1, _r2, _nano -> [meta.id, meta] }

    def pick_real_or_placeholder = { id, files ->
        def real = files.find { f -> !(f.name.startsWith('NO_FILE')) }
        [id, real ?: files[0]]
    }

    ch_rs_illu_full = ch_meta_keyed.map { id, meta -> [id, no_rs_illu] }
        .mix(ch_rs_illu)
        .groupTuple()
        .map(pick_real_or_placeholder)

    ch_rs_nano_full = ch_meta_keyed.map { id, meta -> [id, no_rs_nano] }
        .mix(ch_rs_nano)
        .groupTuple()
        .map(pick_real_or_placeholder)

    ch_auri_full = ch_meta_keyed.map { id, meta -> [id, no_auri] }
        .mix(ch_auri_keyed)
        .groupTuple()
        .map(pick_real_or_placeholder)

    ch_species_input = ch_meta_keyed
        .join(ch_sylph_keyed)
        .join(ch_rs_illu_full)
        .join(ch_rs_nano_full)
        .join(ch_auri_full)
        .map { id, meta, sylph_tsv, rs_illu_tsv, rs_nano_tsv, auri_tsv ->
            [meta, sylph_tsv, rs_illu_tsv, rs_nano_tsv, auri_tsv]
        }

    SPECIES_CALL(ch_species_input, manifest)

    emit:
    sylph_profile   = SYLPH_PROFILE.out.profile_out
    species_call    = SPECIES_CALL.out.tsv
    auriclass       = AURICLASS.out.tsv
}
