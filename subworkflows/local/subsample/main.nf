// SUBSAMPLE_SWF — asymmetric subsampling.
//
//   illumina-only  → seqtk to params.target_coverage (for SPAdes)
//   nanopore-only  → filtlong --target_bases (genome_size × target_coverage)
//   hybrid         → NO illumina subsample (Polypolish + PyPolca need full
//                    coverage); nanopore goes through filtlong like above
//
// Emits per-sample genome size (gsize_bp) for downstream Flye consumption,
// resolved from references_manifest.tsv by meta.reference_tag.

include { SUBSAMPLE as SUBSAMPLE_ILLU } from '../../../modules/local/subsample/main'
include { FILTLONG_SNP as FILTLONG_ASM } from '../../../modules/local/filtlong/main'

workflow SUBSAMPLE_SWF {

    take:
    ch_clean_illu    // [meta, [R1, R2]]            (PASS samples; illumina + hybrid)
    ch_clean_nano    // [meta, clean.fastq.gz]      (PASS samples; nanopore + hybrid)
    ch_species_call  // [meta, species_call.tsv]
    refs_manifest    // path to references_manifest.tsv

    main:

    // ── 1. Build tag → genome_size_bp lookup from the manifest ──
    def ch_gsize_by_tag = channel
        .fromPath(refs_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.tag, row.genome_size_bp as Long] }

    // ── 2. Subsample Illumina ONLY for illumina-only samples ──
    // For hybrid, illumina passes through downstream as the full clean reads —
    // SUBSAMPLE_SWF doesn't emit them here.
    def ch_illu_for_spades = ch_clean_illu
        .filter { meta, _files -> meta.assembly_type == 'illumina' }
        .map { meta, files -> [meta.id, meta, files] }
        .combine(
            ch_species_call.map { meta, tsv -> [meta.id, tsv] }, by: 0
        )
        .map { _id, meta, files, tsv -> [meta, files, tsv] }

    SUBSAMPLE_ILLU(ch_illu_for_spades, file(params.ref_manifest, checkIfExists: true))

    def ch_sub_illu = SUBSAMPLE_ILLU.out.reads.map { meta, files ->
        [meta, files instanceof List ? files : [files]]
    }

    // ── 3. Filtlong on Nanopore (any sample with nano: nanopore-only or hybrid) ──
    // target_bases = genome_size_bp × target_coverage  — dynamic per sample
    def ch_nano_in = ch_clean_nano
        .map { meta, file_ -> [meta.reference_tag, meta, file_] }
        .combine(ch_gsize_by_tag, by: 0)
        .map { _tag, meta, file_, gsize ->
            def target_bases = gsize * (params.target_coverage as Integer)
            [meta + [target_bases: target_bases, genome_size_bp: gsize], file_]
        }

    FILTLONG_ASM(ch_nano_in)

    def ch_sub_nano = FILTLONG_ASM.out.reads.map { meta, file_ ->
        [meta, file_ instanceof List ? file_[0] : file_]
    }

    // ── 4. Emit gsize for ALL samples (Flye -g, etc.) from the manifest ──
    def ch_all_meta = ch_clean_illu.mix(ch_clean_nano)
        .map { meta, _f -> [meta.id, meta] }
        .unique { it[0] }
    def ch_gsize = ch_all_meta
        .map { id, meta -> [meta.reference_tag, id] }
        .combine(ch_gsize_by_tag, by: 0)
        .map { _tag, id, gsize -> [id, gsize] }

    emit:
    sub_illu = ch_sub_illu       // [meta, [R1, R2]]   — illumina-only samples only
    sub_nano = ch_sub_nano       // [meta, fastq.gz]   — all samples with nano
    gsize    = ch_gsize          // [meta.id, gsize_bp_long]
}
