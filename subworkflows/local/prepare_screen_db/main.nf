// Build bowtie2 + minimap2 indices for the 27-reference screen panel once per run.
// Skipped entirely if params.screen_bt2_index AND params.screen_mmi are both set.

include { SCREEN_FASTA   } from '../../../modules/local/screen/fasta/main'
include { BOWTIE2_BUILD  } from '../../../modules/nf-core/bowtie2/build/main'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'

workflow PREPARE_SCREEN_DB {

    take:
    refs_dir       // value channel of Path
    refs_manifest  // value channel of Path to references_manifest.tsv

    main:
    // refs_dir arrives as a value channel wrapping a Path. Nextflow 25.10
    // rejects channel.value([..., <channel>]) — wrap with .map instead so we
    // emit a single value-tuple downstream. Combine with refs_manifest so
    // SCREEN_FASTA can resolve tag -> sylph_slug for header prefixes.
    def ch_input = refs_dir
        .combine(refs_manifest)
        .map { dir, mani -> [[id: 'epicandi_screen'], dir, mani] }
    SCREEN_FASTA(ch_input)
    BOWTIE2_BUILD(SCREEN_FASTA.out.fasta)
    MINIMAP2_INDEX(SCREEN_FASTA.out.fasta)

    emit:
    fasta = SCREEN_FASTA.out.fasta                    // tuple val(meta), path(screen.fa.gz)
    bt2   = BOWTIE2_BUILD.out.index                   // tuple val(meta), path(bowtie2/)
    mmi   = MINIMAP2_INDEX.out.index                  // tuple val(meta), path(*.mmi)
}
