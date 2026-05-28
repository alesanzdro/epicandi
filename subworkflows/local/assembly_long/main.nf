// ASSEMBLY_LONG — Nanopore-only: FLYE → MEDAKA → ASM_FINALIZE.

include { FLYE          } from '../../../modules/nf-core/flye/main'
include { MEDAKA        } from '../../../modules/nf-core/medaka/main'
include { ASM_FINALIZE  } from '../../../modules/local/asm/finalize/main'

workflow ASSEMBLY_LONG {

    take:
    ch_sub_nano   // [meta, sub_nano.fq.gz]
    ch_gsize      // [meta.id, gsize_bp.txt]

    main:

    // FLYE -g <bp> comes via ext.args; gsize file is also used by MEDAKA optionally.
    def ch_flye_in = ch_sub_nano.map { meta, nano -> [meta + [single_end: true], nano] }
    FLYE(ch_flye_in, '--nano-hq')

    def ch_medaka_in = FLYE.out.fasta
        .map { meta, fa -> [meta.id, meta, fa] }
        .join( ch_sub_nano.map { meta, nano -> [meta.id, nano] } )
        .map { id, meta, fa, nano -> [meta + [single_end: true], nano, fa] }
    MEDAKA(ch_medaka_in)
    ASM_FINALIZE(MEDAKA.out.assembly)

    emit:
    assembly = ASM_FINALIZE.out.fasta
    finalize = ASM_FINALIZE.out.stats
    flye_log = FLYE.out.log
}
