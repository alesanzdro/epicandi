// ASSEMBLY_SHORT — Illumina-only:
//   SPAdes  ← subsampled reads (params.target_coverage, faster assembly)
//   Polish  ← FULL cleaned reads (Polypolish needs multi-mapping in repeats,
//             PyPolca needs full local depth for its Bayesian model — Wick 2024)
//
//   SPAdes → POLYPOLISH (full) → PYPOLCA (full) → ASM_FINALIZE.

include { SPADES        } from '../../../modules/nf-core/spades/main'
include { POLYPOLISH    } from '../../../modules/local/polypolish/main'
include { PYPOLCA       } from '../../../modules/local/pypolca/main'
include { ASM_FINALIZE  } from '../../../modules/local/asm/finalize/main'

workflow ASSEMBLY_SHORT {

    take:
    ch_sub_illu     // [meta, [sub_R1, sub_R2]]   — feeds SPAdes
    ch_clean_illu   // [meta, [R1, R2]]           — feeds Polypolish + PyPolca (FULL)

    main:

    def ch_spades_in = ch_sub_illu.map { meta, files -> [meta + [single_end: false], files, [], []] }
    SPADES(ch_spades_in, [], [])

    def ch_polypolish_in = SPADES.out.contigs
        .map { meta, fa -> [meta.id, meta, fa] }
        .join( ch_clean_illu.map { meta, files -> [meta.id, files] } )
        .map { _id, meta, fa, illu -> [meta, fa, illu] }
    POLYPOLISH(ch_polypolish_in)

    def ch_pypolca_in = POLYPOLISH.out.fasta
        .map { meta, fa -> [meta.id, meta, fa] }
        .join( ch_clean_illu.map { meta, files -> [meta.id, files] } )
        .map { _id, meta, fa, illu -> [meta, fa, illu] }
    PYPOLCA(ch_pypolca_in)

    ASM_FINALIZE(PYPOLCA.out.fasta)

    emit:
    assembly       = ASM_FINALIZE.out.fasta
    finalize       = ASM_FINALIZE.out.stats
    spades_log     = SPADES.out.log
    polypolish_log = POLYPOLISH.out.log
    pypolca_log    = PYPOLCA.out.log
}
