// ASSEMBLY_HYBRID — Illumina + Nanopore.
//   Flye   ← filtlong-subsampled Nanopore (params.target_coverage × genome_size)
//   Medaka ← same nano reads (×1 pass)
//   Polish ← FULL cleaned Illumina (NO subsample; needed by Polypolish + PyPolca)
//
//   Flye → MEDAKA → POLYPOLISH (full illu) → PYPOLCA (full illu) → ASM_FINALIZE.

include { FLYE          } from '../../../modules/nf-core/flye/main'
include { MEDAKA        } from '../../../modules/nf-core/medaka/main'
include { POLYPOLISH    } from '../../../modules/local/polypolish/main'
include { PYPOLCA       } from '../../../modules/local/pypolca/main'
include { ASM_FINALIZE  } from '../../../modules/local/asm/finalize/main'

workflow ASSEMBLY_HYBRID {

    take:
    ch_clean_illu   // [meta, [R1, R2]]           FULL (no subsample) — for polish
    ch_sub_nano     // [meta, sub_nano.fq.gz]     filtlong subsampled — for Flye/Medaka

    main:

    def ch_flye_in = ch_sub_nano.map { meta, nano -> [meta + [single_end: true], nano] }
    FLYE(ch_flye_in, '--nano-hq')

    def ch_medaka_in = FLYE.out.fasta
        .map { meta, fa -> [meta.id, fa] }
        .join( ch_sub_nano.map { meta, nano -> [meta.id, meta, nano] } )
        .map { _id, fa, meta, nano -> [meta + [single_end: true], nano, fa] }
    MEDAKA(ch_medaka_in)

    def ch_polypolish_in = MEDAKA.out.assembly
        .map { meta, fa -> [meta.id, fa] }
        .join( ch_clean_illu.map { meta, illu -> [meta.id, meta, illu] } )
        .map { _id, fa, meta, illu -> [meta, fa, illu] }
    POLYPOLISH(ch_polypolish_in)

    def ch_pypolca_in = POLYPOLISH.out.fasta
        .map { meta, fa -> [meta.id, meta, fa] }
        .join( ch_clean_illu.map { meta, illu -> [meta.id, illu] } )
        .map { _id, meta, fa, illu -> [meta, fa, illu] }
    PYPOLCA(ch_pypolca_in)

    ASM_FINALIZE(PYPOLCA.out.fasta)

    emit:
    assembly       = ASM_FINALIZE.out.fasta
    finalize       = ASM_FINALIZE.out.stats
    polypolish_log = POLYPOLISH.out.log
    pypolca_log    = PYPOLCA.out.log
}
