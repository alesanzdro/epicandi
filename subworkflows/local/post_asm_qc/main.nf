// POST_ASM_QC — QUAST (per-sample reference resolved via manifest) + BUSCO (saccharomycetes_odb12).

include { QUAST       } from '../../../modules/nf-core/quast/main'
include { BUSCO_BUSCO } from '../../../modules/nf-core/busco/busco/main'

workflow POST_ASM_QC {

    take:
    ch_assembly      // [meta, final.fasta]
    ch_species_call  // [meta, species_call.tsv]
    refs_manifest    // path to references_manifest.tsv
    busco_lineage    // value (string)
    busco_downloads  // value channel of path (optional)

    main:

    // Per-sample reference resolved from the manifest by meta.reference_tag,
    // not from a slug-based legacy layout. Falls back gracefully if the
    // sample's tag is missing.
    def ch_manifest_fasta = channel
        .fromPath(refs_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.tag, file(row.fasta_path, checkIfExists: false)] }

    def ch_ref = ch_assembly
        .map { meta, fa -> [meta.reference_tag, meta, fa] }
        .combine(ch_manifest_fasta, by: 0)
        .map { _tag, meta, fa, ref ->
            if (!ref.exists()) {
                error("POST_ASM_QC: reference fasta not found for tag '${meta.reference_tag}': ${ref}")
            }
            [meta, fa, ref]
        }

    def ch_quast_consensus = ch_ref.map { meta, fa, ref -> [meta, [fa]] }
    def ch_quast_ref       = ch_ref.map { meta, fa, ref -> [[id: meta.id], ref] }
    def ch_quast_gff       = channel.value([[:], []])

    QUAST(ch_quast_consensus, ch_quast_ref, ch_quast_gff)

    BUSCO_BUSCO(
        ch_assembly,                  // tuple val(meta), path(fasta)
        'genome',
        busco_lineage,
        busco_downloads,
        [],                           // config_file
        true                          // clean_intermediates
    )

    emit:
    quast_tsv  = QUAST.out.tsv
    busco_txt  = BUSCO_BUSCO.out.short_summaries_txt
    busco_log  = BUSCO_BUSCO.out.log
}
