// AMR_REPORT — ChroQueTaS scan + chroquetas_join 5-state report.

include { CHROQUETAS       } from '../../../modules/local/chroquetas/run/main'
include { CHROQUETAS_JOIN  } from '../../../modules/local/chroquetas/join/main'

workflow AMR_REPORT {

    take:
    ch_assembly      // [meta, final.fasta]
    ch_species_call  // [meta, species_call.tsv]
    manifest         // path
    amr_panel        // path
    rosetta          // path

    main:

    // Per-sample manifest lookup: pipeline_id (chroquetas species arg) + has_chroquetas + species_clinical.
    def ch_per_sample = ch_assembly
        .map { meta, fa -> [meta.id, meta, fa] }
        .join( ch_species_call.map { meta, tsv -> [meta.id, tsv] } )
        .map { id, meta, fa, tsv ->
            def species_cols  = tsv.readLines()[1].split('\t')
            def species_clin  = species_cols[2]
            def slug          = species_cols[10]
            def row = manifest.toFile().readLines().drop(1).find { line ->
                def m = line.split('\t')
                m.size() >= 5 && m[0] == slug
            }
            // Manifest column layout:
            //   1 slug · 2 pipeline_id (carries clade suffix for C. auris) ·
            //   3 chroquetas_species (clade-stripped, what ChroQueTaS expects) ·
            //   5 has_chroquetas (true/false)
            def species_arg = 'NA'
            if (row) {
                def m = row.split('\t')
                species_arg = (m[4] == 'true' && m.size() >= 3 && m[2]) ? m[2] : 'NA'
            }
            [meta, fa, species_arg, species_clin]
        }

    def ch_chroquetas_in = ch_per_sample.map { meta, fa, sp_arg, sp_clin -> [meta, fa, sp_arg] }
    CHROQUETAS(ch_chroquetas_in)

    def no_summary = file("${projectDir}/assets/NO_FILE_CHROQUETAS", checkIfExists: true)
    // Pre-fill summary with placeholder per sample, then mix+groupTuple to take
    // the real ChroQueTaS summary when present (same trick as IDENTIFICATION).
    def ch_summary_keyed = ch_per_sample.map { meta, fa, sp_arg, sp_clin -> [meta.id, no_summary] }
        .mix( CHROQUETAS.out.summary.map { meta, summary -> [meta.id, summary] } )
        .groupTuple()
        .map { id, files ->
            def real = files.find { f -> !f.name.startsWith('NO_FILE') }
            [id, real ?: files[0]]
        }

    def ch_join_in = ch_per_sample
        .map { meta, fa, sp_arg, sp_clin -> [meta.id, meta, sp_clin] }
        .join(ch_summary_keyed)
        .map { id, meta, sp_clin, summary ->
            [meta, summary, sp_clin]
        }
    CHROQUETAS_JOIN(ch_join_in, amr_panel, rosetta)

    emit:
    chroquetas_summary = CHROQUETAS.out.summary
    resistance_report  = CHROQUETAS_JOIN.out.report
}
