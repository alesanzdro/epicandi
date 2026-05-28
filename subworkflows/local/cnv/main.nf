// CNV — per-sample CNV detection (mosdepth-targeted + CNVkit WGS) +
//       cohort aggregation with mosdepth↔CNVkit cross-validation and
//       optional AMR↔CNV modulation via curated rules.
//
// Inputs:
//   ch_dedup_bam     [meta, bam, bai]  from SNP_CALLING (GATK4_MARKDUPLICATES).
//                                       Reused — we DO NOT realign for CNV.
//   refs_manifest    path TSV          assets/references_manifest.tsv
//   panel_tsv        path              assets/cnv_loci/panel.tsv
//   coords_dir       path              assets/cnv_loci/coords_pre/
//   cnv_amr_rules    path              assets/cnv_loci/cnv_amr_rules.tsv (or NO_FILE)
//   amr_call_matrix  path              from BUILD_MASTER_TABLE.out.call_matrix (or NO_FILE)
//
// Per-sample dataflow:
//   1. join (meta, bam, bai) ↔ refs_manifest by meta.reference_tag
//      → enrich with fasta + mask
//   2. resolve coords_pre/coords_<tag>.tsv; skip samples without pre-computed
//      coords (log.warn + filter) — no runtime fallback in v2.0.
//   3. MOSDEPTH on dedup_bam + panel.bed derived from coords
//   4. CNV_DETECT on mosdepth output → events + contig_depth
//   5. CNVKIT_REFERENCE per unique (ref_tag, fasta, mask)  [cached]
//   6. CNVKIT_BATCH on (sample × ref_tag) for has_illumina samples
//   7. Cohort: COVERAGE_TRACKS + CNV_AGGREGATE  [collected]
//
// Nanopore-only samples: skip CNVKIT_BATCH (CNVkit WGS needs Illumina depth),
// but MOSDEPTH + CNV_DETECT still run for coarse aneuploidy detection (Chr5x2).
//
include { MOSDEPTH               } from '../../../modules/local/mosdepth/main'
include { CNV_DETECT             } from '../../../modules/local/cnv_detect/main'
include { CNVKIT_REFERENCE       } from '../../../modules/local/cnvkit_reference/main'
include { CNVKIT_BATCH           } from '../../../modules/local/cnvkit_batch/main'
include { COVERAGE_TRACKS        } from '../../../modules/local/coverage_tracks/main'
include { CNV_AGGREGATE          } from '../../../modules/local/cnv_aggregate/main'
include { EXTRACT_GENE_PROTEINS  } from '../../../modules/local/extract_gene_proteins/main'
include { CNV_AMR_ATLAS          } from '../../../modules/local/cnv_amr_atlas/main'


workflow CNV {

    take:
    ch_dedup_bam      // [meta, bam, bai]
    ch_assembly       // [meta, polished.fasta]   for EXTRACT_GENE_PROTEINS
    ch_resistance     // [meta, resistance_report.tsv]   from AMR_REPORT
    refs_manifest     // path
    panel_tsv         // path
    coords_dir        // path (dir)
    canonical_proteins_dir  // path (dir)
    mutation_catalog  // path
    cnv_amr_rules     // path (may be NO_FILE)
    amr_call_matrix   // path (may be NO_FILE)

    main:

    // ── 1. Manifest map: tag → {fasta, mask} ──
    ch_manifest = channel
        .fromPath(refs_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            [ row.tag,
              [ fasta: file(row.fasta_path,    checkIfExists: false),
                mask:  file(row.mask_bed_path, checkIfExists: false) ] ]
        }

    // ── 2. Enrich each sample with ref info + coords lookup ──
    def coords_dir_val = coords_dir
    ch_per_sample = ch_dedup_bam
        .map { meta, bam, bai -> [meta.reference_tag, meta, bam, bai] }
        .combine(ch_manifest, by: 0)
        .map { tag, meta, bam, bai, info ->
            def coords_file = file("${coords_dir_val}/coords_${tag}.tsv")
            if (!coords_file.exists()) {
                log.warn "[CNV] No pre-computed coords for sample '${meta.id}' " +
                         "(reference_tag=${tag}); skipping CNV path for this sample."
                return null
            }
            if (!info.fasta.exists()) {
                log.warn "[CNV] Reference fasta missing for tag '${tag}'; " +
                         "skipping sample '${meta.id}'."
                return null
            }
            def mask = info.mask.exists() ? info.mask : file("${projectDir}/assets/NO_FILE")
            [meta, bam, bai, coords_file, tag, info.fasta, mask]
        }
        .filter { it != null }

    // ── 3. MOSDEPTH (all platforms) ──
    MOSDEPTH(
        ch_per_sample.map { meta, bam, bai, coords, _t, _f, _m -> [meta, bam, bai, coords] }
    )

    // ── 4. CNV_DETECT (consumes mosdepth + coords) ──
    ch_cnv_detect_in = MOSDEPTH.out.regions
        .join(MOSDEPTH.out.summary)
        .join( ch_per_sample.map { meta, _b, _bi, coords, _t, _f, _m -> [meta, coords] } )
    CNV_DETECT(ch_cnv_detect_in, panel_tsv)

    // ── 5. CNVKIT_REFERENCE per unique ref_tag (cached) ──
    ch_unique_refs = ch_per_sample
        .map { _meta, _b, _bi, _c, tag, fasta, mask -> [tag, fasta, mask] }
        .unique { tup -> tup[0] }
    CNVKIT_REFERENCE(ch_unique_refs)
    ch_access_by_tag = CNVKIT_REFERENCE.out.access   // [tag, access.bed]

    // ── 6. CNVKIT_BATCH (Illumina/hybrid only) ──
    ch_cnvkit_in = ch_per_sample
        .filter { meta, _b, _bi, _c, _t, _f, _m -> meta.has_illumina == true }
        .map { meta, bam, bai, coords, tag, fasta, _mask ->
            [tag, meta, bam, bai, coords, fasta]
        }
        .combine(ch_access_by_tag, by: 0)
        .map { tag, meta, bam, bai, coords, fasta, access ->
            tuple(meta, bam, bai, tag, fasta, access, coords)
        }
    CNVKIT_BATCH(ch_cnvkit_in)

    // ── 7. Cohort: coverage tracks ──
    ch_all_regions = MOSDEPTH.out.regions.map { _m, regions -> regions }.collect()
    COVERAGE_TRACKS(ch_all_regions, coords_dir, panel_tsv)

    // ── 8. Cohort: aggregate matrices + AMR↔CNV modulation ──
    def no_file = file("${projectDir}/assets/NO_FILE")
    ch_events_all       = CNV_DETECT.out.events.map      { _m, e -> e }.collect()
    ch_cnvkit_panel_all = CNVKIT_BATCH.out.panel_calls
        .map { _m, p -> p }.collect()
        .ifEmpty([no_file])
    ch_contig_depth_all = CNV_DETECT.out.contig_depth.map{ _m, c -> c }.collect()

    CNV_AGGREGATE(
        ch_events_all,
        ch_cnvkit_panel_all,
        ch_contig_depth_all,
        panel_tsv,
        cnv_amr_rules    ?: no_file,
        amr_call_matrix  ?: no_file,
    )

    // ── 9. Per-sample: extract gene proteins from assembly + classify ──
    // Use meta.id as join key because the assembly path's meta is the original
    // (pre-SNP-calling) one without 'ploidy', while ch_per_sample's meta was
    // enriched downstream — `join` on the whole meta object would fail.
    ch_assembly_filtered = ch_assembly
        .map { meta, fa -> [meta.id, meta, fa] }
        .join( ch_per_sample.map { meta, _b, _bi, _c, _t, _f, _m -> [meta.id, true] }, by: 0 )
        .map { _id, meta, fa, _ok -> [meta, fa] }
    EXTRACT_GENE_PROTEINS(
        ch_assembly_filtered,
        canonical_proteins_dir,
        panel_tsv,
        mutation_catalog,
    )

    // ── 10. Per-sample: integrated CNV + AMR atlas (genome + chrom + gene PNGs) ──
    // Same reason: chain joins via meta.id (string), not the full meta map.
    ch_atlas_in = ch_per_sample
        .map { meta, _b, _bi, coords, _t, _f, _m -> [meta.id, meta, coords] }
        .join( CNVKIT_BATCH.out.bins.map           { m, p -> [m.id, p] }, by: 0 )  // cnr
        .join( CNVKIT_BATCH.out.segments_call.map  { m, p -> [m.id, p] }, by: 0 )  // cns
        .join( CNV_DETECT.out.events.map           { m, p -> [m.id, p] }, by: 0 )
        .join( CNV_DETECT.out.contig_depth.map     { m, p -> [m.id, p] }, by: 0 )
        .join( EXTRACT_GENE_PROTEINS.out.long.map  { m, p -> [m.id, p] }, by: 0 )
        .join( EXTRACT_GENE_PROTEINS.out.alignments.map { m, p -> [m.id, p] }, by: 0 )
        .join( MOSDEPTH.out.regions.map            { m, p -> [m.id, p] }, by: 0 )
        .join( ch_resistance.map                   { m, p -> [m.id, p] }, by: 0 )
        .map { _id, meta, coords, cnr, cns, ev, cd, gp, aln, mreg, rr ->
            // Module input tuple order:
            // meta, cnr, cns, events, contig_depth, gene_proteins, aln_dir,
            // mosdepth_regions, resistance, coords
            [meta, cnr, cns, ev, cd, gp, aln, mreg, rr, coords]
        }
    CNV_AMR_ATLAS(ch_atlas_in, panel_tsv, mutation_catalog)

    emit:
    events            = CNV_DETECT.out.events
    contig_depth      = CNV_DETECT.out.contig_depth
    cnvkit_panel      = CNVKIT_BATCH.out.panel_calls
    cnvkit_segments   = CNVKIT_BATCH.out.segments_call
    cnvkit_scatter    = CNVKIT_BATCH.out.scatter_png
    cnvkit_diagram    = CNVKIT_BATCH.out.diagram_pdf
    consensus_matrix  = CNV_AGGREGATE.out.matrix_consensus
    mosdepth_matrix   = CNV_AGGREGATE.out.matrix_mosdepth
    cnvkit_matrix     = CNV_AGGREGATE.out.matrix_cnvkit
    log2_matrix       = CNV_AGGREGATE.out.log2_matrix
    events_summary    = CNV_AGGREGATE.out.events_summary
    contig_combined   = CNV_AGGREGATE.out.contig_depth_combined
    amr_with_cnv      = CNV_AGGREGATE.out.amr_with_cnv
    amr_modulation    = CNV_AGGREGATE.out.amr_modulation_log
    tracks            = COVERAGE_TRACKS.out.tracks_dir
    gene_proteins     = EXTRACT_GENE_PROTEINS.out.long
    gene_proteins_summary = EXTRACT_GENE_PROTEINS.out.summary
    atlas             = CNV_AMR_ATLAS.out.atlas_dir
}
