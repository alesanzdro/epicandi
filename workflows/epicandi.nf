/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EPICANDI — full workflow (CHECKPOINT D).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC RAW → CLEAN → QC CLEAN → QC GATE → IDENTIFICATION → SUBSAMPLE
        → ASSEMBLY (short / long / hybrid) → POST_ASM_QC → AMR
        → SNP_CALLING (GATK4 + Clair3, cohort-grouped, ploidy+tier aware)
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_epicandi_pipeline'

include { FASTQ_QC_RAW           } from '../subworkflows/local/fastq_qc_raw'
include { FASTQ_CLEAN            } from '../subworkflows/local/fastq_clean'
include { FASTQ_QC_CLEAN         } from '../subworkflows/local/fastq_qc_clean'
include { QC_GATE_SWF            } from '../subworkflows/local/qc_gate'
include { IDENTIFICATION         } from '../subworkflows/local/identification'
include { SUBSAMPLE_SWF          } from '../subworkflows/local/subsample'
include { ASSEMBLY_SHORT         } from '../subworkflows/local/assembly_short'
include { ASSEMBLY_LONG          } from '../subworkflows/local/assembly_long'
include { ASSEMBLY_HYBRID        } from '../subworkflows/local/assembly_hybrid'
include { POST_ASM_QC            } from '../subworkflows/local/post_asm_qc'
include { AMR_REPORT             } from '../subworkflows/local/amr_report'
include { SNP_CALLING            } from '../subworkflows/local/snp_calling'
include { BUILD_MASTER_TABLE     } from '../modules/local/build_master_table/main'
include { PLOT_AMR_HEATMAP       } from '../modules/local/plot_amr_heatmap/main'
include { CNV                    } from '../subworkflows/local/cnv'
include { CNV_VISUALIZATION      } from '../modules/local/cnv_visualization/main'

workflow EPICANDI {

    take:
    ch_samplesheet // [meta, illu_r1, illu_r2, nano] from PIPELINE_INITIALISATION
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // ── §1 QC RAW ──
    FASTQ_QC_RAW(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_RAW.out.fastqc_zip.map { _meta, file_ -> file_ })

    // ── §2 CLEAN ──
    FASTQ_CLEAN(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_CLEAN.out.fastp_json.map { _meta, file_ -> file_ })

    // ── §3 QC CLEAN ──
    FASTQ_QC_CLEAN(FASTQ_CLEAN.out.clean_illu, FASTQ_CLEAN.out.clean_nano)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_CLEAN.out.fastqc_zip.map { _meta, file_ -> file_ })

    // ── §4 IDENTIFICATION (runs on all samples to feed the QC GATE) ──
    def manifest_ch     = channel.value(file(params.ref_manifest,   checkIfExists: true))
    def refs_dir_ch     = channel.value(file(params.refs_dir,       checkIfExists: true))
    def sylph_db_ch     = channel.value(file(params.sylph_db,       checkIfExists: true))
    def amr_panel_ch    = channel.value(file(params.amr_panel,      checkIfExists: true))
    def rosetta_ch      = channel.value(file(params.rosetta,        checkIfExists: true))
    def refs_manifest_f = file(params.refs_manifest, checkIfExists: true)
    def busco_dl_ch     = params.busco_downloads ? channel.value(file(params.busco_downloads, checkIfExists: false)) : channel.value([])

    IDENTIFICATION(
        ch_samplesheet,
        FASTQ_CLEAN.out.clean_illu,
        FASTQ_CLEAN.out.clean_nano,
        sylph_db_ch,
        refs_dir_ch,
        manifest_ch
    )

    // ── §4b Auto-detect reference_tag from species_call when samplesheet leaves it blank ──
    // species_call.tsv emits `reference_slug` (e.g. cladeI_B8441). The manifest now carries
    // a `sylph_slug` column that maps to the canonical tag (e.g. CaurisI_B8441). If
    // meta.reference_tag is empty/NA, we substitute the manifest tag for that slug.
    def ch_slug_to_tag = channel
        .fromPath(refs_manifest_f, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.sylph_slug, row.tag] }
        .collectFile(newLine: true) { slug, tag -> [ "slug_to_tag.tsv", "${slug}\t${tag}" ] }
        .map { f -> f.text.readLines().findAll { it }.collectEntries { def p = it.split('\t'); [(p[0]): p[1]] } }
    ch_samplesheet_resolved = ch_samplesheet
        .map { meta, r1, r2, nano -> [meta.id, meta, r1, r2, nano] }
        .join(IDENTIFICATION.out.species_call.map { meta, tsv -> [meta.id, tsv] })
        .combine(ch_slug_to_tag)
        .map { _id, meta, r1, r2, nano, sp_tsv, slug2tag ->
            def needs_auto = !meta.reference_tag || meta.reference_tag in ['', 'NA']
            def resolved = meta.reference_tag
            if (needs_auto) {
                def slug = sp_tsv.readLines()[1].split('\t')[10]  // reference_slug column
                resolved = slug2tag[slug]
                if (!resolved) {
                    log.warn "[auto-ref] Sample '${meta.id}': species_call slug '${slug}' not in manifest sylph_slug column; QC will FAIL."
                    resolved = ""
                }
                log.info "[auto-ref] Sample '${meta.id}': reference_tag auto-detected → '${resolved}' (slug=${slug})"
            }
            [meta + [reference_tag: resolved], r1, r2, nano]
        }

    // ── §5 QC GATE v2 (5-criteria AND, post-IDENT) ──
    QC_GATE_SWF(
        ch_samplesheet_resolved,
        FASTQ_CLEAN.out.clean_illu,
        FASTQ_CLEAN.out.clean_nano,
        IDENTIFICATION.out.species_call,
        refs_manifest_f
    )

    // Filter clean reads to PASS samples only so all downstream steps skip FAIL.
    def ch_pass_ids = QC_GATE_SWF.out.samplesheet_pass.map { meta, _r1, _r2, _nano -> [meta.id, meta] }
    def ch_clean_illu_pass = FASTQ_CLEAN.out.clean_illu
        .map { meta, files -> [meta.id, files] }
        .join(ch_pass_ids)
        .map { _id, files, meta -> [meta, files] }
    def ch_clean_nano_pass = FASTQ_CLEAN.out.clean_nano
        .map { meta, file_ -> [meta.id, file_] }
        .join(ch_pass_ids)
        .map { _id, file_, meta -> [meta, file_] }

    // ── §6 SUBSAMPLE (PASS samples only) ──
    // Asymmetric: illumina-only → seqtk for SPAdes; nano (any route) → filtlong
    // with target_bases = genome_size × target_coverage; hybrid-illumina → no
    // subsample (polish gets the full file directly from clean reads).
    SUBSAMPLE_SWF(
        ch_clean_illu_pass,
        ch_clean_nano_pass,
        IDENTIFICATION.out.species_call,
        refs_manifest_f
    )

    // ── §7 ASSEMBLY routes ──
    def ch_sub_illu = SUBSAMPLE_SWF.out.sub_illu       // illumina-only samples only
    def ch_sub_nano = SUBSAMPLE_SWF.out.sub_nano       // nanopore + hybrid (filtlong)

    // illumina-only: SPAdes ← subsampled, polish ← FULL
    ASSEMBLY_SHORT(
        ch_sub_illu        .filter { meta, _f -> meta.assembly_type == 'illumina' },
        ch_clean_illu_pass .filter { meta, _f -> meta.assembly_type == 'illumina' }
    )

    // nanopore-only: Flye + Medaka on subsampled nano, no illumina polish
    ASSEMBLY_LONG(
        ch_sub_nano.filter { meta, _f -> meta.assembly_type == 'nanopore' },
        SUBSAMPLE_SWF.out.gsize
    )

    // hybrid: Flye on sub nano, polish on FULL illumina
    ASSEMBLY_HYBRID(
        ch_clean_illu_pass.filter { meta, _f -> meta.assembly_type == 'hybrid' },
        ch_sub_nano       .filter { meta, _f -> meta.assembly_type == 'hybrid' }
    )

    def ch_assembly = ASSEMBLY_SHORT.out.assembly
        .mix(ASSEMBLY_LONG.out.assembly)
        .mix(ASSEMBLY_HYBRID.out.assembly)

    // ── §8 POST_ASM_QC ──
    POST_ASM_QC(
        ch_assembly,
        IDENTIFICATION.out.species_call,
        refs_manifest_f,
        params.busco_lineage,
        busco_dl_ch
    )
    ch_multiqc_files = ch_multiqc_files.mix(POST_ASM_QC.out.quast_tsv.map { _meta, file_ -> file_ })
    ch_multiqc_files = ch_multiqc_files.mix(POST_ASM_QC.out.busco_txt.map { _meta, f -> f })

    // ── §9 AMR ──
    AMR_REPORT(
        ch_assembly,
        IDENTIFICATION.out.species_call,
        file(params.ref_manifest, checkIfExists: true),
        amr_panel_ch,
        rosetta_ch
    )

    // ── §10 SNP_CALLING (multi-platform: GATK4 + Clair3) — PASS samples only ──
    if (!params.skip_snp_calling) {
        SNP_CALLING(
            ch_clean_illu_pass,
            ch_clean_nano_pass,
            refs_manifest_f,
            file(params.input, checkIfExists: true)
        )
    }

    // ── §11 CNV layer — depends on SNP_CALLING's dedup BAM ──
    if (!params.skip_cnv && !params.skip_snp_calling) {
        // Build cohort-wide AMR matrix from per-sample resistance_reports.
        // qc_flags optional; pass placeholder if QC_GATE_SWF has no flag_tsv yet.
        def no_file_qc = file("${projectDir}/assets/NO_FILE")
        def ch_qc_flags = QC_GATE_SWF.out.flag_tsv
            .map { _meta, tsv -> tsv }
            .collect()
            .ifEmpty([no_file_qc])
        def ch_resistance_reports = AMR_REPORT.out.resistance_report
            .map { _meta, tsv -> tsv }
            .collect()

        BUILD_MASTER_TABLE(ch_resistance_reports, ch_qc_flags)

        CNV(
            SNP_CALLING.out.dedup_bam,
            ch_assembly,
            AMR_REPORT.out.resistance_report,
            refs_manifest_f,
            file(params.cnv_panel,       checkIfExists: true),
            file(params.cnv_coords_dir,  checkIfExists: true),
            file(params.cnv_canonical_proteins_dir, checkIfExists: true),
            file(params.cnv_mutation_catalog,       checkIfExists: true),
            file(params.cnv_amr_rules,   checkIfExists: true),
            BUILD_MASTER_TABLE.out.call_matrix
        )

        // §12 Cohort visualisations: clinical AMR heatmap + integrated CNV plots
        PLOT_AMR_HEATMAP(BUILD_MASTER_TABLE.out.call_matrix)

        CNV_VISUALIZATION(
            CNV.out.events.map           { _meta, f -> f }.collect(),
            CNV.out.contig_depth.map     { _meta, f -> f }.collect(),
            BUILD_MASTER_TABLE.out.call_matrix,
            file(params.cnv_panel, checkIfExists: true),
        )
    }

    //
    // Collate and save software versions via topic channels
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'epicandi_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'epicandi'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    versions       = ch_versions
    species_call   = IDENTIFICATION.out.species_call
    qc_flag        = QC_GATE_SWF.out.flag_tsv
    assembly       = ch_assembly
    amr_report     = AMR_REPORT.out.resistance_report
}
