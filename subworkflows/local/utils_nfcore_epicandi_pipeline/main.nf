//
// Subworkflow with functionality specific to the epicandi/epicandi pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <conda/docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,           // parameters_schema (defaults to nextflow_schema.json)
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null            // cli_typecast — use plugin default
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from samplesheet (EpiCandi 7-column schema).
    // Shape per row: [ meta, illumina_r1, illumina_r2, nanopore ]
    //   meta tags from schema: id, dorado_model, batch, is_external
    //
    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { row ->
            def (meta, illu1, illu2, nano) = row
            def has_illu = illu1 && illu1 != 'NA' && illu2 && illu2 != 'NA'
            def has_nano = nano && nano != 'NA'
            if (!has_illu && !has_nano) {
                error("Sample '${meta.id}' has neither Illumina nor Nanopore reads.")
            }
            def assembly_type = (has_illu && has_nano) ? 'hybrid' : (has_illu ? 'illumina' : 'nanopore')
            // tier for SNP-calling layer (drives cohort grouping + caller choice):
            //   short_read = Illumina/hybrid → BWA-MEM2 + GATK4 HaplotypeCaller
            //   long_read  = Nanopore-only   → minimap2 + Clair3
            def tier = (assembly_type == 'nanopore') ? 'long_read' : 'short_read'
            // Normalise samplesheet overrides: empty / 'NA' → null so manifest takes precedence
            def ploidy_ovr_raw = meta.ploidy_override
            def ploidy_ovr = (ploidy_ovr_raw != null && ploidy_ovr_raw.toString() != '' && ploidy_ovr_raw.toString() != 'NA')
                ? (ploidy_ovr_raw.toString() as Integer) : null
            def clade_raw = (meta.clade ?: 'NA').toString()
            def clade = (clade_raw == 'NA' || clade_raw == '') ? '-' : clade_raw
            def enriched = meta + [
                assembly_type   : assembly_type,
                tier            : tier,
                single_end      : false,
                has_illumina    : has_illu,
                has_nanopore    : has_nano,
                is_external_bool: (meta.is_external == true || meta.is_external?.toString()?.toLowerCase() in ['true', 't', '1']),
                ploidy_override : ploidy_ovr,
                clade           : clade
            ]
            return [
                enriched,
                has_illu ? file(illu1, checkIfExists: true) : [],
                has_illu ? file(illu2, checkIfExists: true) : [],
                has_nano ? file(nano,  checkIfExists: true) : []
            ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Refer to docs/usage.md for troubleshooting."
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Custom validation of pipeline parameters
//
def validateInputParameters() {
    if (!params.skip_readscreen && !params.screen_bt2_index && !params.screen_mmi) {
        log.warn("Read screening enabled but no precomputed bowtie2/minimap2 indices set. " +
                 "Indices will be built from refs_dir at runtime by PREPARE_SCREEN_DB.")
    }
    if (params.refs_dir == null || params.ref_manifest == null) {
        log.warn("params.refs_dir / params.ref_manifest are unset. ID/SUBSAMPLE/QUAST/CHROQUETAS/SNIPPY will fail unless provided.")
    }
}

//
// Methods description used by MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used: nf-schema, nf-core/tools, FastQC, NanoPlot, fastp, porechop_abi, chopper,",
            "filtlong, bowtie2, minimap2, samtools, seqkit, seqtk, sylph, AuriClass, SPAdes, Flye,",
            "medaka, Polypolish, PyPolca, BWA-MEM2, QUAST, BUSCO, ChroQueTaS, GATK4 (HaplotypeCaller,",
            "MarkDuplicates, CombineGVCFs, GenotypeGVCFs, VariantFiltration), Clair3, bcftools,",
            "vcf2phylip, snp-sites, snp-dists, IQ-TREE2, NetworkX (UPGMA), MultiQC."
        ].join(' ').trim()
    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>",
            "<li>Ewels P et al. (2016) MultiQC. Bioinformatics 32(19):3047–3048. doi:10.1093/bioinformatics/btw354.</li>"
        ].join(' ').trim()
    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>Cite the pipeline version DOI when available.</li>"

    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    def methods_text = mqc_methods_yaml.text
    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)
    return description_html.toString()
}
