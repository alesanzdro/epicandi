// CNV_AGGREGATE — cohort-wide CNV matrices and AMR modulation.
//
// Consumes per-sample artefacts from CNV_DETECT and CNVKIT_BATCH:
//   - *.cnv_events.tsv     (mosdepth events, 1 per sample)
//   - *.cnvkit_panel.tsv   (cnvkit events, 1 per sample)
//   - *.contig_depth.tsv   (per-contig depth, 1 per sample)
//
// Plus optionally the ChroQueTas AMR call_matrix (cohort-wide) and the
// curated cnv_amr_rules.tsv (16 rules with PMID evidence) to emit a
// call_matrix_with_cnv.tsv where CDR1 deletion → S, ERG11 amp → reinforce R,
// Chr5x2 → reinforce R, FKS1 del → R, etc. Trace of every rule that fired
// is preserved in cnv_amr_modulation_log.tsv.
//
// Outputs (in publishDir aggregated/):
//   - cnv_call_matrix_mosdepth.tsv
//   - cnv_call_matrix_cnvkit.tsv
//   - cnv_call_matrix_consensus.tsv   ← KEY (LOW_CONFIDENCE_<m>_<c> on disagreement)
//   - cnv_log2_matrix.tsv
//   - cnv_events_summary.tsv          (long format, non-normal only)
//   - cnv_contig_depth_combined.tsv   (all samples × all contigs)
//   - call_matrix_with_cnv.tsv        (only if --call_matrix provided)
//   - cnv_amr_modulation_log.tsv      (only if --call_matrix provided)
//
process CNV_AGGREGATE {
    tag 'cohort'
    label 'process_single'

    conda 'conda-forge::python=3.12'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    path events_files            // collected: list of *.cnv_events.tsv
    path cnvkit_files            // collected: list of *.cnvkit_panel.tsv (may be empty for nanopore-only cohorts)
    path contig_depth_files      // collected: list of *.contig_depth.tsv
    path panel_tsv
    path cnv_amr_rules           // optional rules file (use NO_FILE_CNV_RULES placeholder if absent)
    path call_matrix             // optional AMR call_matrix.tsv  (use NO_FILE_AMR_MATRIX placeholder)

    output:
    path "aggregated/cnv_call_matrix_mosdepth.tsv",       emit: matrix_mosdepth
    path "aggregated/cnv_call_matrix_cnvkit.tsv",         emit: matrix_cnvkit
    path "aggregated/cnv_call_matrix_consensus.tsv",      emit: matrix_consensus
    path "aggregated/cnv_log2_matrix.tsv",                emit: log2_matrix
    path "aggregated/cnv_events_summary.tsv",             emit: events_summary
    path "aggregated/cnv_contig_depth_combined.tsv",      emit: contig_depth_combined
    path "aggregated/call_matrix_with_cnv.tsv",           emit: amr_with_cnv,        optional: true
    path "aggregated/cnv_amr_modulation_log.tsv",         emit: amr_modulation_log,  optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def have_cnvkit  = cnvkit_files instanceof List ? !cnvkit_files.isEmpty() : (cnvkit_files != null && !cnvkit_files.name.startsWith('NO_FILE'))
    def have_rules   = cnv_amr_rules != null && !cnv_amr_rules.name.startsWith('NO_FILE')
    def have_amr     = call_matrix   != null && !call_matrix.name.startsWith('NO_FILE')
    def cnvkit_arg   = have_cnvkit ? "--cnvkit_panel ${cnvkit_files}" : ""
    def rules_arg    = (have_rules && have_amr) ? "--cnv_amr_rules ${cnv_amr_rules}" : ""
    def amr_arg      = have_amr ? "--call_matrix ${call_matrix}" : ""
    """
    set -euo pipefail
    mkdir -p aggregated

    aggregate_cnv.py \\
        --events ${events_files} \\
        ${cnvkit_arg} \\
        --contig_depth ${contig_depth_files} \\
        ${amr_arg} \\
        ${rules_arg} \\
        --panel ${panel_tsv} \\
        --outdir aggregated/
    """

    stub:
    """
    mkdir -p aggregated
    for f in cnv_call_matrix_mosdepth cnv_call_matrix_cnvkit cnv_call_matrix_consensus cnv_log2_matrix cnv_events_summary cnv_contig_depth_combined; do
        printf "sample_id\\tERG11\\nstub\\tnormal\\n" > aggregated/\${f}.tsv
    done
    """
}
