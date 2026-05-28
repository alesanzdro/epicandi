// BUILD_MASTER_TABLE — cohort-wide AMR matrix + master sample table.
//
// Inputs come from AMR_REPORT (resistance_report.tsv per sample) and optionally
// from QC_GATE_SWF (qc_flag.tsv per sample). The script pivots ChroQueTas (pos/neg)
// evidence tuples to R/I/S calls per drug (pos<=4 → R, pos>4 → I, NA → S) and
// emits two outputs:
//
//   call_matrix.tsv   sample × drug → {R, I, S}
//                     This is the input that aggregate_cnv.py consumes for
//                     AMR↔CNV modulation (e.g. CDR1 deletion → itraconazole S).
//
//   master_table.tsv  one row per sample with key fields (species, qc_flag,
//                     resistance mutations, resistance drugs affected, plus
//                     all the drug call columns). Extensible — add assembly
//                     stats / contig depth / etc. when those streams come.
//
process BUILD_MASTER_TABLE {
    tag 'cohort'
    label 'process_single'

    conda 'conda-forge::python=3.12'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    path resistance_reports    // collected: list of *.resistance_report.tsv
    path qc_flags              // collected: list of *.qc_flag.tsv (optional, may be NO_FILE)

    output:
    path 'call_matrix.tsv',   emit: call_matrix
    path 'master_table.tsv',  emit: master_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def qc_arg = (qc_flags && !(qc_flags instanceof List && qc_flags.isEmpty()))
        ? "--qc_flags ${qc_flags}" : ""
    """
    set -euo pipefail
    build_master_table.py \\
        --resistance_reports ${resistance_reports} \\
        ${qc_arg} \\
        --outdir .
    """

    stub:
    """
    printf "sample_id\\tfluconazole\\titraconazole\\tamphotericin_b\\techinocandins\\nstub\\tS\\tS\\tS\\tS\\n" > call_matrix.tsv
    printf "sample_id\\tspecies\\tqc_flag\\tqc_notes\\tn_resistance_mutations\\tresistance_mutations\\tresistance_drugs_affected\\tfluconazole\\titraconazole\\tamphotericin_b\\techinocandins\\nstub\\tstub\\tPASS\\tstub\\t0\\t\\t\\tS\\tS\\tS\\tS\\n" > master_table.tsv
    """
}
