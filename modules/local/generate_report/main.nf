// GENERATE_REPORT — run-level HTML report (EPIMOL branding).
//
// The Python script reads from the same results tree the pipeline writes to
// (qc_flag, species_call, resistance_report, cnv_events, snpdists, master_table).
// We pass the results directory and the input samplesheet so the report can
// honour batch / collection_date / is_external columns when present.

process GENERATE_REPORT {
    tag "cohort"
    label 'process_single'

    conda 'conda-forge::python>=3.10 conda-forge::pandas conda-forge::numpy conda-forge::openpyxl'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-pandas_openpyxl:abc' :
        'community.wave.seqera.io/library/pandas_openpyxl:abc' }"

    input:
    path results_dir
    path samplesheet
    path branding_dir

    output:
    path "${run_name}_epicandi_report.html", emit: html
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    run_name = task.ext.run_name ?: 'epicandi'
    """
    generate_epicandi_report.py \\
        --results-dir ${results_dir} \\
        --run-name ${run_name} \\
        --samplesheet ${samplesheet} \\
        --branding ${branding_dir} \\
        --output ${run_name}_epicandi_report.html
    """

    stub:
    run_name = task.ext.run_name ?: 'epicandi'
    """
    echo '<html><body>stub</body></html>' > ${run_name}_epicandi_report.html
    """
}
