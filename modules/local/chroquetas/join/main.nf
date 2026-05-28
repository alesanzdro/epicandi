// Wraps bin/chroquetas_join.py — parses ChroQueTaS summary, joins with EpiCandi AMR panel, emits 5-state report.
process CHROQUETAS_JOIN {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::pandas conda-forge::numpy"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(chroquetas_summary), val(species_clinical)
    path amr_panel
    path rosetta

    output:
    tuple val(meta), path("${meta.id}.resistance_report.tsv"), emit: report
    tuple val("${task.process}"), val('chroquetas_join'), eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_chroquetas_join

    when:
    task.ext.when == null || task.ext.when

    script:
    def summary_arg = !chroquetas_summary.name.startsWith('NO_FILE') ? "--chroquetas ${chroquetas_summary}" : ""
    """
    chroquetas_join.py \\
        --sample_id ${meta.id} \\
        --species "${species_clinical}" \\
        ${summary_arg} \\
        --panel ${amr_panel} \\
        --rosetta ${rosetta} \\
        --output ${meta.id}.resistance_report.tsv \\
        --min_depth 10
    """

    stub:
    """
    printf "sample\\tspecies\\tprotein\\tposition\\taa_reference\\taa_query\\tmutation_id\\n%s\\t%s\\tNA\\tNA\\tNA\\tNA\\tNA\\n" "${meta.id}" "${species_clinical}" > ${meta.id}.resistance_report.tsv
    """
}
