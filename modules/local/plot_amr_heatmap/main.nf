// PLOT_AMR_HEATMAP — clinical heatmap samples × drugs (R/I/S coloured).
// Wraps bin/plot_amr_heatmap.py on the cohort-wide call_matrix.tsv.

process PLOT_AMR_HEATMAP {
    tag 'cohort'
    label 'process_single'

    conda 'conda-forge::matplotlib conda-forge::numpy'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-matplotlib_numpy:abc' :
        'community.wave.seqera.io/library/matplotlib_numpy:abc' }"

    input:
    path call_matrix

    output:
    path 'amr_heatmap.png', emit: png
    path 'amr_heatmap.svg', emit: svg

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_amr_heatmap.py \\
        --call_matrix ${call_matrix} \\
        --output_prefix amr_heatmap
    """

    stub:
    """
    touch amr_heatmap.png amr_heatmap.svg
    """
}
