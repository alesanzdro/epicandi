// CNV_VISUALIZATION — cohort-wide 4 plots:
//   plot_cnv_heatmap.png            depth log2 (samples × genes)
//   plot_amr_cnv_combined.png       AMR × CNV overlay
//   plot_chr5_aneuploidy.png        aneuploidy detection across cohort
//   plot_carolus_signature.png      ERG11+TAC1B+FKS1+CIS2(+PEA2) joint pattern

process CNV_VISUALIZATION {
    tag 'cohort'
    label 'process_single'

    conda 'conda-forge::matplotlib conda-forge::pandas conda-forge::numpy'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-matplotlib_pandas_numpy:abc' :
        'community.wave.seqera.io/library/matplotlib_pandas_numpy:abc' }"

    input:
    path cnv_events_files        // collected: list of *.cnv_events.tsv
    path contig_depth_files      // collected: list of *.contig_depth.tsv
    path call_matrix             // AMR cohort matrix
    path panel_tsv

    output:
    path 'plot_cnv_heatmap.png',         emit: heatmap,   optional: true
    path 'plot_amr_cnv_combined.png',    emit: combined,  optional: true
    path 'plot_chr5_aneuploidy.png',     emit: aneuploidy,optional: true
    path 'plot_carolus_signature.png',   emit: carolus,   optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cnv_visualization.py \\
        --cnv_events ${cnv_events_files} \\
        --contig_depth ${contig_depth_files} \\
        --amr_call ${call_matrix} \\
        --panel ${panel_tsv} \\
        --outdir . || echo "[cnv_visualization] non-fatal: some plots may be missing"
    """

    stub:
    """
    touch plot_cnv_heatmap.png plot_amr_cnv_combined.png plot_chr5_aneuploidy.png plot_carolus_signature.png
    """
}
