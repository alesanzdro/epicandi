// COVERAGE_TRACKS — per-(sample, gene) coverage track gallery.
//
// Cohort-wide visualization: takes all MOSDEPTH regions outputs and renders
// one PNG per (sample, gene) with depth profile and event annotation.
// By default only events != normal are plotted (--only_events) to keep the
// gallery focused on what matters clinically.
//
// Plots are uniform style; depth around the locus with ±padding_bp context.
//
process COVERAGE_TRACKS {
    tag 'cohort'
    label 'process_low'

    conda 'conda-forge::matplotlib conda-forge::pandas conda-forge::numpy'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-12d4ce04c9e7cd1f3fb6e9d5cf0b0b8d1f6c9e8d:abc' :
        'community.wave.seqera.io/library/matplotlib_pandas_numpy:b6e9d5cf0b0b8d1' }"

    input:
    path regions_files     // collected: list of *.regions.bed.gz from all samples
    path coords_dir        // assets/cnv_loci/coords_pre/  (all coords_<tag>.tsv)
    path panel_tsv         // assets/cnv_loci/panel.tsv

    output:
    path "coverage_tracks/", emit: tracks_dir
    path "coverage_tracks/_index.tsv", emit: index, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    // Multi-reference: coverage_tracks.py expects a single --coords TSV, so we
    // merge all coords_<tag>.tsv from coords_dir into one cohort-wide file.
    // The script matches by (chrom, start, end) so the union works.
    def args = task.ext.args ?: '--only_events --padding_bp 15000'
    """
    set -euo pipefail
    mkdir -p coverage_tracks

    # Build a cohort-wide coords TSV: header from the first file, body from all
    head -1 \$(ls ${coords_dir}/coords_*.tsv | head -1) > all_coords.tsv
    for f in ${coords_dir}/coords_*.tsv; do tail -n +2 "\$f" >> all_coords.tsv; done

    coverage_tracks.py \\
        --mosdepth_regions ${regions_files} \\
        --coords all_coords.tsv \\
        --panel ${panel_tsv} \\
        --outdir coverage_tracks/ \\
        ${args}

    rm -f all_coords.tsv
    """

    stub:
    """
    mkdir -p coverage_tracks
    touch coverage_tracks/stub__ERG11.coverage_track.png
    printf "sample_id\\tgene\\tpath\\nstub\\tERG11\\tstub__ERG11.coverage_track.png\\n" > coverage_tracks/_index.tsv
    """
}
