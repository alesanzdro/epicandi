// CNV_DETECT — Round B.2: per-sample CNV event calling from mosdepth output.
//
// Wraps bin/cnv_detect.py. Five-category classifier on normalised depth
// (amplification/dup_high/dup_low/normal/deletion/deletion_total) plus
// scan-by-contig for aneuploidy detection. Emits an extra Chr5_full row
// when the contig hosting ERG11 reaches ratio >= 1.7 vs genome mean
// (Chr5x2, Li 2024).
//
// Inputs:
//   - regions.bed.gz + summary.txt   from MOSDEPTH
//   - coords_<reference_tag>.tsv     from assets/cnv_loci/coords_pre/
//   - panel.tsv                      from assets/cnv_loci/
//
// Outputs (consumed by AGGREGATE_CNV):
//   - <sample>.cnv_events.tsv   one row per panel gene (+ Chr5_full if any)
//   - <sample>.contig_depth.tsv one row per contig
//
process CNV_DETECT {
    tag "${meta.id}"
    label 'process_single'

    conda 'conda-forge::python=3.12'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(regions_bed_gz), path(summary_txt), path(coords_tsv)
    path  panel_tsv

    output:
    tuple val(meta), path("${meta.id}.cnv_events.tsv"),   emit: events
    tuple val(meta), path("${meta.id}.contig_depth.tsv"), emit: contig_depth
    tuple val("${task.process}"), val('python'),
          eval('python3 --version 2>&1 | awk "{print \\$2}"'),
          topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    cnv_detect.py \\
        --regions_bed_gz ${regions_bed_gz} \\
        --summary        ${summary_txt} \\
        --coords         ${coords_tsv} \\
        --panel          ${panel_tsv} \\
        --sample_id      ${meta.id} \\
        --output_prefix  ${meta.id}
    """

    stub:
    """
    printf "sample_id\\tgene\\ttier\\tchrom\\tstart\\tend\\tmean_depth\\tgenome_median\\tnorm\\tlog2_ratio\\tevent\\tconfidence\\texpected_event\\tagreement_with_expected\\tnotes\\n%s\\tERG11\\t1\\tchr1\\t100\\t200\\t30.0\\t30.0\\t1.0\\t0.0\\tnormal\\thigh\\tduplication\\tNA\\tstub\\n" "${meta.id}" > ${meta.id}.cnv_events.tsv
    printf "sample_id\\tcontig\\tlength_bp\\tmean_depth\\tratio_vs_genome\\tlog2_ratio\\tevent\\n%s\\tchr1\\t1000\\t30.0\\t1.0\\t0.0\\tnormal\\n" "${meta.id}" > ${meta.id}.contig_depth.tsv
    """
}
