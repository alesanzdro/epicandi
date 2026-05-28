// MOSDEPTH — per-sample depth across the CNV panel.
//
// Consumes the dedup BAM produced by SNP_CALLING (GATK4_MARKDUPLICATES) and
// the pre-computed coords_<reference_tag>.tsv from assets/cnv_loci/coords_pre/.
// The panel BED is derived on the fly from coords (no need to commit one BED
// per strain).
//
// Outputs (consumed by CNV_DETECT):
//   - <sample>.regions.bed.gz  (5 cols: chrom start end gene mean_depth)
//   - <sample>.mosdepth.summary.txt  (per-contig length/mean/min/max)
//
// Resource profile: fast (seconds on Illumina 30x). label process_low.
//
process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::mosdepth=0.3.8'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(coords_tsv)

    output:
    tuple val(meta), path("${meta.id}.regions.bed.gz"),       emit: regions
    tuple val(meta), path("${meta.id}.regions.bed.gz.csi"),   emit: regions_csi
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"), emit: summary
    tuple val(meta), path("${meta.id}.mosdepth.global.dist.txt"), emit: global_dist, optional: true
    tuple val("${task.process}"), val('mosdepth'),
          eval('mosdepth --version 2>&1 | awk "{print \\$2}"'),
          topic: versions, emit: versions_mosdepth

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--no-per-base --fast-mode --mapq 30'
    """
    set -euo pipefail
    # Derive 4-column BED (chrom, start, end, gene) from coords TSV.
    awk -v OFS='\\t' 'NR>1 && \$2!="" {print \$2, \$3, \$4, \$1}' \\
        ${coords_tsv} | sort -k1,1 -k2,2n > panel.bed

    mosdepth \\
        --threads ${task.cpus} \\
        --by panel.bed \\
        ${args} \\
        ${meta.id} \\
        ${bam}
    """

    stub:
    """
    printf "chr1\\t100\\t200\\tGENE_X\\t30.0\\n" | gzip > ${meta.id}.regions.bed.gz
    touch ${meta.id}.regions.bed.gz.csi
    printf "chrom\\tlength\\tbases\\tmean\\tmin\\tmax\\nchr1\\t1000\\t30000\\t30.0\\t0\\t200\\ntotal\\t1000\\t30000\\t30.0\\t0\\t200\\n" > ${meta.id}.mosdepth.summary.txt
    touch ${meta.id}.mosdepth.global.dist.txt
    """
}
