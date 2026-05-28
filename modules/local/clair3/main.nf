// CLAIR3: Nanopore variant calling. Mask is applied via --bed_fn pointing to
// the COMPLEMENT of the mask (i.e. accessible regions only), not the mask itself.

process CLAIR3 {
    tag { "${meta.id} (ploidy=${meta.ploidy}, model=${meta.clair3_model})" }
    label 'process_high'

    conda 'bioconda::clair3 bioconda::bedtools'
    container 'community.wave.seqera.io/library/clair3_bedtools:f7f7a3e95c8e7c8a'

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(mask_bed), path(model_dir)

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def haploid_arg = (meta.ploidy as Integer) == 1 ? '--haploid_sensitive' : ''
    def all_ctgs    = params.clair3_include_all_ctgs ? '--include_all_ctgs' : ''
    """
    if [ -s ${mask_bed} ]; then
        cut -f1,2 ${fai} > genome.sizes
        bedtools complement -i ${mask_bed} -g genome.sizes | sort -k1,1 -k2,2n > accessible.bed
        BED_ARG="--bed_fn=accessible.bed"
    else
        BED_ARG=""
    fi

    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${fasta} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=${model_dir} \\
        --output=clair3_out \\
        --gvcf \\
        --min_coverage=${params.clair3_min_coverage} \\
        \${BED_ARG} \\
        ${haploid_arg} ${all_ctgs}

    cp clair3_out/merge_output.gvcf.gz     ${meta.id}.g.vcf.gz
    cp clair3_out/merge_output.gvcf.gz.tbi ${meta.id}.g.vcf.gz.tbi
    """

    stub:
    """
    echo | gzip > ${meta.id}.g.vcf.gz
    touch ${meta.id}.g.vcf.gz.tbi
    """
}
