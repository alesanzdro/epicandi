// GATK_FILTER_SNPS: ploidy + tier aware variant filtration.
// Replaces the chain of nf-core gatk4/selectvariants → gatk4/variantfiltration
// → bcftools view with a single process so the filter switch logic stays
// in one place (matching docs/gatk4/main.nf:367-434).

process GATK_FILTER_SNPS {
    tag { meta.id }
    label 'process_medium'

    conda 'bioconda::gatk4 bioconda::bcftools'
    container 'community.wave.seqera.io/library/gatk4_bcftools:b8e7c8a3f7f7a3e9'

    input:
    tuple val(meta), path(joint_vcf), path(joint_tbi), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.snps.pass.vcf.gz"), path("${meta.id}.snps.pass.vcf.gz.tbi"), emit: snps
    tuple val(meta), path("${meta.id}.filter.summary.tsv"),                                       emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def tier = meta.tier
    def ploidy = meta.ploidy as Integer
    def snp_filter, indel_filter, gt_filter, max_miss
    if (tier == 'MEDIUM') {
        snp_filter   = "QUAL < ${params.clair3_min_qual}"
        indel_filter = "QUAL < ${params.clair3_min_qual}"
        gt_filter    = "DP < ${params.clair3_min_dp}"
        max_miss     = params.clair3_max_missing_site
    } else if (ploidy == 1) {
        snp_filter   = 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0'
        indel_filter = 'QD < 2.0 || FS > 200.0 || SOR > 10.0'
        gt_filter    = "DP < ${params.gatk_min_dp} || GQ < ${params.gatk_min_gq}"
        max_miss     = params.gatk_max_missing_site
    } else {
        snp_filter   = 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
        indel_filter = 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0'
        gt_filter    = "DP < ${params.gatk_min_dp} || GQ < ${params.gatk_min_gq}"
        max_miss     = params.gatk_max_missing_site
    }
    def xmx = (task.memory.giga - 1) as Integer
    """
    gatk --java-options "-Xmx${xmx}g" SelectVariants \\
        -R ${fasta} -V ${joint_vcf} --select-type-to-include SNP -O snps.raw.vcf.gz
    gatk --java-options "-Xmx${xmx}g" SelectVariants \\
        -R ${fasta} -V ${joint_vcf} --select-type-to-include INDEL -O indels.raw.vcf.gz

    gatk --java-options "-Xmx${xmx}g" VariantFiltration \\
        -R ${fasta} -V snps.raw.vcf.gz \\
        --filter-expression "${snp_filter}" --filter-name "FAIL_SNP" \\
        -O snps.flagged.vcf.gz
    gatk --java-options "-Xmx${xmx}g" VariantFiltration \\
        -R ${fasta} -V indels.raw.vcf.gz \\
        --filter-expression "${indel_filter}" --filter-name "FAIL_INDEL" \\
        -O indels.flagged.vcf.gz

    gatk --java-options "-Xmx${xmx}g" VariantFiltration \\
        -R ${fasta} -V snps.flagged.vcf.gz \\
        --genotype-filter-expression "${gt_filter}" \\
        --genotype-filter-name "lowQual" \\
        --set-filtered-genotype-to-no-call true \\
        -O snps.genofiltered.vcf.gz

    bcftools view -f PASS snps.genofiltered.vcf.gz \\
      | bcftools view -e 'F_MISSING > ${max_miss}' -Oz \\
        -o ${meta.id}.snps.pass.vcf.gz
    bcftools index -t ${meta.id}.snps.pass.vcf.gz

    n_raw=\$(bcftools view -H snps.raw.vcf.gz | wc -l)
    n_pass=\$(bcftools view -H ${meta.id}.snps.pass.vcf.gz | wc -l)
    {
      printf 'cohort\\ttier\\tploidy\\tn_snps_raw\\tn_snps_pass\\n'
      printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${tier}" "${ploidy}" "\$n_raw" "\$n_pass"
    } > ${meta.id}.filter.summary.tsv
    """

    stub:
    """
    echo | gzip > ${meta.id}.snps.pass.vcf.gz
    touch ${meta.id}.snps.pass.vcf.gz.tbi
    printf 'cohort\\ttier\\tploidy\\tn_snps_raw\\tn_snps_pass\\n${meta.id}\\t${meta.tier}\\t${meta.ploidy}\\t0\\t0\\n' > ${meta.id}.filter.summary.tsv
    """
}
