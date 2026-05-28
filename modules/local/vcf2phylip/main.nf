// VCF2PHYLIP: produce full.aln (REF/ALT/N per sample x site) and snps_only.aln
// (snp-sites -c). full.aln feeds snp-dists; snps_only.aln feeds IQ-TREE.

process VCF2PHYLIP {
    tag { meta.id }
    label 'process_low'

    // vcf2phylip is a single-file Python script distributed only via GitHub
    // (edgardomortiz/vcf2phylip — not on PyPI or bioconda). Pin to v2.8 and
    // fetch on demand into the task workdir.
    conda 'bioconda::snp-sites conda-forge::python>=3.8 conda-forge::curl'
    container 'community.wave.seqera.io/library/vcf2phylip_snp-sites:6c8f5717a1aaf80f'

    input:
    tuple val(meta), path(snps_vcf), path(snps_tbi)

    output:
    tuple val(meta), path("${meta.id}.full.aln"), path("${meta.id}.snps_only.aln"), emit: aln

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Fetch vcf2phylip.py v2.8 from GitHub if not already on PATH.
    if ! command -v vcf2phylip.py >/dev/null 2>&1; then
        curl -fsSL https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/v2.8/vcf2phylip.py \\
            -o vcf2phylip.py
        chmod +x vcf2phylip.py
        export PATH="\$PWD:\$PATH"
    fi
    vcf2phylip.py -i ${snps_vcf} --fasta --output-prefix ${meta.id} --min-samples-locus 1
    mv ${meta.id}.min1.fasta ${meta.id}.full.aln
    snp-sites -c -o ${meta.id}.snps_only.aln ${meta.id}.full.aln
    """

    stub:
    """
    touch ${meta.id}.full.aln ${meta.id}.snps_only.aln
    """
}
