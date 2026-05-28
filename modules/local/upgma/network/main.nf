// UPGMA_NETWORK: epidemiological network from snp-dists matrix.
// Renders SVG (with tier visible, colour by batch) + GEXF (Gephi) + pairs TSV.

process UPGMA_NETWORK {
    tag { meta.id }
    label 'process_low'

    conda 'conda-forge::networkx conda-forge::matplotlib conda-forge::pandas conda-forge::scipy'
    container 'community.wave.seqera.io/library/networkx_matplotlib_pandas_scipy:a3f7f7a3e95c8e7c'

    input:
    tuple val(meta), path(matrix_tsv), path(molten_tsv)
    path samplesheet

    output:
    tuple val(meta), path("${meta.id}.upgma.svg"),                       emit: svg
    tuple val(meta), path("${meta.id}.upgma.gexf"),                      emit: gexf
    tuple val(meta), path("${meta.id}.transmission_pairs.tsv"),          emit: pairs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    upgma_network.py \\
        --matrix ${matrix_tsv} \\
        --molten ${molten_tsv} \\
        --samplesheet ${samplesheet} \\
        --cohort ${meta.id} \\
        --tier ${meta.tier} \\
        --threshold ${params.transmission_threshold_snps} \\
        --out-svg ${meta.id}.upgma.svg \\
        --out-gexf ${meta.id}.upgma.gexf \\
        --out-pairs ${meta.id}.transmission_pairs.tsv
    """

    stub:
    """
    touch ${meta.id}.upgma.svg ${meta.id}.upgma.gexf ${meta.id}.transmission_pairs.tsv
    """
}
