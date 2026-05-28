// CNVKIT_BATCH — per-sample CNVkit WGS pipeline against an auto-generated
// flat reference.
//
// KEY DESIGN: reuse the dedup BAM produced by SNP_CALLING (GATK4_MARKDUPLICATES)
// instead of re-aligning. CNVkit just needs a sorted+indexed BAM. Re-aligning
// would burn 5-10 min per sample on cohort runs for no gain — alignment params
// match exactly what the SNP path already validated.
//
// Calling step uses `cnvkit.py call --method clonal --ploidy 1` (calibrated
// F1.5, 2026-05-24). The threshold method does NOT respect --ploidy 1 for
// haploids and produces cn=2 for normal segments. clonal applies the exact
// formula round(ploidy * 2^log2), which is correct for clonal-pure haploid
// isolates (auris).
//
// cnvkit_extract_panel.py extracts panel-locus calls in the same format as
// cnv_detect.py output, so aggregate_cnv.py can cross-validate cell-to-cell.
//
// Outputs (consumed by AGGREGATE_CNV and exposed as artefacts):
//   - <sample>.call.cns       segmented copy-number calls (haploid)
//   - <sample>.cnr            bin-level log2 ratios
//   - <sample>.cnvkit_panel.tsv   per-panel-gene cn / log2 / event
//   - <sample>-scatter.png    genome-wide log2 scatter
//   - <sample>-diagram.pdf    chromosome ideogram with CNV calls
//
process CNVKIT_BATCH {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::cnvkit=0.9.13 bioconda::samtools'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.13--pyhdfd78af_0' :
        'quay.io/biocontainers/cnvkit:0.9.13--pyhdfd78af_0' }"

    input:
    // Single combined tuple so the channel is keyed by reference_tag upstream
    // and we don't depend on positional channel alignment across samples/refs.
    tuple val(meta), path(bam), path(bai), val(ref_tag), path(fasta), path(access_bed), path(coords_tsv)

    output:
    tuple val(meta), path("${meta.id}.call.cns"),         emit: segments_call
    tuple val(meta), path("${meta.id}.cnr"),              emit: bins
    tuple val(meta), path("${meta.id}.cnvkit_panel.tsv"), emit: panel_calls
    tuple val(meta), path("${meta.id}-scatter.png"),      emit: scatter_png,        optional: true
    tuple val(meta), path("${meta.id}-diagram.pdf"),      emit: diagram_pdf,        optional: true
    tuple val(meta), path("scatter_by_contig"),           emit: scatter_by_contig,  optional: true
    tuple val("${task.process}"), val('cnvkit'),
          eval('cnvkit.py version 2>&1 | awk "{print \\$NF}"'),
          topic: versions, emit: versions_cnvkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def cat_cmd = fasta.toString().endsWith('.gz') ? 'zcat' : 'cat'
    def ploidy  = params.cnvkit_ploidy ?: 1
    """
    set -euo pipefail
    ${cat_cmd} ${fasta} > ref.fasta

    # CNVkit looks for <bam>.bai; the SNP path produces <bam>-stem.bai
    if [[ ! -f ${bam}.bai ]]; then
        ln -sf ${bai} ${bam}.bai
    fi

    # batch without --scatter/--diagram (we'll re-generate them tuned post-call)
    cnvkit.py batch ${bam} \\
        --method wgs \\
        --normal \\
        --fasta ref.fasta \\
        --access ${access_bed} \\
        --output-reference ${ref_tag}.flat.cnn \\
        --output-dir cnvkit_out \\
        --processes ${task.cpus}

    # Re-call with haploid clonal model (overrides batch's diploid threshold call)
    cnvkit.py call \\
        --method clonal \\
        --ploidy ${ploidy} \\
        cnvkit_out/${bam.baseName}.cns \\
        -o ${meta.id}.call.cns

    cp cnvkit_out/${bam.baseName}.cnr ${meta.id}.cnr

    # Tuned genome-wide scatter (y-limits trimmed; haploid baseline = log2 0)
    cnvkit.py scatter ${meta.id}.cnr -s ${meta.id}.call.cns \\
        --y-min -3 --y-max 3 \\
        --title "${meta.id} — genome-wide (ploidy=${ploidy})" \\
        -o ${meta.id}-scatter.png || true

    # Diagram (chromosome ideogram + CNV calls)
    cnvkit.py diagram -s ${meta.id}.call.cns ${meta.id}.cnr \\
        --title "${meta.id} — CNV diagram" \\
        -o ${meta.id}-diagram.pdf || true

    # Per-contig zoom scatters — useful for spotting segmental dups (e.g. Chr5x2)
    mkdir -p scatter_by_contig
    for chrom in \$(awk '{print \$1}' ${meta.id}.cnr | tail -n +2 | sort -u); do
        cnvkit.py scatter ${meta.id}.cnr -s ${meta.id}.call.cns \\
            -c \${chrom} \\
            --y-min -3 --y-max 3 \\
            --title "${meta.id} — \${chrom}" \\
            -o scatter_by_contig/${meta.id}-scatter-\${chrom}.png || true
    done

    # Extract panel-locus calls for cross-validation with mosdepth
    cnvkit_extract_panel.py \\
        --calls ${meta.id}.call.cns \\
        --coords ${coords_tsv} \\
        --sample_id ${meta.id} \\
        --ploidy ${ploidy} \\
        --output ${meta.id}.cnvkit_panel.tsv

    rm -rf cnvkit_out ref.fasta ${ref_tag}.flat.cnn
    """

    stub:
    """
    printf "chromosome\\tstart\\tend\\tgene\\tlog2\\tcn\\tdepth\\tprobes\\tweight\\nchr1\\t0\\t1000\\t-\\t0.0\\t1\\t30\\t10\\t10\\n" > ${meta.id}.call.cns
    printf "chromosome\\tstart\\tend\\tgene\\tlog2\\tdepth\\tweight\\nchr1\\t0\\t1000\\t-\\t0.0\\t30\\t10\\n" > ${meta.id}.cnr
    printf "sample_id\\tgene\\tchrom\\tstart\\tend\\tcnvkit_cn\\tcnvkit_log2\\tcnvkit_event\\tcnvkit_depth\\n%s\\tERG11\\tchr1\\t100\\t200\\t1\\t0.0\\tnormal\\t30.0\\n" "${meta.id}" > ${meta.id}.cnvkit_panel.tsv
    """
}
