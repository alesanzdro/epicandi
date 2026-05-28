// CNVKIT_REFERENCE — build the access BED for one reference (cached per
// reference_tag). The flat .cnn reference is auto-generated inside CNVKIT_BATCH
// when --normal is passed without samples, so this module only needs to
// produce the access.bed that excludes masked regions.
//
// Reference fasta and mask.bed come from references_manifest.tsv. Decompresses
// the fasta if needed and passes mask via `cnvkit.py access -x mask.bed`.
//
// Calibrated 2026-05-24 against B11220: access.bed shrinks to a handful of
// long mappable intervals after applying the funannotate-derived mask.
// See envs/nf-cnvkit.yml for the version pin justification.
//
process CNVKIT_REFERENCE {
    tag "${ref_tag}"
    label 'process_low'

    // NOTE: cnvkit=0.9.13 requires biopython>=1.80 which conflicts with the
    // 'defaults' channel under conda's strict channel priority. Run once on
    // the host:   conda config --set channel_priority flexible
    // (See HANDOFF / setup docs.) Inline conda string matches the repo's
    // pattern (polypolish, clair3) and avoids yml-only quirks.
    conda 'bioconda::cnvkit=0.9.13'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.13--pyhdfd78af_0' :
        'quay.io/biocontainers/cnvkit:0.9.13--pyhdfd78af_0' }"

    input:
    tuple val(ref_tag), path(fasta), path(mask_bed)

    output:
    tuple val(ref_tag), path("${ref_tag}.access.bed"), emit: access
    tuple val("${task.process}"), val('cnvkit'),
          eval('cnvkit.py version 2>&1 | awk "{print \\$NF}"'),
          topic: versions, emit: versions_cnvkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def cat_cmd = fasta.toString().endsWith('.gz') ? 'zcat' : 'cat'
    """
    set -euo pipefail
    ${cat_cmd} ${fasta} > ref.fasta

    cnvkit.py access ref.fasta \\
        -x ${mask_bed} \\
        -o ${ref_tag}.access.bed

    rm -f ref.fasta
    """

    stub:
    """
    printf "chr1\\t0\\t1000\\n" > ${ref_tag}.access.bed
    """
}
