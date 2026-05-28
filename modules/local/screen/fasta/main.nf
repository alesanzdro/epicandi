// Concatenate every reference genome under refs_dir into a single screen.fa.gz
// with contig headers reheadered as `${sylph_slug}__${original_header}` so the
// downstream SCREEN_TALLY emits per-slug counts that line up with the slug
// nomenclature used by Sylph and by SPECIES_CALL's manifest lookup.
//
// The mapping from on-disk directory name (manifest column `tag`) to the
// Sylph nomenclature (`sylph_slug`) is read from references_manifest.tsv.
// If a reference directory has no entry in the manifest, the dirname is used
// as a fallback prefix.
process SCREEN_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.13"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.13.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(refs_dir), path(refs_manifest)

    output:
    tuple val(meta), path('screen.fa.gz'), emit: fasta
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -uo pipefail

    # Build tag -> sylph_slug map from the manifest (header row + tab-delimited).
    awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++)col[\$i]=i; next} \\
                 col["tag"] && col["sylph_slug"] && \$col["tag"]!="" \\
                     {print \$col["tag"] "\\t" \$col["sylph_slug"]}' \\
        ${refs_manifest} > tag_to_slug.tsv

    : > screen.fa
    n_added=0
    for slug_dir in ${refs_dir}/*/; do
        [ -d "\$slug_dir" ] || continue
        tag=\$(basename "\$slug_dir")
        # Resolve sylph_slug for this tag; fall back to tag if not in manifest.
        sylph_slug=\$(awk -F'\\t' -v t="\$tag" '\$1==t {print \$2; exit}' tag_to_slug.tsv)
        [ -z "\$sylph_slug" ] && sylph_slug="\$tag"

        ref=""
        for cand in "\$slug_dir/reference.fasta" "\$slug_dir/reference.fa" "\$slug_dir/genome.fna.gz" "\$slug_dir/genome.fna" "\$slug_dir/genome.fa.gz"; do
            if [ -f "\$cand" ]; then
                ref="\$cand"
                break
            fi
        done
        if [ -z "\$ref" ]; then
            echo "[SCREEN_FASTA] no reference fasta under \$slug_dir - skipping" >&2
            continue
        fi
        seqkit replace -p '^(\\S+)' -r "\${sylph_slug}__\\\${1}" "\$ref" >> screen.fa
        n_added=\$(( n_added + 1 ))
    done
    if [ "\$n_added" -eq 0 ]; then
        echo "[SCREEN_FASTA] ERROR: no reference fastas found under ${refs_dir}/*/" >&2
        exit 1
    fi
    gzip -f screen.fa
    """

    stub:
    """
    echo '>cladeI_B8441__stub' | gzip > screen.fa.gz
    """
}
