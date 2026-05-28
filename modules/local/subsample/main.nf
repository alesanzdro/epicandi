// Per-sample subsampling with hysteresis. Mirrors steps/10_subsample.sh:
//   - Resolves genome size (MB) from manifest by reference_slug; falls back to 12.5 MB.
//   - Computes current Gb of cleaned reads via seqkit stats.
//   - If current > 1.1 × target → seqtk sample to FRAC = target/current; else symlink through.
//   - Writes ${meta.id}.gsize_bp (used by FLYE -g downstream).
//
// Handles Illumina (paired or single) AND Nanopore in the same process via meta.single_end.
process SUBSAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.13 bioconda::seqtk=1.5 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c0d80b8eea0a8e2c9a3fc31d61b22e0a87b78c5cf69a4f2acf9b7b9adf2c8a7/data' :
        'community.wave.seqera.io/library/seqkit_seqtk_pigz:b97a3b1bf60d28a4' }"

    input:
    tuple val(meta), path(reads), path(species_call_tsv)
    path manifest

    output:
    tuple val(meta), path("sub_*"),           emit: reads
    tuple val(meta), path("${meta.id}.gsize_bp"), emit: gsize
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def target_cov = params.target_coverage
    def seed       = params.subsample_seed
    """
    set -euo pipefail

    # ── 1. Resolve genome size MB from manifest by slug ──────────────────
    SLUG=\$(awk -F'\\t' 'NR==2{print \$11}' "${species_call_tsv}")
    # genome_size_mb may live in an extension column (manifest currently has 15
    # cols; bash defaults to 12.5 MB when absent). We probe columns 16/17/18/19
    # to remain forward-compatible if the manifest grows.
    GENOME_MB=""
    for col in 16 17 18 19; do
        GENOME_MB=\$(awk -F'\\t' -v s="\$SLUG" -v c=\$col 'NR>1 && \$1==s{print \$c; exit}' "${manifest}")
        case "\$GENOME_MB" in ''|[!0-9.]*) ;; *) break;; esac
    done
    case "\$GENOME_MB" in ''|[!0-9.]*) GENOME_MB="12.5";; esac

    GSIZE_BP=\$(awk -v g="\$GENOME_MB" 'BEGIN{printf "%d", g*1e6}')
    echo "\$GSIZE_BP" > ${meta.id}.gsize_bp

    # ── 2. Current Gb of cleaned reads ───────────────────────────────────
    CUR_GB=\$(seqkit stats -T --threads ${task.cpus} ${reads.join(' ')} | awk 'NR>1{s+=\$5} END{printf "%.6f", s/1e9}')
    TARGET_GB=\$(awk -v c=${target_cov} -v g="\$GENOME_MB" 'BEGIN{printf "%.6f", c*g/1000}')

    # ── 3. Decide ────────────────────────────────────────────────────────
    NEEDS=\$(awk -v c="\$CUR_GB" -v t="\$TARGET_GB" 'BEGIN{print (c > 1.1*t) ? 1 : 0}')
    if [ "\$NEEDS" = "1" ]; then
        FRAC=\$(awk -v t="\$TARGET_GB" -v c="\$CUR_GB" 'BEGIN{printf "%.6f", t/c}')
        for f in ${reads.join(' ')}; do
            base=\$(basename "\$f")
            seqtk sample -s ${seed} "\$f" "\$FRAC" | pigz -p ${task.cpus} -c > "sub_\${base}"
        done
    else
        for f in ${reads.join(' ')}; do
            base=\$(basename "\$f")
            ln -sf "\$f" "sub_\${base}"
        done
    fi
    """

    stub:
    """
    for f in ${reads.join(' ')}; do
        base=\$(basename "\$f")
        echo "" | gzip > "sub_\${base}"
    done
    echo 12500000 > ${meta.id}.gsize_bp
    """
}
