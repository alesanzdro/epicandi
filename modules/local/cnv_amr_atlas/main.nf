// CNV_AMR_ATLAS — per-sample clinical atlas combining CNV (mosdepth +
// CNVkit), assembly-derived AMR mutations (extract_gene_proteins),
// ChroQueTaS AMR calls and FungAMR-known mutation evidence.
//
// One process per sample produces an atlas/ directory:
//   atlas/<sample>_atlas_genome.{png,svg}       always
//   atlas/<sample>_atlas_<chrom>.{png,svg}      one per contig that hosts
//                                                a panel gene with known
//                                                FungAMR R/S evidence
//   atlas/<sample>_atlas_<gene>.{png,svg}       one per panel gene that
//                                                carries known FungAMR R/S
//                                                evidence (default) — can
//                                                be extended to include
//                                                CNV-event genes via
//                                                params.atlas_emit_cnv_genes.
//
// Selection logic intentionally narrow so we don't write 21 panels per
// sample when only one or two are clinically meaningful; the full TSVs
// (gene_proteins, cnv_events) still carry the rest.
//
process CNV_AMR_ATLAS {
    tag "${meta.id}"
    label 'process_low'

    conda 'conda-forge::matplotlib=3.8 conda-forge::numpy'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-matplotlib_numpy:abc' :
        'community.wave.seqera.io/library/matplotlib_numpy:abc' }"

    input:
    tuple val(meta),
          path(cnr),
          path(cns),
          path(cnv_events),
          path(contig_depth),
          path(gene_proteins),
          path(aln_dir),
          path(mosdepth_regions),
          path(resistance_report),
          path(coords)
    path  panel_tsv
    path  mutation_catalog

    output:
    tuple val(meta), path("atlas/"), emit: atlas_dir
    tuple val(meta), path("atlas/${meta.id}_atlas_genome.png"), emit: genome_png

    when:
    task.ext.when == null || task.ext.when

    script:
    def species = meta.species_clinical ?: 'C. auris'
    def qc      = meta.qc_flag           ?: 'PASS'
    def prefix  = "${meta.id}_atlas"
    """
    set -euo pipefail
    mkdir -p atlas

    # 1) Genome-wide atlas — always
    cnv_amr_atlas.py \\
        --mode genome \\
        --sample_id ${meta.id} \\
        --species "${species}" \\
        --qc_flag "${qc}" \\
        --cnr ${cnr} \\
        --cns ${cns} \\
        --cnv_events ${cnv_events} \\
        --contig_depth ${contig_depth} \\
        --gene_proteins ${gene_proteins} \\
        --resistance ${resistance_report} \\
        --coords ${coords} \\
        --panel ${panel_tsv} \\
        --mutation_catalog ${mutation_catalog} \\
        --output atlas/${prefix}_genome

    # 2) Selection: genes with FungAMR known_resistance OR known_sensitivity
    GENES=\$(awk -F'\\t' 'NR>1 && (\$14=="known_resistance" || \$14=="known_sensitivity") {print \$2}' \\
             ${gene_proteins} | sort -u)
    echo "[atlas] genes with known evidence: \${GENES:-(none)}"

    # 3) Per-chromosome panels for contigs hosting any such gene
    if [[ -n "\${GENES:-}" ]]; then
        CHROMS=\$(for g in \$GENES; do
            awk -v g="\$g" -F'\\t' 'NR>1 && \$1==g {print \$2}' ${coords}
        done | sort -u)
        for chrom in \$CHROMS; do
            cnv_amr_atlas.py \\
                --mode chrom --chrom "\$chrom" \\
                --sample_id ${meta.id} \\
                --species "${species}" \\
                --qc_flag "${qc}" \\
                --cnr ${cnr} --cns ${cns} \\
                --cnv_events ${cnv_events} \\
                --contig_depth ${contig_depth} \\
                --gene_proteins ${gene_proteins} \\
                --resistance ${resistance_report} \\
                --coords ${coords} \\
                --panel ${panel_tsv} \\
                --mutation_catalog ${mutation_catalog} \\
                --output atlas/${prefix}
        done
    fi

    # 4) Per-gene panels with protein alignment
    if [[ -n "\${GENES:-}" ]]; then
        for gene in \$GENES; do
            cnv_amr_atlas.py \\
                --mode gene --gene "\$gene" \\
                --sample_id ${meta.id} \\
                --cnr ${cnr} --cns ${cns} \\
                --gene_proteins ${gene_proteins} \\
                --aln_dir ${aln_dir} \\
                --coords ${coords} \\
                --panel ${panel_tsv} \\
                --mosdepth_regions_bed ${mosdepth_regions} \\
                --padding_bp 20000 \\
                --output atlas/${prefix}
        done
    fi

    ls -la atlas/
    """

    stub:
    """
    mkdir -p atlas
    touch atlas/${meta.id}_atlas_genome.png
    touch atlas/${meta.id}_atlas_genome.svg
    """
}
