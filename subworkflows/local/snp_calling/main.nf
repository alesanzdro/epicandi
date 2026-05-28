// SNP_CALLING — multi-platform ploidy-aware variant calling layer.
//
//   illumina/hybrid (tier HIGH) → BWA-MEM2 → MARKDUPLICATES → HAPLOTYPECALLER
//   nanopore         (tier MEDIUM) → FILTLONG → MINIMAP2 → CLAIR3
//
// gVCFs are grouped by cohort = (reference_tag, tier) and joint-genotyped per
// cohort with CombineGVCFs + GenotypeGVCFs. Filters are ploidy + tier aware.
// Downstream: vcf2phylip → snp-dists (full.aln) + IQ-TREE (snps_only.aln) +
// UPGMA network coloured by batch.

include { PREPARE_REFERENCE        } from '../../../modules/local/refs/prepare/main'
include { FILTLONG_SNP             } from '../../../modules/local/filtlong/main'
include { MINIMAP2_ONT             } from '../../../modules/local/minimap2_ont/main'
include { CLAIR3                   } from '../../../modules/local/clair3/main'
include { GATK_FILTER_SNPS         } from '../../../modules/local/gatk_filter_snps/main'
include { VCF2PHYLIP               } from '../../../modules/local/vcf2phylip/main'
include { UPGMA_NETWORK            } from '../../../modules/local/upgma/network/main'

include { BWAMEM2_MEM              } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MARKDUPLICATES     } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_HAPLOTYPECALLER    } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_COMBINEGVCFS       } from '../../../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_GENOTYPEGVCFS      } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { SNPDISTS as SNPDISTS_MATRIX } from '../../../modules/nf-core/snpdists/main'
include { SNPDISTS as SNPDISTS_MOLTEN } from '../../../modules/nf-core/snpdists/main'
include { IQTREE                   } from '../../../modules/nf-core/iqtree/main'


// Normalise Dorado model name to Clair3 model dir: v5.2.0 → v520
def clair3_model_from_dorado(dorado_model) {
    if (!dorado_model || dorado_model == 'NA') return params.clair3_default_model
    return dorado_model.replaceAll(/v(\d+)\.(\d+)\.(\d+)/, 'v$1$2$3')
}


workflow SNP_CALLING {

    take:
    ch_clean_illu      // [meta, [r1, r2]]      from FASTQ_CLEAN
    ch_clean_nano      // [meta, fq]            from FASTQ_CLEAN
    refs_manifest      // path to the consolidated TSV
    samplesheet_file   // raw samplesheet path (for UPGMA_NETWORK colouring)

    main:

    // ── 1. Cargar manifest → mapa tag → [fasta, mask, ploidy, clade] ──
    ch_manifest_rows = channel
        .fromPath(refs_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            [
                row.tag,
                [
                    fasta:    file(row.fasta_path,    checkIfExists: false),
                    mask:     file(row.mask_bed_path, checkIfExists: false),
                    ploidy:   row.ploidy as Integer,
                    clade:    row.clade,
                    is_auris: row.is_auris == 'true',
                ]
            ]
        }

    // ── 2. Solo preparamos las refs que aparecen en el samplesheet ──
    ch_used_tags = ch_clean_illu.mix(ch_clean_nano)
        .map { meta, _files -> meta.reference_tag }
        .unique()

    ch_to_prepare = ch_used_tags
        .map { tag -> [tag, tag] }   // key = tag, value = tag (used for combine)
        .combine(ch_manifest_rows, by: 0)
        .map { tag, _t, info ->
            if (!info.fasta.exists()) error "Reference fasta not found for tag '${tag}': ${info.fasta}"
            [[id: tag], info.fasta, info.mask.exists() ? info.mask : file("${projectDir}/assets/NO_FILE")]
        }

    PREPARE_REFERENCE(ch_to_prepare)
    // [meta(id=tag), fasta, fai, dict, bwamem2dir, mmi, mask] — keyed by tag
    ch_prepared = PREPARE_REFERENCE.out.prepared
        .map { m, fa, fai, dict, bwa_dir, mmi, mask -> [m.id, [fasta: fa, fai: fai, dict: dict, bwa: bwa_dir, mmi: mmi, mask: mask]] }

    // ── 3. Enriquecer cada muestra con ploidy + ref info ──
    ch_manifest_info = ch_manifest_rows.map { tag, info -> [tag, info] }

    ch_illu_ready = ch_clean_illu
        .map { meta, files -> [meta.reference_tag, meta, files] }
        .combine(ch_manifest_info, by: 0)
        .combine(ch_prepared,      by: 0)
        .map { tag, meta, files, info, prep ->
            def ploidy = meta.ploidy_override ?: info.ploidy ?: 1
            def meta2 = meta + [ploidy: ploidy]
            [meta2, files, prep]
        }

    // Nanopore SNP path runs ONLY on nanopore-only samples. Hybrid samples
    // already get HAPLOTYPECALLER via the Illumina path; routing them through
    // Clair3 too would emit a second gVCF with the same meta.id and collide
    // in GATK4_COMBINEGVCFS (e.g. 'B11103_hybrid.g.vcf.gz' produced twice).
    ch_nano_ready = ch_clean_nano
        .filter { meta, _file -> meta.assembly_type == 'nanopore' }
        .map { meta, file_ -> [meta.reference_tag, meta, file_] }
        .combine(ch_manifest_info, by: 0)
        .combine(ch_prepared,      by: 0)
        .map { tag, meta, file_, info, prep ->
            def ploidy = meta.ploidy_override ?: info.ploidy ?: 1
            def clair3_model = clair3_model_from_dorado(meta.dorado_model)
            def meta2 = meta + [ploidy: ploidy, clair3_model: clair3_model]
            [meta2, file_, prep]
        }

    // ── 4a. Ruta Illumina/Hybrid (tier HIGH) ──
    ch_bwa_in = ch_illu_ready.map { meta, files, prep ->
        [meta, files instanceof List ? files : [files], prep]
    }

    BWAMEM2_MEM(
        ch_bwa_in.map { meta, reads, prep -> [meta, reads] },
        ch_bwa_in.map { meta, _r, prep    -> [meta, prep.bwa] },
        ch_bwa_in.map { meta, _r, prep    -> [meta, prep.fasta] },
        true   // sort_bam
    )

    GATK4_MARKDUPLICATES(
        BWAMEM2_MEM.out.bam,
        ch_illu_ready.map { meta, _f, prep -> prep.fasta },
        ch_illu_ready.map { meta, _f, prep -> prep.fai   }
    )

    ch_hc_in = GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai)
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .combine(
            ch_illu_ready.map { meta, _f, prep -> [meta.id, prep] },
            by: 0
        )
        .map { _id, meta, bam, bai, prep ->
            [
                [meta, bam, bai, prep.mask, []],     // input (intervals=mask, dragstr=none)
                [meta, prep.fasta],
                [meta, prep.fai],
                [meta, prep.dict],
                [[id: 'none'], []],                   // dbsnp
                [[id: 'none'], []]                    // dbsnp_tbi
            ]
        }

    GATK4_HAPLOTYPECALLER(
        ch_hc_in.map { it -> it[0] },
        ch_hc_in.map { it -> it[1] },
        ch_hc_in.map { it -> it[2] },
        ch_hc_in.map { it -> it[3] },
        ch_hc_in.map { it -> it[4] },
        ch_hc_in.map { it -> it[5] }
    )
    ch_illu_gvcf = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi)

    // ── 4b. Ruta Nanopore (tier MEDIUM) ──
    FILTLONG_SNP(ch_nano_ready.map { meta, file_, _p -> [meta, file_] })

    ch_mm2_in = FILTLONG_SNP.out.reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .combine(ch_nano_ready.map { meta, _f, prep -> [meta.id, prep] }, by: 0)
        .map { _id, meta, reads, prep -> [meta, reads, prep.fasta, prep.fai, prep.mmi] }

    MINIMAP2_ONT(ch_mm2_in)

    ch_clair3_in = MINIMAP2_ONT.out.bam
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .combine(ch_nano_ready.map { meta, _f, prep -> [meta.id, prep] }, by: 0)
        .map { _id, meta, bam, bai, prep ->
            def model_dir = file("${params.clair3_models_dir}/${meta.clair3_model}", checkIfExists: false)
            if (!model_dir.exists()) {
                error "Clair3 model dir not found: ${model_dir}. Download from https://github.com/nanoporetech/rerio/tree/master/clair3_models"
            }
            [meta, bam, bai, prep.fasta, prep.fai, prep.mask, model_dir]
        }

    CLAIR3(ch_clair3_in)
    ch_nano_gvcf = CLAIR3.out.gvcf

    // ── 5. Mix + agrupar por cohort_id ──
    ch_all_gvcf = ch_illu_gvcf.mix(ch_nano_gvcf)
        .map { meta, gvcf, tbi ->
            def cohort_id = "${meta.reference_tag}_${meta.tier}"
            [
                [id: cohort_id, reference_tag: meta.reference_tag, tier: meta.tier, ploidy: meta.ploidy],
                gvcf, tbi
            ]
        }

    ch_cohorts = ch_all_gvcf
        .map { cmeta, gvcf, tbi -> [cmeta.id, cmeta, gvcf, tbi] }
        .groupTuple(by: 0)
        .map { cid, cmetas, gvcfs, tbis ->
            def ploidies = cmetas.collect { it.ploidy }.unique()
            def tiers    = cmetas.collect { it.tier    }.unique()
            def refs     = cmetas.collect { it.reference_tag }.unique()
            if (ploidies.size() > 1 || tiers.size() > 1 || refs.size() > 1) {
                error "Cohort ${cid} has inconsistent ploidy/tier/reference: ${cmetas}"
            }
            [
                [id: cid, reference_tag: refs[0], tier: tiers[0], ploidy: ploidies[0], n: gvcfs.size()],
                gvcfs, tbis
            ]
        }
        .filter { cmeta, _g, _t ->
            if (cmeta.n < params.cohort_min_samples) {
                log.warn "Cohort ${cmeta.id}: ${cmeta.n} samples (<${params.cohort_min_samples}), skipping joint genotyping"
                return false
            }
            return true
        }

    // ── 6. CombineGVCFs + GenotypeGVCFs ──
    ch_combine_in = ch_cohorts
        .map { cmeta, gvcfs, tbis -> [cmeta.reference_tag, cmeta, gvcfs, tbis] }
        .combine(ch_prepared, by: 0)
        .map { _tag, cmeta, gvcfs, tbis, prep ->
            [[cmeta, gvcfs, tbis], prep]
        }

    GATK4_COMBINEGVCFS(
        ch_combine_in.map { it -> it[0] },
        ch_combine_in.map { it -> it[1].fasta },
        ch_combine_in.map { it -> it[1].fai   },
        ch_combine_in.map { it -> it[1].dict  }
    )

    ch_gt_in = GATK4_COMBINEGVCFS.out.combined_gvcf
        .join(GATK4_COMBINEGVCFS.out.combined_tbi)
        .map { cmeta, gvcf, tbi -> [cmeta.id, cmeta, gvcf, tbi] }
        .combine(
            ch_cohorts.map { cmeta, _g, _t -> [cmeta.id, cmeta] }, by: 0
        )
        .map { _id, cmeta, gvcf, tbi, _cm2 -> [cmeta.reference_tag, cmeta, gvcf, tbi] }
        .combine(ch_prepared, by: 0)
        .map { _tag, cmeta, gvcf, tbi, prep ->
            [
                [cmeta, gvcf, tbi, [], []],        // input + tbi + intervals + intervals_tbi
                [cmeta, prep.fasta],
                [cmeta, prep.fai],
                [cmeta, prep.dict],
                [[id: 'none'], []],
                [[id: 'none'], []]
            ]
        }

    GATK4_GENOTYPEGVCFS(
        ch_gt_in.map { it -> it[0] },
        ch_gt_in.map { it -> it[1] },
        ch_gt_in.map { it -> it[2] },
        ch_gt_in.map { it -> it[3] },
        ch_gt_in.map { it -> it[4] },
        ch_gt_in.map { it -> it[5] }
    )

    // ── 7. Filtros ploidy + tier aware ──
    ch_filter_in = GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
        .map { cmeta, vcf, tbi -> [cmeta.reference_tag, cmeta, vcf, tbi] }
        .combine(ch_prepared, by: 0)
        .map { _tag, cmeta, vcf, tbi, prep -> [cmeta, vcf, tbi, prep.fasta, prep.fai, prep.dict] }

    GATK_FILTER_SNPS(ch_filter_in)

    // ── 8. Alignment + distance + tree ──
    VCF2PHYLIP(GATK_FILTER_SNPS.out.snps)

    ch_full_aln = VCF2PHYLIP.out.aln.map { cmeta, full, _snps -> [cmeta, full] }
    ch_snps_aln = VCF2PHYLIP.out.aln.map { cmeta, _full, snps -> [cmeta, snps] }

    SNPDISTS_MATRIX(ch_full_aln)
    SNPDISTS_MOLTEN(ch_full_aln)

    ch_iqtree_in = ch_snps_aln.map { cmeta, snps -> [cmeta, [snps], []] }
    IQTREE(
        ch_iqtree_in,
        [],[],[],[],[],[],[],[],[],[],[],[]
    )

    // ── 9. UPGMA network ──
    ch_upgma_in = SNPDISTS_MATRIX.out.tsv
        .join(SNPDISTS_MOLTEN.out.tsv)
        .map { cmeta, matrix, molten -> [cmeta, matrix, molten] }

    UPGMA_NETWORK(ch_upgma_in, samplesheet_file)

    emit:
    cohort_vcf     = GATK_FILTER_SNPS.out.snps
    cohort_dists   = SNPDISTS_MATRIX.out.tsv
    cohort_tree    = IQTREE.out.phylogeny
    cohort_network = UPGMA_NETWORK.out.svg
    // Per-sample dedup BAM (Illumina/hybrid HIGH tier) for downstream CNV reuse.
    dedup_bam      = GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai)
}
