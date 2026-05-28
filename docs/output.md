# epicandi/epicandi: Output

## Introduction

Per-sample outputs are published under `<outdir>/<sample_id>/`. Cohort-level outputs land in `<outdir>/cohort/`. MultiQC and pipeline-info reports live at the top of `<outdir>/`.

> The skeleton iteration does not yet wire processes — the directory layout below is the contract that CHECKPOINT C/D fulfils.

## Per-sample outputs

```
<outdir>/<sample_id>/
├── qc/
│   ├── fastqc_raw/             *.html + *.zip
│   ├── fastp/                  trim_R1/R2.fq.gz, fastp.json, fastp.html
│   ├── nanoplot_raw/           NanoStats.txt + plots
│   ├── nanoplot_clean/         NanoStats.txt
│   ├── fastqc_clean/
│   └── qc_flag.tsv             PASS/WARN/FAIL + counts
├── decont/                     decont_R1/R2.fq.gz (only when DECONT is on)
├── id/
│   ├── sylph.profile.tsv
│   ├── readscreen.illumina.tsv
│   ├── readscreen.nanopore.tsv
│   ├── auriclass.tsv
│   └── species_call.tsv        ← THE summary line per sample
├── subsample/                  sub_R1/R2.fq.gz, sub_nano.fq.gz, .gsize_bp
├── assembly/
│   ├── <sample>.fasta          final assembly (post-FINALIZE)
│   ├── spades.log              (illumina)
│   ├── flye.log + assembly_info.txt   (nanopore/hybrid)
│   └── polypolish.log
├── post_asm_qc/
│   ├── quast/                  report.tsv, NGA50, misassemblies
│   └── busco/                  short_summary.specific.*.txt
├── amr/
│   ├── chroquetas/<sample>.ChroQueTaS.AMR_summary.txt
│   └── resistance_report.tsv   ← 5-state matrix
└── (BAMs + gVCFs published under the SNP-calling cohort dirs — see below)
```

## SNP-calling outputs (multi-platform GATK4 + Clair3)

Per-sample artefacts are grouped by `cohort_id = "${reference_tag}_${tier}"`. Tier is
`HIGH` for Illumina/hybrid samples (GATK HaplotypeCaller) and `MEDIUM` for
Nanopore-only samples (Clair3). Cross-tier comparison is intentionally not
reported because the two callers calibrate confidence differently.

```
<outdir>/
├── 01_refs/<reference_tag>/             prepared reference (faidx + dict + bwa-mem2 index + minimap2 mmi + mask.bed)
├── 03_bams/<cohort_id>/<sample>.dedup.bam[.bai]    MarkDuplicates output
├── 04_gvcf/<cohort_id>/<sample>.g.vcf.gz[.tbi]     per-sample gVCF (HC) or Clair3 gVCF
├── 05_joint/<cohort_id>/<cohort>.joint.vcf.gz      joint-genotyped multi-sample VCF
├── 06_filtered/<cohort_id>/                         ploidy + tier aware filtered SNP VCF + summary.tsv
├── 07_alignment/<cohort_id>/
│   ├── <cohort>.full.aln                          ALL sites (REF/ALT/N), feeds snp-dists
│   └── <cohort>.snps_only.aln                     variable sites only, feeds IQ-TREE
├── 08_snpdist/<cohort_id>/
│   ├── <cohort>.snpdist.tsv                       absolute pairwise SNP matrix
│   └── <cohort>.snpdist.molten.tsv                long form (i, j, dist)
├── 09_phylogeny/<cohort_id>/<cohort>.treefile     IQ-TREE2 Newick with bootstrap + SH-aLRT
└── 10_network/<cohort_id>/
    ├── <cohort>.upgma.svg                         UPGMA epi network, coloured by samplesheet `batch`
    ├── <cohort>.upgma.gexf                        Gephi-importable
    └── <cohort>.transmission_pairs.tsv            pairs ≤ params.transmission_threshold_snps
```

### Cohort gating

Cohorts with fewer than `params.cohort_min_samples` samples (default 3) skip joint
genotyping with a warn — per-sample gVCFs are still emitted to `04_gvcf/` but
nothing downstream is produced for that cohort.

### Why two alignments

snp-dists runs over `full.aln` (REF/ALT/N at every site) because that gives
**absolute pairwise SNP distances** interpretable for transmission analysis.
IQ-TREE runs over `snps_only.aln` (snp-sites -c) because the ascertainment-bias
correction (`+ASC`) needs the input to be variable sites only. The two outputs
are complementary, not redundant — the UPGMA network consumes the snp-dists
matrix (absolute SNPs), the topological tree is the IQ-TREE output.

## Top-level outputs

```
<outdir>/
├── multiqc/
│   ├── multiqc_report.html
│   ├── multiqc_data/
│   └── multiqc_plots/
└── pipeline_info/
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    ├── pipeline_dag_*.html
    ├── epicandi_software_mqc_versions.yml
    └── params_*.json
```

MultiQC aggregates FastQC, fastp, NanoPlot, QUAST, BUSCO, samtools stats, and the per-process software versions collected via topic channels.

## Pipeline information

Nextflow's standard `execution_report`, `execution_timeline`, `execution_trace` and `pipeline_dag` are written to `pipeline_info/` for every run. Use them to diagnose runtime/memory issues. `params_*.json` records the exact parameter set that was applied.
