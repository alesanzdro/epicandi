# epicandi/epicandi

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

## Introduction

**epicandi/epicandi** is the Nextflow DSL2 rewrite of the EpiCandi v2 standalone (bash) pipeline.
It runs *Candida* spp. sequencing data — Illumina short-read, Oxford Nanopore long-read or hybrid — through QC, cleaning, identification, assembly, post-assembly QC, antifungal resistance scanning, per-sample SNP calling and (when ≥3 samples per species/clade) cohort-level phylogeny.

Per-sample topology, identical to the bash version:

1. **QC RAW** — FastQC (Illumina) / NanoPlot (Nanopore)
2. **CLEAN** — fastp · porechop_abi + chopper
3. **DECONT** *(optional, default OFF for SRA cohorts)* — bowtie2/minimap2 vs PhiX + Lambda
4. **QC CLEAN** — FastQC / NanoPlot rerun
5. **QC GATE** — PASS/WARN/FAIL on read count + contamination fraction
6. **IDENTIFICATION** — Sylph (sketch + profile) · per-read screen (bowtie2 / minimap2) · AuriClass (*C. auris* clade) · `species_call_v2`
7. **SUBSAMPLE** — target 100× coverage with 1.1× hysteresis
8. **ASSEMBLY** — SPAdes (Illumina) · Flye + medaka (Nanopore) · + Polypolish (Hybrid / Illumina polish)
9. **FINALIZE** — minlen ≥500 bp + rename to `${id}_contig{nr}`
10. **POST_ASM_QC** — QUAST `-r --fungus` · BUSCO `saccharomycetes_odb12`
11. **AMR** — ChroQueTaS (FungAMR panel) · `chroquetas_join` (5-state matrix)
12. **SNP_CALLING** — per-platform variant calling (multi-platform v2 architecture):
    * Illumina / hybrid (tier HIGH): BWA-MEM2 → GATK4 MarkDuplicates → HaplotypeCaller (ploidy-aware, `--exclude-intervals mask.bed`) → per-sample gVCF
    * Nanopore-only (tier MEDIUM): filtlong → minimap2 → Clair3 (`--bed_fn` accessible regions)
    * Cohorts grouped by `reference_tag × tier`; samples below `params.cohort_min_samples` skip joint genotyping with a warn

Cohort step (added in v3):

13. **COHORT GENOTYPING + PHYLOGENY + NETWORK** (per cohort with ≥`params.cohort_min_samples`):
    * GATK4 CombineGVCFs → GenotypeGVCFs → VariantFiltration (ploidy + tier aware)
    * vcf2phylip → `full.aln` + `snps_only.aln` (snp-sites -c)
    * snp-dists over `full.aln` → absolute pairwise SNP distance matrix
    * IQ-TREE2 over `snps_only.aln` (`-m GTR+G+ASC -B 1000 -alrt 1000`)
    * UPGMA epidemiological network (NetworkX) coloured by samplesheet `batch`, with transmission threshold `params.transmission_threshold_snps` (default 6 SNPs)

## Usage

> First-time Nextflow / nf-core? Read the [setup guide](https://nf-co.re/docs/get_started/environment_setup/overview) and test with `-profile test` before running on real data.

Samplesheet (`samplesheet.csv`):

```csv
id,illumina_r1,illumina_r2,nanopore,dorado_model,batch,is_external
Cauris_01,/path/R1.fq.gz,/path/R2.fq.gz,NA,NA,clinical_2026,false
EPI001,NA,NA,/path/nano.fq.gz,r1041_e82_400bps_sup_v5.2.0,clinical_2026,false
HYB001,/path/R1.fq.gz,/path/R2.fq.gz,/path/nano.fq.gz,r1041_e82_400bps_sup_v5.2.0,clinical_2026,false
```

- `illumina_r1/r2 = NA` → nanopore-only sample
- `nanopore = NA` → illumina-only sample
- both set → hybrid

Run:

```bash
nextflow run epicandi/epicandi \
   -profile <conda/docker/singularity/hpc> \
   --input samplesheet.csv \
   --outdir results/ \
   --sylph_db /path/epicandi.syldb \
   --refs_dir /path/epicandi_refs/refs \
   --ref_manifest /path/ref_candida_info.tsv \
   --amr_panel /path/epicandi_AMR_panel_unique.tsv \
   --rosetta /path/candida_rosetta_v2.tsv \
   --busco_downloads /path/busco_downloads
```

Smoke test (single *C. auris* B11220 Illumina-only sample):

```bash
nextflow run epicandi/epicandi -profile test,conda --outdir results_smoke
```

> Provide pipeline parameters via the CLI or a `-params-file`. Custom `-c` configs must not redefine parameters — see [nf-core docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

epicandi/epicandi was originally written by Alejandro Sanz-Carbonell (FISABIO / IIS La Fe), as the Nextflow rewrite of the `epicandi_standalone_v2` bash pipeline.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in [`CITATIONS.md`](CITATIONS.md).

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
