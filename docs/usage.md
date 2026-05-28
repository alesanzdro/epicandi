# epicandi/epicandi: Usage

> _Per-parameter documentation lives in `nextflow_schema.json`; see `nextflow run … --help_full`._

## Introduction

epicandi/epicandi ingests *Candida* spp. sequencing data (Illumina, Nanopore, or hybrid), then runs QC → identification → assembly → post-assembly QC → AMR → SNP calling → cohort phylogeny. The topology mirrors the `epicandi_standalone_v2` bash pipeline 1:1 and produces directly comparable outputs.

## Samplesheet input

Provide a CSV with the 7 columns below. Use the literal string `NA` for absent reads/models.

```csv title="samplesheet.csv"
id,illumina_r1,illumina_r2,nanopore,dorado_model,batch,is_external
Cauris_01,/path/R1.fq.gz,/path/R2.fq.gz,NA,NA,clinical_2026,false
EPI001,NA,NA,/path/nano.fq.gz,r1041_e82_400bps_sup_v5.2.0,clinical_2026,false
HYB001,/path/R1.fq.gz,/path/R2.fq.gz,/path/nano.fq.gz,r1041_e82_400bps_sup_v5.2.0,clinical_2026,false
```

| Column | Required | Description |
| --- | --- | --- |
| `id` | yes | Unique sample identifier (no spaces) |
| `illumina_r1` | yes (or `NA`) | Paired-end R1 (`.fastq.gz`/`.fq.gz`) |
| `illumina_r2` | yes (or `NA`) | Paired-end R2 (`.fastq.gz`/`.fq.gz`) |
| `nanopore` | yes (or `NA`) | Long reads (`.fastq.gz`/`.fq.gz`) |
| `dorado_model` | yes (or `NA`) | Dorado basecaller model (e.g. `r1041_e82_400bps_sup_v5.2.0`) — required for medaka consensus when `nanopore` is set |
| `batch` | no | Free-text cohort label; cohort-level phylogeny groups by `species_assigned`, not by `batch` |
| `is_external` | no | `true`/`false`; informational only (carried through to the master table) |

Platform is auto-detected from the row: `illumina-only`, `nanopore-only`, or `hybrid`. Validation is performed by `nf-schema` against `assets/schema_input.json`.

## Running the pipeline

```bash
nextflow run epicandi/epicandi \
    --input samplesheet.csv \
    --outdir results/ \
    -profile conda
```

The pipeline writes:

```
work/             Nextflow scratch (delete after the run)
results/          Pipeline outputs (see docs/output.md)
.nextflow_log     Nextflow log
```

Repeatable runs: pass settings via `-params-file params.yaml`:

```yaml title="params.yaml"
input: ./samplesheet.csv
outdir: ./results/
sylph_db: /home/asanzc/epicandi_databases/sylph/sketches/epicandi.syldb
busco_downloads: /home/asanzc/epicandi_databases/busco_downloads
refs_dir: /path/to/epicandi_refs/refs
ref_manifest: /path/to/ref_candida_info.tsv
amr_panel: /path/to/epicandi_AMR_panel_unique.tsv
rosetta: /path/to/candida_rosetta_v2.tsv
skip_decont: true
skip_snp_calling: false
target_coverage: 100
```

## Profiles

The pipeline bundles:

- `conda` / `mamba` — per-module conda envs (default for our HPC nodes; reuses the `/home/asanzc/epicandi-nf/.conda` cache via `hpc.config`)
- `docker` / `singularity` / `apptainer` / `podman` / `charliecloud` / `shifter` — container engines
- `test` — single Illumina-only smoke sample (*C. auris* B11220 / SRR30847274) reaching POST_ASM_QC + AMR
- `test_full` — same shape; override `--input` for a full cohort
- `hpc` — Heimdal/Gru fat node defaults (16 cpu / 128 GB / 48 h resource limits + `conda.cacheDir`)
- `wave` / `gpu` / `arm64` / `debug` — nf-core standard profiles

Multiple profiles compose: `-profile test,conda,hpc`.

## Required pre-built databases

| Param | Default | Notes |
| --- | --- | --- |
| `--sylph_db` | none | Sylph sketch DB built from the 27 reference genomes. Build with `sylph sketch -g <ref>.fna.gz` for each reference. |
| `--busco_downloads` | none | Pre-downloaded BUSCO lineage cache (`saccharomycetes_odb12`). |
| `--refs_dir` | none | Directory with one subdir per slug, each holding `genome.fna.gz` + `annotation.gff.gz`. |
| `--ref_manifest` | none | `ref_candida_info.tsv` (27 refs metadata). |
| `--amr_panel` | none | `epicandi_AMR_panel_unique.tsv`. |
| `--rosetta` | none | `candida_rosetta_v2.tsv` (species-name normalization). |
| `--screen_bt2_index` / `--screen_mmi` | optional | Pre-built read_screen indices; if absent the `PREPARE_SCREEN_DB` subworkflow builds them once at run time from `refs_dir`. |

## Core Nextflow arguments

`-profile`, `-resume`, `-c`, `-params-file`, `-r <release>`. See the [nf-core usage notes](https://nf-co.re/docs/running/run-pipelines).

## Resource customisation

Per-process resources come from labels in `conf/base.config` (`process_low/medium/high/long/high_memory`) and per-module `withName:` overrides in `conf/modules.config`. The smoke `test` profile caps the whole run at 16 cpu / 64 GB / 12 h via `process.resourceLimits`; `hpc` raises it to 16 / 128 GB / 48 h.

## Memory for Nextflow itself

Long runs benefit from:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```
