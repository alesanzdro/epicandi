# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline is

`epicandi/epicandi` is the Nextflow DSL2 rewrite of the `epicandi_standalone_v2` bash pipeline for *Candida* spp. sequencing. Behavioural parity with v2 is the design goal — per-module `ext.args`/`ext.prefix` in `conf/modules.config` mirror the v2 step scripts so a behavioural diff is one config diff away. When something in v2 changes, that file is the first thing to update.

Pipeline name on disk is `epicandi-nf3` (with dash) to disambiguate from prior `epicandi-nf` / `epicandi-nf2` attempts on the same HPC. Nextflow include aliases inside `main.nf` use underscores (`EPICANDI`, `EPICANDI_EPICANDI`) so the parser never sees the dash. See `docs/lint_warnings.md` for the five accepted lint warnings — do not try to "fix" them.

Not an nf-core/* pipeline: `manifest.name` is `epicandi/epicandi`. Several lint checks are intentionally skipped via `.nf-core.yml` for that reason.

## Common commands

```bash
# Smoke test (single C. auris B11220 Illumina-only sample, ends at AMR; PHYLOGENY needs ≥3 samples)
nextflow run . -profile test,conda --outdir results_smoke

# Full real run (requires all six external DB params)
nextflow run . -profile <conda|docker|singularity|hpc> \
   --input samplesheet.csv --outdir results/ \
   --sylph_db ... --refs_dir ... --ref_manifest ... \
   --amr_panel ... --rosetta ... --busco_downloads ...

# Resume after edits (caches per-process)
nextflow run . -profile test,conda --outdir results_smoke -resume

# nf-test (uses the `test` profile defined in nextflow.config)
nf-test test --tag pipeline --profile +conda --verbose
nf-test test --tag pipeline --profile +conda --update-snapshots   # after intentional output changes

# Lint
nf-core pipelines lint .
```

`nf-test` automatically ignores `modules/nf-core/**` and `subworkflows/nf-core/**` tests (see `nf-test.config`). Plugins are disabled there because nf-test does not honour `$https_proxy` on the HPC.

For long runs: `export NXF_OPTS='-Xms1g -Xmx4g'`.

## Architecture

### Top-level flow

`main.nf` → `PIPELINE_INITIALISATION` (samplesheet validation + auto-detect platform) → `EPICANDI` (the analysis) → `PIPELINE_COMPLETION`. The initialisation subworkflow lives at `subworkflows/local/utils_nfcore_epicandi_pipeline/`; it enriches each sample's meta map with `assembly_type` (`illumina` / `nanopore` / `hybrid`), `has_illumina`, `has_nanopore`, `is_external_bool`. Downstream branching depends on `meta.assembly_type`.

### Per-sample topology (`workflows/epicandi.nf`)

11 ordered subworkflows, all under `subworkflows/local/`:

1. `FASTQ_QC_RAW` — FastQC (Illumina) / NanoPlot (Nanopore)
2. `FASTQ_CLEAN` — fastp · porechop_abi + chopper
3. `FASTQ_QC_CLEAN` — same QC tools rerun on cleaned reads
4. `QC_GATE_SWF` — wraps `QC_GATE` module; computes PASS/WARN/FAIL on read count + contamination and **filters the samplesheet** (downstream only sees `samplesheet_pass`)
5. `IDENTIFICATION` — Sylph + per-read screen (bowtie2/minimap2) + AuriClass → `SPECIES_CALL` (combined classification)
6. `SUBSAMPLE_SWF` — target 100× coverage (with 1.1× hysteresis); emits `sub_illu`, `sub_nano`, `gsize`
7. Assembly fan-out by `meta.assembly_type`: `ASSEMBLY_SHORT` (SPAdes), `ASSEMBLY_LONG` (Flye + medaka), `ASSEMBLY_HYBRID` (+ Polypolish). Outputs are `.mix()`'d into one `ch_assembly`.
8. `POST_ASM_QC` — QUAST `--fungus` + BUSCO `saccharomycetes_odb12`
9. `AMR_REPORT` — ChroQueTaS (FungAMR) + `chroquetas_join.py` 5-state matrix
10. `SNP_CALLING` — multi-platform variant calling layer (v2 architecture, replaced Snippy in commit `1707c1c`):
    - Illumina/hybrid (tier HIGH): BWA-MEM2 → `GATK4_MARKDUPLICATES` → `GATK4_HAPLOTYPECALLER` (ploidy from `references_manifest.tsv`, mask via `--exclude-intervals`) → per-sample gVCF
    - Nanopore-only (tier MEDIUM): `FILTLONG_SNP` → `MINIMAP2_ONT` → `CLAIR3` (mask complement via `bedtools complement` + `--bed_fn`)
    - Cohorts grouped by `(reference_tag, tier)`; with `auris` subdivided by clade. Below `params.cohort_min_samples` (default 3) the joint-genotype + downstream chain is skipped with a warn.
11. **Cohort path** (per qualifying cohort): `GATK4_COMBINEGVCFS` → `GATK4_GENOTYPEGVCFS` → `GATK_FILTER_SNPS` (ploidy + tier aware filters) → `VCF2PHYLIP` (full.aln + snps_only.aln) → `SNPDISTS` (over full.aln, absolute pairwise SNPs) + `IQTREE` (snps_only.aln, +ASC) + `UPGMA_NETWORK` (SVG/GEXF coloured by batch, transmission pairs ≤ `params.transmission_threshold_snps`).

### Channel patterns to know

- **Platform-conditional channels.** Many subworkflows take all three of `ch_samplesheet`, `ch_clean_illu`, `ch_clean_nano` and branch internally on `meta.has_illumina` / `meta.has_nanopore` / `meta.assembly_type`. Don't try to pre-split at the workflow level.
- **`NO_FILE_*` placeholders** live in `assets/` (`NO_FILE`, `NO_FILE_AURI`, `NO_FILE_R1/R2`, `NO_FILE_NANO`, `NO_FILE_RS_ILLU/NANO`, `NO_FILE_CHROQUETAS`). They exist because Nextflow's `.join(remainder:true)` pads unmatched rows with nulls in unexpected positions. The convention is: pre-fill every per-sample channel with a placeholder, then `.mix()` real outputs in and `.groupTuple()` → pick-real-or-placeholder. See `IDENTIFICATION` (`subworkflows/local/identification/main.nf:108`) for the canonical example. **Reuse this pattern; don't introduce `remainder:true`.**
- **Topic versions channel.** Software versions are collected via the Nextflow `versions` topic channel (see end of `workflows/epicandi.nf`); local modules must emit into it. Do not hand-mix into `ch_versions`.

### Module organisation

- `modules/nf-core/` — managed by `nf-core modules install`, tracked in `modules.json`. Touch only via the CLI.
- `modules/local/` — pipeline-specific processes (one dir per tool): `asm/finalize`, `auriclass`, `chroquetas`, `clair3`, `filtlong`, `gatk_filter_snps`, `minimap2_ont`, `polypolish`, `pypolca`, `qc/gate`, `refs/prepare`, `screen/{fasta,tally}`, `species/call`, `subsample`, `upgma/network`, `vcf2phylip`. Each should follow nf-core conventions (tag, label, conda+container, topic versions, stub, `meta.yml`).
- `bin/` — Python helpers (`chroquetas_join.py`, `nj_tree.py`, `qc_gate.py`, `species_call_v2.py`) called from local modules via the shared `nf-python-helpers` conda env (Python 3.12 + pandas + numpy + biopython).

### Where to put which kind of change

| Want to change... | Edit... |
|---|---|
| Per-process flags / output paths | `conf/modules.config` (`withName:` blocks — **never inline `ext.args` in the process call**) |
| A new tunable | `nextflow.config` (`params { ... }`) AND `nextflow_schema.json` (run `nf-core pipelines schema build`) |
| Resource defaults | `conf/base.config` (`withLabel:` blocks) |
| HPC overrides | `conf/hpc.config` (Heimdal/Gru fat node defaults; conda cache at `/home/asanzc/epicandi-nf/.conda`) |
| Add a new pipeline step | New local subworkflow under `subworkflows/local/<step>/`, wire into `workflows/epicandi.nf` |
| Tool not on Bioconda | Conda env may list the `nmquijada` channel (already enabled — ChroQueTaS lives there) |

### Things that look weird but are intentional

- `qc_gate` and `subsample` subworkflows are named `*_SWF` because they wrap same-named modules. Documented in `docs/lint_warnings.md` #3, #4 — don't rename to make lint happy.
- `qc_gate/` only includes one module; the ≥2-modules-per-subworkflow lint warning is intentionally accepted.
- iGenomes config (`conf/igenomes*.config`) has been deleted — references come from `--refs_dir` / `--ref_manifest` (curated *Candida* refs), not iGenomes.

## Samplesheet

Validated by nf-schema against `assets/schema_input.json`. 7 columns: `id, illumina_r1, illumina_r2, nanopore, dorado_model, batch, is_external`. Use the literal string `NA` for absent reads/models — the init subworkflow's auto-detect treats `NA` as "missing" and sets `assembly_type` accordingly.

## Required external databases

`--sylph_db`, `--busco_downloads`, `--refs_dir`, `--ref_manifest`, `--amr_panel`, `--rosetta`. Optional `--screen_bt2_index` / `--screen_mmi` are built on the fly by `PREPARE_SCREEN_DB` from `refs_dir` if not provided. The smoke `test` profile points all of these at local paths under `/home/asanzc/`; expect those to be unavailable on other machines.

## Contributing

`docs/CONTRIBUTING.md` is the source of truth for the contribution workflow. Branch off `dev`, not `main`. Use nf-core modules upstream whenever available; only write a local module if no upstream module exists or the flags aren't exposed.
