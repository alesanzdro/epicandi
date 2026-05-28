# epicandi/epicandi: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - 2026-05-28

First clean self-contained release.  Replicates the EpiCandi v2 standalone
bash pipeline (per-sample QC → identification → assembly → AMR → SNP
calling → cohort phylogeny + CNV).

### `Added`

- Top-level workflow EPICANDI (renamed from EPICANDI_NF3) over 11 ordered
  subworkflows:  `FASTQ_QC_RAW` → `FASTQ_CLEAN` → `FASTQ_QC_CLEAN`
  → `QC_GATE_SWF` → `IDENTIFICATION` → `SUBSAMPLE_SWF` → assembly fan-out
  (`ASSEMBLY_SHORT` / `ASSEMBLY_LONG` / `ASSEMBLY_HYBRID`) → `POST_ASM_QC`
  → `AMR_REPORT` → `SNP_CALLING` → `CNV` → `GENERATE_REPORT`.
- Multi-platform SNP-calling layer (BWA-MEM2 + GATK4 HaplotypeCaller for
  Illumina/hybrid; minimap2-ont + Clair3 for nanopore-only).
- Cohort-level joint genotyping + phylogeny + UPGMA epidemiological
  network, with transmission cut-off `params.transmission_threshold_snps`
  (default 12 SNPs, CDC *C. auris* standard).
- CNV detection (`mosdepth` + `CNVkit`) with consensus per-gene calling
  and an AMR × CNV combined heatmap in the HTML report.
- EPIMOL/FISABIO-branded HTML run report (`generate_epicandi_report.py`)
  with Plotly + matplotlib figures, batch-coloured MSN + UPGMA, and an
  Excel companion (`*.xlsx`, 6 sheets).
- Self-contained `resources/` and `assets/`:
  - 13 Candida species/clade references under `resources/references/`
    with associated GFF3 + TRF/Dust masks.
  - `assets/amr_db/` with FungAMR panel, rosetta drug-name table, and tier
    summary (extracted from `epicandi_standalone_v2`).
  - `resources/databases/fungamr/` with the FungAMR CSV catalog.

### `Fixed`

- All in-repo absolute paths replaced with `${projectDir}/...` so the
  pipeline runs from any location.
- Test profile no longer depends on the `epicandi_standalone_v2` sibling
  repo for AMR/reference assets.
- Drug-alias collision in `chroquetas_join.py` (`Amphotericin_B` was being
  split by the FungAMR regex on the `_`) — fixed via explicit alias map.

### `Dependencies`

- Nextflow ≥ 25.10.4
- Conda 23+ (or Docker / Singularity / Apptainer)
- Heavy external databases (assumed present at fixed Heimdal locations,
  overridable via CLI):
  - Sylph DB (`/home/asanzc/epicandi_databases/sylph/sketches/epicandi.syldb`)
  - BUSCO downloads (`/home/asanzc/epicandi/resources/busco_downloads/`)
  - Optional: Clair3 model directory for nanopore samples
