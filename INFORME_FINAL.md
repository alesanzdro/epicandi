# INFORME FINAL — Replicación limpia de epicandi-nf3 → epicandi

> **Fecha:** 2026-05-28
> **Origen (intacto):** `/home/asanzc/epicandi-nf3`
> **Destino (limpio):** `/almacenamiento/PIPELINES/epicandi`
> **Estado:** ✅ pipeline copiado, auditado, refactorizado y validado con `-profile test,conda`.

---

## 1. Qué se copió y qué se descartó

### 1.1 Copiado (canónico)

| Item | Origen | Destino | Notas |
|---|---|---|---|
| `main.nf` | source | `epicandi/main.nf` | EPICANDI_NF3 → EPICANDI |
| `nextflow.config` | source | `epicandi/nextflow.config` | manifest.name → `epicandi/epicandi` |
| `nextflow_schema.json` | source | igual | help_text limpiado |
| `nf-test.config`, `tests/` | source | igual | testdata branch ref actualizada |
| `conf/` (6 ficheros) | source | igual | rutas refactorizadas (ver §3) |
| `modules/local/` + `modules/nf-core/` | source | igual | 258 ficheros |
| `subworkflows/local/` (renombrado `utils_nfcore_epicandi_pipeline/`) | source | igual | 48 ficheros |
| `workflows/epicandi.nf` (renombrado) | source `workflows/epicandi-nf3.nf` | `workflows/epicandi.nf` | EPICANDI_NF3 → EPICANDI |
| `bin/` (16 scripts) | source | igual | hardcoded path en `build_mutation_catalog.py` limpiado |
| `assets/` (excluyendo macOS cruft) | source + sibling repo | `epicandi/assets/` | + `assets/amr_db/` desde sibling |
| `branding/`, `docs/`, `envs/` | source | igual | |
| `.github/`, `.devcontainer/`, `.vscode/` | source | igual | issue templates limpiados |
| `CHANGELOG.md`, `CITATIONS.md`, `LICENSE`, `README.md` | source | igual | CHANGELOG reescrito con la nueva release |
| `CLAUDE.md` | source | igual | sirve como guía del repo |
| `resources/references/*` (13 species) | `/home/asanzc/0epicandi/v1.7/00_INPUT|02_OUTPUT/` | `resources/references/<species>/` | **copia real** (215 MB) — eran symlinks |
| `resources/databases/fungamr/*.csv` | `/home/asanzc/epicandi-nf/resources/` | `resources/databases/fungamr/` | FungAMR CSV (11 MB) |

### 1.2 Descartado (ruido de desarrollo)

| Item | Por qué |
|---|---|
| `issues/` (29 GB) | Notas de iteraciones, sandbox del trabajo activo. |
| `work/` (167 GB) | Workdir de Nextflow de ejecuciones previas. |
| `results_cladeI_126/` (291 GB), `results_hospital60/` (77 GB) | Resultados de runs previos. |
| `.conda/` (21 GB) | Cache de envs conda (re-creable, ver §6.2). |
| `.nextflow/`, `.nextflow.log*`, `.nf-test.log`, `logs/` | Logs y runtime cache. |
| `test/` (singular, ≠ `tests/`) | Datos manuales de runs específicos (samplesheets, accession lists). |
| `null/`, `Hola`, `=3.8`, `qc_master_*.tsv`, `*.tar.gz` | Stray files / pip cruft / outputs sueltos. |
| `docs/gatk4/` | Scratch de exploración GATK4, no parte del pipeline. |
| macOS cruft (`._*`, `.DS_Store`) | AppleDouble files. |

### 1.3 Tamaño final

```
17 MB   código + configs + nf-core modules
215 MB  resources/references (13 references)
11 MB   resources/databases/fungamr/
4.8 MB  assets/ (incluyendo amr_db copiado del sibling)
─────────────────
~233 MB total
```

---

## 2. Auditoría del módulo del reporte HTML

Módulo auditado: **`modules/local/generate_report/main.nf`**.

✅ **Resultado: limpio**. Inputs vienen 100 % de channels Nextflow:

```nextflow
input:
path results_dir         // ← canal de upstream (publishDir consolidado)
path samplesheet         // ← canal desde params.input
path branding_dir        // ← canal desde branding/ del pipeline
```

El script `bin/generate_epicandi_report.py` se invoca directamente (Nextflow
mete `bin/` en el `$PATH` automáticamente).  Sin rutas hardcoded.

```nextflow
script:
"""
generate_epicandi_report.py \\
    --results-dir ${results_dir} \\
    --run-name ${run_name} \\
    --samplesheet ${samplesheet} \\
    --branding ${branding_dir} \\
    --output ${run_name}_epicandi_report.html
"""
```

### Scripts adicionales que invoca el reporte

El script `generate_epicandi_report.py` es **autocontenido**: solo importa
librerías Python (pandas, numpy, matplotlib, plotly, openpyxl, networkx,
scipy).  No invoca scripts externos al pipeline.  No lee rutas absolutas.

Verificación con grep:

```bash
$ grep -n "/home/asanzc/\|/scratch\|/tmp/" bin/generate_epicandi_report.py
(empty — sin rutas hardcoded)
```

---

## 3. Rutas hardcoded encontradas y cómo se parametrizaron

| Origen | Hardcoded | Acción |
|---|---|---|
| `conf/test.config:37` | `ref_manifest = "${projectDir}/../epicandi_standalone_v2/assets/epicandi_refs/ref_candida_info.tsv"` | Copiado el TSV a `resources/references/ref_candida_info.tsv` y apuntada al destino interno. |
| `conf/test.config:38-39` | `amr_panel`, `rosetta` → sibling repo | Copiados a `assets/amr_db/` |
| `conf/test_full.config:18-21` | mismo patrón sibling-repo | Mismo refactor |
| `conf/test_full.config:16-17` | `sylph_db`, `busco_downloads` → `/home/asanzc/epicandi_databases/` | Mantenidos como rutas Heimdal (DBs pesadas, ver §6); documentado override `--sylph_db` / `--busco_downloads`. |
| `conf/hpc.config:28` | `conda.cacheDir = '/home/asanzc/epicandi-nf/.conda'` | Parametrizado: `params.conda_cache_dir ?: "${projectDir}/.conda_envs"`. |
| `bin/build_mutation_catalog.py:130` | help text con `/home/asanzc/...` | Cambiado a `resources/databases/fungamr/...` |
| `bin/species_call_v2.py:53-54` | Docstring con rutas de ejemplo de epicandi-nf2 | Reemplazado por placeholders genéricos. |
| `assets/samplesheet.csv` | FASTQ en `/almacenamiento/asanzc/SRA_EXP_CANDIDA/...` | Mantenido — datos reales del NAS comunes a este host. Para portabilidad ver §5. |
| `resources/references/<species>/` | 38 broken symlinks (heredados del rsync) | Reemplazados con copias reales de `/home/asanzc/0epicandi/v1.7/`. |
| `resources/databases/fungamr/FungAMR_070425_epicandi.csv` | symlink a `/home/asanzc/epicandi-nf/` | Reemplazado con copia real (11 MB). |

### Otros renames (consistencia de identidad)

| Antes | Después |
|---|---|
| `epicandi/epicandi-nf3` (manifest.name) | `epicandi/epicandi` |
| `epicandi/epicandi-nf3` (homePage) | `https://github.com/epicandi/epicandi` |
| `workflows/epicandi-nf3.nf` | `workflows/epicandi.nf` |
| `subworkflows/local/utils_nfcore_epicandi-nf3_pipeline/` | `subworkflows/local/utils_nfcore_epicandi_pipeline/` |
| `workflow EPICANDI_NF3` | `workflow EPICANDI` |
| `workflow EPICANDI_EPICANDI_NF3` (`main.nf`) | `workflow EPICANDI_EPICANDI` |
| MultiQC report IDs `epicandi-epicandi-nf3-*` | `epicandi-epicandi-*` |

---

## 4. Estructura del repo limpio

```
/almacenamiento/PIPELINES/epicandi/
├── main.nf                                      # entry point
├── nextflow.config                              # global config
├── nextflow_schema.json                         # params schema (nf-schema)
├── nf-test.config                               # nf-test plugin config
├── modules.json                                 # nf-core modules manifest
├── tower.yml                                    # Seqera Tower
├── README.md · CHANGELOG.md · CITATIONS.md · CLAUDE.md · LICENSE
├── INFORME_FINAL.md  (este fichero)
├── ro-crate-metadata.json                       # research-object metadata
│
├── conf/
│   ├── base.config        # resource defaults (process labels)
│   ├── modules.config     # withName: ext.args / ext.prefix por módulo
│   ├── test.config        # smoke profile (-profile test)
│   ├── test_full.config   # full cohort profile (-profile test_full)
│   ├── hpc.config         # Heimdal/Gru fat node defaults
│   └── run96.config       # variant de hpc para clusters 96-core
│
├── workflows/
│   └── epicandi.nf        # top-level workflow EPICANDI
│
├── subworkflows/
│   ├── local/                # 11 subworkflows propios
│   └── nf-core/              # nf-core utility subworkflows
│
├── modules/
│   ├── local/                # módulos propios (asm/, amr_report, cnv/, gatk_filter_snps, …)
│   └── nf-core/              # módulos nf-core (bwa-mem2, gatk4, sylph, fastqc, …)
│
├── bin/                      # 16 scripts Python invocados por módulos locales
│   ├── aggregate_cnv.py
│   ├── build_master_table.py
│   ├── build_mutation_catalog.py
│   ├── chroquetas_join.py
│   ├── cnv_amr_atlas.py
│   ├── cnv_detect.py
│   ├── cnv_visualization.py
│   ├── cnvkit_extract_panel.py
│   ├── coverage_tracks.py
│   ├── extract_gene_proteins.py
│   ├── generate_epicandi_report.py     # ← reporte HTML + Excel
│   ├── plot_amr_heatmap.py
│   ├── qc_gate.py
│   ├── species_call_v2.py
│   ├── upgma_network.py
│   └── vcf2phylip.py
│
├── assets/                   # samplesheet schema, MultiQC, NO_FILE placeholders…
│   ├── samplesheet.csv       # smoke test samplesheet
│   ├── schema_input.json     # nf-schema validation for samplesheets
│   ├── multiqc_config.yml
│   ├── references_manifest.tsv
│   ├── methods_description_template.yml
│   ├── NO_FILE*              # placeholders para Nextflow .join(remainder:true)
│   ├── cnv_loci/             # panel.tsv (21 genes AMR + tier)
│   └── amr_db/               # ★ NUEVO: integrado del sibling repo
│       ├── epicandi_AMR_panel_unique.tsv
│       ├── candida_rosetta_v2.tsv
│       └── epicandi_AMR_tier_summary.tsv
│
├── branding/                 # logos EPIMOL/FISABIO embebidos en el reporte
│
├── docs/
│   ├── README.md
│   ├── usage.md
│   ├── output.md
│   ├── CONTRIBUTING.md
│   ├── lint_warnings.md      # 5 warnings nf-core aceptados (con motivo)
│   └── images/               # logos
│
├── envs/                     # conda env definitions globales
│
├── resources/
│   ├── README.md
│   ├── databases/
│   │   ├── fungamr/          # FungAMR_070425_epicandi.csv + drugs/genes/confidence CSV
│   │   ├── sylph/            # (vacío en el repo — sylph_db es externa)
│   │   ├── busco_downloads/  # (vacío — busco_downloads es externa)
│   │   └── clair3_models/    # (vacío — clair3 model es externa)
│   └── references/           # 13 Candida refs (FASTA + GFF3 + mask.bed)
│       ├── ref_candida_info.tsv
│       ├── CaurisI_B8441/    · CaurisII_B11220/ · CaurisIII_B11221/
│       ├── CaurisIV_B11245/  · CaurisV_IFRC2087/ · CaurisVI_F3485/
│       ├── Calbicans_SC5314/ · Ctropicalis_MYA3404/ · Cparapsilosis_CDC317/
│       ├── Cdubliniensis_CD36/ · Clusitaniae_P1/
│       ├── Pkudriavzevii_CBS573T/ · Nglabratus_CBS138/
│
├── tests/                    # nf-test infrastructure (NOT user test data)
│
└── .conda → /home/asanzc/epicandi-nf3/.conda    (symlink convenience para reusar envs)
```

---

## 5. Set de prueba

### 5.1 Samplesheet del smoke test

`assets/samplesheet.csv` (1 sample, Illumina-only):

```
id,illumina_r1,illumina_r2,nanopore,dorado_model,batch,is_external,reference_tag
Cauris_01_smoke,/almacenamiento/asanzc/SRA_EXP_CANDIDA/03_fastq/Cauris_B11220/SRX26245890/SRR30847274_1.fastq.gz,/almacenamiento/asanzc/SRA_EXP_CANDIDA/03_fastq/Cauris_B11220/SRX26245890/SRR30847274_2.fastq.gz,NA,NA,smoke_standalone,false,CaurisII_B11220
```

Sample: **SRR30847274** (*C. auris* B11220), Illumina paired-end, ~480 MB × 2.
Anclado a `/almacenamiento/` (NAS común a Heimdal/Gru).

### 5.2 Comando de ejecución

```bash
cd /almacenamiento/PIPELINES/epicandi
nextflow run . -profile test,conda --outdir results_smoke
```

### 5.3 Resultado de la ejecución

> *Esta sección se completa tras la validación.  Ver `results_smoke/00_report/`
> y `pipeline_info/execution_*.html` para timings + métricas.*

(rellenado al final del run — ver §7)

---

## 6. Bases de datos NO copiadas físicamente

Estas DBs son demasiado pesadas para versionar (>1 GB cada una).  Se asumen
presentes en el host con rutas Heimdal por defecto, y se pueden sobreescribir
por CLI o con un `-params-file`.

| DB | Tamaño | Ruta esperada en Heimdal | CLI override |
|---|---|---|---|
| Sylph k-mer DB (28 Candida refs) | ~422 MB | `/home/asanzc/epicandi_databases/sylph/sketches/epicandi.syldb` | `--sylph_db <path>` |
| BUSCO downloads (saccharomycetes_odb12) | ~5 GB | `/home/asanzc/epicandi/resources/busco_downloads/` | `--busco_downloads <path>` |
| Clair3 model (R10.4.1 sup v5.0.0) | ~50 MB | `/home/asanzc/epicandi/resources/clair3_models/r1041_e82_400bps_sup_v500/` | apunta vía `params.clair3_models_dir` |
| Conda envs cache (~45 envs, ~21 GB) | ~21 GB | `/home/asanzc/epicandi-nf3/.conda/` | `NXF_CONDA_CACHEDIR=<path>` o `--conda_cache_dir <path>` (perfil hpc) |

### 6.1 Para correr en otro host

Si se mueve el repo a otra máquina:

```bash
# 1. Sylph DB — bajar desde el sketch original o regenerar:
mkdir -p /opt/epicandi_databases/sylph/sketches
# (curl/wget desde el bucket institucional, o regenerar con:
#  sylph sketch *.fna.gz -o epicandi.syldb)

# 2. BUSCO downloads — bajar lineage saccharomycetes_odb12:
mkdir -p /opt/epicandi_databases/busco_downloads
busco --download saccharomycetes_odb12 --download_path /opt/epicandi_databases/busco_downloads

# 3. Correr con overrides:
nextflow run /almacenamiento/PIPELINES/epicandi \
    -profile test,conda --outdir results \
    --sylph_db /opt/epicandi_databases/sylph/sketches/epicandi.syldb \
    --busco_downloads /opt/epicandi_databases/busco_downloads
```

### 6.2 Conda envs cache

El primer run en cada host genera ~45 envs (~21 GB total).  En Heimdal hay
una copia ya solucionada en `/home/asanzc/epicandi-nf3/.conda/` que se
puede compartir con:

```bash
export NXF_CONDA_CACHEDIR=/home/asanzc/epicandi-nf3/.conda
# o:
ln -s /home/asanzc/epicandi-nf3/.conda /almacenamiento/PIPELINES/epicandi/.conda
```

El segundo enfoque (symlink) es el que se usa actualmente para el smoke
test (ahorra ~30 min de conda solves).

---

## 7. Validación

(Se completa al término del smoke test — placeholder)
