# EpiCandi – Pipeline para Nanopore (Snakemake)

Este pipeline automatiza QC, limpieza, ensamblaje, pulido, taxonomía, evaluación de calidad y reportes MultiQC para lecturas largas de Nanopore, diseñado específicamente para análisis de Candida.

## Características principales

- Procesamiento completo de lecturas Nanopore
- Soporte para genomas de referencia de diferentes clados de _C. auris_ 
- Detección de mutaciones asociadas a resistencia a antifúngicos
- Reportes MultiQC integrados
- Compatibilidad con entornos Conda/Mamba

## Estructura de carpetas

```
├── Snakefile
├── config.yaml
├── envs
│   ├── epinano_annotation.yml
│   ├── epinano_checkm2.yml
│   └── epinano.yml
├── input/
├── INSTALL_CONDA.md
├── output/
│   ├── 01_data/
│   │   ├── 01.1_fastq_raw/
│   │   ├── 01.2_nanoplot_raw/
│   │   ├── 01.3_fastq_filtered/
│   │   └── 01.4_nanoplot_filtered/
│   ├── 02_assembly/
│   │   ├── 02.1_flye/
│   │   ├── 02.2_medaka/
│   │   ├── 02.3_quast/
│   │   ├── 02.4_checkm2/
│   │   └── 02.5_prokka/
│   ├── 03_taxonomy/
│   │   └── 03.1_kraken2/
│   └── 04_report/
│       ├── 04.1_collection_renamed/
│       └── 04.2_multiqc/
├── logs/
│   ├── 01.1_copy_raw_data/
│   ├── 01.2_nanoplot_raw/
│   ├── 01.3.1_porechop/
│   ├── 01.3.2_filtlong/
│   ├── 01.4_nanoplot_filtered/
│   ├── 02.1_flye_assembly/
│   ├── 02.2_medaka_consensus/
│   ├── 02.3_quast/
│   ├── 02.4_checkm2/
│   ├── 02.5_prokka_annotation/
│   ├── 03.1_kraken2_classification/
│   ├── 03.2_krona/
│   ├── 04.1_collection_renamed/
│   ├── 04.1_multiqc_report/
│   └── 04.2_multiqc_report/
└── resources/
```

Notas:
- Coloca los FASTQ de entrada en `input/` (se aceptan rutas en subcarpetas). El pipeline excluye automáticamente controles que contengan etiquetas como: `negative`, `control`, `blanco`, `unclassified`, `negativo` (mayúsculas/minúsculas indiferentes).
- Las bases de datos/recursos se resuelven desde `resources/` (p. ej., CheckM2 DB, Krona taxonomy, QUAST setup).

## Flujo del pipeline (resumen)

Por muestra:
- 01_data
  - Copia a `01.1_fastq_raw/`
  - QC inicial con NanoPlot → `01.2_nanoplot_raw/`
  - Limpieza: Porechop → `01.3_fastq_filtered/*_porechop.fastq.gz`
  - Filtrado: Filtlong → `01.3_fastq_filtered/*_porechop_filtlong.fastq.gz`
  - QC post-filtrado con NanoPlot → `01.4_nanoplot_filtered/`
- 02_assembly
  - Ensamblaje (Flye) → `02.1_flye/{sample}/assembly.fasta`
  - Pulido (Medaka) → `02.2_medaka/{sample}/consensus.fasta`
  - QUAST → `02.3_quast/{sample}/report.html`
  - CheckM2 → `02.4_checkm2/{sample}/quality_report.tsv`
  - Anotación (Prokka) → `02.5_prokka/{sample}/annotation.*`
- 03_taxonomy
  - Kraken2 → `03.1_kraken2/{sample}_report.txt`
  - (Opcional) Krona → `03.2_krona/{sample}_krona.html`
- 04_report
  - `collect_reports`: consolida archivos clave en `04.1_collection_renamed/`
    - NanoStats (filtradas): `nanostats/{sample}_NanoStats.txt`
    - CheckM2 (re-etiquetado del nombre de muestra): `checkm2/{sample}_quality_report.tsv`
    - Prokka: GFF y logs
    - Kraken2: `kraken2/{sample}.txt` (sin sufijo `_report`)
  - `multiqc_report`: genera `04.2_multiqc/multiqc_report.html` buscando en `04.1_collection_renamed/` y en `02.3_quast/`

## Requisitos
- Conda/Mamba (ver `INSTALL_CONDA.md`).
- Snakemake con soporte conda. Ejemplo: `conda install -c conda-forge -c bioconda snakemake`.

## Ejecución básica

Lanzar el pipeline completo (respetando outputs existentes):
```bash
snakemake --cores 16 --use-conda --rerun-triggers mtime
```
Esto realizará, por muestra: copia → QC → trimming/filtrado → ensamblaje → pulido → QUAST → CheckM2 → anotación → consolidación → MultiQC.

Sólo CheckM2 en todas las muestras (asumiendo consensos listos):
```bash
snakemake --cores 8 --use-conda checkm2_quality_all
```

Sólo consolidación + MultiQC (rápido, sin tocar cómputo pesado):
```bash
snakemake --cores 4 --use-conda multiqc_report
```
Opcional (forzar que sólo se ejecuten esas dos reglas):
```bash
snakemake multiqc_report --allowed-rules collect_reports multiqc_report --cores 4 --use-conda
```

Desbloquear si quedó un bloqueo de Snakemake:
## Ejecución recomendada (64 hilos, Kraken2=1, memoria 256+124 GB)

En una sesión screen para dejarlo corriendo en segundo plano:

```bash
screen -S epinanox -dm bash -lc 'cd $(pwd) && bash scripts/run_epinanox.sh'
screen -ls   # ver la sesión
screen -r epinanox  # adjuntarse cuando quieras
```

El script `scripts/run_epinanox.sh` lanza Snakemake con:
- `--cores 64`
- `--resources mem_mb=389120 kraken_db=1` (1 Kraken2 a la vez, memoria total 256+124 GB)
- `--rerun-triggers mtime` y `--printshellcmds`

```bash
snakemake --unlock
```

## Esquemas del pipeline (DAG / rulegraph)

Requiere Graphviz para exportar a imagen. Instálalo si no lo tienes:
```bash
sudo apt-get update && sudo apt-get install -y graphviz
```

Generar el DAG de trabajos (dot + png):
```bash
snakemake --dag | awk 'f||/^digraph/{f=1;print}' | tee output/04_report/pipeline_dag.dot | dot -Tpng -o output/04_report/pipeline_dag.png
```

Generar el grafo de reglas (dot + png):
```bash
snakemake --rulegraph | awk 'f||/^digraph/{f=1;print}' | tee output/04_report/pipeline_rulegraph.dot | dot -Tpng -o output/04_report/pipeline_rulegraph.png
```

Notas:
- Si quieres limitar el DAG a un objetivo (p.ej. MultiQC): `snakemake -n multiqc_report --dag | awk 'f||/^digraph/{f=1;print}' | dot -Tpng -o output/04_report/dag_multiqc.png`.
- Para DAGs grandes, puedes generar sólo `.dot` y visualizarlo con herramientas que soporten zoom.

## Limpieza

- Limpiar temporales de filtrado:
```bash
snakemake clean
```
- Eliminar todos los outputs generados:
```bash
snakemake clean_all
```
- Limpiar metadatos/procedencia si se tocaron scripts pero no quieres re-ejecutar:
```bash
snakemake --cleanup-metadata output/04_report/04.1_collection_renamed/collection_collected.flag
```

## Configuración (`config.yaml`)

Parámetros principales (ejemplos):
- checkm2_db: ruta a la DB de CheckM2 (puede ser relativa, se resuelve a absoluta en el Snakefile).
  - Ej.: `resources/CheckM2_database/uniref100.KO.1.dmnd`
- checkm2_extension: extensión que CheckM2 buscará en la carpeta de entrada. Por defecto: `.fasta`.
- kraken_db: ruta a la base de datos de Kraken2.
- kraken_conf: umbral de confianza para Kraken2 (ej. `0.1`).
- quast_params: flags extra para QUAST (opcional).
- resources: recursos por regla (threads, mem_mb, walltime) y valores por defecto.

Ejemplo mínimo de bloque de recursos:
```yaml
resources:
  default:
    threads: 4
    mem: 8000
    walltime: 120
  nanoplot:
    threads: 4
    mem: 8000
    walltime: 120
  flye:
    threads: 16
    mem: 64000
    walltime: 720
  # ... (porechop, filtlong, medaka, quast, checkm2, prokka, kraken2, multiqc)
```

Notas de bases de datos:
- CheckM2: si no existe, `setup_checkm2` intentará descargar/registrar la DB bajo `resources/`. El pipeline exporta `CHECKM2DB` automáticamente.
- Krona: `setup_krona` crea/enlaza la taxonomía en `resources/krona/taxonomy` y ejecuta `ktUpdateTaxonomy.sh`.
- QUAST: `setup_quast` descarga SILVA/BUSCO y marca `resources/quast/.setup_done`.

## Consejos y solución de problemas
- Si Snakemake quiere re-ejecutar muchos pasos sólo por cambios de código/entorno, añade `--rerun-triggers mtime` o limpia metadatos específicos con `--cleanup-metadata`.
- Si ya tienes outputs correctos pero faltan marcas de tiempo, puedes usar `--touch` sobre ficheros concretos para marcarlos como hechos.
- Si usas sólo los reportes, recuerda que `collect_reports` espera, como mínimo, los `quality_report.tsv` de CheckM2.

## Créditos
- Entornos Conda en `envs/`.
- Lógica principal en `Snakefile` y utilidades en `scripts/`.
