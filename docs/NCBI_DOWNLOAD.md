# Guía de Descarga con Herramientas NCBI

## Descripción General

Este documento describe cómo usar las herramientas oficiales de NCBI para descargar genomas de referencia de *Candida auris*. El pipeline EpiCandi utiliza `datasets` y `dataformat` del NCBI Datasets toolkit.

## Herramientas NCBI Requeridas

### 1. NCBI Datasets
El NCBI Datasets es una herramienta de línea de comandos para descargar datos genómicos estructurados.

#### Instalación
```bash
# Opción 1: A través de conda (recomendado)
conda install -c conda-forge ncbi-datasets-cli

# Opción 2: Descarga directa
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets
sudo mv datasets /usr/local/bin/

# Opción 3: A través del ambiente EpiCandi
conda activate epicandi  # Ya incluye ncbi-datasets-cli
```

#### Verificación
```bash
datasets --version
# datasets version: 16.27.0
```

### 2. NCBI Dataformat
Herramienta complementaria para extraer metadatos de las descargas.

```bash
# Verificar instalación
dataformat --version
# dataformat version: 16.27.0
```

## Comandos Básicos de Descarga

### Descarga Individual por Accession

```bash
# Ejemplo básico: descargar genoma completo
datasets download genome accession GCA_016772135.1

# Incluir archivos específicos
datasets download genome accession GCA_016772135.1 --include genome,gff3,protein

# Especificar formato de descarga
datasets download genome accession GCA_016772135.1 --include genome,gff3 --filename cauris_b11205.zip
```

### Tipos de Archivos Disponibles

| Parámetro | Descripción | Archivos incluidos |
|-----------|-------------|-------------------|
| `genome` | Secuencia del genoma | `*.fna` (FASTA nucleótidos) |
| `gff3` | Anotación génica | `*.gff` (anotaciones) |
| `protein` | Secuencias proteicas | `*.faa` (FASTA aminoácidos) |
| `rna` | Secuencias de ARN | `*.frn` (FASTA ARN) |
| `cds` | Secuencias codificantes | `*.fna` (CDS) |
| `seq-report` | Reporte de secuencias | `sequence_report.jsonl` |

### Descarga por Taxonomía

```bash
# Descargar todos los genomas de C. auris
datasets download genome taxon "Candida auris"

# Limitar número de genomas
datasets download genome taxon "Candida auris" --limit 5

# Filtrar por nivel de ensamblado
datasets download genome taxon "Candida auris" --assembly-level complete,chromosome
```

### Filtros Avanzados

```bash
# Filtrar por fuente de datos
datasets download genome taxon "Candida auris" --assembly-source refseq

# Filtrar por fecha de publicación
datasets download genome taxon "Candida auris" --released-after 2020-01-01

# Excluir genomas específicos
datasets download genome taxon "Candida auris" --exclude-atypical
```

## Extracción y Organización

### Estructura de Archivos Descargados

```
ncbi_dataset.zip
├── README.md
├── md5sum.txt
└── ncbi_dataset/
    ├── data/
    │   ├── assembly_data_report.jsonl
    │   └── GCA_XXXXXXX.X/
    │       ├── GCA_XXXXXXX.X_STRAIN_genomic.fna
    │       ├── genomic.gff
    │       └── protein.faa
    └── metadata/
        └── assembly_data_report.jsonl
```

### Extracción Automatizada

```bash
# Descomprimir y extraer
unzip ncbi_dataset.zip

# Localizar archivos específicos
find ncbi_dataset -name "*.fna" -type f
find ncbi_dataset -name "*.gff" -type f

# Renombrar con patrón consistente
mv ncbi_dataset/data/*/*.fna genome.fasta
mv ncbi_dataset/data/*/*.gff annotation.gff3
```

### Script de Procesamiento Completo

```bash
#!/bin/bash
# Script para descarga y procesamiento de C. auris

ACCESSION="GCA_016772135.1"
OUTPUT_NAME="CladeI_B11205"

# Descargar
echo "Descargando $ACCESSION..."
datasets download genome accession $ACCESSION --include genome,gff3

# Extraer
echo "Extrayendo archivos..."
unzip -q ncbi_dataset.zip

# Renombrar y organizar
echo "Organizando archivos..."
mv ncbi_dataset/data/*/*.fna ${OUTPUT_NAME}.fasta
mv ncbi_dataset/data/*/*.gff ${OUTPUT_NAME}.gff3

# Obtener metadatos
echo "Extrayendo metadatos..."
dataformat tsv genome --package ncbi_dataset.zip > ${OUTPUT_NAME}_metadata.tsv

# Limpiar
rm -rf ncbi_dataset ncbi_dataset.zip README.md md5sum.txt

echo "Descarga completada: ${OUTPUT_NAME}.fasta y ${OUTPUT_NAME}.gff3"
```

## Extracción de Metadatos

### Usar dataformat para metadatos

```bash
# Extraer información en formato TSV
dataformat tsv genome --package ncbi_dataset.zip

# Campos específicos
dataformat tsv genome --fields accession,organism-name,asm-stats-contig-n50 --package ncbi_dataset.zip

# Formato JSON para procesamiento automatizado
dataformat json genome --package ncbi_dataset.zip | jq '.assemblies[].assembly.assembly_accession'
```

### Campos de Metadatos Útiles

| Campo | Descripción |
|-------|-------------|
| `accession` | Número de acceso |
| `organism-name` | Nombre del organismo |
| `asm-stats-total-sequence-len` | Longitud total del genoma |
| `asm-stats-contig-n50` | N50 de contigs |
| `asm-stats-number-of-contigs` | Número de contigs |
| `assembly-level` | Nivel de ensamblado |
| `assembly-type` | Tipo de ensamblado |
| `submission-date` | Fecha de envío |

## Automatización en el Pipeline

### Integración con Snakemake

El archivo `rules/references.smk` implementa la descarga automatizada:

```python
# Ejemplo de regla Snakemake
rule download_reference:
    output:
        fasta="resources/references/{reference_id}.fasta",
        gff3="resources/references/{reference_id}.gff3"
    params:
        accession=lambda wildcards: get_accession(wildcards.reference_id)
    shell:
        """
        datasets download genome accession {params.accession} --include genome,gff3
        unzip -q ncbi_dataset.zip
        mv ncbi_dataset/data/*/*.fna {output.fasta}
        mv ncbi_dataset/data/*/*.gff {output.gff3}
        rm -rf ncbi_dataset ncbi_dataset.zip README.md md5sum.txt
        """
```

### Uso en el Pipeline EpiCandi

```bash
# Descargar referencia específica
snakemake --cores 4 --use-conda resources/references/CladeI_B11205.fasta

# Descargar todas las referencias
snakemake --cores 8 --use-conda download_all_references

# Verificar descargas
snakemake --cores 1 resources/references/references_summary.txt
```

## Solución de Problemas

### Errores Comunes

#### Error de conexión
```bash
# Error: Failed to download
# Solución: Verificar conectividad y proxy
curl -I https://api.ncbi.nlm.nih.gov/datasets/v2alpha/

# Configurar proxy si es necesario
export https_proxy=http://proxy.servidor:puerto
export http_proxy=http://proxy.servidor:puerto
```

#### Accession no encontrado
```bash
# Error: Accession not found
# Verificar accession en NCBI
datasets summary genome accession GCA_XXXXXXX.X

# Buscar accessions alternativos
datasets summary genome taxon "Candida auris" | head -10
```

#### Archivos corruptos
```bash
# Verificar integridad con checksums
md5sum -c md5sum.txt

# Re-descargar si es necesario
rm ncbi_dataset.zip
datasets download genome accession GCA_XXXXXXX.X --include genome,gff3
```

### Limitaciones y Consideraciones

#### Límites de descarga
- NCBI puede limitar descargas masivas
- Implementar pausas entre descargas: `sleep 5`
- Usar `--limit` para controlar número de genomas

#### Tamaño de archivos
- Genomas completos pueden ser grandes (>50MB)
- Planificar espacio en disco
- Considerar compresión: `gzip *.fasta`

#### Versiones de accessions
- Usar versión más reciente: `.1`, `.2`, etc.
- Verificar actualizaciones periódicamente
- Documentar versiones utilizadas

## Recursos Adicionales

### Documentación Oficial
- [NCBI Datasets Documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/)
- [Command Line Tools User Guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/)

### APIs y Servicios Web
- [NCBI Datasets REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)

### Scripts de Ejemplo
Ver directorio `scripts/` del pipeline para ejemplos adicionales de procesamiento automatizado.

---

*Documento actualizado: $(date)*