# Guía de Uso del Sistema de Referencias

## Comandos Principales

### Descargar Referencias Individuales

```bash
# Descargar una referencia específica del Clado I (recomendada)
snakemake --cores 4 --use-conda resources/references/CladeI_B11205.fasta

# Descargar referencia del Clado II
snakemake --cores 4 --use-conda resources/references/CladeII_B11220.fasta

# Descargar referencia del Clado IV
snakemake --cores 4 --use-conda resources/references/CladeIV_B11245.fasta
```

### Descargar Todas las Referencias

```bash
# Setup completo: descarga + índices + lista FastANI (RECOMENDADO)
snakemake --cores 8 --use-conda setup_references

# Alternativamente, solo descargar referencias (sin índices)
snakemake --cores 8 --use-conda download_all_references
```

## Referencias Disponibles

| ID | Cepa | Clado | Calidad | Recomendado para |
|----|------|-------|---------|------------------|
| `CladeI_B11205` | B11205 | I | ⭐⭐⭐ Complete | Análisis general |
| `CladeI_B8441` | B8441 | I | ⭐⭐⭐ Chromosome | Anotación, snpEff |
| `CladeI_6684` | 6684 | I | ⭐⭐ Scaffold | Análisis histórico |
| `CladeII_B11220` | B11220 | II | ⭐⭐⭐ Complete | Análisis general |
| `CladeII_JCM15448` | JCM15448 | II | ⭐⭐ Contig | Análisis específico |
| `CladeIII_B11221` | B11221 | III | ⭐⭐ Scaffold | Única ref. del clado |
| `CladeIV_B11245` | B11245 | IV | ⭐⭐⭐ Complete | Análisis general |
| `CladeIV_B11243` | B11243 | IV | ⭐⭐ Scaffold | Análisis alternativo |
| `CladeV_B18474` | B18474 | V | ⭐⭐⭐ Complete | Clado más reciente |

## Archivos Generados

Para cada referencia se crean los siguientes archivos:

```
resources/references/
├── CladeI_B11205.fasta          # Secuencia del genoma
├── CladeI_B11205.gff3           # Anotación génica
├── CladeI_B11205_info.txt       # Metadatos de la cepa
├── CladeI_B11205.fasta.fai      # Índice samtools
├── CladeI_B11205.fasta.amb      # Índices BWA
├── CladeI_B11205.fasta.ann      # (múltiples archivos)
├── CladeI_B11205.fasta.bwt
├── CladeI_B11205.fasta.pac
├── CladeI_B11205.fasta.sa
└── references_summary.txt        # Resumen de todas las referencias
```

## Ejemplos de Uso Específicos

### 1. Análisis Filogenético Completo
```bash
# Descargar todas las referencias de alta calidad
snakemake --cores 8 --use-conda \
  resources/references/CladeI_B11205.fasta \
  resources/references/CladeII_B11220.fasta \
  resources/references/CladeIII_B11221.fasta \
  resources/references/CladeIV_B11245.fasta \
  resources/references/CladeV_B18474.fasta
```

### 2. Análisis de Clado Específico
```bash
# Solo referencias del Clado I para comparación detallada
snakemake --cores 4 --use-conda \
  resources/references/CladeI_B11205.fasta \
  resources/references/CladeI_B8441.fasta \
  resources/references/CladeI_6684.fasta
```

### 3. Preparación para Mapeo
```bash
# Descargar referencia e indexar para BWA
snakemake --cores 4 --use-conda resources/references/CladeI_B11205.fasta.bwt

# Indexar para samtools
snakemake --cores 4 --use-conda resources/references/CladeI_B11205.fasta.fai
```

## Verificación de Descargas

```bash
# Verificar que los archivos se descargaron correctamente
ls -la resources/references/

# Verificar contenido del resumen
cat resources/references/references_summary.txt

# Verificar metadatos de una referencia específica
cat resources/references/CladeI_B11205_info.txt
```

## Integración con Análisis

### Usando referencias en reglas personalizadas

```python
# Ejemplo de regla que usa una referencia
rule align_to_reference:
    input:
        reads="data/{sample}.fastq.gz",
        ref="resources/references/CladeI_B11205.fasta",
        ref_index="resources/references/CladeI_B11205.fasta.bwt"
    output:
        "alignments/{sample}_aligned.bam"
    shell:
        "bwa mem {input.ref} {input.reads} | samtools sort > {output}"
```

### Selección automática de referencia

```python
# Función para seleccionar referencia por clado
def get_best_reference_for_clade(clade):
    best_refs = {
        "I": "CladeI_B11205",
        "II": "CladeII_B11220",
        "III": "CladeIII_B11221",
        "IV": "CladeIV_B11245",
        "V": "CladeV_B18474"
    }
    return f"resources/references/{best_refs[clade]}.fasta"
```

## Troubleshooting

### Error: "datasets command not found"
```bash
# Verificar instalación de ncbi-datasets-cli
conda list | grep datasets

# Si no está instalado, activar ambiente correcto
conda activate epicandi
```

### Error: "Permission denied" o "Connection failed"
```bash
# Verificar conectividad
curl -I https://api.ncbi.nlm.nih.gov/datasets/v2alpha/

# Si hay proxy, configurar variables de entorno
export https_proxy=http://tu-proxy:puerto
export http_proxy=http://tu-proxy:puerto
```

### Error: "Accession not found"
```bash
# Verificar que el accession existe
datasets summary genome accession GCA_016772135.1

# Listar todas las referencias disponibles de C. auris
datasets summary genome taxon "Candida auris"
```

### Archivos incompletos o corruptos
```bash
# Limpiar y volver a descargar
rm -rf resources/references/CladeI_B11205.*
snakemake --cores 4 --use-conda resources/references/CladeI_B11205.fasta

# Forzar re-descarga
snakemake --cores 4 --use-conda --forcerun download_reference \
  --wildcards reference_id=CladeI_B11205
```

## Consideraciones de Espacio y Tiempo

### Tamaños aproximados
- Cada genoma FASTA: ~12-15 MB
- Cada archivo GFF3: ~3-5 MB
- Índices BWA: ~50-60 MB por genoma
- **Total para todas las referencias**: ~650 MB

### Tiempos de descarga
- Referencia individual: 30-60 segundos
- Todas las referencias: 8-15 minutos
- Indexado completo: 5-10 minutos

### Recomendaciones
```bash
# Para uso ocasional: descargar solo referencias necesarias
snakemake --cores 4 --use-conda resources/references/CladeI_B11205.fasta

# Para análisis completo: descargar todo
snakemake --cores 8 --use-conda download_all_references

# Para análisis de producción: setup completo en un comando
snakemake --cores 8 --use-conda setup_references
```

---

*Para más detalles técnicos, consultar `docs/REFERENCES.md` y `docs/NCBI_DOWNLOAD.md`*