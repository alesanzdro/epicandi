# Genomas de Referencia para C. auris

## Descripción General

Este documento describe los genomas de referencia de *Candida auris* utilizados en el pipeline EpiCandi. Los genomas están organizados por clados filogenéticos y representan cepas de alta calidad de diferentes regiones geográficas.

## Clados de C. auris

*C. auris* se clasifica en cinco clados filogenéticos principales basados en análisis genómicos:

- **Clado I**: Sudeste asiático (India, Pakistán)
- **Clado II**: Este asiático (Japón, Corea)
- **Clado III**: África
- **Clado IV**: Sudamérica (Venezuela)
- **Clado V**: Irán

## Genomas de Referencia Disponibles

### Clado I (Sudeste Asiático)

| Cepa | Accession | Calidad | Tecnología | Origen | Año | Descripción |
|------|-----------|---------|------------|--------|-----|-------------|
| **B11205** | GCA_016772135.1 | ⭐⭐⭐ Complete | PacBio | India | N/A | Genoma completo, alta calidad |
| **B8441** | GCA_002759435.3 | ⭐⭐⭐ Chromosome | PacBio HGAP v3 | Pakistán | 2008 | Nivel cromosómico, referencia snpEff |
| 6684 | GCA_001189475.1 | ⭐⭐ Scaffold | Illumina | India, Bengaluru | 2013 | Referencia histórica snpEff |

### Clado II (Este Asiático)

| Cepa | Accession | Calidad | Tecnología | Origen | Año | Descripción |
|------|-----------|---------|------------|--------|-----|-------------|
| **B11220** | GCA_003013715.2 | ⭐⭐⭐ Complete | Oxford Nanopore + Flye v2.4.2 | Japón | 2009 | Genoma completo, referencia snpEff |
| JCM15448 | GCA_007168705.1 | ⭐⭐ Contig | MiSeq + múltiples herramientas | N/A | 2005 | Ensamblado híbrido complejo |

### Clado III (África)

| Cepa | Accession | Calidad | Tecnología | Origen | Año | Descripción |
|------|-----------|---------|------------|--------|-----|-------------|
| **B11221** | GCA_002775015.1 | ⭐⭐ Scaffold | PacBio HGAP v3 | Sudáfrica | 2012 | Única referencia del clado III |

### Clado IV (Sudamérica)

| Cepa | Accession | Calidad | Tecnología | Origen | Año | Descripción |
|------|-----------|---------|------------|--------|-----|-------------|
| **B11245** | GCA_008275145.1 | ⭐⭐⭐ Complete | Oxford Nanopore + Canu v1.5 | Venezuela | 2012 | Genoma completo, referencia snpEff |
| B11243 | GCA_003014415.1 | ⭐⭐ Scaffold | Illumina + SPAdes v3.1.1 | Venezuela | 2013 | Referencia snpEff |

### Clado V (Irán)

| Cepa | Accession | Calidad | Tecnología | Origen | Año | Descripción |
|------|-----------|---------|------------|--------|-----|-------------|
| **B18474** | GCA_016809505.1 | ⭐⭐⭐ Complete | Oxford Nanopore + Canu v1.6 | Irán, Sari | 2018 | Genoma más reciente, clado emergente |

## Clasificación de Calidad

- **⭐⭐⭐ Complete/Chromosome**: Genomas de máxima calidad, secuencias continuas o cromosómicas
- **⭐⭐ Scaffold**: Genomas de buena calidad con algunos gaps
- **⭐ Contig**: Genomas fragmentados pero utilizables

## Recomendaciones de Uso

### Para Análisis Filogenético
- **Recomendado**: Usar todas las referencias de alta calidad (⭐⭐⭐)
- **Mínimo**: Al menos una referencia por clado

### Para Alineamiento y Mapeo
- **Clado I**: B11205 (completo) o B8441 (cromosómico)
- **Clado II**: B11220 (completo)
- **Clado III**: B11221 (única opción)
- **Clado IV**: B11245 (completo)
- **Clado V**: B18474 (completo)

### Para Anotación Génica
Las cepas marcadas con "snpEff" tienen anotaciones funcionales validadas:
- B8441 (Clado I)
- B11220 (Clado II)
- B11221 (Clado III)
- B11245 (Clado IV)

## Uso en el Pipeline

### Descarga de Referencias

```bash
# Descargar una referencia específica
snakemake --cores 4 --use-conda download_reference --wildcards reference_id=CladeI_B11205

# Descargar todas las referencias
snakemake --cores 8 --use-conda download_all_references

# Crear índices para alineamiento
snakemake --cores 4 --use-conda setup_reference_indexes
```

### Archivos Generados

Para cada referencia se generan:
- `{reference_id}.fasta`: Secuencia del genoma
- `{reference_id}.gff3`: Anotación génica
- `{reference_id}_info.txt`: Metadatos de la cepa
- `{reference_id}.fasta.fai`: Índice samtools
- `{reference_id}.fasta.bwt`: Índice BWA (y archivos asociados)

### Ubicación de Archivos
```
resources/references/
├── CladeI_B11205.fasta
├── CladeI_B11205.gff3
├── CladeI_B11205_info.txt
├── CladeI_B11205.fasta.fai
├── CladeI_B11205.fasta.bwt
├── ...
└── references_summary.txt
```

## Consideraciones Importantes

### Diversidad Genética
- La divergencia entre clados es significativa (~1-3% SNPs)
- Importante usar referencia del clado correcto para análisis precisos
- Para muestras de origen desconocido, usar múltiples referencias

### Características del Genoma
- **Tamaño**: ~12-13 Mb
- **GC Content**: ~45%
- **Cromosomas**: 7 (en genomas completos)
- **Genes**: ~5,400-5,800
- **Naturaleza diploide**: Considerar heterocigosidad

### Limitaciones Tecnológicas
- Genomas Illumina: Pueden tener gaps en regiones repetitivas
- Genomas PacBio/Nanopore: Mayor continuidad pero posibles errores de homopolímeros
- Genomas híbridos: Mejor balance entre continuidad y precisión

## Referencias Bibliográficas

1. Lockhart SR, et al. (2017). *Simultaneous Emergence of Multidrug-Resistant Candida auris on 3 Continents Confirmed by Whole-Genome Sequencing and Epidemiological Analyses.* Clin Infect Dis. 64(2):134-140.

2. Chow NA, et al. (2020). *Multiple introductions and subsequent transmission of multidrug-resistant Candida auris in the USA.* Nat Commun. 11(1):4577.

3. Muñoz JF, et al. (2018). *Genomic insights into multidrug-resistance, mating and virulence in Candida auris and related emerging species.* Nat Commun. 9(1):5346.

---

*Documento actualizado: $(date)*