# Implementación de Checkpoint Dinámico - EpiCandi Pipeline

## Resumen

Se ha implementado exitosamente un sistema de filtrado dinámico de muestras usando **checkpoints** de Snakemake en el pipeline `Snakefile_robust`. Esta implementación permite al pipeline:

1. **Iniciar con todas las muestras** disponibles en `samplesinfo.csv`
2. **Ejecutar validación de calidad** para determinar qué muestras pasan los umbrales
3. **Filtrar dinámicamente** las muestras para ensamblado basándose en los resultados de QC
4. **Adaptarse automáticamente** el DAG a las muestras validadas

## Arquitectura del Checkpoint

### Checkpoint Principal: `samples_validation`

**Ubicación**: `rules/qc_validation.smk:321`

**Función**: Analiza la tabla de métricas de QC y determina qué muestras pueden proceder al ensamblado.

**Salidas clave**:
```
output/01_data/01.5_samples_pass/passed_samples.txt          # Lista de todas las muestras validadas
output/01_data/01.5_samples_pass/hybrid_samples.txt         # Muestras para ensamblado híbrido
output/01_data/01.5_samples_pass/short_reads_samples.txt    # Muestras para ensamblado Illumina
output/01_data/01.5_samples_pass/long_reads_samples.txt     # Muestras para ensamblado Nanopore
output/01_data/01.5_samples_pass/checkpoint_validation_report.txt  # Reporte detallado
```

### Funciones de Lectura Dinámica

**Ubicación**: `rules/common.smk:211`

**Funciones principales**:

1. **`get_validated_samples(wildcards)`**
   - Lee `passed_samples.txt` del checkpoint
   - Devuelve lista de todas las muestras que pasaron validación

2. **`get_validated_samples_by_strategy(wildcards, strategy)`**
   - Lee archivos específicos por estrategia
   - Strategies: `'hybrid'`, `'short_reads'`, `'long_reads'`

3. **`get_assembly_outputs(wildcards)`**
   - Genera dinámicamente las salidas de ensamblado
   - Solo para muestras validadas

### Reglas Agregadas Dinámicas

**Ubicación**: `rules/assembly.smk:305`

**Reglas implementadas**:
- `assemblies_hybrid`: Ensambla solo muestras híbridas validadas
- `assemblies_short_reads`: Ensambla solo muestras Illumina validadas
- `assemblies_long_reads`: Ensambla solo muestras Nanopore validadas
- `all_assemblies`: Consolida todos los ensamblados

## Flujo de Ejecución

### Fase 1: Validación Inicial (Checkpoint)
```bash
# Ejecutar hasta la validación de muestras
snakemake -s Snakefile_robust output/01_data/01.5_samples_pass/passed_samples.txt --cores 8
```

### Fase 2: Ensamblado Dinámico
```bash
# Ejecutar solo ensamblados validados
snakemake -s Snakefile_robust validated_assemblies --cores 16

# O ejecutar tipos específicos
snakemake -s Snakefile_robust assemblies_hybrid --cores 16
snakemake -s Snakefile_robust assemblies_short_reads --cores 16
snakemake -s Snakefile_robust assemblies_long_reads --cores 16
```

### Fase 3: Pipeline Completo
```bash
# Ejecutar pipeline completo con filtrado dinámico
snakemake -s Snakefile_robust --cores 16
```

## Criterios de Filtrado

Las muestras se filtran basándose en:

1. **Umbrales de QC**:
   - Illumina: ≥ 50,000 lecturas limpias (`illumina_min_reads`)
   - Nanopore: ≥ 5,000 lecturas limpias (`nanopore_min_reads`)

2. **Estrategia de Ensamblado**:
   - `hybrid`: Tanto Illumina como Nanopore pasan QC
   - `short_reads`: Solo Illumina pasa QC
   - `long_reads`: Solo Nanopore pasa QC
   - `fail`: Ninguna plataforma pasa QC (excluida del ensamblado)

## Ventajas de la Implementación

### ✅ Robustez
- **Manejo seguro de errores**: Funciones con try/catch y fallbacks
- **Validación de archivos**: Verificación de existencia de archivos antes de leer
- **Logging detallado**: Información completa sobre decisiones de filtrado

### ✅ Flexibilidad
- **Adaptación automática**: DAG se ajusta dinámicamente a muestras disponibles
- **Estrategias múltiples**: Soporte para ensamblados híbridos, Illumina-only, Nanopore-only
- **Configuración fácil**: Umbrales ajustables en `config.yaml`

### ✅ Transparencia
- **Reportes detallados**: Información completa sobre muestras filtradas y razones
- **Trazabilidad**: Archivos de checkpoint permiten auditar decisiones
- **Logging completo**: Registro de todas las decisiones de filtrado

### ✅ Eficiencia
- **Sin desperdicios**: Solo se procesan muestras que pasan QC
- **Paralelización**: Diferentes estrategias pueden ejecutarse en paralelo
- **Reutilización**: Checkpoint evita re-ejecución innecesaria

## Archivos de Salida del Checkpoint

```
output/01_data/01.5_samples_pass/
├── passed_samples.txt                    # ARCHIVO PRINCIPAL DEL CHECKPOINT
├── hybrid_samples.txt                    # Muestras para ensamblado híbrido
├── short_reads_samples.txt              # Muestras para ensamblado Illumina
├── long_reads_samples.txt               # Muestras para ensamblado Nanopore
├── checkpoint_validation_report.txt     # Reporte detallado de validación
├── samples_metrics_table.tsv           # Tabla completa de métricas
└── validation_report.txt               # Reporte de validación original
```

## Debugging y Troubleshooting

### Verificar muestras validadas:
```bash
cat output/01_data/01.5_samples_pass/passed_samples.txt
```

### Ver reporte de validación:
```bash
cat output/01_data/01.5_samples_pass/checkpoint_validation_report.txt
```

### Ejecutar solo validación:
```bash
snakemake -s Snakefile_robust samples_validation --cores 4
```

### Verificar DAG dinámico:
```bash
snakemake -s Snakefile_robust --dry-run --quiet --cores 1
```

## Compatibilidad

- **Snakemake**: Versión 5.0+ (checkpoints)
- **Python**: 3.6+
- **Pandas**: Para procesamiento de tablas de métricas

## Comparación con Implementación Anterior

| Aspecto | Anterior | Nueva Implementación |
|---------|----------|---------------------|
| **Filtrado** | Estático, basado en listas predefinidas | Dinámico, basado en checkpoint |
| **Robustez** | Fallbacks a "todas las muestras" | Manejo de errores específico |
| **Transparencia** | Limitada | Reportes detallados en cada paso |
| **Eficiencia** | Procesamiento de muestras fallidas | Solo muestras validadas |
| **Adaptabilidad** | Manual | Automática basada en QC |

## Conclusión

La implementación del checkpoint dinámico convierte el pipeline EpiCandi en un sistema verdaderamente robusto y eficiente, que:

1. **Se adapta automáticamente** a la calidad de los datos
2. **Evita desperdiciar recursos** en muestras de baja calidad
3. **Proporciona transparencia completa** sobre las decisiones de filtrado
4. **Mantiene compatibilidad** con el workflow existente

El sistema está listo para uso en producción y cumple con las mejores prácticas de Snakemake para checkpoints dinámicos.