#!/usr/bin/env python3
"""
Pipeline Snakemake robusto para análisis de Candida auris
Sistema avanzado de gestión de muestras con validaciones de integridad
Jerarquía: samplesinfo.csv ID > normalización de archivos (fallback)
"""

import os
import pandas as pd
import re
import sys
from pathlib import Path

# ==============================
# CONFIGURACIÓN
# ==============================

configfile: "config.yaml"

# ==============================
# IMPORTAR REGLAS Y FUNCIONES COMUNES
# ==============================

include: "rules/common.smk"

# ==============================
# CARGA DE MUESTRAS CON VALIDACIONES
# ==============================

try:
    SAMPLES, SAMPLE_TYPES = load_samples()
    SAMPLE_NAMES = list(SAMPLES.keys())

    # Clasificar muestras por tipo
    ILLUMINA_SAMPLES = [s for s, t in SAMPLE_TYPES.items() if t in ['illumina', 'hybrid']]
    NANOPORE_SAMPLES = [s for s, t in SAMPLE_TYPES.items() if t in ['nanopore', 'hybrid']]
    HYBRID_SAMPLES = [s for s, t in SAMPLE_TYPES.items() if t == 'hybrid']

    print("=" * 50)
    print(f"RESUMEN DE MUESTRAS:")
    print(f"Total: {len(SAMPLE_NAMES)}")
    print(f"Illumina: {len(ILLUMINA_SAMPLES)}")
    print(f"Nanopore: {len(NANOPORE_SAMPLES)}")
    print(f"Hibridas: {len(HYBRID_SAMPLES)}")
    print("=" * 50)

except Exception as e:
    print(f"ERROR CRITICO en carga de muestras: {e}")
    sys.exit(1)

# ==============================
# PREPARACIÓN DE DATOS CON ENLACES SIMBÓLICOS
# ==============================

rule stage_raw_data:
    output:
        # Marcador para indicar que todos los enlaces se han creado
        touch("output/01_data/00_staged_data/.staging_complete")
    log:
        "output/logs/00_staged_data/staging.log"
    run:
        import os
        import pandas as pd
        
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        with open(log[0], "w") as log_file:
            log_file.write("=== ORGANIZANDO DATOS CRUDOS ===\n")
            
            # Crear directorios para enlaces simbólicos
            illdir = "output/01_data/00_staged_data/illumina"
            nanodir = "output/01_data/00_staged_data/nanopore"
            os.makedirs(illdir, exist_ok=True)
            os.makedirs(nanodir, exist_ok=True)
            
            log_file.write(f"Creados directorios de datos: {illdir}, {nanodir}\n")
            
            # Crear enlaces simbólicos con nombres estandarizados
            for sample_id, data in SAMPLES.items():
                # Enlaces para datos Illumina
                if data['illumina_r1'] and data['illumina_r2']:
                    r1_link = os.path.join(illdir, f"{sample_id}_R1.fastq.gz")
                    r2_link = os.path.join(illdir, f"{sample_id}_R2.fastq.gz")
                    
                    # Crear enlaces simbólicos relativos (mejor para portabilidad)
                    try:
                        if os.path.exists(r1_link):
                            os.unlink(r1_link)
                        if os.path.exists(r2_link):
                            os.unlink(r2_link)
                            
                        os.symlink(os.path.abspath(data['illumina_r1']), r1_link)
                        os.symlink(os.path.abspath(data['illumina_r2']), r2_link)
                        log_file.write(f"Enlaces para Illumina {sample_id}: {r1_link}, {r2_link}\n")
                    except Exception as e:
                        log_file.write(f"ERROR creando enlaces Illumina para {sample_id}: {e}\n")
                
                # Enlaces para datos Nanopore
                if data['nanopore']:
                    nano_link = os.path.join(nanodir, f"{sample_id}.nanopore.fastq.gz")
                    
                    try:
                        if os.path.exists(nano_link):
                            os.unlink(nano_link)
                            
                        os.symlink(os.path.abspath(data['nanopore']), nano_link)
                        log_file.write(f"Enlace para Nanopore {sample_id}: {nano_link}\n")
                    except Exception as e:
                        log_file.write(f"ERROR creando enlace Nanopore para {sample_id}: {e}\n")
            
            log_file.write("=== STAGING COMPLETADO ===\n")

# ==============================
# REGLA PRINCIPAL
# ==============================

def get_all_outputs(wildcards=None):
    """Generar dinámicamente las salidas según el tipo de muestras"""
    outputs = []

    # Primero asegurarnos que los datos están correctamente organizados
    outputs.append("output/01_data/00_staged_data/.staging_complete")

    # Outputs para muestras con Illumina
    if ILLUMINA_SAMPLES:
        # Utilizamos extend para cada expand individual para evitar listas anidadas
        outputs.extend(expand("output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R1_fastqc.html", sample=ILLUMINA_SAMPLES))
        outputs.extend(expand("output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R2_fastqc.html", sample=ILLUMINA_SAMPLES))
        outputs.extend(expand("output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz", sample=ILLUMINA_SAMPLES))
        outputs.extend(expand("output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz", sample=ILLUMINA_SAMPLES))
        outputs.extend(expand("output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R1_fastqc.html", sample=ILLUMINA_SAMPLES))
        outputs.extend(expand("output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R2_fastqc.html", sample=ILLUMINA_SAMPLES))
        outputs.extend(expand("output/01_data/01.4_kraken2/illumina/{sample}_report.txt", sample=ILLUMINA_SAMPLES))

    # Outputs para muestras con Nanopore
    if NANOPORE_SAMPLES:
        outputs.extend(expand("output/01_data/01.1_fastq_raw_qc/nanoplot/{sample}/NanoPlot-report.html", sample=NANOPORE_SAMPLES))
        outputs.extend(expand("output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz", sample=NANOPORE_SAMPLES))
        outputs.extend(expand("output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}/NanoPlot-report.html", sample=NANOPORE_SAMPLES))
        outputs.extend(expand("output/01_data/01.4_kraken2/nanopore/{sample}_report.txt", sample=NANOPORE_SAMPLES))

    # Outputs comunes para todas las muestras
    if SAMPLE_NAMES:
        outputs.extend([
            "output/01_data/01.5_samples_pass/validation_report.txt",
            "output/01_data/01.5_samples_pass/samples_metrics_table.tsv",
            # CHECKPOINT: Archivos de validación dinámicos
            "output/01_data/01.5_samples_pass/passed_samples.txt",
            "output/01_data/01.5_samples_pass/checkpoint_validation_report.txt",
            "output/04_report/multiqc_report.html"
        ])

        # CHECKPOINT DINÁMICO: Añadir ensamblados basados en validación
        # Solo intentar obtener ensamblados si el checkpoint ya existe
        checkpoint_file = "output/01_data/01.5_samples_pass/passed_samples.txt"
        if os.path.exists(checkpoint_file):
            try:
                # Leer muestras validadas directamente desde el archivo
                with open(checkpoint_file, 'r') as f:
                    validated_samples = [line.strip() for line in f if line.strip()]

                # Añadir las salidas del ensamblado
                assembly_outputs = []
                assembly_outputs.extend(expand("output/02_assembly/02.3_consensus/{sample}.fasta", sample=validated_samples))
                assembly_outputs.extend([
                    "output/02_assembly/hybrid_assemblies_summary.txt",
                    "output/02_assembly/short_reads_assemblies_summary.txt",
                    "output/02_assembly/long_reads_assemblies_summary.txt",
                    "output/02_assembly/all_assemblies_summary.txt"
                ])
                outputs.extend(assembly_outputs)
                print(f"CHECKPOINT: Se añadieron {len(assembly_outputs)} ensamblados")

                # Añadir las evaluaciones de ensamblados
                evaluation_outputs = [
                    "output/02_assembly/02.4_evaluation/evaluation_summary.txt"
                ]
                outputs.extend(evaluation_outputs)
                print(f"CHECKPOINT: Se añadieron {len(evaluation_outputs)} evaluaciones")

                # Añadir análisis de cobertura
                coverage_outputs = []
                coverage_outputs.extend(expand("output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam", sample=validated_samples))
                coverage_outputs.extend(expand("output/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt", sample=validated_samples))
                outputs.extend(coverage_outputs)
                print(f"CHECKPOINT: Se añadieron {len(coverage_outputs)} análisis de cobertura")

                # Añadir caracterización completa
                characterization_outputs = [
                    "output/03_characterization/phylogenetic_characterization_report.txt",
                    "output/03_characterization/resistance_characterization_report.txt"
                ]
                outputs.extend(characterization_outputs)
                print(f"CHECKPOINT: Se añadieron {len(characterization_outputs)} caracterizaciones")
            except Exception as e:
                print(f"CHECKPOINT: Error obteniendo ensamblados/evaluaciones: {e}")
        else:
            print("CHECKPOINT: Archivo de validación no existe aún, se ejecutará en fase posterior")

    # DEBUG: Mostrar cuántos outputs se van a generar
    print(f"DEBUG: get_all_outputs() retornando {len(outputs)} archivos objetivo")
    return outputs

# ==============================
# IMPORTAR RESTO DE REGLAS
# ==============================

include: "rules/illumina_processing.smk"
include: "rules/nanopore_processing.smk"
include: "rules/qc_validation.smk"
include: "rules/assembly.smk"
include: "rules/assembly_evaluation.smk"
include: "rules/coverage_analysis.smk"
include: "rules/characterization_phylo.smk"
include: "rules/characterization_resistance.smk"
include: "rules/characterization_phylogenetic.smk"
include: "rules/reports.smk"



# ==============================
# REGLA ALL
# ==============================

rule all:
    input:
        lambda wildcards: get_all_outputs(wildcards)

rule assemblies:
    input:
        lambda wildcards: get_assembly_outputs(wildcards)

# Nueva regla para ejecutar solo los ensamblados validados
rule validated_assemblies:
    input:
        "output/02_assembly/all_assemblies_summary.txt"

# ==============================
# REGLAS DE LIMPIEZA
# ==============================

rule clean:
    shell:
        """
        rm -rf output/01_data/01.2_porechop_filtlong/
        rm -rf output/01_data/01.2_fastp/
        rm -rf output/01_data/01.4_kraken2/
        """

rule clean_all:
    shell:
        """
        rm -rf output/
        rm -rf logs/
        """