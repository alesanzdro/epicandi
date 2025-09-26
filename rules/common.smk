# ==============================
# FUNCIONES COMUNES Y UTILIDADES
# ==============================

import os
import pandas as pd
import re
from pathlib import Path

# Parámetros principales
GENOME_SIZE = config.get("general", {}).get("genome_size", "13m")
MIN_LENGTH = config.get("general", {}).get("min_length", 1000)
MIN_MEAN_Q = config.get("general", {}).get("min_mean_q", 12)
OUTPUT_DIR = config.get("general", {}).get("output_dir", "output")

# Umbrales de QC
QC_MIN_READS_ILLUMINA = config.get("qc_thresholds", {}).get("illumina_min_reads", 1400000)
QC_MIN_READS_NANOPORE = config.get("qc_thresholds", {}).get("nanopore_min_reads", 5000)

# ==============================
# FUNCIONES DE VALIDACIÓN Y NORMALIZACIÓN
# ==============================

def validate_samplesinfo_integrity(csv_path="samplesinfo.csv"):
    """
    Validar integridad del archivo samplesinfo.csv
    - Verificar que existe
    - Comprobar columnas requeridas
    - Detectar IDs duplicados
    - Validar rutas de archivos
    """
    print("=== VALIDANDO INTEGRIDAD DE SAMPLESINFO.CSV ===")

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"ERROR: No se encuentra {csv_path}")

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        raise ValueError(f"ERROR: No se puede leer {csv_path}: {e}")

    # Verificar columnas requeridas
    required_cols = ['id', 'illumina_r1', 'illumina_r2', 'nanopore']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"ERROR: Columnas faltantes en CSV: {missing_cols}")

    # Verificar IDs duplicados
    duplicated_ids = df[df['id'].duplicated()]['id'].tolist()
    if duplicated_ids:
        raise ValueError(f"ERROR: IDs duplicados encontrados: {duplicated_ids}")

    # Verificar IDs vacíos o nulos
    empty_ids = df[df['id'].isna() | (df['id'] == '')].index.tolist()
    if empty_ids:
        raise ValueError(f"ERROR: IDs vacíos en filas: {[i+2 for i in empty_ids]}")  # +2 por header y 0-index

    print(f"VALIDACION OK: {len(df)} muestras, sin duplicados")
    return True

def normalize_filename(filepath, suffix_patterns):
    """
    Normalizar nombres de archivo eliminando sufijos de secuenciador
    IMPORTANTE: Esta es una función auxiliar, NO reemplaza el ID del CSV
    """
    if not filepath or pd.isna(filepath):
        return None

    filename = os.path.basename(str(filepath))
    normalized = filename

    for pattern in suffix_patterns:
        normalized = re.sub(pattern, '', normalized)

    return normalized

def sanitize_sample_name(name):
    """Sanitizar nombres de muestras para compatibilidad con Snakemake"""
    sanitized = re.sub(r'[^\w\-.]', '_', str(name))
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_')
    return sanitized

def is_valid_path(path_str):
    """Verificar si una ruta es válida (no NA, vacía o None)"""
    if pd.isna(path_str) or str(path_str).upper() in ['NA', 'NULL', '']:
        return False
    return True

def validate_file_paths(samples_dict):
    """Validar que los archivos especificados existen"""
    print("=== VALIDANDO RUTAS DE ARCHIVOS ===")

    missing_files = []
    for sample_id, data in samples_dict.items():
        for key, path in data.items():
            if key in ['illumina_r1', 'illumina_r2', 'nanopore'] and path:
                if not os.path.exists(path):
                    missing_files.append(f"{sample_id}: {key} -> {path}")

    if missing_files:
        print("ADVERTENCIA: Archivos faltantes:")
        for missing in missing_files:
            print(f"  - {missing}")
        print("Continuando con muestras que tienen archivos válidos...")
    else:
        print("VALIDACION OK: Todos los archivos especificados existen")

def load_samples():
    """
    Cargar y validar muestras desde samplesinfo.csv
    Implementa la jerarquía: CSV id > normalización (fallback)
    """
    # Validar integridad del CSV primero
    validate_samplesinfo_integrity()

    df = pd.read_csv("samplesinfo.csv")
    suffix_patterns = config.get("illumina_suffix_patterns", [])

    # Excluir controles
    exclude_tags = ['negative', 'control', 'blanco', 'unclassified', 'negativo']
    original_count = len(df)
    df = df[~df['id'].str.lower().str.contains('|'.join(exclude_tags), na=False)]
    excluded_count = original_count - len(df)

    if excluded_count > 0:
        print(f"EXCLUIDAS {excluded_count} muestras de control")

    samples_dict = {}
    sample_types = {}

    for _, row in df.iterrows():
        # FUENTE DE VERDAD: El ID del CSV es el nombre final de la muestra
        sample_id = sanitize_sample_name(row['id'])

        # Verificar qué datos tiene cada muestra
        has_illumina = is_valid_path(row['illumina_r1']) and is_valid_path(row['illumina_r2'])
        has_nanopore = is_valid_path(row['nanopore'])

        if not has_illumina and not has_nanopore:
            print(f"SKIP {sample_id}: Sin datos válidos")
            continue

        # Clasificar tipo de muestra
        if has_illumina and has_nanopore:
            sample_type = "hybrid"
        elif has_illumina:
            sample_type = "illumina"
        elif has_nanopore:
            sample_type = "nanopore"

        # Guardar información de la muestra
        # IMPORTANTE: Usar las rutas exactas del CSV, NO nombres normalizados
        samples_dict[sample_id] = {
            'illumina_r1': row['illumina_r1'] if has_illumina else None,
            'illumina_r2': row['illumina_r2'] if has_illumina else None,
            'nanopore': row['nanopore'] if has_nanopore else None,
            'dorado_model': row.get('dorado_model', 'unknown'),
            'type': sample_type,
            'original_id': row['id']  # Preservar ID original para trazabilidad
        }
        sample_types[sample_id] = sample_type

        # Mostrar mapeo de archivos (para debugging)
        print(f"ADDED {sample_id} ({sample_type.upper()})")
        if has_illumina:
            print(f"  Illumina R1: {row['illumina_r1']}")
            print(f"  Illumina R2: {row['illumina_r2']}")
        if has_nanopore:
            print(f"  Nanopore: {row['nanopore']}")

    # Validar que los archivos existen
    validate_file_paths(samples_dict)

    return samples_dict, sample_types

def get_multiqc_inputs(ILLUMINA_SAMPLES, NANOPORE_SAMPLES):
    """Obtener todos los archivos para MultiQC según las muestras disponibles"""
    inputs = []

    # FastQC Illumina - usar extend para aplanar las listas
    if ILLUMINA_SAMPLES:
        inputs.extend(expand("output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R{read}_fastqc.zip",
                           sample=ILLUMINA_SAMPLES, read=[1,2]))
        inputs.extend(expand("output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R{read}_fastqc.zip",
                           sample=ILLUMINA_SAMPLES, read=[1,2]))
        inputs.extend(expand("output/01_data/01.2_filtered/fastp/{sample}_fastp.json",
                           sample=ILLUMINA_SAMPLES))
        inputs.extend(expand("output/01_data/01.4_kraken2/illumina/{sample}_report.txt",
                           sample=ILLUMINA_SAMPLES))

    # NanoPlot Nanopore - usar extend para aplanar las listas
    if NANOPORE_SAMPLES:
        inputs.extend(expand("output/01_data/01.1_fastq_raw_qc/nanoplot/{sample}/NanoStats.txt",
                           sample=NANOPORE_SAMPLES))
        inputs.extend(expand("output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}/NanoStats.txt",
                           sample=NANOPORE_SAMPLES))
        inputs.extend(expand("output/01_data/01.4_kraken2/nanopore/{sample}_report.txt",
                           sample=NANOPORE_SAMPLES))

    return inputs

def get_multiqc_inputs_with_evaluations(ILLUMINA_SAMPLES, NANOPORE_SAMPLES):
    """Obtener todos los archivos para MultiQC incluyendo evaluaciones de ensamblajes"""
    inputs = get_multiqc_inputs(ILLUMINA_SAMPLES, NANOPORE_SAMPLES)

    # Añadir reportes de validación y métricas
    inputs.extend([
        "output/01_data/01.5_samples_pass/validation_report.txt",
        "output/01_data/01.5_samples_pass/samples_metrics_table.tsv"
    ])

    # Si existe el checkpoint, añadir evaluaciones solo si existen
    checkpoint_file = "output/01_data/01.5_samples_pass/passed_samples.txt"
    if os.path.exists(checkpoint_file):
        # Añadir resúmenes de evaluaciones para MultiQC solo si existen
        evaluation_files = [
            "output/02_assembly/02.4_evaluation/quast_summary.tsv",
            "output/02_assembly/02.4_evaluation/busco_summary.tsv",
            "output/02_assembly/02.4_evaluation/checkm2_summary.tsv"
        ]

        # Solo añadir archivos que realmente existen
        for eval_file in evaluation_files:
            if os.path.exists(eval_file):
                inputs.append(eval_file)

    return inputs

# ==============================
# FUNCIONES PARA CHECKPOINT DINÁMICO
# ==============================

def get_validated_samples(wildcards):
    """
    Función crítica para leer las muestras validadas desde el checkpoint.
    Esta función lee el archivo creado por el checkpoint samples_validation
    y devuelve la lista de muestras que han pasado el control de calidad.
    """
    try:
        # Obtener el output del checkpoint
        checkpoint_output = checkpoints.samples_validation.get(**wildcards).output[0]

        # Leer el archivo de muestras validadas
        passed_samples_file = "output/01_data/01.5_samples_pass/passed_samples.txt"

        if not os.path.exists(passed_samples_file):
            print(f"ADVERTENCIA: No se encontró el archivo {passed_samples_file}")
            return []

        # Leer las muestras validadas
        with open(passed_samples_file, 'r') as f:
            validated_samples = [line.strip() for line in f if line.strip()]

        print(f"CHECKPOINT: Encontradas {len(validated_samples)} muestras validadas")
        return validated_samples

    except Exception as e:
        print(f"ERROR en get_validated_samples: {e}")
        return []

def get_validated_samples_by_strategy(wildcards, strategy):
    """
    Obtener muestras validadas filtradas por estrategia de ensamblado
    strategy: 'hybrid', 'short_reads', 'long_reads'
    """
    try:
        # Obtener el output del checkpoint usando la forma correcta
        checkpoint_output = checkpoints.samples_validation.get().output[0]

        # Mapear estrategia a archivo
        strategy_files = {
            'hybrid': 'hybrid_samples.txt',
            'short_reads': 'short_reads_samples.txt',
            'long_reads': 'long_reads_samples.txt'
        }

        if strategy not in strategy_files:
            print(f"ERROR: Estrategia {strategy} no reconocida")
            return []

        strategy_file = f"output/01_data/01.5_samples_pass/{strategy_files[strategy]}"

        if not os.path.exists(strategy_file):
            print(f"ADVERTENCIA: No se encontró el archivo {strategy_file}")
            return []

        # Leer las muestras de esta estrategia
        with open(strategy_file, 'r') as f:
            strategy_samples = [line.strip() for line in f if line.strip()]

        print(f"CHECKPOINT: Encontradas {len(strategy_samples)} muestras para estrategia {strategy}")
        return strategy_samples

    except Exception as e:
        print(f"ERROR en get_validated_samples_by_strategy: {e}")
        return []

# Función para obtener las salidas específicas del ensamblado
def get_assembly_outputs(wildcards):
    """
    Función principal para generar las salidas de ensamblado basadas en el checkpoint.
    Esta función lee las muestras validadas y genera los archivos de salida correspondientes.
    """
    try:
        # Obtener las muestras validadas desde el checkpoint
        validated_samples = get_validated_samples(wildcards)

        if not validated_samples:
            print("CHECKPOINT: No hay muestras validadas para ensamblado")
            return []

        # Generar salidas para todas las muestras validadas
        assembly_outputs = expand("output/02_assembly/02.3_consensus/{sample}.fasta",
                                sample=validated_samples)

        print(f"CHECKPOINT: Generando {len(assembly_outputs)} salidas de ensamblado")
        return assembly_outputs

    except Exception as e:
        print(f"ERROR en get_assembly_outputs: {e}")
        return []

# Función para compatibilidad con código existente
def get_valid_samples(wildcards):
    """
    Función de compatibilidad que utiliza el nuevo sistema de checkpoint.
    Redirige a get_validated_samples para mantener compatibilidad.
    """
    return get_validated_samples(wildcards)

# Function removed - duplicated above


# Función para obtener recursos de la configuración de manera segura
def get_resource(rule, resource, default_value=None):
    """Obtener de manera segura un recurso desde la configuración"""
    try:
        return config["resources"][rule][resource]
    except KeyError:
        try:
            return config["resources"]["default"][resource]
        except KeyError:
            return default_value