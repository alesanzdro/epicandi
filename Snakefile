#!/usr/bin/env python3
"""
Pipeline Snakemake para ensamblaje de datos Nanopore
Incluye: QC, trimming, filtrado, ensamblaje, polishing, taxonomía y calidad
"""

import os
import glob
from pathlib import Path

# ==============================
# CONFIGURACIÓN
# ==============================

# Cargar configuración
configfile: "config.yaml"

# Parámetros principales
GENOME_SIZE = "13m"
THREADS = 8
MIN_LENGTH = 1000
MIN_MEAN_Q = 12

# Directorios
INPUT_DIR = "input"
OUTPUT_DIR = "output"

# CheckM2: parámetros por configuración
CHECKM2_DB = config.get("checkm2_db", "resources/CheckM2_database/uniref100.KO.1.dmnd")
# Resolver a ruta absoluta si es relativa (basado en la raíz del workflow)
if not os.path.isabs(CHECKM2_DB):
    CHECKM2_DB = str(Path(workflow.basedir) / CHECKM2_DB)
CHECKM2_EXTENSION = config.get("checkm2_extension", ".fasta")

# Modelo de Medaka
MEDAKA_MODEL = "r1041_e82_400bps_sup_v5.2.0"

# ==============================
# FUNCIONES AUXILIARES
# ==============================

def get_samples():
    """Obtener lista de muestras desde el directorio de entrada, excluyendo controles"""
    samples = []
    
    # Tags a excluir (controles y blancos)
    EXCLUDE_TAGS = [
        'negative', 'Negative', 'NEGATIVE',
        'control', 'Control', 'CONTROL', 
        'blanco', 'Blanco', 'BLANCO',
        'unclassified', 'Unclassified', 'UNCLASSIFIED',
        'negativo', 'Negativo', 'NEGATIVO'
    ]
    
    # Buscar en subdirectorios también
    patterns = [
        f"{INPUT_DIR}/*.fastq.gz", 
        f"{INPUT_DIR}/*.fq.gz",
        f"{INPUT_DIR}/*/*.fastq.gz", 
        f"{INPUT_DIR}/*/*.fq.gz"
    ]
    
    for pattern in patterns:
        for file_path in glob.glob(pattern):
            basename = os.path.basename(file_path)
            # Remover extensiones
            sample = basename.replace('.fastq.gz', '').replace('.fq.gz', '')
            
            # Excluir muestras de control
            if any(tag in sample for tag in EXCLUDE_TAGS):
                print(f"⚠️  Excluyendo muestra de control: {sample}")
                continue
            
            if sample not in samples:
                samples.append(sample)
                print(f"✓ Muestra añadida: {sample}")
    
    return samples

# Obtener muestras
SAMPLES = get_samples()

def get_sample_path(sample_name):
    """Encontrar la ruta completa de un archivo de muestra.

    Busca primero en input/ y en subdirectorios inmediatos (dos patrones):
      input/{sample}.fastq.gz
      input/{sample}.fq.gz
      input/*/{sample}.fastq.gz
      input/*/{sample}.fq.gz
    Devuelve la primera coincidencia absoluta. Lanza FileNotFoundError si no lo encuentra.
    """
    candidate_patterns = [
        f"{INPUT_DIR}/{sample_name}.fastq.gz",
        f"{INPUT_DIR}/{sample_name}.fq.gz",
        f"{INPUT_DIR}/*/{sample_name}.fastq.gz",
        f"{INPUT_DIR}/*/{sample_name}.fq.gz",
    ]
    for pattern in candidate_patterns:
        matches = glob.glob(pattern)
        if matches:
            # Preferir la coincidencia más corta (menos ambigua) si varias
            matches.sort(key=len)
            return matches[0]
    raise FileNotFoundError(f"No se encontró archivo FASTQ/FQ para la muestra: {sample_name}")

# ==============================
# REGLA PRINCIPAL
# ==============================

rule all:
    input:
        # Ejecutar NanoPlot RAW por muestra (no incluido en MultiQC)
        expand(
            OUTPUT_DIR + "/01_data/01.2_nanoplot_raw/{sample}/NanoPlot-report.html",
            sample=SAMPLES,
        ),

        # Ejecutar Krona plots por muestra (no incluido en MultiQC)
        expand(
            OUTPUT_DIR + "/03_taxonomy/03.2_krona/{sample}_krona.html",
            sample=SAMPLES,
        ),

        # Mantener el pipeline anidado a través de MultiQC
        OUTPUT_DIR + "/04_report/04.2_multiqc/multiqc_report.html",

        # Asegurar cobertura por muestra (BAMs ordenados e índices)
        expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.stats.txt", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.coverage.txt", sample=SAMPLES),

        # Asegurar salidas de mosdepth
        expand(OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.global.dist.txt", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.per-base.bed.gz", sample=SAMPLES),

        # Consolidado de métricas
        OUTPUT_DIR + "/02_assembly/02.8_coverage_metrics/all_samples_coverage.tsv",
        OUTPUT_DIR + "/02_assembly/02.8_coverage_metrics/coverage_multiqc.json",

# ==============================
# REGLAS DE CONFIGURACIÓN DE BASES DE DATOS
# ==============================

rule setup_krona:
    output:
        "resources/krona/taxonomy/taxonomy.tab"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        # Create the directory structure and symbolic link for Krona taxonomy
        mkdir -p resources/krona/taxonomy
        
        # Get conda environment path dynamically
        CONDA_PREFIX=$(python -c "import os; print(os.environ['CONDA_PREFIX'])")
        
        # Remove existing krona taxonomy directory in conda env if it exists
        rm -rf $CONDA_PREFIX/opt/krona/taxonomy
        
        # Create symbolic link from conda env to our resources directory
        ln -s $(pwd)/resources/krona/taxonomy $CONDA_PREFIX/opt/krona/taxonomy
        
        # Update taxonomy database using ktUpdateTaxonomy.sh
        ktUpdateTaxonomy.sh
        """

rule setup_checkm2:
    output:
        CHECKM2_DB
    conda:
        "envs/epicandi_checkm2.yml"
    shell:
        """
    # Descargar base de datos si no existe
    BASEDIR="{workflow.basedir}"
    mkdir -p "$BASEDIR"/resources
    checkm2 database --download --path "$BASEDIR"/resources || echo "CheckM2 DB: download step skipped/ok"

        # Validar que el archivo esperado exista
        if [ ! -f {output} ]; then
            echo "ERROR: No se encontró la DB esperada: {output}" >&2
            echo "Contenido de resources/CheckM2_database:" >&2
            ls -lah resources/CheckM2_database || true
            exit 1
        fi

    # Registrar la ubicación en la configuración interna (best-effort)
    checkm2 database --setdblocation $(readlink -f {output}) || true
        """

rule setup_quast:
    output:
        "resources/quast/.setup_done"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        # Create resources directory
        mkdir -p resources/quast
        
        # Download SILVA database for rRNA gene finding
        quast-download-silva
        
        # Download BUSCO database for assessment  
        quast-download-busco
        
        # Mark setup as complete
        touch {output}
        """

# -------------------- RULE FOR BUSCO SETUP (Robust Version) -------------------- #

rule setup_busco:
    """
    Downloads the specific BUSCO lineage dataset defined in the config file
    into a local project directory for reproducible, offline use.
    """
    output:
        # We track a specific file inside the downloaded dataset to ensure it's complete
        touch("resources/busco/lineages/{lineage}/dataset.cfg".format(lineage=config["busco_lineage"]))
    log:
        OUTPUT_DIR + "/logs/00_setup/busco.log"
    params:
        # Define a local path within our project for BUSCO downloads. Usamos lambda para evitar expansión de wildcards.
        download_path = lambda wildcards: f"{workflow.basedir}/resources/busco",
        lineage = config["busco_lineage"]
    conda:
        "envs/busco.yml"
    shell:
        """
        mkdir -p $(dirname {log}) {params.download_path}

        echo "Downloading BUSCO lineage '{params.lineage}' to '{params.download_path}'..." &> {log}
        
        # Use --download and --download_path to control exactly where the data goes
    busco --download {params.lineage} --download_path {params.download_path} &>> {log}

        # Verify that the expected output was created successfully
        if [ ! -f {output} ]; then
            echo "ERROR: BUSCO download failed. Expected file not found: {output}" >&2
            echo "Please check the log for details: {log}" >&2
            exit 1
        fi
        """

# ==============================
# REGLAS DE PROCESAMIENTO DE DATOS
# ==============================

rule copy_raw_data:
    input:
        lambda wildcards: get_sample_path(wildcards.sample)
    output:
        OUTPUT_DIR + "/01_data/01.1_fastq_raw/{sample}.fastq.gz"
    log:
        OUTPUT_DIR + "/logs/01.1_copy_raw_data/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/01.1_copy_raw_data/{sample}_benchmark.txt"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        walltime=config["resources"]["default"]["walltime"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/01_data/01.1_fastq_raw {OUTPUT_DIR}/logs/01.1_copy_raw_data
        cp {input} {output} 2>&1 | tee {log}
        """

rule nanoplot_raw:
    input:
        OUTPUT_DIR + "/01_data/01.1_fastq_raw/{sample}.fastq.gz"
    output:
        OUTPUT_DIR + "/01_data/01.2_nanoplot_raw/{sample}/NanoPlot-report.html"
    log:
        OUTPUT_DIR + "/logs/01.2_nanoplot_raw/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/01.2_nanoplot_raw/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/01_data/01.2_nanoplot_raw/{sample}"
    threads: config["resources"]["nanoplot"]["threads"]
    resources:
        mem_mb=config["resources"]["nanoplot"]["mem"],
        walltime=config["resources"]["nanoplot"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/01.2_nanoplot_raw
        NanoPlot --fastq {input} \
                 --threads {threads} \
                 --outdir {params.outdir} \
                 --plots dot \
                 --N50 2>&1 | tee {log} || echo "NanoPlot falló, continuando..." | tee -a {log}
        """

rule porechop:
    input:
        OUTPUT_DIR + "/01_data/01.1_fastq_raw/{sample}.fastq.gz"
    output:
        OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop.fastq.gz"
    log:
        OUTPUT_DIR + "/logs/01.3.1_porechop/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/01.3.1_porechop/{sample}_benchmark.txt"
    threads: config["resources"]["porechop"]["threads"]
    resources:
        mem_mb=config["resources"]["porechop"]["mem"],
        walltime=config["resources"]["porechop"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/01_data/01.3_fastq_filtered {OUTPUT_DIR}/logs/01.3.1_porechop
        porechop \
            -i {input} \
            -o {output} \
            --threads {threads} \
            --discard_middle 2>&1 | tee {log}
        """

rule filtlong:
    input:
        OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop.fastq.gz"
    output:
        OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop_filtlong.fastq.gz"
    log:
        OUTPUT_DIR + "/logs/01.3.2_filtlong/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/01.3.2_filtlong/{sample}_benchmark.txt"
    params:
        min_length = MIN_LENGTH,
        min_mean_q = MIN_MEAN_Q
    threads: config["resources"]["filtlong"]["threads"]
    resources:
        mem_mb=config["resources"]["filtlong"]["mem"],
        walltime=config["resources"]["filtlong"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/logs/01.3.2_filtlong
        filtlong \
            --min_length {params.min_length} \
            --min_mean_q {params.min_mean_q} \
            {input} 2> {log} | gzip > {output}
        """

rule nanoplot_filtered:
    input:
        OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop_filtlong.fastq.gz"
    output:
        OUTPUT_DIR + "/01_data/01.4_nanoplot_filtered/{sample}/NanoPlot-report.html"
    log:
        OUTPUT_DIR + "/logs/01.4_nanoplot_filtered/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/01.4_nanoplot_filtered/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/01_data/01.4_nanoplot_filtered/{sample}"
    threads: config["resources"]["nanoplot"]["threads"]
    resources:
        mem_mb=config["resources"]["nanoplot"]["mem"],
        walltime=config["resources"]["nanoplot"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/01.4_nanoplot_filtered
        NanoPlot --fastq {input} \
                 --threads {threads} \
                 --outdir {params.outdir} \
                 --plots dot \
                 --N50 2>&1 | tee {log} || echo "NanoPlot falló, continuando..." | tee -a {log}
        """

# ==============================
# COMPARATIVA NANOSTATS RAW vs TRIM
# ==============================

rule nanostats_trim_compare:
    input:
        raw_stats = expand(OUTPUT_DIR + "/01_data/01.2_nanoplot_raw/{sample}/NanoStats.txt", sample=SAMPLES),
        trim_stats = expand(OUTPUT_DIR + "/01_data/01.4_nanoplot_filtered/{sample}/NanoStats.txt", sample=SAMPLES)
    output:
        OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_comparison.tsv"
    log:
        OUTPUT_DIR + "/logs/04.1_collection_renamed/nanostats_trim_compare.log"
    params:
        collection_dir = OUTPUT_DIR + "/04_report/04.1_collection_renamed",
        samples = SAMPLES,
        output_dir = OUTPUT_DIR
    conda:
        "envs/epicandi.yml"
    script:
        "scripts/generate_nanostats_trim_compare.py"

# Genera los JSONs para los scatters RAW vs TRIM que consumirá MultiQC
rule nanostats_trim_scatters:
    input:
        tsv = OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_comparison.tsv"
    output:
        reads    = OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_reads.json",
        mean_len = OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_mean_len.json",
        mean_q   = OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_mean_q.json",
        n50      = OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_n50.json"
    log:
        OUTPUT_DIR + "/logs/04.1_collection_renamed/nanostats_trim_scatters.log"
    benchmark:
        OUTPUT_DIR + "/logs/04.1_collection_renamed/nanostats_trim_scatters_benchmark.txt"
    threads: 1
    resources:
        mem_mb=2000,
        walltime="00:10:00"
    conda:
        "envs/epicandi.yml"
    script:
        "scripts/generate_nanostats_trim_scatter.py"

# ==============================
# REGLAS DE TAXONOMÍA
# ==============================

rule kraken2_classification:
    input:
        OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop_filtlong.fastq.gz"
    output:
        report = OUTPUT_DIR + "/03_taxonomy/03.1_kraken2/{sample}_report.txt",
        output = OUTPUT_DIR + "/03_taxonomy/03.1_kraken2/{sample}_output.txt"
    log:
        OUTPUT_DIR + "/logs/03.1_kraken2_classification/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/03.1_kraken2_classification/{sample}_benchmark.txt"
    params:
        db = config["kraken_db"],
        confidence = config["kraken_conf"]
    threads: config["resources"]["kraken2"]["threads"]
    resources:
        mem_mb=config["resources"]["kraken2"]["mem"],
        walltime=config["resources"]["kraken2"]["walltime"],
        kraken_db=1  # Solo 1 proceso Kraken2 a la vez
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/03_taxonomy/03.1_kraken2 {OUTPUT_DIR}/logs/03.1_kraken2_classification
        kraken2 \\
            --db {params.db} \\
            --threads {threads} \\
            --confidence {params.confidence} \\
            --report {output.report} \\
            --output {output.output} \\
            --use-names \\
            {input} 2>&1 | tee {log}
        """

rule krona_plot:
    input:
        kraken_output = OUTPUT_DIR + "/03_taxonomy/03.1_kraken2/{sample}_output.txt",
        krona_db = "resources/krona/taxonomy/taxonomy.tab"
    output:
        OUTPUT_DIR + "/03_taxonomy/03.2_krona/{sample}_krona.html"
    log:
        OUTPUT_DIR + "/logs/03.2_krona/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/03.2_krona/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/03_taxonomy/03.2_krona"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        walltime=config["resources"]["default"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/03.2_krona
        ktImportTaxonomy -t 5 -m 3 -tax {workflow.basedir}/resources/krona/taxonomy -o {output} {input.kraken_output} 2>&1 | tee {log}
    """

# ==============================
# REGLAS DE ENSAMBLAJE
# ==============================

rule flye_assembly:
    input:
        OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop_filtlong.fastq.gz"
    output:
        OUTPUT_DIR + "/02_assembly/02.1_flye/{sample}/assembly.fasta"
    log:
        OUTPUT_DIR + "/logs/02.1_flye_assembly/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.1_flye_assembly/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/02_assembly/02.1_flye/{sample}",
        genome_size = GENOME_SIZE
    threads: config["resources"]["flye"]["threads"]
    resources:
        mem_mb=config["resources"]["flye"]["mem"],
        walltime=config["resources"]["flye"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/02.1_flye_assembly
        flye \
            --nano-hq {input} \
            --genome-size {params.genome_size} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --asm-coverage 80 \
            --iterations 3 \
            --scaffold 2>&1 | tee {log}
        """

rule medaka_consensus:
    input:
        reads = OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop_filtlong.fastq.gz",
        assembly = OUTPUT_DIR + "/02_assembly/02.1_flye/{sample}/assembly.fasta"
    output:
        OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}/consensus.fasta"
    log:
        OUTPUT_DIR + "/logs/02.2_medaka_consensus/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.2_medaka_consensus/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}",
        model = MEDAKA_MODEL
    threads: config["resources"]["medaka"]["threads"]
    resources:
        mem_mb=config["resources"]["medaka"]["mem"],
        walltime=config["resources"]["medaka"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/02.2_medaka_consensus
        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.outdir} \
            -t {threads} \
            -m {params.model} 2>&1 | tee {log}
        """

rule quast:
    input:
        assembly = OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}/consensus.fasta",
        setup_flag = "resources/quast/.setup_done"
    output:
        OUTPUT_DIR + "/02_assembly/02.3_quast/{sample}/report.html"
    log:
        OUTPUT_DIR + "/logs/02.3_quast/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.3_quast/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/02_assembly/02.3_quast/{sample}",
        extra_params = config.get("quast_params", "")
    threads: config["resources"]["quast"]["threads"]
    resources:
        mem_mb=config["resources"]["quast"]["mem"],
        walltime=config["resources"]["quast"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/02.3_quast
        quast.py {input.assembly} \
            -o {params.outdir} \
            --eukaryote \
            {params.extra_params} \
            -t {threads} \
            -l {wildcards.sample} 2>&1 | tee {log}
        """

# ==============================
# REGLAS DE EVALUACIÓN DE CALIDAD
# ==============================

rule checkm2_quality:
    input:
        fasta = OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}/consensus.fasta",
        checkm2_db = CHECKM2_DB
    output:
        OUTPUT_DIR + "/02_assembly/02.4_checkm2/{sample}/quality_report.tsv"
    log:
        OUTPUT_DIR + "/logs/02.4_checkm2/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.4_checkm2/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/02_assembly/02.4_checkm2/{sample}",
        input_dir = OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}"
    threads: config["resources"]["checkm2"]["threads"]
    resources:
        mem_mb=config["resources"]["checkm2"]["mem"],
        walltime=config["resources"]["checkm2"]["walltime"]
    conda:
        "envs/epicandi_checkm2.yml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/logs/02.4_checkm2
        
        # Remove output directory if it exists to avoid conflicts
        rm -rf {params.outdir}
        
    # Export database path from config/input (already absolute) y mostrar
    export CHECKM2DB=$(readlink -f {input.checkm2_db})
    echo "Using CHECKM2DB=$CHECKM2DB" | tee -a {log}
        
        # Run CheckM2 with the working configuration
        checkm2 predict \
            --threads {threads} \
            --force \
            --extension {CHECKM2_EXTENSION} \
            --input {params.input_dir} \
            --output-directory {params.outdir} 2>&1 | tee {log}
        """

# Target agregado para facilitar ejecutar solo CheckM2 en todas las muestras
rule checkm2_quality_all:
    input:
        expand(OUTPUT_DIR + "/02_assembly/02.4_checkm2/{sample}/quality_report.tsv", sample=SAMPLES)

rule prokka_annotation:
    input:
        OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}/consensus.fasta"
    output:
        gff = OUTPUT_DIR + "/02_assembly/02.5_prokka/{sample}/annotation.gff",
        faa = OUTPUT_DIR + "/02_assembly/02.5_prokka/{sample}/annotation.faa",
        ffn = OUTPUT_DIR + "/02_assembly/02.5_prokka/{sample}/annotation.ffn"
    log:
        OUTPUT_DIR + "/logs/02.5_prokka_annotation/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.5_prokka_annotation/{sample}_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/02_assembly/02.5_prokka/{sample}"
    threads: config["resources"]["prokka"]["threads"]
    resources:
        mem_mb=config["resources"]["prokka"]["mem"],
        walltime=config["resources"]["prokka"]["walltime"]
    conda:
        "envs/epicandi_annotation.yml"
    shell:
        """
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/02.5_prokka_annotation
        prokka \\
            --outdir {params.outdir} \\
            --prefix annotation \\
            --cpus {threads} \\
            --force \\
            {input} 2>&1 | tee {log}
        """

# ==============================
# REGLAS DE COBERTURA
# ==============================

rule alignment_minimap2:
    input:
        reads = OUTPUT_DIR + "/01_data/01.3_fastq_filtered/{sample}_porechop_filtlong.fastq.gz",
        assembly = OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}/consensus.fasta"
    output:
        bam = temp(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.unsorted.bam")
    log:
        OUTPUT_DIR + "/logs/02.6_coverage/{sample}_minimap2.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.6_coverage/{sample}_minimap2_benchmark.txt"
    threads: config["resources"]["minimap2"]["threads"] if "minimap2" in config["resources"] else 12
    resources:
        mem_mb=config["resources"]["minimap2"]["mem"] if "minimap2" in config["resources"] else 16000,
        walltime=config["resources"]["minimap2"]["walltime"] if "minimap2" in config["resources"] else "04:00:00"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        
        # Log versions and parameters
        echo "Running minimap2 alignment for {wildcards.sample}" | tee {log}
        minimap2 --version | tee -a {log}
        echo "Assembly: {input.assembly}" | tee -a {log}
        echo "Reads: {input.reads}" | tee -a {log}
        
        # Run alignment with ONT-specific parameters
        minimap2 -ax map-ont \
            -t {threads} \
            --secondary=no \
            -K 500M \
            {input.assembly} \
            {input.reads} 2>> {log} | \
        samtools view -b -@ 2 -F 4 -F 256 -F 2048 - > {output.bam} 2>> {log}
        
        # Report alignment statistics
        echo -e "\nAlignment completed. Quick stats:" | tee -a {log}
        samtools flagstat {output.bam} 2>/dev/null | tee -a {log}
        """

rule samtools_processing:
    input:
        bam = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.unsorted.bam"
    output:
        sorted_bam = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam",
        sorted_bai = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam.bai",
        stats = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.stats.txt",
        coverage = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.coverage.txt",
        flagstat = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.flagstat.txt",
        idxstats = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.idxstats.txt"
    log:
        OUTPUT_DIR + "/logs/02.6_coverage/{sample}_samtools.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.6_coverage/{sample}_samtools_benchmark.txt"
    threads: config["resources"]["samtools"]["threads"] if "samtools" in config["resources"] else 8
    resources:
        mem_mb=config["resources"]["samtools"]["mem"] if "samtools" in config["resources"] else 8000,
        walltime=config["resources"]["samtools"]["walltime"] if "samtools" in config["resources"] else "02:00:00"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        echo "Processing BAM file for {wildcards.sample}" | tee {log}
        
        # Step 1: Sort BAM file
        echo "Sorting BAM file..." | tee -a {log}
        samtools sort \
            -@ {threads} \
            -m 1G \
            -o {output.sorted_bam} \
            {input.bam} 2>> {log}
        
        # Step 2: Index sorted BAM
        echo "Indexing sorted BAM..." | tee -a {log}
        samtools index -@ {threads} {output.sorted_bam} 2>> {log}
        
        # Step 3: Generate comprehensive statistics
        echo "Generating statistics..." | tee -a {log}
        samtools stats -@ {threads} {output.sorted_bam} > {output.stats} 2>> {log}
        
        # Step 4: Generate alignment statistics
        echo "Generating flagstat..." | tee -a {log}
        samtools flagstat -@ {threads} {output.sorted_bam} > {output.flagstat} 2>> {log}
        
        # Step 5: Generate index statistics
        echo "Generating idxstats..." | tee -a {log}
        samtools idxstats {output.sorted_bam} > {output.idxstats} 2>> {log}
        
        # Step 6: Calculate per-contig coverage
        echo "Calculating coverage..." | tee -a {log}
        samtools coverage \
            --min-MQ 10 \
            --min-BQ 10 \
            {output.sorted_bam} > {output.coverage} 2>> {log}
        
        # Summary
        echo -e "\nProcessing complete for {wildcards.sample}" | tee -a {log}
        echo "Coverage summary:" | tee -a {log}
        tail -n +2 {output.coverage} | \
            awk '{{sum+=$7*$3; total+=$3}} END {{if(total>0) printf "Average coverage: %.2fX\\n", sum/total}}' | \
            tee -a {log}
        """

rule mosdepth_coverage:
    input:
        sorted_bam = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam",
        sorted_bai = OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam.bai"
    output:
        dist = OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.global.dist.txt",
        summary = OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt",
        per_base = OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.per-base.bed.gz",
        regions = OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.regions.bed.gz"
    log:
        OUTPUT_DIR + "/logs/02.7_mosdepth/{sample}_mosdepth.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.7_mosdepth/{sample}_mosdepth_benchmark.txt"
    params:
        prefix = OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}",
        window_size = 500
    threads: config["resources"]["mosdepth"]["threads"] if "mosdepth" in config["resources"] else 4
    resources:
        mem_mb=config["resources"]["mosdepth"]["mem"] if "mosdepth" in config["resources"] else 4000,
        walltime=config["resources"]["mosdepth"]["walltime"] if "mosdepth" in config["resources"] else "01:00:00"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {params.prefix}) $(dirname {log})
        
        echo "Running mosdepth for {wildcards.sample}" | tee {log}
        mosdepth --version | tee -a {log}
        
        # Run mosdepth with comprehensive options
        mosdepth \
            --threads {threads} \
            --fast-mode \
            --by {params.window_size} \
            --mapq 10 \
            {params.prefix} \
            {input.sorted_bam} 2>> {log}
        
        # Report summary
        echo -e "\nMosdepth summary:" | tee -a {log}
        cat {output.summary} | tee -a {log}
        """

rule extract_coverage_metrics:
    input:
        coverages = expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.coverage.txt", sample=SAMPLES),
        stats = expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.stats.txt", sample=SAMPLES),
        mos_summaries = expand(OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt", sample=SAMPLES)
    output:
        tsv = OUTPUT_DIR + "/02_assembly/02.8_coverage_metrics/all_samples_coverage.tsv",
        json = OUTPUT_DIR + "/02_assembly/02.8_coverage_metrics/coverage_multiqc.json"
    log:
        OUTPUT_DIR + "/logs/02.8_coverage_metrics/consolidate_metrics.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.8_coverage_metrics/consolidate_metrics_benchmark.txt"
    params:
        samples = SAMPLES,
        output_dir = OUTPUT_DIR
    threads: 1
    resources:
        mem_mb=2000,
        walltime="00:30:00"
    conda:
        "envs/epicandi.yml"
    script:
        "scripts/extract_coverage_metrics.py"

# -------------------- RULE FOR RUNNING BUSCO (Modified) -------------------- #

rule busco:
    """
    Runs BUSCO to assess genome assembly completeness against a specific lineage.
    """
    input:
        assembly = OUTPUT_DIR + "/02_assembly/02.2_medaka/{sample}/consensus.fasta",
        # This now depends on the specific dataset file, making it more robust
        setup_flag = "resources/busco/lineages/{lineage}/dataset.cfg".format(lineage=config["busco_lineage"])
    output:
        summary = OUTPUT_DIR + "/02_assembly/02.9_busco/{sample}/short_summary.txt"
    log:
        OUTPUT_DIR + "/logs/02.9_busco/{sample}.log"
    benchmark:
        OUTPUT_DIR + "/logs/02.9_busco/{sample}_benchmark.txt"
    params:
        # Directorio padre donde BUSCO creará la carpeta de la muestra (-o <sample>)
        outdir = OUTPUT_DIR + "/02_assembly/02.9_busco",
        lineage = config["busco_lineage"],
        download_path = lambda wildcards: f"{workflow.basedir}/resources/busco",
        extra = config.get("busco_params", "")
    threads: config["resources"]["busco"]["threads"]
    resources:
        mem_mb=config["resources"]["busco"]["mem"],
        walltime=config["resources"]["busco"]["walltime"]
    conda:
        "envs/busco.yml" # o la que uses para busco
    shell:
        """
    mkdir -p {params.outdir} {OUTPUT_DIR}/logs/02.9_busco

        echo "Running BUSCO for sample {wildcards.sample}" > {log}
        echo "Assembly: {input.assembly}" >> {log}
        echo "Lineage: {params.lineage}" >> {log}
    echo "Download path: {params.download_path}" >> {log}
    echo "Parent outdir: {params.outdir}" >> {log}

        # Ejecutar BUSCO: -o define el nombre y --out_path el directorio destino estable que controlamos.
        busco -i {input.assembly} \
            -o {wildcards.sample} \
            -l {params.lineage} \
            -m genome \
            --cpu {threads} \
            --download_path {params.download_path} \
            --out_path {params.outdir} \
            {params.extra} &>> {log}

        SPECIFIC_SUMMARY="{params.outdir}/{wildcards.sample}/short_summary.specific.{params.lineage}.{wildcards.sample}.txt"
        BASIC_SUMMARY="{params.outdir}/{wildcards.sample}/run_{params.lineage}/short_summary.txt"
        if [ -f "$SPECIFIC_SUMMARY" ]; then
            cp "$SPECIFIC_SUMMARY" {output.summary}
            echo "Copiado resumen específico -> {output.summary}" >> {log}
        elif [ -f "$BASIC_SUMMARY" ]; then
            cp "$BASIC_SUMMARY" {output.summary}
            echo "Copiado resumen básico run_{params.lineage} -> {output.summary}" >> {log}
        else
            echo "ERROR: No se encontró ni $SPECIFIC_SUMMARY ni $BASIC_SUMMARY" >> {log}
            ls -R {params.outdir}/{wildcards.sample} >> {log} 2>&1 || true
            exit 1
        fi
        """

# ==============================
# REGLAS DE RECOLECCIÓN Y REPORTE
# ==============================

rule collect_reports:
    input:
        # Asegurar dependencias explícitas para el grafo de reglas
        expand(OUTPUT_DIR + "/02_assembly/02.4_checkm2/{sample}/quality_report.tsv", sample=SAMPLES),
        expand(OUTPUT_DIR + "/03_taxonomy/03.1_kraken2/{sample}_report.txt", sample=SAMPLES),
        # NanoPlot filtrado (usamos el HTML como proxy; genera NanoStats.txt en el mismo dir)
        expand(OUTPUT_DIR + "/01_data/01.4_nanoplot_filtered/{sample}/NanoPlot-report.html", sample=SAMPLES),
        # Comparativa NanoStats trim vs raw
        OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_comparison.tsv",
        # Scatters de NanoStats (JSON) para MultiQC
        OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_reads.json",
        OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_mean_len.json",
        OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_mean_q.json",
        OUTPUT_DIR + "/04_report/04.1_collection_renamed/nanostats/nanostats_trim_n50.json",
        # Prokka
        expand(OUTPUT_DIR + "/02_assembly/02.5_prokka/{sample}/annotation.gff", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.3_quast/{sample}/report.html", sample=SAMPLES),
        # Cobertura e idxstats por muestra (para renombrado en colección)
        expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.coverage.txt", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.6_coverage/{sample}/{sample}.idxstats.txt", sample=SAMPLES),
        # Métricas consolidadas de cobertura (para copiar a colección)
        OUTPUT_DIR + "/02_assembly/02.8_coverage_metrics/coverage_multiqc.json",
        OUTPUT_DIR + "/02_assembly/02.8_coverage_metrics/all_samples_coverage.tsv",
        # Mosdepth
        expand(OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.global.dist.txt", sample=SAMPLES),
        expand(OUTPUT_DIR + "/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt", sample=SAMPLES),
        # BUSCO
        expand(OUTPUT_DIR + "/02_assembly/02.9_busco/{sample}/short_summary.txt", sample=SAMPLES)
    output:
        flag = OUTPUT_DIR + "/04_report/04.1_collection_renamed/collection_collected.flag"
    log:
        OUTPUT_DIR + "/logs/04.1_collection_renamed/collect_reports.log"
    benchmark:
        OUTPUT_DIR + "/logs/04.1_collection_renamed/collect_reports_benchmark.txt"
    params:
        outdir = OUTPUT_DIR + "/04_report/04.1_collection_renamed",
        samples = SAMPLES,
        output_dir = OUTPUT_DIR
    script:
        "scripts/collect_reports.py"

rule render_multiqc_config:
    input:
        template = "conf/multiqc/multiqc_config.template.yaml",
        config_yaml = "config.yaml",
        script = "scripts/render_multiqc_config.py"
    output:
        rendered = "conf/multiqc/multiqc_config.rendered.yaml"
    log:
        OUTPUT_DIR + "/logs/04.2_multiqc_report/render_multiqc_config.log"
    threads: 1
    resources:
        mem_mb=1000,
        walltime="00:05:00"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        python {input.script} 2>&1 | tee {log}
        test -s {output.rendered}
        """

###############################
# Resumen BUSCO para MultiQC (TSV)
###############################

rule busco_summary:
    """Genera un TSV resumido de BUSCO con cabecera comentada para MultiQC.
    Lee busco/*.busco_summary.txt ya copiados por collect_reports."""
    input:
        results = expand(OUTPUT_DIR + "/02_assembly/02.9_busco/{sample}/short_summary.txt", sample=SAMPLES),    
        collection_flag = OUTPUT_DIR + "/04_report/04.1_collection_renamed/collection_collected.flag",
        script = "scripts/generate_busco_metrics_json.py"
    output:
        tsv = OUTPUT_DIR + "/04_report/04.1_collection_renamed/busco/busco_summary_mqc.tsv"
    log:
        OUTPUT_DIR + "/logs/04.1_collection_renamed/busco_summary.log"
    threads: 1
    resources:
        mem_mb=1000,
        walltime="00:05:00"
    conda:
        "envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        python {input.script} --tsv 2>&1 | tee {log}
        test -s {output.tsv}
        """

rule multiqc_report:
    input:
        collection_flag = OUTPUT_DIR + "/04_report/04.1_collection_renamed/collection_collected.flag",
        busco_tsv = OUTPUT_DIR + "/04_report/04.1_collection_renamed/busco/busco_summary_mqc.tsv",
        rendered_config = "conf/multiqc/multiqc_config.rendered.yaml",
        custom_logo = "conf/multiqc/logo.png"
    output:
        OUTPUT_DIR + "/04_report/04.2_multiqc/multiqc_report.html"
    params:
        outdir = OUTPUT_DIR + "/04_report/04.2_multiqc",
        search_dirs_str = " ".join([
            OUTPUT_DIR + "/04_report/04.1_collection_renamed",
            OUTPUT_DIR + "/02_assembly/02.3_quast"
        ])
    log:
        OUTPUT_DIR + "/logs/04.2_multiqc_report/multiqc_report.log"
    benchmark:
        OUTPUT_DIR + "/logs/04.2_multiqc_report/multiqc_report_benchmark.txt"
    threads: config["resources"]["multiqc"]["threads"]
    resources:
        mem_mb=config["resources"]["multiqc"]["mem"],
        walltime=config["resources"]["multiqc"]["walltime"]
    conda:
        "envs/epicandi.yml"
    shell:
        r"""
        mkdir -p {params.outdir} {OUTPUT_DIR}/logs/04.2_multiqc_report
        
        # Generar el reporte MultiQC
        multiqc \
            --outdir {params.outdir} \
            --filename multiqc_report \
            --config {input.rendered_config} \
            --force \
            --verbose \
            --no-ai \
            --no-megaqc-upload \
            --ignore '*flagstat*' \
            {params.search_dirs_str} 2>&1 | tee {log}
        
        # Pequeña pausa para asegurar que el archivo está completamente escrito
        sleep 1
        
        sed -i '\|<a href="http://multiqc.info"|,\|</a>|d' {params.outdir}/multiqc_report.html
        sed -i '\|id="mqc_welcome"|,\|</div>|d' {params.outdir}/multiqc_report.html
        """

# ==============================
# REGLAS DE LIMPIEZA
# ==============================

rule clean:
    shell:
        """
        rm -rf {OUTPUT_DIR}/01_data/01.3_fastq_filtered/*_temp*
        """

rule clean_all:
    shell:
        """
        rm -rf {OUTPUT_DIR}
        """

# ==============================
# CONFIGURACIÓN ADICIONAL
# ==============================

# Configuración de recursos por defecto
def get_mem_mb(wildcards, attempt):
    return attempt * 8000

# Configuración de reintentos
__default__ = {
    "retries": 2,
    "resources": {
        "mem_mb": get_mem_mb
    }
}
