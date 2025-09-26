# ==============================
# REGLAS DE ENSAMBLADO
# ==============================

# Ensamblado con lecturas largas (Flye + Medaka)
rule flye:
    input:
        reads="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    output:
        assembly="output/02_assembly/02.1_draft_assembly/flye/{sample}/assembly.fasta",
        info="output/02_assembly/02.1_draft_assembly/flye/{sample}/assembly_info.txt"
    params:
        outdir="output/02_assembly/02.1_draft_assembly/flye/{sample}",
        genome_size=config.get("assembly", {}).get("flye", {}).get("genome_size", "13m"),
        extra=config.get("assembly", {}).get("flye", {}).get("params", "--min-overlap 1000")
    threads: get_resource("flye", "threads", 16)
    resources:
        mem_mb=get_resource("flye", "mem", 32000),
        walltime=get_resource("flye", "walltime", "24:00:00")
    log:
        "output/logs/02_assembly/02.1_draft_assembly/flye/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        
        flye --nano-raw {input.reads} \
             --out-dir {params.outdir} \
             --genome-size {params.genome_size} \
             --threads {threads} \
             {params.extra} > {log} 2>&1
        """

rule medaka:
    input:
        assembly="output/02_assembly/02.1_draft_assembly/flye/{sample}/assembly.fasta",
        reads="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    output:
        polished="output/02_assembly/02.3_consensus/{sample}.fasta"
    params:
        outdir="output/02_assembly/02.2_polish/medaka/{sample}",
        model=config.get("assembly", {}).get("medaka", {}).get("model", "r1041_e82_400bps_sup_v4.2.0")
    threads: config.get("resources", {}).get("medaka", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("medaka", {}).get("mem", 32000),
        walltime=config.get("resources", {}).get("medaka", {}).get("walltime", "12:00:00")
    log:
        "output/logs/02_assembly/02.2_polish/medaka/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        mkdir -p $(dirname {output.polished})

        medaka_consensus -i {input.reads} \
                        -d {input.assembly} \
                        -o {params.outdir} \
                        -t {threads} \
                        -m {params.model} > {log} 2>&1

        # Copiar resultado a ubicación final
        cp {params.outdir}/consensus.fasta {output.polished}
        """

# Ensamblado con lecturas cortas (SPAdes + Pilon)
rule spades:
    input:
        r1="output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz",
        r2="output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz"
    output:
        contigs="output/02_assembly/02.1_draft_assembly/spades/{sample}/contigs.fasta",
        scaffolds="output/02_assembly/02.1_draft_assembly/spades/{sample}/scaffolds.fasta"
    params:
        outdir="output/02_assembly/02.1_draft_assembly/spades/{sample}",
        extra=config.get("assembly", {}).get("spades", {}).get("params", "--careful")
    threads: config.get("resources", {}).get("spades", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("spades", {}).get("mem", 64000),
        walltime=config.get("resources", {}).get("spades", {}).get("walltime", "24:00:00")
    log:
        "output/logs/02_assembly/02.1_draft_assembly/spades/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        
        spades.py -1 {input.r1} \
                 -2 {input.r2} \
                 -o {params.outdir} \
                 --threads {threads} \
                 {params.extra} > {log} 2>&1
        """

rule pilon:
    input:
        assembly="output/02_assembly/02.1_draft_assembly/spades/{sample}/scaffolds.fasta",
        r1="output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz",
        r2="output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz"
    output:
        polished="output/02_assembly/02.3_consensus/{sample}.fasta"
    params:
        outdir="output/02_assembly/02.2_polish/pilon/{sample}",
        prefix="pilon",
        extra=config.get("assembly", {}).get("pilon", {}).get("params", "--changes --fix all"),
        polishers=config.get("assembly", {}).get("pilon", {}).get("polishers", 3)
    threads: config.get("resources", {}).get("pilon", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("pilon", {}).get("mem", 32000),
        walltime=config.get("resources", {}).get("pilon", {}).get("walltime", "24:00:00")
    log:
        "output/logs/02_assembly/02.2_polish/pilon/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        
        # Alineamos lecturas al ensamblado
        bwa index {input.assembly}
        bwa mem -t {threads} {input.assembly} {input.r1} {input.r2} | samtools sort -@ {threads} -o {params.outdir}/aligned.bam -
        samtools index {params.outdir}/aligned.bam
        
        # Ejecutamos Pilon
        pilon --genome {input.assembly} \
              --bam {params.outdir}/aligned.bam \
              --output {params.prefix} \
              --outdir {params.outdir} \
              --threads {threads} \
              {params.extra} > {log} 2>&1

        # Copiar resultado a ubicación final
        mkdir -p $(dirname {output.polished})
        cp {params.outdir}/{params.prefix}.fasta {output.polished}
        """

# Ensamblado híbrido (Flye + Medaka + Pypolca)
rule hybrid_flye:
    input:
        reads="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    output:
        assembly="output/02_assembly/02.1_draft_assembly/hybrid_flye/{sample}/assembly.fasta",
        info="output/02_assembly/02.1_draft_assembly/hybrid_flye/{sample}/assembly_info.txt"
    params:
        outdir="output/02_assembly/02.1_draft_assembly/hybrid_flye/{sample}",
        genome_size=config.get("assembly", {}).get("flye", {}).get("genome_size", "13m"),
        extra=config.get("assembly", {}).get("flye", {}).get("params", "--min-overlap 1000")
    threads: config.get("resources", {}).get("flye", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("flye", {}).get("mem", 32000),
        walltime=config.get("resources", {}).get("flye", {}).get("walltime", "24:00:00")
    log:
        "output/logs/02_assembly/02.1_draft_assembly/hybrid_flye/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        
        flye --nano-raw {input.reads} \
             --out-dir {params.outdir} \
             --genome-size {params.genome_size} \
             --threads {threads} \
             {params.extra} > {log} 2>&1
        """

rule hybrid_medaka:
    input:
        assembly="output/02_assembly/02.1_draft_assembly/hybrid_flye/{sample}/assembly.fasta",
        reads="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    output:
        polished="output/02_assembly/02.2_polish/hybrid_medaka/{sample}/consensus.fasta"
    params:
        outdir="output/02_assembly/02.2_polish/hybrid_medaka/{sample}",
        model=config.get("assembly", {}).get("medaka", {}).get("model", "r1041_e82_400bps_sup_v4.2.0")
    threads: config.get("resources", {}).get("medaka", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("medaka", {}).get("mem", 32000),
        walltime=config.get("resources", {}).get("medaka", {}).get("walltime", "12:00:00")
    log:
        "output/logs/02_assembly/02.2_polish/hybrid_medaka/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        
        medaka_consensus -i {input.reads} \
                        -d {input.assembly} \
                        -o {params.outdir} \
                        -t {threads} \
                        -m {params.model} > {log} 2>&1
        """

rule pypolca:
    input:
        assembly="output/02_assembly/02.2_polish/hybrid_medaka/{sample}/consensus.fasta",
        r1="output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz",
        r2="output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz"
    output:
        polished="output/02_assembly/02.3_consensus/{sample}.fasta"
    params:
        outdir="output/02_assembly/02.2_polish/hybrid_pypolca/{sample}",
        prefix="{sample}",
        extra=config.get("assembly", {}).get("pypolca", {}).get("params", "--min_alt 4 --min_ratio 2")
    threads: config.get("resources", {}).get("pypolca", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("pypolca", {}).get("mem", 32000),
        walltime=config.get("resources", {}).get("pypolca", {}).get("walltime", "12:00:00")
    log:
        "output/logs/02_assembly/02.2_polish/hybrid_pypolca/{sample}.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}

        # Ejecutamos pypolca con rutas absolutas
        pypolca run -a {input.assembly} \
                   -1 {input.r1} \
                   -2 {input.r2} \
                   -o {params.outdir} \
                   -t {threads} \
                   --force \
                   {params.extra} \
                   -p {params.prefix} > {log} 2>&1

        # Verificamos que el output existe y lo copiamos a la ubicación final
        mkdir -p $(dirname {output.polished})
        if [ -f {params.outdir}/{params.prefix}_corrected.fasta ]; then
            cp {params.outdir}/{params.prefix}_corrected.fasta {output.polished}
        elif [ -f {params.outdir}/pypolca.fasta ]; then
            cp {params.outdir}/pypolca.fasta {output.polished}
        else
            echo "Error: No se encontró el archivo de salida de pypolca"
            exit 1
        fi
        """

# Establecer prioridades para resolver ambigüedades de reglas
ruleorder: pypolca > pilon > medaka

# ==============================
# REGLAS AGREGADAS DINÁMICAS USANDO CHECKPOINT
# ==============================

def get_hybrid_assembly_inputs(wildcards):
    """Función de entrada para regla aggregate de ensamblados híbridos"""
    hybrid_samples = []
    strategy_file = "output/01_data/01.5_samples_pass/hybrid_samples.txt"
    if os.path.exists(strategy_file):
        with open(strategy_file, 'r') as f:
            hybrid_samples = [line.strip() for line in f if line.strip()]
    return expand("output/02_assembly/02.3_consensus/{sample}.fasta", sample=hybrid_samples)

def get_short_reads_assembly_inputs(wildcards):
    """Función de entrada para regla aggregate de ensamblados de lecturas cortas"""
    sr_samples = []
    strategy_file = "output/01_data/01.5_samples_pass/short_reads_samples.txt"
    if os.path.exists(strategy_file):
        with open(strategy_file, 'r') as f:
            sr_samples = [line.strip() for line in f if line.strip()]
    return expand("output/02_assembly/02.3_consensus/{sample}.fasta", sample=sr_samples)

def get_long_reads_assembly_inputs(wildcards):
    """Función de entrada para regla aggregate de ensamblados de lecturas largas"""
    lr_samples = []
    strategy_file = "output/01_data/01.5_samples_pass/long_reads_samples.txt"
    if os.path.exists(strategy_file):
        with open(strategy_file, 'r') as f:
            lr_samples = [line.strip() for line in f if line.strip()]
    return expand("output/02_assembly/02.3_consensus/{sample}.fasta", sample=lr_samples)

# Reglas agregadas que utilizan el checkpoint para determinar qué muestras ensamblar
rule assemblies_hybrid:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        assemblies=get_hybrid_assembly_inputs
    output:
        summary="output/02_assembly/hybrid_assemblies_summary.txt"
    log:
        "output/logs/02_assembly/hybrid_assemblies_summary.log"
    run:
        shell("mkdir -p $(dirname {log})")

        with open(log[0], "w") as log_file, open(output.summary, "w") as summary_file:
            # Leer muestras híbridas directamente del archivo
            hybrid_samples = []
            strategy_file = "output/01_data/01.5_samples_pass/hybrid_samples.txt"
            if os.path.exists(strategy_file):
                with open(strategy_file, 'r') as f:
                    hybrid_samples = [line.strip() for line in f if line.strip()]

            log_file.write(f"=== RESUMEN DE ENSAMBLADOS HÍBRIDOS ===\n")
            log_file.write(f"Total de muestras híbridas: {len(hybrid_samples)}\n")

            summary_file.write(f"=== ENSAMBLADOS HÍBRIDOS COMPLETADOS ===\n")
            summary_file.write(f"Fecha: $(date)\n")
            summary_file.write(f"Total de muestras: {len(hybrid_samples)}\n\n")

            for sample in sorted(hybrid_samples):
                log_file.write(f"Ensamblado híbrido completado: {sample}\n")
                summary_file.write(f"{sample}: Flye + Medaka + Pypolca\n")

rule assemblies_short_reads:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        assemblies=get_short_reads_assembly_inputs
    output:
        summary="output/02_assembly/short_reads_assemblies_summary.txt"
    log:
        "output/logs/02_assembly/short_reads_assemblies_summary.log"
    run:
        shell("mkdir -p $(dirname {log})")

        with open(log[0], "w") as log_file, open(output.summary, "w") as summary_file:
            # Leer muestras de lecturas cortas directamente del archivo
            sr_samples = []
            strategy_file = "output/01_data/01.5_samples_pass/short_reads_samples.txt"
            if os.path.exists(strategy_file):
                with open(strategy_file, 'r') as f:
                    sr_samples = [line.strip() for line in f if line.strip()]

            log_file.write(f"=== RESUMEN DE ENSAMBLADOS DE LECTURAS CORTAS ===\n")
            log_file.write(f"Total de muestras de lecturas cortas: {len(sr_samples)}\n")

            summary_file.write(f"=== ENSAMBLADOS DE LECTURAS CORTAS COMPLETADOS ===\n")
            summary_file.write(f"Fecha: $(date)\n")
            summary_file.write(f"Total de muestras: {len(sr_samples)}\n\n")

            for sample in sorted(sr_samples):
                log_file.write(f"Ensamblado de lecturas cortas completado: {sample}\n")
                summary_file.write(f"{sample}: SPAdes + Pilon\n")

rule assemblies_long_reads:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        assemblies=get_long_reads_assembly_inputs
    output:
        summary="output/02_assembly/long_reads_assemblies_summary.txt"
    log:
        "output/logs/02_assembly/long_reads_assemblies_summary.log"
    run:
        shell("mkdir -p $(dirname {log})")

        with open(log[0], "w") as log_file, open(output.summary, "w") as summary_file:
            # Leer muestras de lecturas largas directamente del archivo
            lr_samples = []
            strategy_file = "output/01_data/01.5_samples_pass/long_reads_samples.txt"
            if os.path.exists(strategy_file):
                with open(strategy_file, 'r') as f:
                    lr_samples = [line.strip() for line in f if line.strip()]

            log_file.write(f"=== RESUMEN DE ENSAMBLADOS DE LECTURAS LARGAS ===\n")
            log_file.write(f"Total de muestras de lecturas largas: {len(lr_samples)}\n")

            summary_file.write(f"=== ENSAMBLADOS DE LECTURAS LARGAS COMPLETADOS ===\n")
            summary_file.write(f"Fecha: $(date)\n")
            summary_file.write(f"Total de muestras: {len(lr_samples)}\n\n")

            for sample in sorted(lr_samples):
                log_file.write(f"Ensamblado de lecturas largas completado: {sample}\n")
                summary_file.write(f"{sample}: Flye + Medaka\n")

# Regla maestra que consolida todos los ensamblados
rule all_assemblies:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        hybrid="output/02_assembly/hybrid_assemblies_summary.txt",
        short_reads="output/02_assembly/short_reads_assemblies_summary.txt",
        long_reads="output/02_assembly/long_reads_assemblies_summary.txt"
    output:
        final_summary="output/02_assembly/all_assemblies_summary.txt"
    log:
        "output/logs/02_assembly/all_assemblies_summary.log"
    run:
        shell("mkdir -p $(dirname {log})")

        with open(log[0], "w") as log_file, open(output.final_summary, "w") as summary_file:
            # Leer todas las muestras validadas directamente de los archivos
            with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
                all_validated = [line.strip() for line in f if line.strip()]

            # Leer muestras por estrategia
            hybrid_samples = []
            if os.path.exists("output/01_data/01.5_samples_pass/hybrid_samples.txt"):
                with open("output/01_data/01.5_samples_pass/hybrid_samples.txt", 'r') as f:
                    hybrid_samples = [line.strip() for line in f if line.strip()]

            sr_samples = []
            if os.path.exists("output/01_data/01.5_samples_pass/short_reads_samples.txt"):
                with open("output/01_data/01.5_samples_pass/short_reads_samples.txt", 'r') as f:
                    sr_samples = [line.strip() for line in f if line.strip()]

            lr_samples = []
            if os.path.exists("output/01_data/01.5_samples_pass/long_reads_samples.txt"):
                with open("output/01_data/01.5_samples_pass/long_reads_samples.txt", 'r') as f:
                    lr_samples = [line.strip() for line in f if line.strip()]

            total_assemblies = len(hybrid_samples) + len(sr_samples) + len(lr_samples)

            log_file.write(f"=== RESUMEN FINAL DE TODOS LOS ENSAMBLADOS ===\n")
            log_file.write(f"Total de muestras validadas: {len(all_validated)}\n")
            log_file.write(f"Total de ensamblados completados: {total_assemblies}\n")
            log_file.write(f"Ensamblados híbridos: {len(hybrid_samples)}\n")
            log_file.write(f"Ensamblados lecturas cortas: {len(sr_samples)}\n")
            log_file.write(f"Ensamblados lecturas largas: {len(lr_samples)}\n")

            summary_file.write(f"=== RESUMEN FINAL DE ENSAMBLADOS ===\n")
            summary_file.write(f"Fecha: $(date)\n")
            summary_file.write(f"Pipeline: EpiCandi\n\n")
            summary_file.write(f"ESTADÍSTICAS:\n")
            summary_file.write(f"Total de muestras validadas: {len(all_validated)}\n")
            summary_file.write(f"Total de ensamblados completados: {total_assemblies}\n")
            summary_file.write(f"Ensamblados híbridos: {len(hybrid_samples)}\n")
            summary_file.write(f"Ensamblados lecturas cortas: {len(sr_samples)}\n")
            summary_file.write(f"Ensamblados lecturas largas: {len(lr_samples)}\n\n")

            summary_file.write(f"DETALLE POR ESTRATEGIA:\n")
            if hybrid_samples:
                summary_file.write(f"\nHíbridos ({len(hybrid_samples)}):\n")
                for sample in sorted(hybrid_samples):
                    summary_file.write(f"  {sample}: Flye + Medaka + Pypolca\n")

            if sr_samples:
                summary_file.write(f"\nLecturas cortas ({len(sr_samples)}):\n")
                for sample in sorted(sr_samples):
                    summary_file.write(f"  {sample}: SPAdes + Pilon\n")

            if lr_samples:
                summary_file.write(f"\nLecturas largas ({len(lr_samples)}):\n")
                for sample in sorted(lr_samples):
                    summary_file.write(f"  {sample}: Flye + Medaka\n")