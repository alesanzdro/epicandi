# ==============================
# REGLAS DE SETUP
# ==============================

rule setup_quast:
    output:
        "resources/quast/.setup_done"
    conda: "../envs/epicandi.yml"
    log:
        "output/logs/00_setup/quast_setup.log"
    shell:
        """
        mkdir -p resources/quast
        mkdir -p $(dirname {log})

        # Download SILVA database for rRNA gene finding
        echo "Setting up QUAST databases..." > {log} 2>&1
        quast-download-silva >> {log} 2>&1 || echo "SILVA download failed, but QUAST will still work" >> {log} 2>&1

        # Download BUSCO database for gene finding
        quast-download-busco >> {log} 2>&1 || echo "BUSCO database download failed, but QUAST will still work" >> {log} 2>&1

        touch {output}
        echo "QUAST setup completed" >> {log} 2>&1
        """


rule setup_checkm2:
    output:
        "resources/checkm2/.setup_done"
    conda: "../envs/epicandi.yml"
    log:
        "output/logs/00_setup/checkm2_setup.log"
    shell:
        """
        mkdir -p resources/checkm2
        mkdir -p $(dirname {log})

        # Download CheckM2 database
        checkm2 database --download --path resources/checkm2/ > {log} 2>&1
        touch {output}
        """

# ==============================
# REGLAS DE EVALUACIÓN DE ENSAMBLAJES
# ==============================

rule quast_evaluation:
    input:
        assembly="output/02_assembly/02.3_consensus/{sample}.fasta",
        setup_flag="resources/quast/.setup_done"
    output:
        html="output/02_assembly/02.4_evaluation/quast/{sample}/report.html",
        tsv="output/02_assembly/02.4_evaluation/quast/{sample}/report.tsv"
    params:
        outdir="output/02_assembly/02.4_evaluation/quast/{sample}",
        extra_params=config.get("quast_params", "")
    threads: config.get("resources", {}).get("quast", {}).get("threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("quast", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("quast", {}).get("walltime", "2:00:00")
    log:
        "output/logs/02_assembly/02.4_evaluation/quast/{sample}.log"
    benchmark:
        "output/logs/02_assembly/02.4_evaluation/quast/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        quast.py {input.assembly} \\
            -o {params.outdir} \\
            --eukaryote \\
            {params.extra_params} \\
            -t {threads} \\
            -l {wildcards.sample} > {log} 2>&1
        """

rule busco_evaluation:
    input:
        assembly="output/02_assembly/02.3_consensus/{sample}.fasta",
        setup_flag="resources/busco_downloads/.setup_done"
    output:
        summary="output/02_assembly/02.4_evaluation/busco/{sample}/short_summary.txt"
    params:
        outdir="output/02_assembly/02.4_evaluation/busco/{sample}",
        lineage=config.get("busco_lineage", "saccharomycetes_odb12")
    threads: config.get("resources", {}).get("busco", {}).get("threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("busco", {}).get("mem", 16000),
        walltime=config.get("resources", {}).get("busco", {}).get("walltime", "4:00:00")
    log:
        "output/logs/02_assembly/02.4_evaluation/busco/{sample}.log"
    benchmark:
        "output/logs/02_assembly/02.4_evaluation/busco/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}

        # Set BUSCO config and download path environment variables
        export BUSCO_CONFIG_FILE=resources/busco_downloads/config.ini
        export BUSCO_DATA_PATH=$(pwd)/resources/busco_downloads

        busco -i {input.assembly} \\
              -o {wildcards.sample} \\
              --out_path {params.outdir} \\
              -l {params.lineage} \\
              -m genome \\
              --force \\
              --cpu {threads} > {log} 2>&1

        # Move summary to expected location
        mv {params.outdir}/{wildcards.sample}/short_summary.*.txt {output.summary}
        """

rule checkm2_all_assemblies:
    input:
        assemblies=lambda wildcards: expand("output/02_assembly/02.3_consensus/{sample}.fasta",
                                          sample=get_validated_samples(wildcards)),
        checkm2_db=config.get("checkm2_db", "resources/checkm2/CheckM2_database/uniref100.KO.1.dmnd")
    output:
        summary="output/02_assembly/02.4_evaluation/checkm2_all/quality_report.tsv"
    params:
        outdir="output/02_assembly/02.4_evaluation/checkm2_all",
        input_dir="output/02_assembly/02.3_consensus"
    threads: config["resources"]["checkm2"]["threads"]
    resources:
        mem_mb=config["resources"]["checkm2"]["mem"],
        walltime=config["resources"]["checkm2"]["walltime"]
    log:
        "output/logs/02_assembly/02.4_evaluation/checkm2_all.log"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        rm -rf {params.outdir}

        # Export database path from config/input (already absolute) y mostrar
        export CHECKM2DB=$(readlink -f {input.checkm2_db})
        echo "Using CHECKM2DB=$CHECKM2DB" >> {log}

        checkm2 predict \\
            --threads {threads} \\
            --force \\
            --extension fasta \\
            --input {params.input_dir} \\
            --output-directory {params.outdir} 2>&1 >> {log}
        """

def get_evaluation_inputs(wildcards):
    """Función para obtener inputs de evaluación basados en muestras validadas"""
    checkpoint_output = checkpoints.samples_validation.get(**wildcards).output[0]

    with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
        validated_samples = [line.strip() for line in f if line.strip()]

    inputs = []
    inputs.extend(expand("output/02_assembly/02.4_evaluation/quast/{sample}/report.tsv", sample=validated_samples))
    inputs.extend(expand("output/02_assembly/02.4_evaluation/busco/{sample}/short_summary.txt", sample=validated_samples))
    inputs.append("output/02_assembly/02.4_evaluation/checkm2_all/quality_report.tsv")
    return inputs

rule collect_reports:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        evaluations=get_evaluation_inputs
    output:
        summary="output/02_assembly/02.4_evaluation/evaluation_summary.txt"
    log:
        "output/logs/02_assembly/02.4_evaluation/collect_reports.log"
    run:
        shell("mkdir -p $(dirname {log})")

        with open(log[0], "w") as log_file, open(output.summary, "w") as summary_file:
            # Count validated samples
            with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
                validated_samples = [line.strip() for line in f if line.strip()]

            log_file.write("=== RECOLECTANDO REPORTES DE EVALUACIÓN ===\\n")
            summary_file.write("=== EVALUACIONES DE ENSAMBLAJES COMPLETADAS ===\\n")
            summary_file.write(f"Total de muestras evaluadas: {len(validated_samples)}\\n")
            summary_file.write(f"Herramientas: QUAST, BUSCO, CheckM2\\n\\n")

            for sample in validated_samples:
                summary_file.write(f"{sample}: Evaluación completa\\n")
                log_file.write(f"Evaluación completada para {sample}\\n")

            log_file.write("=== RECOLECCIÓN COMPLETADA ===\\n")

# Regla agregada para toda la evaluación incluyendo CheckM2
rule evaluation_all:
    input:
        collect_reports="output/02_assembly/02.4_evaluation/evaluation_summary.txt",
        checkm2_all="output/02_assembly/02.4_evaluation/checkm2_all/quality_report.tsv"
    output:
        touch("output/02_assembly/02.4_evaluation/.evaluation_complete")
    log:
        "output/logs/02_assembly/02.4_evaluation/evaluation_all.log"
    shell:
        """
        echo "=== EVALUACIÓN COMPLETA ===" > {log}
        echo "QUAST, BUSCO y CheckM2 completados" >> {log}
        echo "$(date)" >> {log}
        """