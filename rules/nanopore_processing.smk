# ==============================
# REGLAS NANOPORE
# ==============================

rule nanoplot_raw:
    input:
        fastq="output/01_data/00_staged_data/nanopore/{sample}.nanopore.fastq.gz",
        staging="output/01_data/00_staged_data/.staging_complete"
    output:
        html="output/01_data/01.1_fastq_raw_qc/nanoplot/{sample}/NanoPlot-report.html",
        stats="output/01_data/01.1_fastq_raw_qc/nanoplot/{sample}/NanoStats.txt"
    params:
        outdir="output/01_data/01.1_fastq_raw_qc/nanoplot/{sample}"
    threads: config.get("resources", {}).get("nanoplot", {}).get("threads", 6)
    resources:
        mem_mb=config.get("resources", {}).get("nanoplot", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("nanoplot", {}).get("walltime", "8:00:00")
    log:
        "output/logs/01.1_fastq_raw_qc/nanoplot/{sample}.log"
    benchmark:
        "output/logs/01.1_fastq_raw_qc/nanoplot/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}
        NanoPlot --fastq {input.fastq} --outdir {params.outdir} --threads {threads} > {log} 2>&1
        """

rule porechop:
    input:
        fastq="output/01_data/00_staged_data/nanopore/{sample}.nanopore.fastq.gz",
        staging="output/01_data/00_staged_data/.staging_complete"
    output:
        fastq="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop.fastq.gz"
    threads: config.get("resources", {}).get("porechop", {}).get("threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("porechop", {}).get("mem", 16000),
        walltime=config.get("resources", {}).get("porechop", {}).get("walltime", "12:00:00")
    log:
        "output/logs/01.2_filtered/porechop_filtlong/{sample}_porechop.log"
    benchmark:
        "output/logs/01.2_filtered/porechop_filtlong/{sample}_porechop_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        porechop -i {input.fastq} -o {output.fastq} --discard_middle --threads {threads} > {log} 2>&1
        """

rule filtlong:
    input:
        fastq="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop.fastq.gz"
    output:
        fastq="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    params:
        min_length=config.get("general", {}).get("min_length", 1000),
        min_mean_q=config.get("general", {}).get("min_mean_q", 12)
    threads: config.get("resources", {}).get("filtlong", {}).get("threads", 2)
    resources:
        mem_mb=config.get("resources", {}).get("filtlong", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("filtlong", {}).get("walltime", "8:00:00")
    log:
        "output/logs/01.2_filtered/porechop_filtlong/{sample}_filtlong.log"
    benchmark:
        "output/logs/01.2_filtered/porechop_filtlong/{sample}_filtlong_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        # Aseguramos que la salida de filtlong esté comprimida correctamente con gzip
        filtlong --min_length {params.min_length} --min_mean_q {params.min_mean_q} {input.fastq} 2> {log} | gzip > {output.fastq}
        """

rule nanoplot_filtered:
    input:
        fastq="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    output:
        html="output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}/NanoPlot-report.html",
        stats="output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}/NanoStats.txt"
    params:
        outdir="output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}"
    threads: config.get("resources", {}).get("nanoplot", {}).get("threads", 6)
    resources:
        mem_mb=config.get("resources", {}).get("nanoplot", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("nanoplot", {}).get("walltime", "8:00:00")
    log:
        "output/logs/01.3_fastq_filtered_qc/nanoplot/{sample}.log"
    benchmark:
        "output/logs/01.3_fastq_filtered_qc/nanoplot/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        NanoPlot --fastq {input.fastq} --outdir {params.outdir} --threads {threads} > {log} 2>&1
        """

rule kraken2_nanopore:
    input:
        fastq="output/01_data/01.2_filtered/porechop_filtlong/{sample}_porechop_filtlong.fastq.gz"
    output:
        report="output/01_data/01.4_kraken2/nanopore/{sample}_report.txt",
        output="output/01_data/01.4_kraken2/nanopore/{sample}_output.txt"
    params:
        db=config["kraken_db"],
        confidence=config["kraken_conf"]
    threads: config.get("resources", {}).get("kraken2", {}).get("threads", 32)
    resources:
        mem_mb=config.get("resources", {}).get("kraken2", {}).get("mem", 64000),
        walltime=config.get("resources", {}).get("kraken2", {}).get("walltime", "12:00:00"),
        kraken_db=1
    log:
        "output/logs/01.4_kraken2/nanopore/{sample}.log"
    benchmark:
        "output/logs/01.4_kraken2/nanopore/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.report})
        
        kraken2 --db {params.db} \\
                --memory-mapping \\
                --threads {threads} \\
                --confidence {params.confidence} \\
                --report {output.report} \\
                --output {output.output} \\
                --use-names \\
                {input.fastq} > {log} 2>&1
        """