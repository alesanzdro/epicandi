# ==============================
# REGLAS ILLUMINA
# ==============================

rule fastqc_raw_illumina:
    input:
        r1="output/01_data/00_staged_data/illumina/{sample}_R1.fastq.gz",
        r2="output/01_data/00_staged_data/illumina/{sample}_R2.fastq.gz",
        staging="output/01_data/00_staged_data/.staging_complete"
    output:
        r1="output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R1_fastqc.html",
        r2="output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R2_fastqc.html",
        r1_zip="output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R1_fastqc.zip",
        r2_zip="output/01_data/01.1_fastq_raw_qc/fastqc/{sample}_R2_fastqc.zip"
    params:
        outdir="output/01_data/01.1_fastq_raw_qc/fastqc"
    threads: get_resource("default", "threads", 4)
    resources:
        mem_mb=get_resource("default", "mem", 4000),
        walltime=get_resource("default", "walltime", "1:00:00")
    log:
        "output/logs/01.1_fastq_raw_qc/fastqc/{sample}.log"
    benchmark:
        "output/logs/01.1_fastq_raw_qc/fastqc/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}

        fastqc {input.r1} {input.r2} --outdir {params.outdir} --threads {threads} > {log} 2>&1
        """

rule fastp:
    input:
        r1="output/01_data/00_staged_data/illumina/{sample}_R1.fastq.gz",
        r2="output/01_data/00_staged_data/illumina/{sample}_R2.fastq.gz",
        staging="output/01_data/00_staged_data/.staging_complete"
    output:
        r1="output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz",
        r2="output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz",
        html="output/01_data/01.2_filtered/fastp/{sample}_fastp.html",
        json="output/01_data/01.2_filtered/fastp/{sample}_fastp.json"
    threads: config.get("resources", {}).get("fastp", {}).get("threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("fastp", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("fastp", {}).get("walltime", "2:00:00")
    log:
        "output/logs/01.2_filtered/fastp/{sample}.log"
    benchmark:
        "output/logs/01.2_filtered/fastp/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        
        fastp -i {input.r1} -I {input.r2} \\
              -o {output.r1} -O {output.r2} \\
              --html {output.html} --json {output.json} \\
              --thread {threads} \\
              --detect_adapter_for_pe \\
              --correction --cut_front --cut_tail \\
              --cut_window_size 4 --cut_mean_quality 15 \\
              --qualified_quality_phred 15 \\
              --unqualified_percent_limit 40 \\
              --length_required 36 > {log} 2>&1
        """

rule fastqc_trim_illumina:
    input:
        r1="output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz",
        r2="output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz"
    output:
        r1="output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R1_fastqc.html",
        r2="output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R2_fastqc.html",
        r1_zip="output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R1_fastqc.zip",
        r2_zip="output/01_data/01.3_fastq_filtered_qc/fastqc/{sample}_fastp_R2_fastqc.zip"
    params:
        outdir="output/01_data/01.3_fastq_filtered_qc/fastqc"
    threads: config.get("resources", {}).get("default", {}).get("threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("default", {}).get("mem", 4000),
        walltime=config.get("resources", {}).get("default", {}).get("walltime", "1:00:00")
    log:
        "output/logs/01.3_fastq_filtered_qc/fastqc/{sample}.log"
    benchmark:
        "output/logs/01.3_fastq_filtered_qc/fastqc/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.outdir}

        fastqc {input.r1} {input.r2} --outdir {params.outdir} --threads {threads} > {log} 2>&1
        """

rule kraken2_illumina:
    input:
        r1="output/01_data/01.2_filtered/fastp/{sample}_fastp_R1.fastq.gz",
        r2="output/01_data/01.2_filtered/fastp/{sample}_fastp_R2.fastq.gz"
    output:
        report="output/01_data/01.4_kraken2/illumina/{sample}_report.txt",
        output="output/01_data/01.4_kraken2/illumina/{sample}_output.txt"
    params:
        db=config["kraken_db"],
        confidence=config["kraken_conf"]
    threads: config.get("resources", {}).get("kraken2", {}).get("threads", 32)
    resources:
        mem_mb=config.get("resources", {}).get("kraken2", {}).get("mem", 64000),
        walltime=config.get("resources", {}).get("kraken2", {}).get("walltime", "12:00:00"),
        kraken_db=1
    log:
        "output/logs/01.4_kraken2/illumina/{sample}.log"
    benchmark:
        "output/logs/01.4_kraken2/illumina/{sample}_benchmark.txt"
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.report})
        
        kraken2 --db {params.db} \\
                --threads {threads} \\
                --memory-mapping \\
                --confidence {params.confidence} \\
                --paired \\
                --report {output.report} \\
                --output {output.output} \\
                --use-names \\
                {input.r1} {input.r2} > {log} 2>&1
        """