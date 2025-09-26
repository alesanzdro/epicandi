# ==============================
# REGLAS DE ANÁLISIS DE COBERTURA
# ==============================

def get_coverage_inputs(wildcards):
    """Obtener inputs correctos según el tipo de muestra para análisis de cobertura"""
    sample_type = SAMPLE_TYPES.get(wildcards.sample)

    # Todos los tipos incluyen assembly
    inputs = {
        "assembly": f"output/02_assembly/02.3_consensus/{wildcards.sample}.fasta"
    }

    if sample_type == "hybrid":
        # Para híbridas, usar lectura largas (nanopore)
        inputs["reads"] = f"output/01_data/01.2_filtered/porechop_filtlong/{wildcards.sample}_porechop_filtlong.fastq.gz"
    elif sample_type == "nanopore":
        # Para nanopore only
        inputs["reads"] = f"output/01_data/01.2_filtered/porechop_filtlong/{wildcards.sample}_porechop_filtlong.fastq.gz"
    elif sample_type == "illumina":
        # Para illumina, usar BWA mem
        inputs["reads_r1"] = f"output/01_data/01.2_filtered/fastp/{wildcards.sample}_fastp_R1.fastq.gz"
        inputs["reads_r2"] = f"output/01_data/01.2_filtered/fastp/{wildcards.sample}_fastp_R2.fastq.gz"

    return inputs

def get_reads_for_shell(wildcards, input):
    """Obtener archivos de reads como string para shell"""
    sample_type = SAMPLE_TYPES.get(wildcards.sample)

    if sample_type == "illumina":
        return f"{input.reads_r1} {input.reads_r2}"
    else:
        return f"{input.reads}"

rule alignment_minimap2:
    input:
        unpack(get_coverage_inputs)
    output:
        bam=temp("output/02_assembly/02.6_coverage/{sample}/{sample}.unsorted.bam")
    params:
        sample_type=lambda wildcards: SAMPLE_TYPES.get(wildcards.sample),
        reads_files=lambda wildcards, input: get_reads_for_shell(wildcards, input)
    log:
        "output/logs/02_assembly/02.6_coverage/{sample}_minimap2.log"
    benchmark:
        "output/logs/02_assembly/02.6_coverage/{sample}_minimap2_benchmark.txt"
    threads: config.get("resources", {}).get("minimap2", {}).get("threads", 12)
    resources:
        mem_mb=config.get("resources", {}).get("minimap2", {}).get("mem", 16000),
        walltime=config.get("resources", {}).get("minimap2", {}).get("walltime", "04:00:00")
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})

        echo "Running alignment for {wildcards.sample} (type: {params.sample_type})" > {log} 2>&1

        if [[ "{params.sample_type}" == "illumina" ]]; then
            # BWA MEM para Illumina
            echo "Using BWA MEM for Illumina reads" >> {log} 2>&1

            # Indexar ensamblado si no existe
            if [ ! -f "{input.assembly}.bwt" ]; then
                echo "Indexing assembly with BWA..." >> {log} 2>&1
                bwa index {input.assembly} >> {log} 2>&1
            fi

            echo "Reads files: {params.reads_files}" >> {log} 2>&1

            bwa mem -t {threads} {input.assembly} {params.reads_files} 2>> {log} |
            samtools view -b -@ 2 -F 4 -F 256 -F 2048 - > {output.bam} 2>> {log}
        else
            # Minimap2 para Nanopore
            echo "Using minimap2 for Nanopore reads" >> {log} 2>&1
            minimap2 --version >> {log} 2>&1

            echo "Nanopore reads: {params.reads_files}" >> {log} 2>&1

            minimap2 -ax map-ont \\
                -t {threads} \\
                --secondary=no \\
                -K 500M \\
                {input.assembly} \\
                {params.reads_files} 2>> {log} |
            samtools view -b -@ 2 -F 4 -F 256 -F 2048 - > {output.bam} 2>> {log}
        fi

        echo "Alignment completed. Quick stats:" >> {log} 2>&1
        samtools flagstat {output.bam} 2>/dev/null >> {log} 2>&1
        """

rule samtools_processing:
    input:
        bam="output/02_assembly/02.6_coverage/{sample}/{sample}.unsorted.bam"
    output:
        sorted_bam="output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam",
        sorted_bai="output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam.bai",
        stats="output/02_assembly/02.6_coverage/{sample}/{sample}.stats.txt",
        coverage="output/02_assembly/02.6_coverage/{sample}/{sample}.coverage.txt",
        flagstat="output/02_assembly/02.6_coverage/{sample}/{sample}.flagstat.txt",
        idxstats="output/02_assembly/02.6_coverage/{sample}/{sample}.idxstats.txt"
    log:
        "output/logs/02_assembly/02.6_coverage/{sample}_samtools.log"
    benchmark:
        "output/logs/02_assembly/02.6_coverage/{sample}_samtools_benchmark.txt"
    threads: config.get("resources", {}).get("samtools", {}).get("threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("samtools", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("samtools", {}).get("walltime", "02:00:00")
    conda: "../envs/epicandi.yml"
    shell:
        """
        echo "Processing BAM file for {wildcards.sample}" > {log} 2>&1

        # Step 1: Sort BAM file
        echo "Sorting BAM file..." >> {log} 2>&1
        samtools sort \\
            -@ {threads} \\
            -m 1G \\
            -o {output.sorted_bam} \\
            {input.bam} >> {log} 2>&1

        # Step 2: Index sorted BAM
        echo "Indexing sorted BAM..." >> {log} 2>&1
        samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1

        # Step 3: Generate comprehensive statistics
        echo "Generating statistics..." >> {log} 2>&1
        samtools stats -@ {threads} {output.sorted_bam} > {output.stats} 2>> {log}

        # Step 4: Generate alignment statistics
        echo "Generating flagstat..." >> {log} 2>&1
        samtools flagstat -@ {threads} {output.sorted_bam} > {output.flagstat} 2>> {log}

        # Step 5: Generate index statistics
        echo "Generating idxstats..." >> {log} 2>&1
        samtools idxstats {output.sorted_bam} > {output.idxstats} 2>> {log}

        # Step 6: Calculate per-contig coverage
        echo "Calculating coverage..." >> {log} 2>&1
        samtools coverage \\
            --min-MQ 10 \\
            --min-BQ 10 \\
            {output.sorted_bam} > {output.coverage} 2>> {log}

        # Summary
        echo "Processing complete for {wildcards.sample}" >> {log} 2>&1
        echo "Coverage summary:" >> {log} 2>&1
        tail -n +2 {output.coverage} | \\
            awk '{{sum+=$7*$3; total+=$3}} END {{if(total>0) printf "Average coverage: %.2fX\\n", sum/total}}' >> {log} 2>&1
        """

rule mosdepth_coverage:
    input:
        sorted_bam="output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam",
        sorted_bai="output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam.bai"
    output:
        dist="output/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.global.dist.txt",
        summary="output/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt",
        per_base="output/02_assembly/02.7_mosdepth/{sample}/{sample}.per-base.bed.gz",
        regions="output/02_assembly/02.7_mosdepth/{sample}/{sample}.regions.bed.gz"
    log:
        "output/logs/02_assembly/02.7_mosdepth/{sample}_mosdepth.log"
    benchmark:
        "output/logs/02_assembly/02.7_mosdepth/{sample}_mosdepth_benchmark.txt"
    params:
        prefix="output/02_assembly/02.7_mosdepth/{sample}/{sample}",
        window_size=500
    threads: config.get("resources", {}).get("mosdepth", {}).get("threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("mosdepth", {}).get("mem", 4000),
        walltime=config.get("resources", {}).get("mosdepth", {}).get("walltime", "01:00:00")
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p $(dirname {params.prefix}) $(dirname {log})

        echo "Running mosdepth for {wildcards.sample}" > {log} 2>&1
        mosdepth --version >> {log} 2>&1

        # Run mosdepth with comprehensive options
        mosdepth \\
            --threads {threads} \\
            --fast-mode \\
            --by {params.window_size} \\
            --mapq 10 \\
            {params.prefix} \\
            {input.sorted_bam} >> {log} 2>&1

        echo "Mosdepth summary:" >> {log} 2>&1
        cat {output.summary} >> {log} 2>&1
        """