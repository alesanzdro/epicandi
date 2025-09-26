# ==============================
# REGLAS DE REPORTE
# ==============================

rule assembly_report:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        validation_report="output/01_data/01.5_samples_pass/checkpoint_validation_report.txt",
        assemblies=lambda wildcards: expand("output/02_assembly/02.3_consensus/{sample}.fasta",
                                       sample=get_validated_samples(wildcards))
    output:
        report="output/02_assembly/assembly_report.txt"
    log:
        "output/logs/02_assembly/assembly_report.log"
    run:
        import pandas as pd
        from datetime import datetime
        import os
        
        shell("mkdir -p $(dirname {log})")
        
        with open(log[0], "w") as log_file:
            log_file.write("=== GENERANDO REPORTE DE ENSAMBLADO ===\n")
            log_file.write(f"Fecha: {datetime.now()}\n\n")
            
            # Recopilar información de los ensamblados
            valid_samples = get_validated_samples(wildcards)
            # Obtener información de estrategia de ensamblado para cada muestra
            sample_strategy = {}
            
            # Cargar métricas para obtener la estrategia de ensamblado
            metrics_file = "output/01_data/01.5_samples_pass/samples_metrics_table.tsv"
            if os.path.exists(metrics_file):
                metrics_df = pd.read_csv(metrics_file, sep='\t')
                for _, row in metrics_df.iterrows():
                    if row['final_assembly'] in ['hybrid', 'short_reads', 'long_reads']:
                        sample_strategy[row['Sample']] = row['final_assembly']
            
            # Conteo de tipos de ensamblado
            hybrid_count = sum(1 for s, strategy in sample_strategy.items() if strategy == 'hybrid')
            sr_count = sum(1 for s, strategy in sample_strategy.items() if strategy == 'short_reads')
            lr_count = sum(1 for s, strategy in sample_strategy.items() if strategy == 'long_reads')
            
            all_assemblies = []
            
            # Procesar todos los ensamblados
            for sample in valid_samples:
                fasta_file = f"output/02_assembly/02.3_consensus/{sample}.fasta"
                if os.path.exists(fasta_file):
                    # Determinar la estrategia de ensamblado
                    strategy = sample_strategy.get(sample, "unknown")
                    
                    # Contar contigs y calcular N50
                    contigs = 0
                    total_length = 0
                    lengths = []
                    
                    with open(fasta_file, 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                contigs += 1
                            else:
                                length = len(line.strip())
                                total_length += length
                                lengths.append(length)
                    
                    # Calcular N50
                    lengths.sort(reverse=True)
                    cumulative = 0
                    n50 = 0
                    for l in lengths:
                        cumulative += l
                        if cumulative >= total_length / 2:
                            n50 = l
                            break
                    
                    all_assemblies.append({
                        'Sample': sample,
                        'Strategy': strategy,
                        'Contigs': contigs,
                        'Total_Length': total_length,
                        'N50': n50
                    })
                    
                    log_file.write(f"Procesado ensamblado {strategy} para {sample}: {contigs} contigs, {total_length} bp, N50: {n50}\n")
            
            # Crear reporte final
            with open(output.report, "w") as f:
                f.write("=== REPORTE DE ENSAMBLADO FINAL ===\n")
                f.write(f"Fecha: {datetime.now()}\n\n")
                
                f.write("RESUMEN DE ENSAMBLADOS:\n")
                f.write(f"Total de ensamblados: {len(all_assemblies)}\n")
                f.write(f"Híbridos: {hybrid_count}\n")
                f.write(f"Lecturas cortas: {sr_count}\n")
                f.write(f"Lecturas largas: {lr_count}\n\n")
                
                f.write("DETALLE POR MUESTRA:\n")
                for assembly in all_assemblies:
                    f.write(f"\n{assembly['Sample']} ({assembly['Strategy']}):\n")
                    f.write(f"  Contigs: {assembly['Contigs']}\n")
                    f.write(f"  Longitud total: {assembly['Total_Length']} bp\n")
                    f.write(f"  N50: {assembly['N50']} bp\n")

rule multiqc:
    input:
        files=lambda wildcards: get_multiqc_inputs_with_evaluations(ILLUMINA_SAMPLES, NANOPORE_SAMPLES)
    output:
        html="output/04_report/multiqc_report.html"
    params:
        outdir="output/04_report",
        search_dirs="output/01_data output/02_assembly/02.4_evaluation",
        custom_data="--data-dir output/04_report/"
    threads: config.get("resources", {}).get("multiqc", {}).get("threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("multiqc", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("multiqc", {}).get("walltime", "1:00:00")
    log:
        "output/logs/04_report/multiqc.log"
    benchmark:
        "output/logs/04_report/multiqc_benchmark.txt"
    conda: "../envs/epicandi_phylo.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        
        multiqc {params.search_dirs} --outdir {params.outdir} {params.custom_data} --force \\
                --title "EpiCandi QC Report - Todas las Muestras" \\
                --comment "Pipeline completo: análisis de ensamblado óptimo" \\
                --filename multiqc_report.html > {log} 2>&1
        """