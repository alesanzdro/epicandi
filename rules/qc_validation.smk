# ==============================
# QC Y VALIDACIÓN DE MUESTRAS
# ==============================

rule samples_metrics_table:
    input:
        # Archivos de validación QC
        illumina_fastp=expand("output/01_data/01.2_filtered/fastp/{sample}_fastp.json", sample=ILLUMINA_SAMPLES) if ILLUMINA_SAMPLES else [],
        nanopore_stats_raw=expand("output/01_data/01.1_fastq_raw_qc/nanoplot/{sample}/NanoStats.txt", sample=NANOPORE_SAMPLES) if NANOPORE_SAMPLES else [],
        nanopore_stats_filtered=expand("output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}/NanoStats.txt", sample=NANOPORE_SAMPLES) if NANOPORE_SAMPLES else [],
        # Taxonomía
        illumina_kraken=expand("output/01_data/01.4_kraken2/illumina/{sample}_report.txt", sample=ILLUMINA_SAMPLES) if ILLUMINA_SAMPLES else [],
        nanopore_kraken=expand("output/01_data/01.4_kraken2/nanopore/{sample}_report.txt", sample=NANOPORE_SAMPLES) if NANOPORE_SAMPLES else []
    output:
        metrics_table="output/01_data/01.5_samples_pass/samples_metrics_table.tsv"
    log:
        "output/logs/01.5_samples_pass/metrics_table.log"
    benchmark:
        "output/logs/01.5_samples_pass/metrics_table_benchmark.txt"
    run:
        import json
        import re
        import pandas as pd
        from datetime import datetime
        
        shell("mkdir -p $(dirname {log})")
        
        with open(log[0], "w") as log_file:
            log_file.write("=== GENERANDO TABLA DE MÉTRICAS DE MUESTRAS ===\n")
            log_file.write(f"Fecha: {datetime.now()}\n\n")
            
            # Definir umbrales para QC (ahora desde config)
            illumina_min_reads = config.get("qc_thresholds", {}).get("illumina_min_reads", QC_MIN_READS_ILLUMINA)
            nanopore_min_reads = config.get("qc_thresholds", {}).get("nanopore_min_reads", QC_MIN_READS_NANOPORE)
            
            # Lista para almacenar los datos de cada muestra
            metrics_data = []
            
            # Procesar cada muestra
            for sample_id in SAMPLE_NAMES:
                data = SAMPLES[sample_id]
                sample_type = data['type']
                
                # Valores por defecto
                metrics = {
                    'Sample': sample_id,
                    'target_assembly': sample_type,  # Por defecto el tipo de muestra
                    'final_assembly': None,  # Se determinará luego
                    'Nanopore_total_reads': 'NA',
                    'Nanopore_len_mean': 'NA',
                    'Nanopore_q15_bases': 'NA',
                    'Nanopore_q15_rate': 'NA',
                    'Nanopore_gc_content': 'NA',
                    'Nanopore_keep_rate': 'NA',
                    'Illumina_total_reads': 'NA',
                    'Illumina_len_mean': 'NA',
                    'Illumina_q30_bases': 'NA',
                    'Illumina_q30_rate': 'NA', 
                    'Illumina_gc_content': 'NA',
                    'Illumina_keep_rate': 'NA',
                    'Nanopore_sample': 'NA',
                    'Illumina_sample': 'NA'
                }
                
                # Variables para rastrear si pasa el QC
                illumina_pass = False
                nanopore_pass = False
                
                # Procesar métricas de Illumina
                if sample_type in ['illumina', 'hybrid']:
                    fastp_file = f"output/01_data/01.2_filtered/fastp/{sample_id}_fastp.json"
                    
                    try:
                        with open(fastp_file, 'r') as f:
                            fastp_data = json.load(f)
                            
                            # Datos antes de filtrado
                            raw_reads = fastp_data['summary']['before_filtering']['total_reads'] // 2
                            
                            # Datos después de filtrado
                            clean_reads = fastp_data['summary']['after_filtering']['total_reads'] // 2
                            total_bases = fastp_data['summary']['after_filtering']['total_bases']
                            q30_bases = fastp_data['summary']['after_filtering']['q30_bases']
                            q30_rate = (q30_bases / total_bases) if total_bases > 0 else 0
                            gc_content = fastp_data['summary']['after_filtering']['gc_content']
                            keep_rate = clean_reads / raw_reads if raw_reads > 0 else 0
                            
                            # Longitud media
                            len_mean = total_bases / clean_reads / 2 if clean_reads > 0 else 0
                            
                            # Actualizar métricas
                            metrics['Illumina_total_reads'] = clean_reads
                            metrics['Illumina_len_mean'] = round(len_mean, 1)
                            metrics['Illumina_q30_bases'] = q30_bases
                            metrics['Illumina_q30_rate'] = round(q30_rate, 6)
                            metrics['Illumina_gc_content'] = round(gc_content, 6)
                            metrics['Illumina_keep_rate'] = round(keep_rate, 6)
                            
                            # Determinar si pasa QC
                            illumina_pass = clean_reads >= illumina_min_reads
                            metrics['Illumina_sample'] = 'pass' if illumina_pass else 'fail'
                            
                            log_file.write(f"{sample_id} (Illumina): {clean_reads:,} reads - {'PASS' if illumina_pass else 'FAIL'}\n")
                    
                    except Exception as e:
                        log_file.write(f"ERROR procesando datos Illumina para {sample_id}: {e}\n")
                
                # Procesar métricas de Nanopore
                if sample_type in ['nanopore', 'hybrid']:
                    nano_stats_file = f"output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample_id}/NanoStats.txt"
                    nano_raw_file = f"output/01_data/01.1_fastq_raw_qc/nanoplot/{sample_id}/NanoStats.txt"
                    
                    try:
                        # Datos de lecturas filtradas
                        with open(nano_stats_file, 'r') as f:
                            content = f.read()
                            
                            # Extraer métricas generales con regex
                            reads_match = re.search(r'Number of reads:\s*([0-9,]+)', content)
                            len_mean_match = re.search(r'Mean read length:\s*([0-9,.]+)', content)
                            q_mean_match = re.search(r'Mean read quality:\s*([0-9,.]+)', content)
                            gc_match = re.search(r'GC content \(%\):\s*([0-9,.]+)', content)
                            
                            # Extraer métricas de calidad específicas
                            total_bases_match = re.search(r'Total bases:\s*([0-9,\.]+)', content)
                            q15_match = re.search(r'>Q15:\s*(\d+)\s+\((\d+\.\d+)%\)\s+([\d\.]+)Mb', content)
                            
                            # Procesar valores extraídos
                            clean_reads = int(reads_match.group(1).replace(',', '')) if reads_match else 0
                            len_mean = float(len_mean_match.group(1).replace(',', '')) if len_mean_match else 0
                            q_mean = float(q_mean_match.group(1)) if q_mean_match else 0
                            
                            # Obtener total de bases
                            total_bases = float(total_bases_match.group(1).replace(',', '')) if total_bases_match else (clean_reads * len_mean)
                            
                            # Obtener información de Q15 directamente del archivo
                            q15_reads = 0
                            q15_rate = 0
                            q15_bases = 0
                            
                            if q15_match:
                                q15_reads = int(q15_match.group(1))
                                q15_rate = float(q15_match.group(2)) / 100  # Convertir porcentaje a decimal
                                q15_bases = float(q15_match.group(3)) * 1_000_000  # Convertir Mb a bases
                            else:
                                # Si no se encuentra el patrón específico, hacer una aproximación
                                q15_rate = max(0, min(1, (q_mean - 7) / 8))  # Aproximación para Q15: 7→0%, 15→100%
                                q15_bases = int(total_bases * q15_rate)
                            
                            # Obtener GC content
                            gc_content = float(gc_match.group(1))/100 if gc_match else 0
                            
                        # Datos de lecturas crudas para calcular keep_rate
                        with open(nano_raw_file, 'r') as f:
                            raw_content = f.read()
                            raw_reads_match = re.search(r'Number of reads:\s*([0-9,]+)', raw_content)
                            raw_reads = int(raw_reads_match.group(1).replace(',', '')) if raw_reads_match else 0
                            
                        keep_rate = clean_reads / raw_reads if raw_reads > 0 else 0
                            
                        # Actualizar métricas
                        metrics['Nanopore_total_reads'] = clean_reads
                        metrics['Nanopore_len_mean'] = round(len_mean, 1)
                        metrics['Nanopore_q15_bases'] = q15_bases
                        metrics['Nanopore_q15_rate'] = round(q15_rate, 6)
                        metrics['Nanopore_gc_content'] = round(gc_content, 6)
                        metrics['Nanopore_keep_rate'] = round(keep_rate, 6)
                        
                        # Determinar si pasa QC
                        nanopore_pass = clean_reads >= nanopore_min_reads
                        metrics['Nanopore_sample'] = 'pass' if nanopore_pass else 'fail'
                        
                        log_file.write(f"{sample_id} (Nanopore): {clean_reads:,} reads - {'PASS' if nanopore_pass else 'FAIL'}\n")
                    
                    except Exception as e:
                        log_file.write(f"ERROR procesando datos Nanopore para {sample_id}: {e}\n")
                
                # Determinar el tipo de ensamblado final basado en QC
                if sample_type == 'hybrid':
                    if illumina_pass and nanopore_pass:
                        metrics['final_assembly'] = 'hybrid'
                    elif illumina_pass:
                        metrics['final_assembly'] = 'short_reads'
                    elif nanopore_pass:
                        metrics['final_assembly'] = 'long_reads'
                    else:
                        metrics['final_assembly'] = 'fail'
                elif sample_type == 'illumina':
                    metrics['final_assembly'] = 'short_reads' if illumina_pass else 'fail'
                elif sample_type == 'nanopore':
                    metrics['final_assembly'] = 'long_reads' if nanopore_pass else 'fail'
                
                metrics_data.append(metrics)
                log_file.write(f"{sample_id}: target={sample_type}, final={metrics['final_assembly']}\n")
            
            # Crear DataFrame y guardar como TSV
            df = pd.DataFrame(metrics_data)
            df.to_csv(output.metrics_table, sep='\t', index=False)
            log_file.write(f"Tabla de métricas guardada en: {output.metrics_table}\n")
            log_file.write("=== TABLA DE MÉTRICAS COMPLETADA ===\n")

rule samples_pass:
    input:
        # Utilizamos la tabla de métricas como entrada
        metrics_table="output/01_data/01.5_samples_pass/samples_metrics_table.tsv",
        # Archivos de validación QC
        illumina_fastp=expand("output/01_data/01.2_filtered/fastp/{sample}_fastp.json", sample=ILLUMINA_SAMPLES) if ILLUMINA_SAMPLES else [],
        nanopore_stats=expand("output/01_data/01.3_fastq_filtered_qc/nanoplot/{sample}/NanoStats.txt", sample=NANOPORE_SAMPLES) if NANOPORE_SAMPLES else [],
        # Taxonomía
        illumina_kraken=expand("output/01_data/01.4_kraken2/illumina/{sample}_report.txt", sample=ILLUMINA_SAMPLES) if ILLUMINA_SAMPLES else [],
        nanopore_kraken=expand("output/01_data/01.4_kraken2/nanopore/{sample}_report.txt", sample=NANOPORE_SAMPLES) if NANOPORE_SAMPLES else []
    output:
        report="output/01_data/01.5_samples_pass/validation_report.txt",
        valid_samples="output/01_data/01.5_samples_pass/samplesinfo_valid.csv"
    log:
        "output/logs/01.5_samples_pass/validation.log"
    benchmark:
        "output/logs/01.5_samples_pass/validation_benchmark.txt"
    run:
        import pandas as pd
        from datetime import datetime
        
        shell("mkdir -p $(dirname {log})")
        
        with open(log[0], "w") as log_file:
            log_file.write("=== VALIDACION DE MUESTRAS ===\n")
            log_file.write(f"Fecha: {datetime.now()}\n\n")
            
            # Cargar la tabla de métricas
            metrics_df = pd.read_csv(input.metrics_table, sep='\t')
            
            # Generar reporte basado en la tabla de métricas
            passed_samples = metrics_df[metrics_df['final_assembly'] != 'fail']
            failed_samples = metrics_df[metrics_df['final_assembly'] == 'fail']
            
            log_file.write(f"Total de muestras: {len(metrics_df)}\n")
            log_file.write(f"Muestras que pasaron: {len(passed_samples)}\n")
            log_file.write(f"Muestras que fallaron: {len(failed_samples)}\n\n")
            
            for _, row in passed_samples.iterrows():
                sample_id = row['Sample']
                target = row['target_assembly']
                final = row['final_assembly']
                log_file.write(f"{sample_id}: {target} → {final}\n")
            
            log_file.write("\nMuestras que fallaron QC:\n")
            for _, row in failed_samples.iterrows():
                sample_id = row['Sample']
                target = row['target_assembly']
                log_file.write(f"{sample_id}: {target} → FAIL\n")
        
        # Crear archivo de reporte
        with open(output.report, "w") as f:
            f.write("=== REPORTE DE VALIDACIÓN DE MUESTRAS ===\n")
            f.write(f"Fecha: {datetime.now()}\n\n")
            
            # Resumen de estrategias de ensamblado
            hybrid_count = len(metrics_df[metrics_df['final_assembly'] == 'hybrid'])
            short_reads_count = len(metrics_df[metrics_df['final_assembly'] == 'short_reads'])
            long_reads_count = len(metrics_df[metrics_df['final_assembly'] == 'long_reads'])
            failed_count = len(metrics_df[metrics_df['final_assembly'] == 'fail'])
            
            f.write("RESUMEN DE ESTRATEGIAS DE ENSAMBLADO:\n")
            f.write(f"Total de muestras: {len(metrics_df)}\n")
            f.write(f"Ensamblado híbrido: {hybrid_count}\n")
            f.write(f"Ensamblado lecturas cortas: {short_reads_count}\n")
            f.write(f"Ensamblado lecturas largas: {long_reads_count}\n")
            f.write(f"Muestras fallidas: {failed_count}\n\n")
            
            f.write("DETALLE POR MUESTRA:\n")
            for _, row in metrics_df.iterrows():
                sample_id = row['Sample']
                target = row['target_assembly']
                final = row['final_assembly']
                
                f.write(f"\n{sample_id}:\n")
                f.write(f"  Estrategia objetivo: {target}\n")
                f.write(f"  Estrategia final: {final}\n")
                
                # Detalles para Illumina si disponibles
                if pd.notna(row['Illumina_sample']) and row['Illumina_sample'] != 'NA':
                    ill_status = str(row['Illumina_sample'])
                    ill_reads = row['Illumina_total_reads']
                    if pd.notna(ill_reads) and ill_reads != 'NA':
                        f.write(f"  Illumina: {ill_status.upper()} ({ill_reads:,} lecturas)\n")
                    else:
                        f.write(f"  Illumina: {ill_status.upper()}\n")

                # Detalles para Nanopore si disponibles
                if pd.notna(row['Nanopore_sample']) and row['Nanopore_sample'] != 'NA':
                    nano_status = str(row['Nanopore_sample'])
                    nano_reads = row['Nanopore_total_reads']
                    if pd.notna(nano_reads) and nano_reads != 'NA':
                        f.write(f"  Nanopore: {nano_status.upper()} ({nano_reads:,} lecturas)\n")
                    else:
                        f.write(f"  Nanopore: {nano_status.upper()}\n")
                        
                # Indicar si hubo cambio de estrategia
                if target != final and final != 'fail':
                    f.write(f"  Nota: Estrategia cambiada de {target} a {final} basado en control de calidad\n")
            
        # Crear samplesinfo_valid.csv con las muestras que pasaron
        valid_assemblies = ['hybrid', 'short_reads', 'long_reads']
        valid_samples_df = metrics_df[metrics_df['final_assembly'].isin(valid_assemblies)]

        # Cargar el samplesinfo.csv original
        original_csv = pd.read_csv("samplesinfo.csv")

        # Filtrar por las muestras válidas
        valid_ids = set(valid_samples_df['Sample'].tolist())
        filtered_csv = original_csv[original_csv['id'].isin(valid_ids)]

        # Guardar el nuevo CSV
        filtered_csv.to_csv(output.valid_samples, index=False)

# Checkpoint para determinar qué muestras han pasado el control de calidad
checkpoint samples_validation:
    input:
        metrics_table="output/01_data/01.5_samples_pass/samples_metrics_table.tsv"
    output:
        # Archivo principal del checkpoint con las muestras que pasan validación
        passed_samples="output/01_data/01.5_samples_pass/passed_samples.txt",
        # Archivos secundarios por estrategia de ensamblado
        hybrid_samples="output/01_data/01.5_samples_pass/hybrid_samples.txt",
        short_reads_samples="output/01_data/01.5_samples_pass/short_reads_samples.txt",
        long_reads_samples="output/01_data/01.5_samples_pass/long_reads_samples.txt",
        # Reporte de validación
        validation_report="output/01_data/01.5_samples_pass/checkpoint_validation_report.txt"
    log:
        "output/logs/01.5_samples_pass/checkpoint_validation.log"
    run:
        import pandas as pd
        from datetime import datetime

        # Crear directorios necesarios
        shell("mkdir -p $(dirname {output.passed_samples})")
        shell("mkdir -p $(dirname {log})")

        with open(log[0], "w") as log_file:
            log_file.write("=== CHECKPOINT: VALIDACIÓN DE MUESTRAS PARA ENSAMBLADO ===\n")
            log_file.write(f"Fecha: {datetime.now()}\n\n")

            # Cargar tabla de métricas
            metrics_df = pd.read_csv(input.metrics_table, sep='\t')

            # Filtrar muestras con estrategia de ensamblado válida (no 'fail')
            valid_assemblies = ['hybrid', 'short_reads', 'long_reads']
            valid_samples = metrics_df[metrics_df['final_assembly'].isin(valid_assemblies)]

            # Clasificar por tipo de ensamblado
            hybrid_samples = valid_samples[valid_samples['final_assembly'] == 'hybrid']
            sr_samples = valid_samples[valid_samples['final_assembly'] == 'short_reads']
            lr_samples = valid_samples[valid_samples['final_assembly'] == 'long_reads']

            # Obtener listas de IDs
            all_valid_ids = sorted(valid_samples['Sample'].tolist())
            hybrid_ids = sorted(hybrid_samples['Sample'].tolist())
            sr_ids = sorted(sr_samples['Sample'].tolist())
            lr_ids = sorted(lr_samples['Sample'].tolist())

            log_file.write(f"Total muestras para ensamblar: {len(all_valid_ids)}\n")
            log_file.write(f"Ensamblado híbrido (Flye + Medaka + Pypolca): {len(hybrid_ids)} muestras\n")
            log_file.write(f"Ensamblado lecturas cortas (SPAdes + Pilon): {len(sr_ids)} muestras\n")
            log_file.write(f"Ensamblado lecturas largas (Flye + Medaka): {len(lr_ids)} muestras\n\n")

            # CHECKPOINT CRÍTICO: Escribir archivo principal con todas las muestras válidas
            with open(output.passed_samples, "w") as f:
                for sample_id in all_valid_ids:
                    f.write(f"{sample_id}\n")
            log_file.write(f"Archivo checkpoint creado: {output.passed_samples} con {len(all_valid_ids)} muestras\n")

            # Escribir archivos por estrategia
            with open(output.hybrid_samples, "w") as f:
                for sample_id in hybrid_ids:
                    f.write(f"{sample_id}\n")
            log_file.write(f"Muestras híbridas: {output.hybrid_samples} con {len(hybrid_ids)} muestras\n")

            with open(output.short_reads_samples, "w") as f:
                for sample_id in sr_ids:
                    f.write(f"{sample_id}\n")
            log_file.write(f"Muestras lecturas cortas: {output.short_reads_samples} con {len(sr_ids)} muestras\n")

            with open(output.long_reads_samples, "w") as f:
                for sample_id in lr_ids:
                    f.write(f"{sample_id}\n")
            log_file.write(f"Muestras lecturas largas: {output.long_reads_samples} con {len(lr_ids)} muestras\n")

            # Escribir informe detallado de validación
            with open(output.validation_report, "w") as f:
                f.write("=== CHECKPOINT: MUESTRAS VÁLIDAS PARA ENSAMBLADO ===\n")
                f.write(f"Fecha: {datetime.now()}\n\n")

                f.write("RESUMEN DE ESTRATEGIAS:\n")
                f.write(f"Total muestras aprobadas: {len(all_valid_ids)}\n")
                f.write(f"Ensamblado híbrido (Flye + Medaka + Pypolca): {len(hybrid_ids)} muestras\n")
                f.write(f"Ensamblado lecturas cortas (SPAdes + Pilon): {len(sr_ids)} muestras\n")
                f.write(f"Ensamblado lecturas largas (Flye + Medaka): {len(lr_ids)} muestras\n\n")

                f.write("DETALLE POR MUESTRA:\n")
                for _, row in valid_samples.iterrows():
                    sample_id = row['Sample']
                    strategy = row['final_assembly']
                    f.write(f"{sample_id}: {strategy}\n")

                # Añadir información sobre muestras fallidas
                failed_samples = metrics_df[metrics_df['final_assembly'] == 'fail']
                if len(failed_samples) > 0:
                    f.write(f"\nMuestras que fallaron QC ({len(failed_samples)}):\n")
                    for _, row in failed_samples.iterrows():
                        sample_id = row['Sample']
                        target = row['target_assembly']
                        f.write(f"{sample_id}: {target} → FAIL\n")

            log_file.write("=== CHECKPOINT COMPLETADO ===\n")