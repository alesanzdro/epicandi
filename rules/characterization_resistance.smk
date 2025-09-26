# ==============================
# CARACTERIZACIÓN DE RESISTENCIA Y VIRULENCIA
# Paquete 2: ChroQueTas + análisis de resistencia Python
# ==============================

# ==============================
# REGLAS DE CHROQUETAS
# ==============================

rule chroquetas_analysis:
    input:
        assembly="output/02_assembly/02.3_consensus/{sample}.fasta"
    output:
        summary="output/03_characterization/03.3_chroquetas/{sample}/{sample}.ChroQueTaS.AMR_summary.txt"
    params:
        outdir="output/03_characterization/03.3_chroquetas/{sample}",
        fungamr_db="/home/asanzc/.conda/envs/amr/FungAMR_db",
        species="Candidozyma_auris"
    log:
        "output/logs/03_characterization/03.3_chroquetas/{sample}_chroquetas.log"
    threads: 24
    conda: "../envs/epicandi_chroquetas.yml"
    shell:
        """
        mkdir -p $(dirname {log})
        rm -rf {params.outdir}

        ChroQueTas.sh \\
            -g {input.assembly} \\
            -o {params.outdir} \\
            -f {params.fungamr_db} \\
            -s {params.species} \\
            -t {threads} > {log} 2>&1
        """

# ==============================
# REGLAS DE ANÁLISIS DE RESISTENCIA PYTHON
# ==============================

rule custom_resistance_analysis:
    input:
        assembly="output/02_assembly/02.3_consensus/{sample}.fasta",
        bam="output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam",
        bai="output/02_assembly/02.6_coverage/{sample}/{sample}.sorted.bam.bai"
    output:
        all_mutations="output/03_characterization/03.4_resistance/{sample}/reports/todas_las_mutaciones.csv",
        resistance_mutations="output/03_characterization/03.4_resistance/{sample}/reports/mutaciones_de_resistencia.csv",
        coverage_report="output/03_characterization/03.4_resistance/{sample}/reports/coverage_report.tsv",
        summary_report="output/03_characterization/03.4_resistance/{sample}/resistance_analysis_summary.txt"
    params:
        outdir="output/03_characterization/03.4_resistance/{sample}",
        ref_proteins=config.get("resistance_analysis", {}).get("ref_proteins", "resources/mutations/ref_protein.faa"),
        mutations_tsv=config.get("resistance_analysis", {}).get("mutations_tsv", "resources/mutations/mutations.tsv"),
        script_path="scripts/analize_resistance.py"
    log:
        "output/logs/03_characterization/03.4_resistance/{sample}_resistance.log"
    benchmark:
        "output/logs/03_characterization/03.4_resistance/{sample}_resistance_benchmark.txt"
    threads: config.get("resources", {}).get("resistance_analysis", {}).get("threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("resistance_analysis", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("resistance_analysis", {}).get("walltime", "02:00:00")
    conda: "../envs/epicandi.yml"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        echo "=== EJECUTANDO ANÁLISIS DE RESISTENCIA PYTHON PARA {wildcards.sample} ===" > {log} 2>&1
        echo "Ensamblado: {input.assembly}" >> {log} 2>&1
        echo "BAM: {input.bam}" >> {log} 2>&1
        echo "Proteínas referencia: {params.ref_proteins}" >> {log} 2>&1
        echo "Mutaciones conocidas: {params.mutations_tsv}" >> {log} 2>&1

        # Verificar archivos de entrada
        if [ ! -f "{params.ref_proteins}" ]; then
            echo "ERROR: Archivo de proteínas de referencia no encontrado: {params.ref_proteins}" >> {log} 2>&1
            exit 1
        fi

        if [ ! -f "{params.mutations_tsv}" ]; then
            echo "ERROR: Archivo de mutaciones conocidas no encontrado: {params.mutations_tsv}" >> {log} 2>&1
            exit 1
        fi

        if [ ! -f "{params.script_path}" ]; then
            echo "ERROR: Script de análisis no encontrado: {params.script_path}" >> {log} 2>&1
            exit 1
        fi

        # Ejecutar análisis de resistencia
        python3 {params.script_path} \\
            --ref_proteins {params.ref_proteins} \\
            --assembly {input.assembly} \\
            --bam {input.bam} \\
            --mutations_tsv {params.mutations_tsv} \\
            --output_dir {params.outdir} >> {log} 2>&1

        # Verificar archivos de salida
        if [ -f "{output.all_mutations}" ] && [ -f "{output.resistance_mutations}" ]; then
            echo "Análisis de resistencia completado exitosamente" >> {log} 2>&1

            # Crear resumen
            echo "# Resumen del análisis de resistencia para {wildcards.sample}" > {output.summary_report}
            echo "Fecha: $(date)" >> {output.summary_report}
            echo "Ensamblado: {input.assembly}" >> {output.summary_report}
            echo "" >> {output.summary_report}

            # Contar mutaciones
            total_mutations=$(tail -n +2 {output.all_mutations} | wc -l)
            resistance_mutations=$(tail -n +2 {output.resistance_mutations} | wc -l)

            echo "Total de mutaciones detectadas: $total_mutations" >> {output.summary_report}
            echo "Mutaciones de resistencia conocidas: $resistance_mutations" >> {output.summary_report}
            echo "" >> {output.summary_report}

            if [ "$resistance_mutations" -gt 0 ]; then
                echo "## Mutaciones de resistencia encontradas:" >> {output.summary_report}
                tail -n +2 {output.resistance_mutations} | \\
                awk -F',' '{{print "- " $1 " " $2 ": " $3 " -> " $4 " (" $5 ")"}}' >> {output.summary_report}
            else
                echo "No se encontraron mutaciones de resistencia conocidas." >> {output.summary_report}
            fi

            echo "" >> {output.summary_report}
            echo "## Cobertura por gen:" >> {output.summary_report}
            if [ -f "{output.coverage_report}" ]; then
                tail -n +2 {output.coverage_report} | \\
                awk -F'\\t' '{{printf "- %s: %.2fx\\n", $1, $6}}' >> {output.summary_report}
            fi

            echo "Archivos de salida generados:" >> {log} 2>&1
            echo "- Todas las mutaciones: {output.all_mutations}" >> {log} 2>&1
            echo "- Mutaciones de resistencia: {output.resistance_mutations}" >> {log} 2>&1
            echo "- Reporte de cobertura: {output.coverage_report}" >> {log} 2>&1
        else
            echo "ERROR: El script de análisis no generó los archivos esperados" >> {log} 2>&1
            exit 1
        fi

        echo "=== ANÁLISIS DE RESISTENCIA COMPLETADO ===" >> {log} 2>&1
        """

# ==============================
# REGLAS DE CONSOLIDACIÓN DE RESISTENCIA
# ==============================

def get_characterization_resistance_inputs(wildcards):
    """Obtener inputs de caracterización de resistencia basados en muestras validadas"""
    try:
        checkpoint_output = checkpoints.samples_validation.get(**wildcards).output[0]

        with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
            validated_samples = [line.strip() for line in f if line.strip()]

        inputs = []
        inputs.extend(expand("output/03_characterization/03.3_chroquetas/{sample}/{sample}.ChroQueTaS.AMR_summary.txt",
                           sample=validated_samples))
        inputs.extend(expand("output/03_characterization/03.4_resistance/{sample}/resistance_analysis_summary.txt",
                           sample=validated_samples))

        return inputs
    except:
        return []

rule consolidate_resistance_characterization:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        characterizations=get_characterization_resistance_inputs
    output:
        chroquetas_summary="output/03_characterization/chroquetas_summary.tsv",
        resistance_summary="output/03_characterization/resistance_analysis_summary.tsv",
        resistance_report="output/03_characterization/resistance_characterization_report.txt"
    log:
        "output/logs/03_characterization/consolidate_resistance.log"
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
            validated_samples = [line.strip() for line in f if line.strip()]

        with open(log[0], "w") as log_file:
            log_file.write("=== CONSOLIDANDO CARACTERIZACIÓN DE RESISTENCIA ===\n")

            # Consolidar resultados de ChroQueTas
            chroquetas_data = []
            for sample in validated_samples:
                chroquetas_file = f"output/03_characterization/03.3_chroquetas/{sample}/{sample}.ChroQueTaS.AMR_summary.txt"
                if os.path.exists(chroquetas_file):
                    try:
                        with open(chroquetas_file, 'r') as f:
                            lines = f.readlines()
                            mutations_count = len([l for l in lines if not l.startswith('#') and l.strip()])
                            chroquetas_data.append({
                                'Sample': sample,
                                'ChroQueTas_Mutations': mutations_count,
                                'Status': 'Completed'
                            })
                            log_file.write(f"ChroQueTas {sample}: {mutations_count} mutaciones\n")
                    except Exception as e:
                        log_file.write(f"Error procesando ChroQueTas para {sample}: {e}\n")
                        chroquetas_data.append({
                            'Sample': sample,
                            'ChroQueTas_Mutations': 'Error',
                            'Status': 'Failed'
                        })
                else:
                    log_file.write(f"Archivo ChroQueTas no encontrado para {sample}\n")
                    chroquetas_data.append({
                        'Sample': sample,
                        'ChroQueTas_Mutations': 'Not_Found',
                        'Status': 'Missing'
                    })

            # Guardar resumen ChroQueTas
            if chroquetas_data:
                chroquetas_df = pd.DataFrame(chroquetas_data)
                chroquetas_df.to_csv(output.chroquetas_summary, sep='\t', index=False)
                log_file.write(f"Resumen ChroQueTas guardado: {len(chroquetas_data)} muestras\n")

            # Consolidar resultados de análisis de resistencia Python
            resistance_data = []
            for sample in validated_samples:
                resistance_file = f"output/03_characterization/03.4_resistance/{sample}/mutaciones_de_resistencia.csv"
                if os.path.exists(resistance_file):
                    try:
                        df = pd.read_csv(resistance_file)
                        resistance_count = len(df) if not df.empty else 0
                        resistance_data.append({
                            'Sample': sample,
                            'Python_Resistance_Mutations': resistance_count,
                            'Status': 'Completed'
                        })
                        log_file.write(f"Análisis Python {sample}: {resistance_count} mutaciones de resistencia\n")
                    except Exception as e:
                        log_file.write(f"Error procesando análisis Python para {sample}: {e}\n")
                        resistance_data.append({
                            'Sample': sample,
                            'Python_Resistance_Mutations': 'Error',
                            'Status': 'Failed'
                        })
                else:
                    log_file.write(f"Archivo de resistencia Python no encontrado para {sample}\n")
                    resistance_data.append({
                        'Sample': sample,
                        'Python_Resistance_Mutations': 'Not_Found',
                        'Status': 'Missing'
                    })

            # Guardar resumen de análisis de resistencia
            if resistance_data:
                resistance_df = pd.DataFrame(resistance_data)
                resistance_df.to_csv(output.resistance_summary, sep='\t', index=False)
                log_file.write(f"Resumen análisis de resistencia guardado: {len(resistance_data)} muestras\n")

            log_file.write("=== CONSOLIDACIÓN COMPLETADA ===\n")

        # Crear reporte de resistencia
        with open(output.resistance_report, "w") as report:
            report.write("# Reporte de Caracterización de Resistencia\n\n")
            report.write(f"Muestras analizadas: {len(validated_samples)}\n")
            report.write(f"ChroQueTas completados: {len([d for d in chroquetas_data if d['Status'] == 'Completed'])}\n")
            report.write(f"Análisis Python completados: {len([d for d in resistance_data if d['Status'] == 'Completed'])}\n\n")

            if chroquetas_data:
                report.write("## Resumen ChroQueTas\n")
                for data in chroquetas_data:
                    if data['Status'] == 'Completed':
                        report.write(f"- {data['Sample']}: {data['ChroQueTas_Mutations']} mutaciones\n")
                report.write("\n")

            if resistance_data:
                report.write("## Resumen Análisis de Resistencia Python\n")
                for data in resistance_data:
                    if data['Status'] == 'Completed':
                        report.write(f"- {data['Sample']}: {data['Python_Resistance_Mutations']} mutaciones de resistencia\n")

