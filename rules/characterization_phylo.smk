# ==============================
# CARACTERIZACIÓN FILOGENÉTICA Y TAXONOMÍA
# Paquete 1: descarga de referencias + fastANI + auriclass
# ==============================

# Lista de referencias C. auris
CAURIS_REF_IDS = [
    "CladeI_6684", "CladeI_B8441", "CladeI_B11205",
    "CladeII_B11220", "CladeII_JCM15448",
    "CladeIII_B11221",
    "CladeIV_B11243", "CladeIV_B11245",
    "CladeV_B18474"
]

# ==============================
# REGLAS DE DESCARGA DE REFERENCIAS
# ==============================

rule download_references:
    input:
        tsv="resources/references/ncbi_info_download.tsv"
    output:
        flag=touch("resources/references/.download_complete"),
        summary="resources/references/references_summary.txt",
        reference_list="resources/references/reference_list.txt"
    conda: "../envs/epicandi_phylo.yml"
    shell:
        """
        cd resources/references

        # Leer TSV y descargar cada referencia
        tail -n +2 ncbi_info_download.tsv | while IFS=$'\t' read -r ref_id accession strain clade origin year tech; do
            echo "Descargando $ref_id ($accession)..."
            if datasets download genome accession $accession --include genome,gff3; then
                if [ -f "ncbi_dataset.zip" ]; then
                    unzip -q ncbi_dataset.zip
                    mv ncbi_dataset/data/*/*.fna $ref_id.fasta
                    mv ncbi_dataset/data/*/*.gff $ref_id.gff3
                    rm -rf ncbi_dataset ncbi_dataset.zip README.md md5sum.txt
                    echo "✓ Descargado exitosamente: $ref_id"
                else
                    echo "✗ Error: No se encontró ncbi_dataset.zip para $ref_id"
                fi
            else
                echo "✗ Error: Fallo en la descarga de $ref_id ($accession)"
            fi
            sleep 2
        done

        # Generar reference_list.txt para FastANI
        ls -1 *.fasta | while read f; do echo "$(pwd)/$f"; done > reference_list.txt

        # Generar resumen
        echo "# Genomas de referencia de C. auris disponibles" > references_summary.txt
        echo "" >> references_summary.txt
        echo "Total de referencias encontradas: $(ls -1 *.fasta | wc -l)" >> references_summary.txt
        echo "Generado: $(date)" >> references_summary.txt
        """

# ==============================
# REGLAS DE FASTANI
# ==============================

rule fastani_analysis:
    input:
        query="output/02_assembly/02.3_consensus/{sample}.fasta",
        reference_list="resources/references/reference_list.txt"
    output:
        results="output/03_characterization/03.1_fastani/{sample}/{sample}_fastani.tsv"
    params:
        outdir="output/03_characterization/03.1_fastani/{sample}"
    log:
        "output/logs/03_characterization/03.1_fastani/{sample}_fastani.log"
    benchmark:
        "output/logs/03_characterization/03.1_fastani/{sample}_fastani_benchmark.txt"
    threads: config.get("resources", {}).get("fastani", {}).get("threads", 16)
    resources:
        mem_mb=config.get("resources", {}).get("fastani", {}).get("mem", 8000),
        walltime=config.get("resources", {}).get("fastani", {}).get("walltime", "02:00:00")
    conda: "../envs/epicandi_phylo.yml"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        echo "=== EJECUTANDO FASTANI PARA {wildcards.sample} ===" > {log} 2>&1
        fastANI --version >> {log} 2>&1
        echo "Query: {input.query}" >> {log} 2>&1
        echo "Referencias: $(wc -l < {input.reference_list}) archivos" >> {log} 2>&1

        # Ejecutar fastANI
        fastANI \\
            -t {threads} \\
            -q {input.query} \\
            --rl {input.reference_list} \\
            -o {output.results} >> {log} 2>&1

        # Verificar resultados
        if [ -f "{output.results}" ]; then
            echo "FastANI completado exitosamente" >> {log} 2>&1
            echo "Resultados encontrados: $(wc -l < {output.results}) comparaciones" >> {log} 2>&1

            # Mostrar las mejores coincidencias
            echo "Mejores coincidencias (ANI > 95%):" >> {log} 2>&1
            awk '$3 > 95 {{print $2, $3}}' {output.results} | sort -k2nr | head -5 >> {log} 2>&1
        else
            echo "ERROR: FastANI no generó resultados" >> {log} 2>&1
            exit 1
        fi

        echo "=== FASTANI COMPLETADO ===" >> {log} 2>&1
        """

# ==============================
# REGLAS DE AURICLASS
# ==============================

def get_auriclass_input(wildcards):
    """Obtener input correcto para auriclass según tipo de muestra"""
    sample_type = SAMPLE_TYPES.get(wildcards.sample)

    # AuriClass funciona mejor con lecturas largas
    if sample_type in ["nanopore", "hybrid"]:
        return f"output/01_data/01.2_filtered/porechop_filtlong/{wildcards.sample}_porechop_filtlong.fastq.gz"
    else:
        # Para Illumina, usar R1 (AuriClass puede trabajar con single-end)
        return f"output/01_data/01.2_filtered/fastp/{wildcards.sample}_fastp_R1.fastq.gz"

rule auriclass_analysis:
    input:
        reads=get_auriclass_input
    output:
        results="output/03_characterization/03.2_auriclass/{sample}/{sample}_auriclass.tsv"
    params:
        outdir="output/03_characterization/03.2_auriclass/{sample}",
        sample_type=lambda wildcards: SAMPLE_TYPES.get(wildcards.sample)
    log:
        "output/logs/03_characterization/03.2_auriclass/{sample}_auriclass.log"
    benchmark:
        "output/logs/03_characterization/03.2_auriclass/{sample}_auriclass_benchmark.txt"
    threads: config.get("resources", {}).get("auriclass", {}).get("threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("auriclass", {}).get("mem", 4000),
        walltime=config.get("resources", {}).get("auriclass", {}).get("walltime", "01:00:00")
    conda: "../envs/epicandi_phylo.yml"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})

        echo "=== EJECUTANDO AURICLASS PARA {wildcards.sample} ===" > {log} 2>&1
        echo "Tipo de muestra: {params.sample_type}" >> {log} 2>&1
        echo "Input: {input.reads}" >> {log} 2>&1

        # Verificar que auriclass esté disponible
        which auriclass >> {log} 2>&1 || (echo "ERROR: auriclass no encontrado" >> {log} 2>&1 && exit 1)

        # Ejecutar auriclass
        auriclass \\
            --fastq \\
            --name {wildcards.sample} \\
            -o {output.results} \\
            {input.reads} >> {log} 2>&1

        # Verificar resultados
        if [ -f "{output.results}" ]; then
            echo "AuriClass completado exitosamente" >> {log} 2>&1
            echo "Resultados:" >> {log} 2>&1
            cat {output.results} >> {log} 2>&1
        else
            echo "ERROR: AuriClass no generó resultados" >> {log} 2>&1
            exit 1
        fi

        echo "=== AURICLASS COMPLETADO ===" >> {log} 2>&1
        """

# ==============================
# REGLAS DE CONSOLIDACIÓN Y RESUMEN
# ==============================

def get_characterization_phylo_inputs(wildcards):
    """Obtener inputs de caracterización filogenética basados en muestras validadas"""
    try:
        checkpoint_output = checkpoints.samples_validation.get(**wildcards).output[0]

        with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
            validated_samples = [line.strip() for line in f if line.strip()]

        inputs = []
        inputs.extend(expand("output/03_characterization/03.1_fastani/{sample}/{sample}_fastani.tsv",
                           sample=validated_samples))
        inputs.extend(expand("output/03_characterization/03.2_auriclass/{sample}/{sample}_auriclass.tsv",
                           sample=validated_samples))

        return inputs
    except:
        return []

rule consolidate_phylo_characterization:
    input:
        checkpoint="output/01_data/01.5_samples_pass/passed_samples.txt",
        characterizations=get_characterization_phylo_inputs
    output:
        fastani_summary="output/03_characterization/fastani_summary.tsv",
        auriclass_summary="output/03_characterization/auriclass_summary.tsv",
        phylo_report="output/03_characterization/phylogenetic_characterization_report.txt"
    log:
        "output/logs/03_characterization/consolidate_phylo.log"
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        with open("output/01_data/01.5_samples_pass/passed_samples.txt", 'r') as f:
            validated_samples = [line.strip() for line in f if line.strip()]

        with open(log[0], "w") as log_file:
            log_file.write("=== CONSOLIDANDO CARACTERIZACIÓN FILOGENÉTICA ===\n")

            # Consolidar resultados de FastANI
            fastani_data = []
            for sample in validated_samples:
                fastani_file = f"output/03_characterization/03.1_fastani/{sample}/{sample}_fastani.tsv"
                if os.path.exists(fastani_file):
                    try:
                        df = pd.read_csv(fastani_file, sep='\t', header=None,
                                       names=['query', 'reference', 'ANI', 'bidirectional_fragments', 'total_fragments'])
                        if not df.empty:
                            # Obtener la mejor coincidencia
                            best_match = df.loc[df['ANI'].idxmax()]
                            ref_name = os.path.basename(best_match['reference']).replace('.fasta', '')
                            fastani_data.append({
                                'Sample': sample,
                                'Best_Reference': ref_name,
                                'ANI': best_match['ANI'],
                                'Bidirectional_Fragments': best_match['bidirectional_fragments'],
                                'Total_Fragments': best_match['total_fragments']
                            })
                            log_file.write(f"FastANI {sample}: {ref_name} (ANI: {best_match['ANI']:.2f}%)\n")
                    except Exception as e:
                        log_file.write(f"Error procesando FastANI para {sample}: {e}\n")
                else:
                    log_file.write(f"Archivo FastANI no encontrado para {sample}\n")

            # Guardar resumen FastANI
            if fastani_data:
                fastani_df = pd.DataFrame(fastani_data)
                fastani_df.to_csv(output.fastani_summary, sep='\t', index=False)
                log_file.write(f"Resumen FastANI guardado: {len(fastani_data)} muestras\n")

            # Consolidar resultados de AuriClass
            auriclass_data = []
            for sample in validated_samples:
                auriclass_file = f"output/03_characterization/03.2_auriclass/{sample}/{sample}_auriclass.tsv"
                if os.path.exists(auriclass_file):
                    try:
                        df = pd.read_csv(auriclass_file, sep='\t')
                        if not df.empty:
                            auriclass_data.append(df.iloc[0].to_dict())
                            log_file.write(f"AuriClass {sample}: {df.iloc[0]['Clade']} (QC: {df.iloc[0]['QC_decision']})\n")
                    except Exception as e:
                        log_file.write(f"Error procesando AuriClass para {sample}: {e}\n")
                else:
                    log_file.write(f"Archivo AuriClass no encontrado para {sample}\n")

            # Guardar resumen AuriClass
            if auriclass_data:
                auriclass_df = pd.DataFrame(auriclass_data)
                auriclass_df.to_csv(output.auriclass_summary, sep='\t', index=False)
                log_file.write(f"Resumen AuriClass guardado: {len(auriclass_data)} muestras\n")

            log_file.write("=== CONSOLIDACIÓN COMPLETADA ===\n")

        # Crear reporte filogenético
        with open(output.phylo_report, "w") as report:
            report.write("# Reporte de Caracterización Filogenética\n\n")
            report.write(f"Muestras analizadas: {len(validated_samples)}\n")
            report.write(f"FastANI completados: {len(fastani_data)}\n")
            report.write(f"AuriClass completados: {len(auriclass_data)}\n\n")

            if fastani_data:
                report.write("## Resumen FastANI\n")
                for data in fastani_data:
                    report.write(f"- {data['Sample']}: {data['Best_Reference']} (ANI: {data['ANI']:.2f}%)\n")
                report.write("\n")

            if auriclass_data:
                report.write("## Resumen AuriClass\n")
                for data in auriclass_data:
                    report.write(f"- {data['Sample']}: {data['Clade']} (QC: {data['QC_decision']})\n")

