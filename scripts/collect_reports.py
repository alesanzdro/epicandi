import os
import shutil
import re

# Snakemake provides a 'snakemake' object when running scripts via the 'script:' directive
# We expect the following to be available:
# - snakemake.params.outdir: destination base directory for collection
# - snakemake.params.samples: list of sample names
# - snakemake.params.output_dir: base 'output' directory in the workflow
# - snakemake.log[0]: log file path
# - snakemake.output.flag: flag file path

base = snakemake.params.outdir
OUTPUT_DIR = snakemake.params.output_dir
SAMPLES = snakemake.params.samples

# Añadimos subdirectorio busco en la colección
subdirs = {
    'nanostats': os.path.join(base, 'nanostats'),
    'checkm2': os.path.join(base, 'checkm2'),
    'prokka': os.path.join(base, 'prokka'),
    'kraken2': os.path.join(base, 'kraken2'),
    'busco': os.path.join(base, 'busco'),  # nuevo
}

# Coverage metrics collection directories
coverage_metrics_base = os.path.join(base, 'coverage_metrics')
coverage_dir = os.path.join(coverage_metrics_base, 'coverage')
mosdepth_dir = os.path.join(coverage_metrics_base, 'mosdepth')

# Ensure directories
for d in [base, coverage_metrics_base, coverage_dir, mosdepth_dir] + list(subdirs.values()) + [os.path.dirname(snakemake.log[0])]:
    os.makedirs(d, exist_ok=True)

with open(snakemake.log[0], 'w') as logfile:
    for sample in SAMPLES:
        # --- BUSCO normalización ---
        busco_dir = f"{OUTPUT_DIR}/02_assembly/02.9_busco/{sample}"
        # Rutas posibles de summary
        specific = os.path.join(busco_dir, f"short_summary.specific.{snakemake.config['busco_lineage']}.{sample}.txt") if 'busco_lineage' in snakemake.config else None
        basic = os.path.join(busco_dir, f"run_{snakemake.config['busco_lineage']}", "short_summary.txt") if 'busco_lineage' in snakemake.config else None
        os.makedirs(subdirs['busco'], exist_ok=True)
        target_summary = os.path.join(subdirs['busco'], f"{sample}.busco_summary.txt")
        chosen = None
        for cand in [specific, basic]:
            if cand and os.path.exists(cand) and os.path.getsize(cand) > 0:
                shutil.copy2(cand, target_summary)
                chosen = cand
                logfile.write(f"BUSCO summary for {sample} -> busco/{sample}.busco_summary.txt (from {os.path.basename(cand)})\n")
                break
        if not chosen:
            logfile.write(f"WARN: No BUSCO summary found for {sample} in {busco_dir}\n")
        # Copiar full_table.tsv si existe (sin modificar) con nombre consistente
        full_table = os.path.join(busco_dir, f"run_{snakemake.config['busco_lineage']}", "full_table.tsv") if 'busco_lineage' in snakemake.config else None
        if full_table and os.path.exists(full_table) and os.path.getsize(full_table) > 0:
            dest_full = os.path.join(subdirs['busco'], f"{sample}.full_table.tsv")
            shutil.copy2(full_table, dest_full)
            logfile.write(f"Copied BUSCO full_table for {sample} -> busco/{sample}.full_table.tsv\n")

        # Coverage per-sample (rename to {sample}.txt)
        cov_src = f"{OUTPUT_DIR}/02_assembly/02.6_coverage/{sample}/{sample}.coverage.txt"
        if os.path.exists(cov_src) and os.path.getsize(cov_src) > 0:
            cov_dest = os.path.join(coverage_dir, f"{sample}.txt")
            shutil.copy2(cov_src, cov_dest)
            logfile.write(f"Copied coverage for {sample} -> coverage_metrics/coverage/{sample}.txt\n")

        # NanoStats (solo filtrados, sin sufijo TRIM)
        filtered_file = f"{OUTPUT_DIR}/01_data/01.4_nanoplot_filtered/{sample}/NanoStats.txt"
        if os.path.exists(filtered_file):
            # Nuevo destino simplificado: {sample}.txt
            dest = f"{subdirs['nanostats']}/{sample}.txt"
            shutil.copy2(filtered_file, dest)
            logfile.write(f"Copied filtered NanoStats for {sample} -> nanostats/{sample}.txt\n")

        # CheckM2: copy and relabel 'consensus' to sample name
        checkm2_file = f"{OUTPUT_DIR}/02_assembly/02.4_checkm2/{sample}/quality_report.tsv"
        if os.path.exists(checkm2_file):
            dest = f"{subdirs['checkm2']}/{sample}_quality_report.tsv"
            try:
                with open(checkm2_file, 'r') as src:
                    content = src.read()
                content = re.sub(r"(?m)^consensus\b", sample, content)
                with open(dest, 'w') as dst:
                    dst.write(content)
                logfile.write(f"Copied and relabeled CheckM2 report for {sample}\n")
            except Exception as e:
                shutil.copy2(checkm2_file, dest)
                logfile.write(f"Copied CheckM2 report (no relabel) for {sample} due to error: {e}\n")

        # Prokka
        prokka_gff = f"{OUTPUT_DIR}/02_assembly/02.5_prokka/{sample}/annotation.gff"
        prokka_log = f"{OUTPUT_DIR}/logs/02.5_prokka_annotation/{sample}.log"
        prokka_txt = f"{OUTPUT_DIR}/02_assembly/02.5_prokka/{sample}/annotation.txt"
        if os.path.exists(prokka_gff):
            dest = f"{subdirs['prokka']}/{sample}_annotation.gff"
            shutil.copy2(prokka_gff, dest)
            logfile.write(f"Copied Prokka GFF for {sample}\n")
        if os.path.exists(prokka_log):
            dest = f"{subdirs['prokka']}/{sample}_prokka.log"
            shutil.copy2(prokka_log, dest)
            logfile.write(f"Copied Prokka log for {sample}\n")
        if os.path.exists(prokka_txt):
            dest = f"{subdirs['prokka']}/{sample}_annotation.txt"
            try:
                # Leer contenido y reconstruir la línea organism conservando Genus species
                with open(prokka_txt, 'r') as src:
                    lines = src.readlines()
                replaced = False
                for i, line in enumerate(lines):
                    m = re.match(r"^(\s*)organism:\s*(.+)$", line)
                    if m:
                        leading, rest = m.groups()
                        # Tomar las dos primeras palabras como Genus species
                        tokens = rest.strip().split()
                        if len(tokens) >= 2:
                            genus_species = " ".join(tokens[:2])
                        else:
                            genus_species = rest.strip()
                        lines[i] = f"{leading}organism: {genus_species} {sample}\n"
                        replaced = True
                        break
                with open(dest, 'w') as dst:
                    dst.writelines(lines)
                if replaced:
                    logfile.write(f"Copied Prokka TXT summary for {sample} (organism genus/species preserved, sample appended)\n")
                else:
                    logfile.write(f"Copied Prokka TXT summary for {sample} (no organism line found)\n")
            except Exception as e:
                shutil.copy2(prokka_txt, dest)
                logfile.write(f"Copied Prokka TXT summary for {sample} (raw copy due to error: {e})\n")

        # QUAST: do not copy TSV; MultiQC will scan the original quast dirs
        # This section intentionally left minimal
        # Copiar report.tsv a la colección para que MultiQC lo encuentre
        #quast_report_tsv = f"{OUTPUT_DIR}/02_assembly/02.3_quast/{sample}/report.tsv"
        #if os.path.exists(quast_report_tsv) and os.path.getsize(quast_report_tsv) > 0:
        #    # 1. Asegurarse de que el directorio base 'quast' existe.
        #    #    Ya no creamos un subdirectorio por muestra.
        #    os.makedirs(subdirs['quast'], exist_ok=True)
        #    # 2. Definir la ruta de destino con el nuevo formato: {sample}.tsv
        #    dest = os.path.join(subdirs['quast'], f"{sample}.tsv")
        #    # Copiar el fichero
        #    shutil.copy2(quast_report_tsv, dest)
        #    # 3. Actualizar el mensaje del log para que refleje la nueva ruta
        #    logfile.write(f"Copied QUAST TSV for {sample} -> quast/{sample}.tsv\n")

        # Kraken2 report: copiar como {sample}.txt (sin sufijo _report y sin subcarpeta)
        kraken_report = f"{OUTPUT_DIR}/03_taxonomy/03.1_kraken2/{sample}_report.txt"
        if os.path.exists(kraken_report):
            dest = os.path.join(subdirs['kraken2'], f"{sample}.txt")
            shutil.copy2(kraken_report, dest)
            logfile.write(f"Copied Kraken2 report for {sample} as {sample}.txt\n")

        # Mosdepth summary: copiar *.mosdepth.summary.txt manteniendo el nombre
        mosdepth_summary = f"{OUTPUT_DIR}/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.summary.txt"
        if os.path.exists(mosdepth_summary) and os.path.getsize(mosdepth_summary) > 0:
            # Renombrar para eliminar '.mosdepth' del nombre
            dest = os.path.join(mosdepth_dir, f"{sample}.summary.txt")
            try:
                with open(mosdepth_summary, 'r') as src, open(dest, 'w') as dst:
                    header = src.readline()
                    dst.write(header)
                    removed = 0
                    for line in src:
                        # Filtrar líneas donde el primer campo acaba en _region
                        first = line.split('\t', 1)[0]
                        if first.endswith('_region'):
                            removed += 1
                            continue
                        dst.write(line)
                logfile.write(f"Copied mosdepth summary for {sample} (filtered {removed} *_region lines) -> coverage_metrics/mosdepth/{sample}.summary.txt (cleaned name)\n")
            except Exception as e:
                shutil.copy2(mosdepth_summary, dest)
                logfile.write(f"Copied mosdepth summary for {sample} (unfiltered due to error: {e}) as {sample}.summary.txt\n")

        # Mosdepth global distribution: copiar *.mosdepth.global.dist.txt manteniendo el nombre
        mosdepth_global = f"{OUTPUT_DIR}/02_assembly/02.7_mosdepth/{sample}/{sample}.mosdepth.global.dist.txt"
        if os.path.exists(mosdepth_global) and os.path.getsize(mosdepth_global) > 0:
            # Renombrar para eliminar '.mosdepth' del nombre, así MultiQC captura 'sample' sin sufijo
            dest_g = os.path.join(mosdepth_dir, f"{sample}.global.dist.txt")
            shutil.copy2(mosdepth_global, dest_g)
            logfile.write(f"Copied mosdepth global.dist for {sample} -> coverage_metrics/mosdepth/{sample}.global.dist.txt (cleaned name)\n")

    # Copy consolidated coverage metrics (TSV and MultiQC JSON) into coverage_metrics root
    cov_tsv = f"{OUTPUT_DIR}/02_assembly/02.8_coverage_metrics/all_samples_coverage.tsv"
    if os.path.exists(cov_tsv) and os.path.getsize(cov_tsv) > 0:
        shutil.copy2(cov_tsv, os.path.join(coverage_metrics_base, os.path.basename(cov_tsv)))
        logfile.write("Copied consolidated TSV: all_samples_coverage.tsv\n")

    cov_json = f"{OUTPUT_DIR}/02_assembly/02.8_coverage_metrics/coverage_multiqc.json"
    if os.path.exists(cov_json) and os.path.getsize(cov_json) > 0:
        shutil.copy2(cov_json, os.path.join(coverage_metrics_base, os.path.basename(cov_json)))
        logfile.write("Copied MultiQC custom JSON: coverage_multiqc.json\n")

# Final flag
with open(snakemake.output.flag, 'w') as flag_file:
    flag_file.write("Reports collection completed\n")
