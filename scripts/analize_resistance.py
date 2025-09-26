#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO

# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================

def run_command(command, shell=False):
    """Ejecuta un comando externo y gestiona errores."""
    try:
        # Usamos shell=True si el comando es una cadena, lo que es más simple
        # para comandos con pipes o redirecciones, aunque aquí lo mantenemos simple.
        if isinstance(command, str):
            shell = True
        
        process = subprocess.run(
            command,
            shell=shell,
            check=True,
            capture_output=True,
            text=True
        )
        return process.stdout
    except subprocess.CalledProcessError as e:
        print(f"ERROR al ejecutar el comando: {' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd}", file=sys.stderr)
        print(f"Código de salida: {e.returncode}", file=sys.stderr)
        print(f"Salida de error (stderr):\n{e.stderr}", file=sys.stderr)
        sys.exit(1) # Detiene el script si un comando falla

def translate_dna(dna_sequence, table_id=12):
    """Traduce una secuencia de ADN a proteína usando la tabla correcta."""
    return dna_sequence.translate(table=table_id)


# =============================================================================
# FUNCIONES PRINCIPALES DEL PIPELINE
# =============================================================================
def process_gene(ref_protein_path, assembly_db_path, assembly_path, gene_output_dir):
    """Procesa un único gen: blast, extracción, traducción y alineamiento."""
    gene_name = ref_protein_path.stem.split('_XP_')[0]
    print(f"  -> Procesando gen: {gene_name}...")

    # Definir rutas de ficheros
    blast_out = gene_output_dir / f"{gene_name}.blast.tsv"
    bed_file = gene_output_dir / f"{gene_name}.bed"
    extracted_dna = gene_output_dir / f"{gene_name}_extracted.fna"
    extracted_protein = gene_output_dir / f"{gene_name}_extracted.faa"
    mafft_aligned = gene_output_dir / f"{gene_name}_mafft_aligned.fasta"
    needle_aligned = gene_output_dir.parent / "visual_alignments" / f"{gene_name}_visual_alignment.txt"

    # 1. tblastn
    tblastn_cmd = f"tblastn -query {ref_protein_path} -db {assembly_db_path} -outfmt '6 sseqid sstart send' -max_target_seqs 1"
    blast_result = run_command(tblastn_cmd)
    
    if not blast_result.strip():
        print(f"     ¡ADVERTENCIA! No se encontró ninguna coincidencia BLAST para {gene_name}.")
        return None

    # 2. Crear fichero BED
    best_hit_line = blast_result.strip().split('\n')[0]
    contig, sstart, send = best_hit_line.split('\t')
    sstart, send = int(sstart), int(send)
    
    if sstart < send:
        strand, start, end = "+", sstart - 1, send
    else:
        strand, start, end = "-", send - 1, sstart
    
    with open(bed_file, 'w') as f:
        f.write(f"{contig}\t{start}\t{end}\t{gene_name}\t0\t{strand}\n")

    # 3. bedtools getfasta
    bedtools_cmd = f"bedtools getfasta -fi {assembly_path} -bed {bed_file} -s -fo {extracted_dna}"
    run_command(bedtools_cmd)

    # 4. Traducir (integrado en Python)
    dna_record = SeqIO.read(extracted_dna, "fasta")
    protein_seq = translate_dna(dna_record.seq)
    protein_record = SeqIO.SeqRecord(protein_seq, id=f"{gene_name}_from_assembly", description="")
    SeqIO.write(protein_record, extracted_protein, "fasta")

    # 5. Alinear (mafft y needle)
    # --- MODIFICACIÓN CLAVE ---
    # Usamos 'cat' y una tubería '|' para pasar los datos a mafft. El '-' final le indica
    # a mafft que lea de la entrada estándar, lo cual es mucho más robusto.
    mafft_cmd = f"cat {ref_protein_path} {extracted_protein} | mafft --auto -"
    mafft_result = run_command(mafft_cmd)
    with open(mafft_aligned, 'w') as f:
        f.write(mafft_result)

    needle_cmd = f"needle -asequence {ref_protein_path} -bsequence {extracted_protein} -outfile {needle_aligned} -gapopen 10.0 -gapextend 0.5"
    run_command(needle_cmd)

    return bed_file


def analyze_mutations(results_dir, reports_dir, mutations_tsv):
    """Analiza los alineamientos para encontrar y reportar mutaciones."""
    print("\nPASO 4: Identificando mutaciones y comparando con base de datos de resistencias...")
    
    try:
        # Usamos r'\s+' para evitar la SyntaxWarning
        known_muts_df = pd.read_csv(mutations_tsv, sep=r'\s+', header=0)
        known_muts_df['Mutation_Code'] = known_muts_df['Ref_AA'] + known_muts_df['Position'].astype(str) + known_muts_df['Mut_AA']
        print(f"  -> Se cargaron {len(known_muts_df)} mutaciones de resistencia desde {mutations_tsv}")
    except Exception as e:
        print(f"  Error al cargar el fichero de mutaciones: {e}")
        known_muts_df = pd.DataFrame()

    all_mutations_found = []
    alignment_files = list(results_dir.glob("*_mafft_aligned.fasta"))

    for align_file in alignment_files:
        gene_name = align_file.stem.split('_mafft_aligned')[0]
        alignment = AlignIO.read(align_file, "fasta")
        ref_seq, query_seq = alignment[0], alignment[1]
        ref_pos = 0

        for i, ref_aa in enumerate(ref_seq.seq):
            if ref_aa != '-':
                ref_pos += 1
            query_aa = query_seq.seq[i]
            
            if ref_aa != query_aa and ref_aa != '-' and query_aa != '-':
                mutation_code = f"{ref_aa}{ref_pos}{query_aa}"
                is_known = known_muts_df[(known_muts_df['Gene'] == gene_name) & (known_muts_df['Mutation_Code'] == mutation_code)]
                significance = "Resistencia Conocida" if not is_known.empty else "Desconocida"

                all_mutations_found.append({
                    'Gene': gene_name, 'Mutation': mutation_code, 'Position': ref_pos,
                    'Ref_AA': ref_aa, 'Query_AA': query_aa, 'Significance': significance
                })

    if not all_mutations_found:
        print("  -> No se encontraron mutaciones en los genes analizados.")
        return
    
    # Reporte completo
    report_df = pd.DataFrame(all_mutations_found)
    report_path = reports_dir / 'todas_las_mutaciones.csv'
    report_df.to_csv(report_path, index=False)
    print(f"  -> Reporte con TODAS las mutaciones guardado en: {report_path}")

    # Reporte de resistencia
    resistance_df = report_df[report_df['Significance'] == 'Resistencia Conocida']
    if not resistance_df.empty:
        resistance_path = reports_dir / 'mutaciones_de_resistencia.csv'
        resistance_df.to_csv(resistance_path, index=False)
        print("  \033[1;31m⚠️  ¡ALERTA! Se encontraron mutaciones de resistencia conocidas.\033[0m")
        print(f"     -> Detalles en: {resistance_path}")
    else:
        print("  -> No se encontraron mutaciones de resistencia conocidas.")

# =============================================================================
# SCRIPT PRINCIPAL
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline completo para identificar genes de resistencia y mutaciones en un ensamblado de C. auris.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-r", "--ref_proteins", required=True, help="Fichero multi-FASTA con las 5 proteínas de referencia.")
    parser.add_argument("-a", "--assembly", required=True, help="Fichero FASTA del ensamblado del genoma a analizar.")
    parser.add_argument("-b", "--bam", required=True, help="Fichero BAM alineado contra el ensamblado (debe estar indexado).")
    parser.add_argument("-m", "--mutations_tsv", required=True, help="Fichero TSV con las mutaciones de resistencia conocidas.")
    parser.add_argument("-o", "--output_dir", default="analisis_resistencia_py", help="Directorio donde se guardarán todos los resultados.")
    
    args = parser.parse_args()

    # --- PASO 0: SETUP Y COMPROBACIONES ---
    print("### INICIANDO PIPELINE DE ANÁLISIS DE RESISTENCIA (Python) ###")
    
    # Crear estructura de directorios
    output_dir = Path(args.output_dir)
    ref_split_dir = output_dir / "ref_split"
    results_dir = output_dir / "results"
    reports_dir = output_dir / "reports"
    visual_dir = output_dir / "visual_alignments"
    for d in [ref_split_dir, results_dir, reports_dir, visual_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    print(f"-> Los resultados se guardarán en: {output_dir.resolve()}")
    
    # --- PASO 1: PREPARACIÓN ---
    print("\nPASO 1: Preparando la base de datos BLAST y las referencias...")
    assembly_db_path = output_dir / "assembly_db"
    run_command(f"makeblastdb -in {args.assembly} -dbtype nucl -out {assembly_db_path}")

    ref_protein_paths = []
    for record in SeqIO.parse(args.ref_proteins, "fasta"):
        sanitized_id = record.id.replace("|", "_").replace(".", "_")
        out_path = ref_split_dir / f"{sanitized_id}.faa"
        SeqIO.write(record, out_path, "fasta")
        ref_protein_paths.append(out_path)

    # --- PASO 2: ANÁLISIS POR GEN ---
    print("\nPASO 2: Procesando cada gen individualmente...")
    all_bed_files = []
    for ref_path in ref_protein_paths:
        bed_file = process_gene(ref_path, assembly_db_path, args.assembly, results_dir)
        if bed_file:
            all_bed_files.append(bed_file)

    # --- PASO 3: ANÁLISIS DE COBERTURA ---
    print("\nPASO 3: Calculando profundidad de cobertura para los genes...")
    if all_bed_files:
        all_genes_bed = results_dir / "all_genes.bed"
        with open(all_genes_bed, 'w') as outfile:
            for fname in all_bed_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
        
        coverage_prefix = reports_dir / "coverage"
        run_command(f"mosdepth --threads 4 -b {all_genes_bed} {coverage_prefix} {args.bam}")
        
        coverage_report_file = reports_dir / "coverage_report.tsv"
        mosdepth_output = coverage_prefix.with_suffix(".regions.bed.gz")
        
        awk_cmd = f"zcat {mosdepth_output} | awk 'BEGIN {{OFS=\"\\t\"}} {{len=$3-$2; print $4, $1, $2, $3, len, $5}}'"
        coverage_data = run_command(awk_cmd)
        
        with open(coverage_report_file, 'w') as f:
            f.write("Gene\tContig\tStart\tEnd\tLength\tMean_Coverage\n")
            f.write(coverage_data)
        
        print(f"  -> Reporte de cobertura generado en: {coverage_report_file}")
        print("  -----------------------------------------------------------------")
        run_command(f"column -t -s $'\t' {coverage_report_file}")
        print("  -----------------------------------------------------------------")
    else:
        print("  No se encontraron genes para calcular la cobertura.")

    # --- PASO 4: ANÁLISIS DE MUTACIONES ---
    analyze_mutations(results_dir, reports_dir, args.mutations_tsv)

    print("\n### ANÁLISIS FINALIZADO ###")

if __name__ == "__main__":
    main()
