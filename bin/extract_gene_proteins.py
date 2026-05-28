#!/usr/bin/env python3
"""
extract_gene_proteins.py — assembly-based per-gene protein extraction +
alignment + FungAMR classification.

For each panel gene, we ask:
  1. What does the protein look like in THIS sample's polished assembly?
  2. How does it differ from the canonical reference protein?
  3. Are those mismatches FungAMR-known resistance mutations?

We deliberately do NOT use the GATK4 VCF here: in C. auris smoke testing,
GATK4 filters out known AMR mutations (e.g. Y132F in Cyp51/ERG11) due to
stringent priors or mask exclusions, while the polished assembly carries
them faithfully — which is also why ChroQueTaS works directly on the
assembly. This is the same logic that drove the original
analize_resistance.py (tblastn → bedtools getfasta → translate → mafft).

Pipeline per (sample, gene):
  1. miniprot --gff canonical_<gene>.faa  vs  sample_assembly.fasta
     → coords (contig, start, end, strand) + CIGAR
  2. samtools faidx sample_assembly.fasta region → DNA
     (reverse-complement on strand '-')
  3. translate (genetic code 12 — Yeast Alternative — for Candida)
  4. Bio.Align.PairwiseAligner global, BLOSUM62, gap_open -10, gap_extend -1
  5. Walk alignment, emit one row per mismatching residue
  6. Cross-reference each mismatch with the FungAMR mutation catalog
     (gene, ref_aa+pos+sample_aa exact match) → annotate drugs, evidence

Outputs (under --outdir):
  <sample>.gene_proteins.tsv          long format, one row per mismatch
  <sample>.gene_proteins_summary.tsv  one row per (sample, gene)
  alignments/<sample>__<gene>.aa_alignment.fasta  2-seq fasta for plotting
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Align import PairwiseAligner, substitution_matrices
except ImportError as e:
    print(f"[extract_gene_proteins] ERROR: biopython is required ({e})", file=sys.stderr)
    sys.exit(2)


GENETIC_CODE_CANDIDA = 12   # Yeast Alternative: CUG → Ser (critical for Candida)


def need(tool: str) -> str:
    p = shutil.which(tool)
    if not p:
        sys.exit(f"[extract_gene_proteins] ERROR: '{tool}' not in PATH")
    return p


def load_panel_genes(panel_path: str) -> set[str]:
    out: set[str] = set()
    with open(panel_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("gene") or "").strip()
            if g and not g.startswith("#") and g != "Chr5_full":
                out.add(g)
    return out


def load_canonical_protein(faa_path: Path) -> str:
    rec = next(SeqIO.parse(str(faa_path), "fasta"))
    return str(rec.seq).rstrip("*")


def load_catalog(path: str) -> dict[tuple[str, str], dict]:
    out: dict[tuple[str, str], dict] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("gene") or "").strip().upper()
            m = (row.get("mutation_id") or "").strip()
            if g and m:
                out[(g, m)] = row
    return out


def miniprot_locate(assembly: str, query_faa: Path, threads: int = 4
                    ) -> tuple[str, int, int, str] | None:
    """Run miniprot, return the best mRNA hit (chrom, start, end, strand)
    in 1-based GFF coords."""
    try:
        res = subprocess.run(
            ["miniprot", "--gff", "-t", str(threads), assembly, str(query_faa)],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"[extract] miniprot failed for {query_faa.stem}: {e.stderr[:200]}",
              file=sys.stderr)
        return None
    best: tuple[str, int, int, str] | None = None
    best_score = -1
    for line in res.stdout.splitlines():
        if "\tmRNA\t" not in line:
            continue
        parts = line.split("\t")
        chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
        # Use 'Identity=' attribute as tie-breaker
        ident = 0.0
        for attr in parts[8].split(";"):
            if attr.startswith("Identity="):
                try:
                    ident = float(attr.split("=", 1)[1])
                except ValueError:
                    pass
        score = ident
        if score > best_score:
            best_score = score
            best = (chrom, start, end, strand)
    return best


def samtools_faidx(assembly: str, chrom: str, start: int, end: int) -> str | None:
    """Return uppercase DNA for a 1-based inclusive region."""
    try:
        res = subprocess.run(
            ["samtools", "faidx", assembly, f"{chrom}:{start}-{end}"],
            capture_output=True, text=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"[extract] samtools faidx failed for {chrom}:{start}-{end}: "
              f"{e.stderr[:200]}", file=sys.stderr)
        return None
    return "".join(ln for ln in res.stdout.splitlines() if not ln.startswith(">")).upper()


def revcomp(dna: str) -> str:
    return str(Seq(dna).reverse_complement())


def translate(dna: str, table: int = GENETIC_CODE_CANDIDA) -> str:
    trim = (len(dna) // 3) * 3
    s = Seq(dna[:trim]).translate(table=table, to_stop=False)
    return str(s).rstrip("*")


def align_proteins(ref_aa: str, sample_aa: str):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    alns = aligner.align(ref_aa, sample_aa)
    return alns[0] if alns else None


def walk_alignment(alignment) -> list[dict]:
    """Bio.Align: alignment[0]/[1] are the aligned strings (with '-' gaps).
    Returns one dict per non-identical position."""
    ref_aligned = str(alignment[0])
    sample_aligned = str(alignment[1])
    rows = []
    ref_pos = 0
    for ra, sa in zip(ref_aligned, sample_aligned):
        if ra != "-":
            ref_pos += 1
        if ra == sa:
            continue
        if ra == "-":
            cls = "insertion"
        elif sa == "-":
            cls = "deletion"
        elif sa == "*" or sa == "X":
            cls = "premature_stop"
        else:
            cls = "missense"
        rows.append({"ref_pos": ref_pos, "ref_aa": ra, "sample_aa": sa, "class": cls})
    return rows


def classify_against_catalog(mm: dict, gene: str,
                             catalog: dict[tuple[str, str], dict]) -> dict:
    if mm["class"] == "missense":
        mut_id = f"{mm['ref_aa']}{mm['ref_pos']}{mm['sample_aa']}"
    elif mm["class"] == "premature_stop":
        mut_id = f"{mm['ref_aa']}{mm['ref_pos']}X"
    else:
        mut_id = ""
    cat = catalog.get((gene.upper(), mut_id)) if mut_id else None
    out = {
        "mutation_id": mut_id,
        "fungamr_match": "yes" if cat else "no",
        "drugs": cat["drugs"] if cat else "",
        "evidence_strength": cat["evidence_strength"] if cat else "",
        "confidence_score": cat["confidence_score"] if cat else "",
        "companion_mutations": cat["companion_mutations"] if cat else "",
        "n_reports_fungamr": cat["n_reports"] if cat else "0",
    }
    if cat and (cat.get("evidence_strength") or "").startswith("R"):
        out["classification"] = "known_resistance"
    elif cat and (cat.get("evidence_strength") or "").startswith("S"):
        out["classification"] = "known_sensitivity"
    elif mm["class"] == "premature_stop":
        out["classification"] = "frameshift_or_LoF"
    elif mm["class"] in ("insertion", "deletion"):
        out["classification"] = "frameshift_or_LoF"
    else:
        out["classification"] = "unknown_missense"
    return out


def write_alignment_fasta(out_path: Path, sample_id: str, gene: str,
                          ref_aa: str, sample_aa: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(f">{gene}__reference\n")
        for i in range(0, len(ref_aa), 60):
            fh.write(ref_aa[i:i+60] + "\n")
        fh.write(f">{gene}__{sample_id}\n")
        for i in range(0, len(sample_aa), 60):
            fh.write(sample_aa[i:i+60] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--assembly", required=True,
                    help="Polished sample assembly fasta (post-Polypolish+PyPolca).")
    ap.add_argument("--canonical_proteins_dir", required=True,
                    help="Directory with <GENE>.faa files (assets/cnv_loci/canonical_proteins/).")
    ap.add_argument("--panel", required=True, help="assets/cnv_loci/panel.tsv")
    ap.add_argument("--mutation_catalog", required=True,
                    help="assets/cnv_loci/mutation_catalog_auris.tsv")
    ap.add_argument("--genetic_code", type=int, default=GENETIC_CODE_CANDIDA)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    need("miniprot")
    need("samtools")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    aln_dir = outdir / "alignments"
    aln_dir.mkdir(parents=True, exist_ok=True)

    panel_genes = load_panel_genes(args.panel)
    catalog = load_catalog(args.mutation_catalog)
    cp_dir = Path(args.canonical_proteins_dir)

    # Pre-index the assembly once
    if not Path(f"{args.assembly}.fai").exists():
        subprocess.run(["samtools", "faidx", args.assembly], check=True)

    long_rows: list[dict] = []
    summary_rows: list[dict] = []

    for gene in sorted(panel_genes):
        canonical_path = cp_dir / f"{gene}.faa"
        if not canonical_path.exists():
            continue
        ref_aa = load_canonical_protein(canonical_path)

        loc = miniprot_locate(args.assembly, canonical_path, args.threads)
        if not loc:
            print(f"[extract] no miniprot hit for {gene}", file=sys.stderr)
            summary_rows.append({
                "sample_id": args.sample_id, "gene": gene,
                "ref_aa_length": len(ref_aa), "sample_aa_length": 0,
                "n_mismatches": "", "n_known_resistance": "",
                "n_known_sensitivity": "", "n_unknown_missense": "",
                "alignment_pct_identity": 0.0,
                "miniprot_status": "no_hit",
            })
            continue
        chrom, start, end, strand = loc

        dna = samtools_faidx(args.assembly, chrom, start, end)
        if dna is None or len(dna) < 30:
            print(f"[extract] no DNA for {gene} ({chrom}:{start}-{end})",
                  file=sys.stderr)
            continue
        if strand == "-":
            dna = revcomp(dna)
        sample_aa = translate(dna, args.genetic_code)

        aln = align_proteins(ref_aa, sample_aa)
        if aln is None:
            continue
        mismatches = walk_alignment(aln)

        n_known_R = n_known_S = n_unknown = 0
        for mm in mismatches:
            cls = classify_against_catalog(mm, gene, catalog)
            row = {
                "sample_id": args.sample_id, "gene": gene,
                "ref_pos": mm["ref_pos"], "ref_aa": mm["ref_aa"],
                "sample_aa": mm["sample_aa"],
                "alignment_class": mm["class"],
                **cls,
            }
            long_rows.append(row)
            if cls["classification"] == "known_resistance":
                n_known_R += 1
            elif cls["classification"] == "known_sensitivity":
                n_known_S += 1
            elif cls["classification"] == "unknown_missense":
                n_unknown += 1

        summary_rows.append({
            "sample_id": args.sample_id, "gene": gene,
            "ref_aa_length": len(ref_aa),
            "sample_aa_length": len(sample_aa),
            "n_mismatches": len(mismatches),
            "n_known_resistance": n_known_R,
            "n_known_sensitivity": n_known_S,
            "n_unknown_missense": n_unknown,
            "alignment_pct_identity": (
                round(100.0 * (len(ref_aa) - len(mismatches)) / len(ref_aa), 2)
                if ref_aa else 0.0),
            "miniprot_status": f"{chrom}:{start}-{end}({strand})",
        })

        write_alignment_fasta(
            aln_dir / f"{args.sample_id}__{gene}.aa_alignment.fasta",
            args.sample_id, gene, ref_aa, sample_aa)

    long_path = outdir / f"{args.sample_id}.gene_proteins.tsv"
    cols_long = ["sample_id", "gene", "ref_pos", "ref_aa", "sample_aa",
                 "alignment_class", "mutation_id", "fungamr_match", "drugs",
                 "evidence_strength", "confidence_score", "companion_mutations",
                 "n_reports_fungamr", "classification"]
    with open(long_path, "w") as fh:
        fh.write("\t".join(cols_long) + "\n")
        for r in long_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols_long) + "\n")

    summary_path = outdir / f"{args.sample_id}.gene_proteins_summary.tsv"
    cols_sum = ["sample_id", "gene", "ref_aa_length", "sample_aa_length",
                "n_mismatches", "n_known_resistance", "n_known_sensitivity",
                "n_unknown_missense", "alignment_pct_identity", "miniprot_status"]
    with open(summary_path, "w") as fh:
        fh.write("\t".join(cols_sum) + "\n")
        for r in summary_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols_sum) + "\n")

    n_R = sum(1 for r in long_rows if r["classification"] == "known_resistance")
    n_S = sum(1 for r in long_rows if r["classification"] == "known_sensitivity")
    print(f"[extract] {args.sample_id}: {len(long_rows)} mismatches across "
          f"{sum(1 for r in summary_rows if r.get('miniprot_status') != 'no_hit')} genes "
          f"(known_R={n_R}, known_S={n_S})", file=sys.stderr)
    print(f"[extract] wrote {long_path}", file=sys.stderr)
    print(f"[extract] wrote {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
