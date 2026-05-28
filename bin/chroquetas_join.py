#!/usr/bin/env python3
"""
Join ChroQueTas (FungAMR) AMR_summary.txt with the curated EpiCandi
cross-reference panel and emit a per-sample resistance_report.tsv.

ChroQueTas summary format (FungAMR / ChroQueTas v1.0.1):
    Protein  Fragment  Position_reference  AA_reference  AA_query  Fungicide_resistance
    Cyp51    1         125                 V             A         Fluconazole(NA/-2),Itraconazole(NA/-2),Posaconazole(NA/-2),Voriconazole(NA/-2)
    Cyp51    1         126                 F             L         Fluconazole(2/NA),Itraconazole(2/NA),Posaconazole(NA/-2),Voriconazole(2/NA)

Each `Drug(pos/neg)` tuple gives:
    pos = positive evidence confidence (1=strongest, 8=weakest, NA = none)
    neg = evidence the mutation does NOT confer resistance
          (always reported as a NEGATIVE integer, e.g. -2)

The output emits two complementary representations:

  Per-mutation rows (long format):
    sample, species, protein, position, aa_reference, aa_query, mutation_id
    + one column per drug class with the (pos/neg) tuple as text

  Aggregate per-sample columns (added to the first row of each sample):
    resistance_mutations           "Cyp51:F126L;Cyp51:V125A"
    resistance_classes             "Fluconazole(2/NA),...;Fluconazole(NA/-2),..."
    n_resistance_mutations         integer
    resistance_evidence_max        best (lowest) positive score; '' if none
    resistance_drugs_affected      "Fluconazole;Itraconazole;Voriconazole"
"""

import argparse
import csv
import re
import sys
from pathlib import Path


# Drug classes we explicitly surface in the long-format table. All other
# antifungals reported by ChroQueTas are still captured in
# resistance_classes (long string) but won't get a dedicated column.
DRUG_COLUMNS = [
    "Fluconazole", "Itraconazole", "Posaconazole", "Voriconazole",
    "Amphotericin", "Echinocandins", "5-Flucytosine", "Caspofungin",
    "Anidulafungin", "Micafungin",
]

# ChroQueTas (FungAMR) names that differ from the canonical column names above.
# Without this map the per-drug column for these aliases stays empty and the
# downstream call_matrix never flips R/I, even though the raw evidence is in
# the AMR_summary.txt.  Keys are the FungAMR names (case-sensitive).
DRUG_ALIASES = {
    "Amphotericin_B":   "Amphotericin",
    "5-fluorocytosine": "5-Flucytosine",
}

# Drug names may contain underscores (e.g. Amphotericin_B); allow them so the
# regex does not split the name and capture the trailing letter as a separate drug.
DRUG_TUPLE_RE = re.compile(r"([A-Za-z0-9_\-]+)\(([^)]*)\)")


def fail(msg, code=1):
    print(f"[chroquetas_join] ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample_id",  required=True)
    p.add_argument("--species",    required=True, help="Clinical name e.g. 'C. auris'")
    p.add_argument("--chroquetas", required=True, help="ChroQueTas AMR_summary.txt")
    p.add_argument("--panel",      required=True, help="epicandi_AMR_panel_unique.tsv (kept for future cross-ref)")
    p.add_argument("--rosetta",    required=True, help="candida_rosetta_v2.tsv (species name normaliser)")
    p.add_argument("--min_depth",  type=int, default=10, help="(reserved for future depth gating)")
    p.add_argument("--output",     required=True)
    return p.parse_args()


def normalize_species(clinical, rosetta_path):
    """Map any species notation to the canonical clinical name in the panel."""
    target = clinical.strip()
    target_no_underscore = target.replace("_", " ")
    candidates = {target, target_no_underscore}
    try:
        with open(rosetta_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if (row.get("clinical_name", "").strip() in candidates
                        or row.get("ncbi_official", "").strip() in candidates):
                    return row["clinical_name"].strip()
    except FileNotFoundError:
        pass
    return target


def parse_drug_field(field):
    """'Fluconazole(2/NA),Itraconazole(NA/-2),...' →
        [{drug, pos, neg, raw}, ...]
    where pos/neg are int or None (NA -> None)."""
    out = []
    for m in DRUG_TUPLE_RE.finditer(field or ""):
        drug = m.group(1).strip()
        body = m.group(2).strip()
        pos_s, sep, neg_s = body.partition("/")
        if not sep:
            neg_s = ""
        def to_int(s):
            s = s.strip()
            if not s or s.upper() == "NA":
                return None
            try:
                return int(s)
            except ValueError:
                return None
        out.append({
            "drug": drug,
            "pos":  to_int(pos_s),
            "neg":  to_int(neg_s),
            "raw":  f"{pos_s}/{neg_s}" if sep else pos_s,
        })
    return out


def parse_chroquetas(path):
    """Parse the ChroQueTas AMR_summary.txt (FungAMR format)."""
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return []
    rows = []
    with p.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            protein = (row.get("Protein") or "").strip()
            position = (row.get("Position_reference") or "").strip()
            aa_ref = (row.get("AA_reference") or "").strip()
            aa_qry = (row.get("AA_query") or "").strip()
            drugs_field = (row.get("Fungicide_resistance") or "").strip()
            if not protein or not position:
                continue
            mutation_id = f"{protein}:{aa_ref}{position}{aa_qry}"
            drugs = parse_drug_field(drugs_field)
            rows.append({
                "protein":      protein,
                "fragment":     (row.get("Fragment") or "").strip(),
                "position":     position,
                "aa_reference": aa_ref,
                "aa_query":     aa_qry,
                "mutation_id":  mutation_id,
                "drugs":        drugs,
                "raw_drugs":    drugs_field,
            })
    return rows


def aggregate(rows):
    """Build per-sample aggregate fields from the per-mutation rows."""
    if not rows:
        return {
            "resistance_mutations":      "",
            "resistance_classes":        "",
            "n_resistance_mutations":    0,
            "resistance_evidence_max":   "",
            "resistance_drugs_affected": "",
        }
    # mutation_ids in the order seen, joined
    muts = [r["mutation_id"] for r in rows]
    classes_strs = [r["raw_drugs"] for r in rows]
    # Best (lowest) positive evidence score across all mutations / drugs
    pos_scores = []
    drugs_with_resistance = set()
    for r in rows:
        for d in r["drugs"]:
            if d["pos"] is not None and d["pos"] > 0:
                pos_scores.append(d["pos"])
                drugs_with_resistance.add(d["drug"])
    return {
        "resistance_mutations":      ";".join(muts),
        "resistance_classes":        ";".join(classes_strs),
        "n_resistance_mutations":    len(rows),
        "resistance_evidence_max":   str(min(pos_scores)) if pos_scores else "",
        "resistance_drugs_affected": ";".join(sorted(drugs_with_resistance)),
    }


def main():
    args = parse_args()
    species_clinical = normalize_species(args.species, args.rosetta)
    rows = parse_chroquetas(args.chroquetas)
    agg  = aggregate(rows)

    print(f"[chroquetas_join] {args.sample_id} ({species_clinical}): "
          f"{len(rows)} mutations, "
          f"{len(agg['resistance_drugs_affected'].split(';')) if agg['resistance_drugs_affected'] else 0} "
          f"drugs affected", file=sys.stderr)

    fixed_cols = [
        "sample", "species", "protein", "position", "aa_reference", "aa_query",
        "mutation_id",
    ]
    drug_cols = list(DRUG_COLUMNS)
    agg_cols = [
        "resistance_mutations", "resistance_classes", "n_resistance_mutations",
        "resistance_evidence_max", "resistance_drugs_affected",
    ]
    all_cols = fixed_cols + drug_cols + agg_cols + ["note"]

    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=all_cols, delimiter="\t")
        writer.writeheader()

        if not rows:
            # Sample with no detected mutations — single row, all blank
            row = {c: "" for c in all_cols}
            row.update({
                "sample":                args.sample_id,
                "species":               species_clinical,
                "note":                  "ChroQueTas ran but reported no AMR-associated mutations",
            })
            row.update(agg)
            writer.writerow(row)
            return

        for i, r in enumerate(rows):
            row = {c: "" for c in all_cols}
            row.update({
                "sample":       args.sample_id,
                "species":      species_clinical,
                "protein":      r["protein"],
                "position":     r["position"],
                "aa_reference": r["aa_reference"],
                "aa_query":     r["aa_query"],
                "mutation_id":  r["mutation_id"],
            })
            # Per-drug tuple cells (only for drugs we explicitly surface). The
            # FungAMR name is normalised through DRUG_ALIASES so synonyms like
            # 'Amphotericin_B' or '5-fluorocytosine' land in their canonical
            # column instead of being silently dropped.
            drug_map = {DRUG_ALIASES.get(d["drug"], d["drug"]): d
                        for d in r["drugs"]}
            for col in DRUG_COLUMNS:
                if col in drug_map:
                    row[col] = drug_map[col]["raw"]
            # Aggregate fields go on the FIRST row only — leaves other rows uncluttered
            if i == 0:
                row.update(agg)
            writer.writerow(row)

    print(f"[chroquetas_join] Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
