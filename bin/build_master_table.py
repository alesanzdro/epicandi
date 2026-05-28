#!/usr/bin/env python3
"""
build_master_table.py — collapse per-sample artefacts into cohort-wide tables.

Two outputs:

  master_table.tsv  one row per sample with key fields:
      sample_id, species, qc_flag, n_resistance_mutations,
      resistance_mutations, resistance_drugs_affected

  call_matrix.tsv   cohort AMR matrix sample × drug → {R, I, S}
      Drug columns are normalised to lowercase / panel-CNV nomenclature
      (fluconazole, itraconazole, amphotericin_b, echinocandins, ...) so
      this matrix plugs directly into aggregate_cnv.py for AMR↔CNV
      modulation.

Mapping ChroQueTas (pos/neg) tuple → R/I/S per drug:
  pos <= 4      → R    (strong FungAMR evidence; 1=strongest, 8=weakest)
  pos > 4       → I    (weak evidence)
  pos == NA / "" → S   (no resistance mutation)
  multiple rows: max call (R > I > S)
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


# ChroQueTas drug column → call_matrix canonical name (lowercase / panel-CNV).
DRUG_NAME_MAP = {
    "Fluconazole":   "fluconazole",
    "Itraconazole":  "itraconazole",
    "Posaconazole":  "posaconazole",
    "Voriconazole":  "voriconazole",
    "Amphotericin":  "amphotericin_b",
    "Echinocandins": "echinocandins",
    "5-Flucytosine": "5fc",
    "Caspofungin":   "caspofungin",
    "Anidulafungin": "anidulafungin",
    "Micafungin":    "micafungin",
}
CALL_ORDER = {"R": 2, "I": 1, "S": 0}


def evidence_to_call(field: str) -> str:
    """Map a 'pos/neg' tuple cell to R / I / S."""
    if field is None or not field.strip():
        return "S"
    pos_s, _, _ = field.partition("/")
    pos_s = pos_s.strip()
    if not pos_s or pos_s.upper() == "NA":
        return "S"
    try:
        pos = int(pos_s)
    except ValueError:
        return "S"
    return "R" if pos <= 4 else "I"


def max_call(a: str, b: str) -> str:
    return a if CALL_ORDER[a] >= CALL_ORDER[b] else b


def aggregate_one_sample(rows: list[dict]) -> dict:
    """For a single sample (multiple rows in resistance_report) compute
    {drug_canonical → R/I/S} taking the max call across mutations."""
    out: dict[str, str] = {d: "S" for d in DRUG_NAME_MAP.values()}
    summary: dict[str, object] = {
        "sample_id": rows[0].get("sample", ""),
        "species":   rows[0].get("species", ""),
        "n_resistance_mutations": 0,
        "resistance_mutations":   "",
        "resistance_drugs_affected": "",
    }
    mut_ids = []
    drugs_affected = set()
    for r in rows:
        mid = r.get("mutation_id", "").strip()
        if mid and mid != "NA":
            mut_ids.append(mid)
        for chroq_col, canonical in DRUG_NAME_MAP.items():
            cell = (r.get(chroq_col) or "").strip()
            call = evidence_to_call(cell)
            out[canonical] = max_call(out[canonical], call)
            if call in ("R", "I"):
                drugs_affected.add(canonical)
    summary["resistance_mutations"]      = ";".join(sorted(set(mut_ids)))
    summary["n_resistance_mutations"]    = len(set(mut_ids))
    summary["resistance_drugs_affected"] = ";".join(sorted(drugs_affected))
    return {"summary": summary, "calls": out}


def load_resistance_report(path: str) -> list[dict]:
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [r for r in reader if (r.get("sample") or "").strip()]


def load_qc_flags(paths: list[str]) -> dict[str, dict]:
    """Optional. {sample_id → {qc_flag, qc_notes}}"""
    out: dict[str, dict] = {}
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = r.get("sample_id", "").strip()
                if sid:
                    out[sid] = {"qc_flag": r.get("qc_flag", ""),
                                "qc_notes": r.get("qc_notes", "")}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resistance_reports", nargs="+", required=True,
                    help="One or more *.resistance_report.tsv files (from AMR_REPORT).")
    ap.add_argument("--qc_flags", nargs="*", default=[],
                    help="Optional: per-sample qc_flag.tsv files (from QC_GATE).")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    qc = load_qc_flags(args.qc_flags)

    summaries: list[dict] = []
    matrix: dict[str, dict] = {}    # sample_id → {drug → call}
    samples: list[str] = []

    for rp in args.resistance_reports:
        rows = load_resistance_report(rp)
        if not rows:
            print(f"[master_table] WARN: {rp} empty or unreadable", file=sys.stderr)
            continue
        agg = aggregate_one_sample(rows)
        sid = agg["summary"]["sample_id"]
        if not sid:
            print(f"[master_table] WARN: {rp} has empty sample_id", file=sys.stderr)
            continue
        samples.append(sid)
        summaries.append(agg["summary"])
        matrix[sid] = agg["calls"]

    samples.sort()
    drug_cols = list(DRUG_NAME_MAP.values())

    # call_matrix.tsv
    cm_path = outdir / "call_matrix.tsv"
    with open(cm_path, "w") as fh:
        fh.write("sample_id\t" + "\t".join(drug_cols) + "\n")
        for s in samples:
            row = [s] + [matrix[s][d] for d in drug_cols]
            fh.write("\t".join(row) + "\n")

    # master_table.tsv
    mt_path = outdir / "master_table.tsv"
    mt_cols = ["sample_id", "species", "qc_flag", "qc_notes",
               "n_resistance_mutations", "resistance_mutations",
               "resistance_drugs_affected"] + drug_cols
    summaries.sort(key=lambda r: r["sample_id"])
    with open(mt_path, "w") as fh:
        fh.write("\t".join(mt_cols) + "\n")
        for sm in summaries:
            sid = sm["sample_id"]
            row = {
                "sample_id":                 sid,
                "species":                   sm.get("species", ""),
                "qc_flag":                   qc.get(sid, {}).get("qc_flag", ""),
                "qc_notes":                  qc.get(sid, {}).get("qc_notes", ""),
                "n_resistance_mutations":    sm.get("n_resistance_mutations", 0),
                "resistance_mutations":      sm.get("resistance_mutations", ""),
                "resistance_drugs_affected": sm.get("resistance_drugs_affected", ""),
            }
            for d in drug_cols:
                row[d] = matrix[sid][d]
            fh.write("\t".join(str(row[c]) for c in mt_cols) + "\n")

    print(f"[master_table] wrote {cm_path} ({len(samples)} samples × {len(drug_cols)} drugs)",
          file=sys.stderr)
    print(f"[master_table] wrote {mt_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
